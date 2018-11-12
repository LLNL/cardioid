#include "stateLoader.hh"
#include "Anatomy.hh"
#include "BucketOfBits.hh"
#include "pio.h"
#include "ioUtils.h"
#include "object_cc.hh"
#include "mpiTpl.hh"
#include "mpiUtils.h"
#include "Vector.hh"
#include "IndexToVector.hh"
#include "IndexToThreeVector.hh"
#include "GridAssignmentObject.h"
#include "readPioFile.hh"
#include <algorithm>
#include <map>
//ddt #include <iostream>

using namespace std;

struct SortMap
{
   unsigned index;
   unsigned dest;
};
bool operator<(const SortMap& a, const SortMap& b)
{
   if (a.dest == b.dest)
      return a.index < b.index;
   return a.dest < b.dest;
}

struct RecordRequest
{
   Long64 gid;
   unsigned dest;
};

bool operator<(const RecordRequest& a, const RecordRequest& b)
{
   if (a.dest == b.dest)
      return a.gid < b.gid;
   return a.dest < b.dest;
}

/** This is an internal self test that makes sure that the number of
 *  records headed for a drop-off site and the number of requests headed
 *  for that site match. */
void testDestinations(const vector<unsigned>& recordDest,
                      const vector<unsigned>& requestDest,
                      MPI_Comm comm)
{
   int nTasks;
   int myRank;
   MPI_Comm_size(comm, &nTasks);
   MPI_Comm_rank(comm, &myRank);

   vector<unsigned> sendBuf(nTasks, 0);
   vector<int> recvCnt(nTasks, 1);

   for (unsigned ii=0; ii<recordDest.size(); ++ii)
      ++sendBuf[recordDest[ii]];

   unsigned nRecDest;
   MPI_Reduce_scatter(&sendBuf[0], &nRecDest, &recvCnt[0], MPI_INT, MPI_SUM, comm);

   sendBuf.assign(nTasks, 0);
   for (unsigned ii=0; ii<requestDest.size(); ++ii)
      ++sendBuf[requestDest[ii]];

   unsigned nReqDest;
   MPI_Reduce_scatter(&sendBuf[0], &nReqDest, &recvCnt[0], MPI_INT, MPI_SUM, comm);

   assert(nRecDest == nReqDest);
//   cout << "Rank " << myRank << " Records " << nRecDest << " Requests " << nReqDest << endl;
//   cout << "Test Passed" << endl;
}



/** Populates the records and gid array from the named pio file.  The
 *  records array will contain the raw record data exactly as it was read
 *  with record ii starting at records[ii*lRec].  This routine does
 *  extract the gid data from the record and store it separately since
 *  the gid data is needed for other steps in the load and distribute
 *  process.
 *
 *  It should be possible to support VARRECORDASCII files by looking
 *  through the BucketOfBits created by reading the file to find the
 *  maximum record length.  Then we just need to pad out shorter records
 *  when they are copied from the BucketOfBits into the records vector
 *  for further processing. */
BucketOfBits* readStateData(const string& filename, MPI_Comm comm,
                            vector<unsigned char>& records,
                            vector<Long64>& gid,
                            unsigned& lRec)
{
   PFILE* file = Popen(filename.c_str(), "r", comm);
   assert(file);

   string recordType;
   objectGet(file->headerObject, "datatype", recordType, "unknown");
   assert(recordType == "FIXRECORDASCII" || recordType == "FIXRECORDBINARY");
   lRec = file->recordLength;
   
   BucketOfBits* bucket = readPioFile(file);

   unsigned gidIndex = bucket->getIndex("gid");
   assert(gidIndex < bucket->nFields());
   unsigned nRecords = bucket->nRecords();
   gid.resize(nRecords);
   records.reserve(nRecords*lRec);
   for (unsigned ii=0; ii<nRecords; ++ii)
   {
      BucketOfBits::Record rr = bucket->getRecord(ii);
      rr.getValue(gidIndex, gid[ii]);
      const char* src = rr.getRawData().data();
      records.push_back('\0'); 
      unsigned char* dest = &records.back();
      records.resize(records.size()+lRec-1);
      copyBytes(dest, src, lRec);
   }

   Pclose(file);

   bucket->clearRecords();
   return bucket;
}

/** There are two sets of destinations that we need to calculate.  Each
 *  task that reads records needs to find the drop off rank to which
 *  each record in the state data should be sent.  This is the
 *  recordDest.  Each task that has cell data in the Anatomy needs to
 *  calculate the drop off rank for each cell that it owns so that it
 *  can send a request for the record data.  This is the requestDest.
 *  The drop off rank depends only on the gid.  This implementation
 *  assigns each gid to the task whose center is closest (voronoi rule).
 *  In the case where two or more centers are equally close, the task
 *  with the lowest rank is chosen.  */
void findDestinations(const Anatomy& anatomy,
                      const vector<Long64>& gid,
                      MPI_Comm comm,
                      vector<unsigned>& recordDest,
                      vector<unsigned>& requestDest)
{
   // calculate domain centers and allGather
   int nTasks;
   int myRank;
   MPI_Comm_size(comm, &nTasks);
   MPI_Comm_rank(comm, &myRank);
   IndexToVector indexToVector(anatomy.nx(), anatomy.ny(), anatomy.nz());
   vector<Vector> centers(nTasks);
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
      centers[myRank] += indexToVector(anatomy.gid(ii));
   centers[myRank] /= anatomy.nLocal();
   allGather(centers, comm);
   
   // calculate destination of each input record (voronoi rule)
   GRID_ASSIGNMENT_OBJECT* gao = gao_init(centers.size(),
                                          (const void*) &(centers[0]),
                                          sizeof(Vector));
   
   IndexToThreeVector indexTo3Vector(anatomy.nx(), anatomy.ny(), anatomy.nz());
   recordDest.resize(gid.size());
   for (unsigned ii=0; ii<gid.size(); ++ii)
   {
      THREE_VECTOR r = indexTo3Vector(gid[ii]);
      recordDest[ii] = gao_nearestCenter(gao, r);
   }

   requestDest.resize(anatomy.nLocal());
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      THREE_VECTOR r = indexTo3Vector(anatomy.gid(ii));
      requestDest[ii] = gao_nearestCenter(gao, r);
   }

   testDestinations(recordDest, requestDest, comm);
   
   gao_destroy(gao);
}

/** Returns the number of records that were received on this task */
unsigned sendRecordsToDropOff(vector<unsigned>& dropTask, MPI_Comm comm,
                              const unsigned lRec,
                              vector<unsigned char>& records,
                              vector<Long64>& gid)
{
   int nTasks;
   int myRank;
   MPI_Comm_size(comm, &nTasks);
   MPI_Comm_rank(comm, &myRank);
   // Build a map that will allow us to sort the records by destination
   vector<SortMap> sortMap(dropTask.size());
   for (unsigned ii=0; ii<dropTask.size(); ++ii)
   {
      sortMap[ii].dest = dropTask[ii];
      sortMap[ii].index = ii;
   }
   sort(sortMap.begin(), sortMap.end());
   sort(dropTask.begin(), dropTask.end());
   
   // Find the number of records each task will receive.
   int nRecv;
   {
      vector<int> buf(nTasks, 0);
      for (unsigned ii=0; ii<sortMap.size(); ++ii)
         ++buf[sortMap[ii].dest];
      vector<int> recvCnt(nTasks, 1);
      MPI_Reduce_scatter(&buf[0], &nRecv, &recvCnt[0], MPI_INT, MPI_SUM, comm);
   }

   // Form virtual global array consisting of gid and record data
   // using sort key to sort records by destination
   int capacity = max(int(gid.size()), nRecv);
   int itemSize = lRec + sizeof(Long64);
   vector<unsigned char> rawData(records.begin(), records.end());
   records.resize(capacity * itemSize);
   for (unsigned ii=0; ii<gid.size(); ++ii)
   {
      unsigned iRec = sortMap[ii].index;
      void* to = &records[ii*itemSize];
      void* from = &gid[iRec];
      int size = sizeof(Long64);
      copyBytes(to, from, size);
      to = &records[ii*itemSize] + size;
      from = &rawData[iRec*lRec];
      copyBytes(to, from, lRec);
   }
   
   // Use assign array to move data
   unsigned nLocal = gid.size();
   assignArray(&records[0],
               &nLocal,
               records.size(),
               itemSize,
               &dropTask[0],
               0,
               comm);
   records.resize(nLocal*itemSize);
   return nLocal;
}

void sendRequestsToDropOff(vector<unsigned>& requestDest,
                           MPI_Comm comm,
                           const Anatomy& anatomy,
                           size_t nRecv,
                           vector<RecordRequest>& requests)
{
   int nTasks;
   int myRank;
   MPI_Comm_size(comm, &nTasks);
   MPI_Comm_rank(comm, &myRank);
   vector<SortMap> sortMap(requestDest.size());
   for (unsigned ii=0; ii<requestDest.size(); ++ii)
   {
      sortMap[ii].dest = requestDest[ii];
      sortMap[ii].index = ii;
   }
   sort(sortMap.begin(), sortMap.end());
   sort(requestDest.begin(), requestDest.end());
   
   unsigned capacity = max(nRecv, requestDest.size());
   requests.resize(capacity);
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      unsigned iRec = sortMap[ii].index;
      requests[ii].gid = anatomy.gid(iRec);
      requests[ii].dest = myRank;
   }

   unsigned nLocal = requestDest.size();
   assignArray((unsigned char*) &requests[0],
               &nLocal,
               requests.size()*sizeof(RecordRequest),
               sizeof(RecordRequest),
               &requestDest[0],
               0,
               comm);
   requests.resize(nLocal);
   assert(nLocal == nRecv);
}

void sendRecordsToRequests(const vector<RecordRequest>& requests,
                           unsigned lRec,
                           vector<unsigned char>& records,
                           const Anatomy& anatomy,
                           MPI_Comm comm)
{
   vector<SortMap> sortMap(requests.size());
   vector<unsigned> finalDest(requests.size());
   for (unsigned ii=0; ii<sortMap.size(); ++ii)
   {
      sortMap[ii].dest = requests[ii].dest;
      sortMap[ii].index = ii;
      finalDest[ii] = requests[ii].dest;
   }
   sort(sortMap.begin(), sortMap.end());
   sort(finalDest.begin(), finalDest.end());


   map<Long64, unsigned> recordMap;
   unsigned itemSize = lRec + sizeof(Long64);
   unsigned nRecords = records.size()/itemSize;
   assert(nRecords == requests.size());
   for (unsigned ii=0; ii<nRecords; ++ii)
   {
      Long64* gidPtr = (Long64*) &(records[ii*itemSize]);
      recordMap[*gidPtr] = ii;
   }

   unsigned capacity = max(records.size(), (size_t)anatomy.nLocal()*itemSize);
   vector<unsigned char> buf(capacity);

   for (unsigned ii=0; ii<requests.size(); ++ii)
   {
      unsigned iRequest = sortMap[ii].index;
      Long64 requestedGid = requests[iRequest].gid;
      unsigned iRecord = recordMap[requestedGid];

      void* from = &records[iRecord*itemSize];
      void* to = &buf[ii*itemSize];
      copyBytes(to, from, itemSize);
   }
   
   assignArray(&buf[0],
               &nRecords,
               capacity*itemSize,
               itemSize,
               &finalDest[0],
               0,
               comm);
   buf.resize(nRecords*itemSize);
   records = buf;
}


BucketOfBits* loadAndDistributeState(const std::string& filename,
                                    const Anatomy& anatomy)
{
   MPI_Comm comm = MPI_COMM_WORLD;

   vector<unsigned char> records;
   vector<Long64> gid;
   unsigned lRec=0;

   BucketOfBits* bucket = readStateData(filename, comm, // inputs
                                        records, gid, lRec); // outputs

   //ddt
//    for (unsigned ii=0; ii<gid.size(); ++ii)
//    {
//       string tmp((char*)&records[ii*lRec], lRec);
//       cout << "gid /"<< gid[ii] <<"/ : \"" << tmp <<"\""<< endl;
//    }
   


   
   vector<unsigned> recordDest;
   vector<unsigned> requestDest;
   findDestinations(anatomy, gid, comm, // inputs
                    recordDest, requestDest); // outputs

   
   unsigned nRecv = sendRecordsToDropOff(recordDest, comm, lRec,
                                         records, gid);

   // From this point on, gid array is no longer in sync with records.
   // It should not be used.  In fact, we ought to dellocate its memory.
   
   //ddt
//    for (unsigned ii=0; ii<nRecv; ++ii)
//    {
//       string tmp((char*)&records[ii*(lRec+8)], lRec+8);
//       cout << tmp << endl;
//    }
   
   vector<RecordRequest> requests;
   sendRequestsToDropOff(requestDest, comm, anatomy, nRecv,
                         requests);

   sendRecordsToRequests(requests, lRec, records, anatomy, comm);

   // At this point we should have all of our records.

   map<Long64, unsigned> recordMap;
   unsigned itemSize = lRec + sizeof(Long64);
   unsigned nRecords = records.size()/itemSize;
   assert(nRecords == anatomy.nLocal());
   for (unsigned ii=0; ii<nRecords; ++ii)
   {
      Long64* gidPtr = (Long64*) &(records[ii*itemSize]);
      recordMap[*gidPtr] = ii;
   }

   for (unsigned ii=0; ii<nRecords; ++ii)
   {
      map<Long64, unsigned>::const_iterator here;
      here = recordMap.find(anatomy.gid(ii));
      assert(here != recordMap.end());
      unsigned iRec = here->second;

      string tmp( (char*) &records[iRec*itemSize+sizeof(Long64)], lRec);
      bucket->addRecord(tmp);
   }
   

   return bucket;   
}
