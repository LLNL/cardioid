#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>

#include "Simulate.hh"
#include "PerformanceTimers.hh"
#include "checkpointIO.hh"
#include "BucketOfBits.hh"
#include "AnatomyReader.hh"
#include "stateLoader.hh"
#include "cellState.hh"
#include "object_cc.hh"
#include "readPioFile.hh"
#include "mpiUtils.h"
#include "pio.h"
#include "units.h"
#include "ioUtils.h"

// compareSnapshots goes through all snapshot subdirectories in two simulation directories and
// computes the max, min, and rms differences in all floating-point fields for plotting and
// comparison.

using namespace std;
namespace
{
    void readAnatomyBucket(Anatomy& anatomy, BucketOfBits* bucketP);
    BucketOfBits* fillDataBucket(Simulate& sim, string dirName, double &time, int &nrecord);
}

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char** argv)
{
   int nTasks, myRank;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  

   if (argc != 3)
   {
      if (myRank == 0)
         cout << "Usage:  compareSnapshots [simulation directory 1] [simulation directory 2]" << endl << endl;
      exit(1);
   }
   
   string runDir1(argv[1]);
   string runDir2(argv[2]);
   const string snapshot("snapshot");  // text that identifies snapshot directory name

   // test if runDir1 exists
   struct stat st1;
   if (!(stat(runDir1.c_str(),&st1)==0 && S_ISDIR( st1.st_mode )) )
   {
      if (myRank == 0)
         cout << "ERROR:  " << runDir1 << " does not exist!" << endl;
      exit(1);
   }

   // go through all subdirectories of runDir1, save snapshot directory names
   vector<string> snapshotDirs1;
   struct dirent *dptr;
   DIR *dir = opendir(runDir1.c_str());
   struct dirent *entry = readdir(dir);
   while (entry != NULL)
   {
      if (entry->d_type == DT_DIR)
      {
         string subdir(entry->d_name);
         size_t skipit = subdir.find(snapshot);
         if (!skipit)
            snapshotDirs1.push_back(subdir);
      }
      entry = readdir(dir);
   }
   
   // Compute the subset of snapshot directories that exist in both run directories.
   // (If runDir2 doesn't exist, this union will be empty and the code will exit without
   // doing anything.)
   vector<string> snapshotUnion;
   struct stat st2;
   for (int ii=0; ii<snapshotDirs1.size(); ii++)
   {
      string testdir2 = runDir2 + '/' + snapshotDirs1[ii];
      if (stat(testdir2.c_str(),&st2)==0 && S_ISDIR(st2.st_mode))
         snapshotUnion.push_back(snapshotDirs1[ii]);
   }

   if (snapshotUnion.size() <= 0)
   {
      if (myRank == 0)
         cout << "ERROR: no common snapshot directories in " << runDir1 << " and " << runDir1 << endl;
      exit(1);
   }

   // sort snapshotUnion
   sort(snapshotUnion.begin(),snapshotUnion.end());
   

   // compare each set of snapshots across run directories
   for (int isnap=0; isnap<snapshotUnion.size(); isnap++)
   {
      if (myRank == 0)
      {
         cout << endl;
         cout << "Comparing snapshot " << snapshotUnion[isnap] << "..." << endl;
      }
      
      string snapshotDir1 = runDir1 + '/' + snapshotUnion[isnap];
      string snapshotDir2 = runDir2 + '/' + snapshotUnion[isnap];
      
      // read state data from first directory
      Simulate sim1;
      int nrecord1;
      double time1;
      BucketOfBits* data1 = fillDataBucket(sim1,snapshotDir1,time1,nrecord1);
      
      // read state data from second directory
      Simulate sim2;
      int nrecord2;
      double time2;
      BucketOfBits* data2 = fillDataBucket(sim2,snapshotDir2,time2,nrecord2);
      

      // determine which fields both directories have in common
      vector<string> fieldUnion;
      vector<int> fieldIndex1;
      vector<int> fieldIndex2;
      int gidindex1 = -1;
      int gidindex2 = -1;
      int nFields1 = data1->nFields();
      int nFields2 = data2->nFields();
      int nFloats = 0;

      for (int ifield=0; ifield<nFields1; ifield++)
      {
         for (int jfield=0; jfield<nFields2; jfield++)
         {
            if (data1->fieldName(ifield) == data2->fieldName(jfield))
            {
               if (data1->fieldName(ifield) == "gid")
               {
                  gidindex1 = ifield;
                  gidindex2 = jfield;
               }
               else
               {
                  fieldUnion.push_back(data1->fieldName(ifield));
                  fieldIndex1.push_back(ifield);
                  fieldIndex2.push_back(jfield);
                  nFloats++;
               }
            }
         }
      }
      assert(nFloats > 0);
      
      assert(gidindex1 >= 0 && gidindex2 >= 0);
      assert(sim1.anatomy_.nx() == sim2.anatomy_.nx());
      assert(sim1.anatomy_.ny() == sim2.anatomy_.ny());
      assert(sim1.anatomy_.nz() == sim2.anatomy_.nz());
      
      const int nGlobal = sim1.anatomy_.nx()*sim1.anatomy_.ny()*sim1.anatomy_.nz();
      const int nLocal1 = sim1.anatomy_.nLocal();
      const int nLocal2 = sim2.anatomy_.nLocal();
      vector<cellState> state1(nLocal1*nFloats);
      vector<cellState> state2(nLocal2*nFloats);
      
      // distribute data across all tasks in a consistent way so we can compare state files
      // generated with different load balance algorithms
      
      // directory 1
      for (unsigned ii=0; ii<nLocal1; ++ii)
      {
         BucketOfBits::Record iRec = data1->getRecord(ii);
         int gid;
         iRec.getValue(gidindex1, gid);
         assert(gid >= 0 && gid < nGlobal);
         
         double gidFrac = (double)gid/(double)nGlobal;
         int dest = gidFrac*nTasks;
         if (dest >= nTasks) dest = nTasks-1;
         
         int floatCnt1 = 0;
         for (unsigned jj=0; jj<nFloats; ++jj)
         {
            int find = fieldIndex1[jj];
            if (data1->dataType(find) == BucketOfBits::floatType)
            {
               double value;
               iRec.getValue(find, value);
               state1[nFloats*ii + floatCnt1].value_ = value;
               state1[nFloats*ii + floatCnt1].gid_ = gid;
               state1[nFloats*ii + floatCnt1].dest_ = dest;
               state1[nFloats*ii + floatCnt1].sortind_ = nFloats*gid + floatCnt1;
               floatCnt1++;
            }
         }
      }
      
      // directory 2
      for (unsigned ii=0; ii<nLocal2; ++ii)
      {
         BucketOfBits::Record iRec = data2->getRecord(ii);
         int gid;
         iRec.getValue(gidindex2, gid);
         assert(gid >= 0 && gid < nGlobal);
         double gidFrac = (double)gid/(double)nGlobal;
         int dest = gidFrac*nTasks;
         if (dest >= nTasks) dest = nTasks-1;
         
         int floatCnt2 = 0;
         for (unsigned jj=0; jj<nFloats; ++jj)
         {
            int find = fieldIndex2[jj];
            if (data2->dataType(find) == BucketOfBits::floatType)
            {
               double value;
               iRec.getValue(find, value);
               state2[nFloats*ii + floatCnt2].value_ = value;
               state2[nFloats*ii + floatCnt2].gid_ = gid;
               state2[nFloats*ii + floatCnt2].dest_ = dest;
               state2[nFloats*ii + floatCnt2].sortind_ = nFloats*gid + floatCnt2;
               floatCnt2++;
            }
         }
      }
      unsigned nSize1 = nLocal1*nFloats;
      unsigned nSize2 = nLocal2*nFloats;
      
      // distribute state values
      
      // directory 1
      sort(state1.begin(),state1.end(),cellState::destLessThan);
      vector<unsigned> dest1(nSize1);
      for (unsigned ii=0; ii<state1.size(); ++ii)
         dest1[ii] = state1[ii].dest_;
      
      int maxLocbuf1 = nSize1;
      int nMax1;
      MPI_Allreduce(&maxLocbuf1, &nMax1, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      
      state1.resize(10*nMax1);
      assignArray((unsigned char*)&(state1[0]), &nSize1, state1.capacity(),
                  sizeof(cellState), &(dest1[0]), 0, MPI_COMM_WORLD);
      
      assert(nSize1 <= (unsigned)nMax1);
      state1.resize(nSize1);
      
      // directory 2
      sort(state2.begin(),state2.end(),cellState::destLessThan);
      vector<unsigned> dest2(nSize2);
      for (unsigned ii=0; ii<state2.size(); ++ii)
         dest2[ii] = state2[ii].dest_;
      
      int maxLocbuf2 = nSize2;
      int nMax2;
      MPI_Allreduce(&maxLocbuf2, &nMax2, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      
      state2.resize(10*nMax2);
      assignArray((unsigned char*)&(state2[0]), &nSize2, state2.capacity(),
                  sizeof(cellState), &(dest2[0]), 0, MPI_COMM_WORLD);
      
      assert(nSize2 <= (unsigned)nMax2);
      state2.resize(nSize2);
      
      // sort data arrays so that local gids are in same order for both
      sort(state1.begin(),state1.end(),cellState::indLessThan);
      sort(state2.begin(),state2.end(),cellState::indLessThan);
      
      assert(nSize1 == nSize2);
      int nLoc = nSize1/nFloats;
      
      int floatCnt1 = 0;
      for (unsigned jj=0; jj<nFloats; ++jj)
      {
         double rmssqLoc = 0.0;
         int rmscnt = 0;
         double maxdiffLoc = 0.0;
         for (unsigned ii=0; ii<nLoc; ++ii)
         {
            Long64 gid1 = state1[nFloats*ii+jj].gid_;
            Long64 gid2 = state2[nFloats*ii+jj].gid_;
            assert(gid1 == gid2);
            
            double val1 = state1[nFloats*ii+jj].value_;
            double val2 = state2[nFloats*ii+jj].value_;
            rmssqLoc += (val1-val2)*(val1-val2);
            rmscnt++;
            if (abs(val1-val2) > maxdiffLoc) maxdiffLoc = abs(val1-val2);
         }
         double rmssq;
         MPI_Allreduce(&rmssqLoc, &rmssq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         double rms = sqrt(rmssq);
         double maxdiff;
         MPI_Allreduce(&maxdiffLoc, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

         //ewd DEBUG:  print out gid of all cells with Vm difference within 0.95 of max
         if (fieldUnion[jj] == "Vm")
         {
            cout.setf(ios::scientific,ios::floatfield);
            for (unsigned ii=0; ii<nLoc; ++ii)
            {
               Long64 gid1 = state1[nFloats*ii+jj].gid_;
               Long64 gid2 = state2[nFloats*ii+jj].gid_;
               assert(gid1 == gid2);
               
               double val1 = state1[nFloats*ii+jj].value_;
               double val2 = state2[nFloats*ii+jj].value_;
               double diff = abs(val1-val2);
               if (diff > 0.99*maxdiff && maxdiff > 1.E-9) 
                  cout << "  gid " << gid1 << ", Vm difference = " << setprecision(10) << diff << ", maximum difference = " << maxdiff << endl;               
            }
         }
         //ewd DEBUG

         
         // save output to file
         if (myRank == 0)
         {
            string fieldFile = "verif." + fieldUnion[jj] + ".dat";
            ofstream ofout;
            if (isnap == 0)
            {
               ofout.open(fieldFile.c_str(),ofstream::out);
               ofout << "#      time            rms avg            rms tot           max diff" << endl;
            }
            else
               ofout.open(fieldFile.c_str(),ofstream::app);
            double rmsavg = (nrecord1 > 0 ? rms/nrecord1 : 0.0);
            ofout.setf(ios::scientific,ios::floatfield);
            ofout << setprecision(10) << time1 << "    " << rmsavg << "    " << rms << "    " << maxdiff << endl;
            ofout.close();
         }
      }

      delete data1;
      delete data2;

   }
   MPI_Finalize();
   return 0;
}

namespace
{
    void readAnatomyBucket(Anatomy& anatomy, BucketOfBits* bucketP)
    {
       vector<AnatomyCell>& cells = anatomy.cellArray();
       
       unsigned nRecords =  bucketP->nRecords();
       unsigned nFields =   bucketP->nFields();
       unsigned gidIndex =  bucketP->getIndex("gid");
       assert( gidIndex != nFields);
       for (unsigned ii=0; ii<nRecords; ++ii)
       {
          AnatomyCell tmp;
          BucketOfBits::Record rr = bucketP->getRecord(ii);
          rr.getValue(gidIndex, tmp.gid_);
          tmp.cellType_ = 100; // hack since state files don't have cell type information
          cells.push_back(tmp);
       }
    }
    
    BucketOfBits* fillDataBucket(Simulate& sim, string dirName, double &time, int &nrecord)
    {
       string stateFile = dirName + "/state#";
       PFILE* file = Popen(stateFile.c_str(), "r", MPI_COMM_WORLD);
       OBJECT* hObj = file->headerObject;

       int nx, ny, nz;
       objectGet(hObj, "nx", nx, "0");
       objectGet(hObj, "ny", ny, "0");
       objectGet(hObj, "nz", nz, "0");
       assert(nx*ny*nz > 0);
       objectGet(hObj, "time", time, "0.0");
       objectGet(hObj, "nrecord", nrecord, "0");
       
       
       //string rTmp;
       //objectGet(hObj, "reactionMethod", rTmp, "reaction");
       //sim.reaction_ = reactionFactory(rTmp, sim.anatomy_, sim.reactionThreads_);

       sim.anatomy_.setGridSize(nx,ny,nz);

       // fill sim.anatomy_ with gids, arbitrary cell type
       BucketOfBits* bucketP = readPioFile(file);
       Pclose(file);
       readAnatomyBucket(sim.anatomy_, bucketP);
       delete bucketP;
       
       //cout << "First state directory, myRank = " << myRank << ", grid = " << nx << " x " << ny << " x " << nz << ", sim.nLocal = " << sim.anatomy_.nLocal() << endl;

       BucketOfBits* data = loadAndDistributeState(stateFile, sim.anatomy_);
       assert(data->nRecords() == sim.anatomy_.nLocal());

       return data;
    }

}
