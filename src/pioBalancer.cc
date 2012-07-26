#include "pioBalancer.hh"

#include <vector>
#include <cassert>
#include <algorithm>

#include "Simulate.hh"
#include "BucketOfBits.hh"
#include "pio.h"
#include "stateLoader.hh"
#include "mpiUtils.h"
#include "readPioFile.hh"


using namespace std;

namespace
{
   void loadAndAssignDomains(
      const string& filename, Simulate& sim, MPI_Comm comm)
   {
      int nTasks;
      MPI_Comm_size(comm, &nTasks);
      
      BucketOfBits* data = loadAndDistributeState(filename, sim.anatomy_);
      unsigned gidIndex = data->getIndex("gid");
      unsigned domainIndex = data->getIndex("domain");
      assert(gidIndex != data->nFields());
      assert(domainIndex != data->nFields());
      
      vector<AnatomyCell>& cells = sim.anatomy_.cellArray();
      for (unsigned ii=0; ii<sim.anatomy_.size(); ++ii)
      {
         BucketOfBits::Record rr = data->getRecord(ii);
         Long64 gidTmp;
         rr.getValue(gidIndex, gidTmp);
         assert(gidTmp == sim.anatomy_.gid(ii));
         rr.getValue(domainIndex, cells[ii].dest_);
         assert(cells[ii].dest_ < nTasks); // need better error handling
      }
      sort(cells.begin(),cells.end(),AnatomyCell::destLessThan);
      unsigned nLocal = cells.size();
      vector<unsigned> dest(nLocal);
      vector<unsigned> nRecv(nTasks, 0);
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         dest[ii] = cells[ii].dest_;
         ++nRecv[dest[ii]];
      }
      vector<int> recvCnts(nTasks, 1);
      int nFinal;
      MPI_Reduce_scatter(&nRecv[0], &nFinal, &recvCnts[0], MPI_UNSIGNED, MPI_SUM, comm);
      
      unsigned capacity = max(vector<AnatomyCell>::size_type(nFinal), cells.size());
      cells.resize(capacity);
      assignArray((unsigned char*)&(cells[0]), &nLocal, cells.capacity(),
                  sizeof(AnatomyCell), &(dest[0]), 0, comm);
      cells.resize(nLocal);
      assert(nFinal == nLocal);
      delete data;
   }
}

class DomainData
{
 public:
   int rank;
   int nD;
   bool operator<(const DomainData& b) const {return this->rank < b.rank;}
};


namespace
{
   int loadAndDistributeDiffusionCores(const string& filename, MPI_Comm comm)
   {
      int nTasks;
      int myRank;
      MPI_Comm_size(comm, &nTasks);
      MPI_Comm_rank(comm, &myRank);
      
      PFILE* pfile = Popen(filename.c_str(), "r", comm);
      BucketOfBits* data = readPioFile(pfile);
      Pclose(pfile);

      // Check that we have the same number of records as tasks.  We
      // should have a better error handling strategy for a mismatch.
      unsigned nGlobal;
      unsigned nLocal = data->nRecords();
      MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_UNSIGNED, MPI_SUM, comm);
      assert(nGlobal == nTasks);
      
      vector<DomainData> records(nLocal);
      unsigned rankIndex = data->getIndex("domain");
      unsigned nDIndex = data->getIndex("nD");
      assert(rankIndex != data->nFields());
      assert(nDIndex != data->nFields());
      
      for (unsigned ii=0; ii<nLocal; ++ii)
      {
         BucketOfBits::Record rr = data->getRecord(ii);
         rr.getValue(rankIndex, records[ii].rank);
         rr.getValue(nDIndex,   records[ii].nD);
      }
      delete data;

      sort(records.begin(), records.end());
      vector<unsigned> dest(records.size());
      for (unsigned ii=0; ii<records.size(); ++ii)
         dest[ii] = records[ii].rank;
      
      unsigned capacity = max(records.size(), vector<int>::size_type(1));
      records.resize(capacity);
      assignArray((unsigned char*)&records[0],
                  &nLocal,
                  capacity,
                  sizeof(DomainData),
                  &dest[0],
                  0,
                  comm);
      assert(nLocal == 1);
      records.resize(nLocal);

      assert(records[0].rank == myRank);
      return records[0].nD;
   }
}

int pioBalancer(const string& domainFile,
                const string& pxyzFile,
                Simulate& sim,
                MPI_Comm comm)
{
   loadAndAssignDomains(domainFile, sim, comm);
   int nD = loadAndDistributeDiffusionCores(pxyzFile, comm);
   return nD;
}
