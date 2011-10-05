#include "KoradiTest.hh"

#include <mpi.h>
#include <cmath>
#include "GridAssignmentObject.h"
#include "mpiUtils.h"

using namespace std;


class AnatomyCellDestSort
{
 public:
   bool operator()(const AnatomyCell& a, const AnatomyCell& b)
   {
      return a._dest < b._dest;
   }
};






KoradiTest::KoradiTest(AnatomyReader& anatomy,
		      int nCentersPerTask)
:_nCentersPerTask(nCentersPerTask),
 _indexToVector(anatomy._nx, anatomy._ny, anatomy._nz),
 _cells(anatomy._anatomy)
{
   MPI_Comm_size(MPI_COMM_WORLD, &_nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &_myRank);

   _localOffset = _myRank*_nCentersPerTask;
   
   pickInitialCenters();
   initialAssignment();
   computeRadii();
}

void KoradiTest::balanceStep()
{
//    updateBias();
//    findInteractingTasks();
//    computeDestinations();
//    assignCells();
//    moveCenters();
//    printStatistics();
}



void KoradiTest::pickInitialCenters()
{
   _centers.resize(_nTasks*_nCentersPerTask);
   for (int ii=0; ii<_nCentersPerTask; ++ii)
   {
      int jj = drand48() * _cells.size();
      _centers[ii+_localOffset] = _indexToVector(_cells[jj]._gid);
   }

   allGatherCenters();
}

void KoradiTest::initialAssignment()
{
   GRID_ASSIGNMENT_OBJECT* gao = gao_init(_centers.size(),
					  (const void*) &(_centers[0]),
					  sizeof(THREE_VECTOR));
   for (unsigned ii=0; ii<_cells.size(); ++ii)
   {
      THREE_VECTOR r = _indexToVector(_cells[ii]._gid);
      _cells[ii]._dest = gao_nearestCenter(gao, r);
   }
   gao_destroy(gao);

   AnatomyCellDestSort sortByDest;
   sort(_cells.begin(), _cells.end(), sortByDest);
   vector<unsigned> dest(_cells.size());
   for (unsigned ii=0; ii<_cells.size(); ++ii)
      dest[ii] = _cells[ii]._dest/_nCentersPerTask;

   unsigned nLocal = _cells.size();
   unsigned capacity = 2*_cells.size();
   _cells.reserve(capacity);
   assignArray((unsigned char*) &(_cells[0]),
	       &nLocal,
	       capacity,
	       sizeof(AnatomyCell),
	       &(dest[0]),
	       0,
	       MPI_COMM_WORLD);
}



void KoradiTest::computeRadii()
{
   for (unsigned ii=0; ii<_radii.size(); ++ii)
      _radii[ii] = 0;
   
   for (int ii=0; ii<_cells.size(); ++ii)
   {
      THREE_VECTOR v = _indexToVector(_cells[ii]._gid);
      int cellOwner = _cells[ii]._dest;
      double r2 = DOT(v, _centers[cellOwner]);
      _radii[cellOwner] = max(_radii[cellOwner], r2);
   }

   for (unsigned ii=0; ii<_radii.size(); ++ii)
      _radii[ii] = sqrt(_radii[ii]);

   allReduceRadii();
}



void KoradiTest::allGatherCenters()
{
   void* sendBuf = (void*)&(_centers[_localOffset]);
   void* recvBuf = (void*)&(_centers[0]);
   int nSend = _nCentersPerTask*sizeof(THREE_VECTOR);
   int nRecv = _nCentersPerTask*sizeof(THREE_VECTOR);
   
   MPI_Allgather(sendBuf, nSend, MPI_CHAR,
		 recvBuf, nRecv, MPI_CHAR, MPI_COMM_WORLD);
}

void KoradiTest::allReduceRadii()
{
   void* sendBuf = (void*)&(_radii[_localOffset]);
   void* recvBuf = (void*)&(_radii[0]);
   int nSend = _nCentersPerTask;
   
   MPI_Allreduce(sendBuf, recvBuf, nSend, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}
