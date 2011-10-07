#include "KoradiTest.hh"

#include <mpi.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "GridAssignmentObject.h"
#include "mpiUtils.h"
#include "Drand48Object.hh"
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

   distributeCellsEvenly();
   pickInitialCenters();
   initialAssignment();
   moveCenters();
   computeRadii();
   printStatistics();
}

void KoradiTest::balanceStep()
{
//    updateBias();
//    findInteractingTasks();
//    computeDestinations();
//    assignCells();
   initialAssignment();
   moveCenters();
   printStatistics();
}

void KoradiTest::distributeCellsEvenly()
{
   Long64 nLocal = _cells.size();
   Long64 nGlobal = 0;
   MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   size_t nWant = nGlobal / _nTasks;
   unsigned nExtra = nGlobal % _nTasks;
   if (_myRank < nExtra)
      ++nWant;

   _cells.resize(max(nWant, _cells.size()));
   distributeArray((unsigned char*)&(_cells[0]),
		   nLocal,
		   nWant,
		   sizeof(AnatomyCell),
		   MPI_COMM_WORLD);
}


void KoradiTest::pickInitialCenters()
{
   assert(_cells.size() >= _nCentersPerTask);
   
   _centers.resize(_nTasks*_nCentersPerTask);
   
   vector<int> indexArray(_cells.size());
   for (unsigned ii=0; ii< indexArray.size(); ++ii)
      indexArray[ii] = ii;
   Drand48Object rand(_myRank);
   random_shuffle(indexArray.begin(), indexArray.end(), rand);
   
   for (int ii=0; ii<_nCentersPerTask; ++ii)
   {
      Long64 gid = _cells[indexArray[ii]]._gid;
      _centers[ii+_localOffset] = _indexToVector(gid);
      _centers[ii+_localOffset].x += .1;
      _centers[ii+_localOffset].y += .2;
      _centers[ii+_localOffset].z += .3;
   }

   allGatherCenters();
   if (_myRank ==0)
   {
//       for (unsigned ii=0; ii<_centers.size(); ++ii)
// 	 cout << "centers " << ii <<": "
// 	      << _centers[ii].x << " "
// 	      << _centers[ii].y << " "
// 	      << _centers[ii].z << " " << endl;
      cout <<"Finished initial centers" <<endl;
   }
   
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

   bruteForceDistanceCheck();

   AnatomyCellDestSort sortByDest;
   sort(_cells.begin(), _cells.end(), sortByDest);
   vector<unsigned> dest(_cells.size());
   for (unsigned ii=0; ii<_cells.size(); ++ii)
      dest[ii] = _cells[ii]._dest/_nCentersPerTask;


   bruteForceDistanceCheck();




   
   unsigned nLocal = _cells.size();
   unsigned capacity = 3*_cells.size();
   _cells.resize(capacity);
   if (_myRank == 0) cout << "Entering assignArray" << endl;
   assignArray((unsigned char*) &(_cells[0]),
	       &nLocal,
	       capacity,
	       sizeof(AnatomyCell),
	       &(dest[0]),
	       0,
	       MPI_COMM_WORLD);
   _cells.resize(nLocal);
   sort(_cells.begin(), _cells.end(), sortByDest);
   if (_myRank ==0) cout <<"Finished initial assignment" <<endl;
}

void KoradiTest::moveCenters()
{
   _centers.assign(_centers.size(), vzero);

   vector<int> nCells(_centers.size(), 0);
   for (unsigned ii=0; ii<_cells.size(); ++ii)
   {
      int dest = _cells[ii]._dest;
      THREE_VECTOR r = _indexToVector(_cells[ii]._gid);
      ++nCells[dest];
      _centers[dest].x += r.x;
      _centers[dest].y += r.y;
      _centers[dest].z += r.z;
   }

   for (unsigned ii=0; ii<_nCentersPerTask; ++ii)
      VSCALE(_centers[ii+_localOffset], 1/double(nCells[ii+_localOffset]));
   allGatherCenters();
}


void KoradiTest::computeRadii()
{
   _radii.assign(_centers.size(), 0);
   
   for (int ii=0; ii<_cells.size(); ++ii)
   {
      THREE_VECTOR v = _indexToVector(_cells[ii]._gid);
      int cellOwner = _cells[ii]._dest;
      assert(cellOwner/_nCentersPerTask == _myRank);
      double r2 = DIFFSQ(v, _centers[cellOwner]);
      _radii[cellOwner] = max(_radii[cellOwner], r2);
   }

   for (unsigned ii=0; ii<_radii.size(); ++ii)
      _radii[ii] = sqrt(_radii[ii]);

   allReduceRadii();
}

void KoradiTest::printStatistics()
{
   computeRadii();

   vector<unsigned> nCellsLocal(_centers.size(), 0);
   for (unsigned ii=0; ii<_cells.size(); ++ii)
   {
      int index = _cells[ii]._dest - _localOffset;
      assert(index >= 0 && index < _centers.size());
      ++nCellsLocal[index];
   }
   
   vector<int> nCells(_nTasks*_nCentersPerTask);
   void* sendBuf = (void*)&(nCellsLocal[0]);
   void* recvBuf = (void*)&(nCells[0]);
   int nSend = _nCentersPerTask;
   int nRecv = _nCentersPerTask;
   
   MPI_Allgather(sendBuf, nSend, MPI_INT,
		 recvBuf, nRecv, MPI_INT, MPI_COMM_WORLD);

   int maxCells = *max_element(nCells.begin(), nCells.end());
   int minCells = *min_element(nCells.begin(), nCells.end());
   double maxRadius = *max_element(_radii.begin(), _radii.end());
   double minRadius = *min_element(_radii.begin(), _radii.end());

   if (_myRank == 0)
      cout << "min/max nCells = " << minCells << " " << maxCells << "\n"
	   << "min/max radius = " << minRadius << " " << maxRadius << endl;

}




void KoradiTest::allGatherCenters()
{
   vector<THREE_VECTOR> tmp(_centers.size());
   void* sendBuf = (void*)&(_centers[_localOffset]);
   void* recvBuf = (void*)&(tmp[0]);
   int nSend = _nCentersPerTask*sizeof(THREE_VECTOR);
   int nRecv = _nCentersPerTask*sizeof(THREE_VECTOR);
   
   MPI_Allgather(sendBuf, nSend, MPI_CHAR,
		 recvBuf, nRecv, MPI_CHAR, MPI_COMM_WORLD);
   _centers = tmp;
}

void KoradiTest::allReduceRadii()
{
   vector<double> tmp(_radii.size());
   void* sendBuf = (void*)&(_radii[0]);
   void* recvBuf = (void*)&(tmp[0]);
   int nSend = _radii.size();
   
   MPI_Allreduce(sendBuf, recvBuf, nSend, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   _radii = tmp;
}

void KoradiTest::bruteForceDistanceCheck()
{
   for (unsigned ii=0; ii<_cells.size(); ++ii)
   {
      double r2Min = 1e30;
      int dest = -1;
      THREE_VECTOR r = _indexToVector(_cells[ii]._gid);
      for (unsigned jj=0; jj<_centers.size(); ++jj)
      {
	 double r2 = DIFFSQ(r, _centers[jj]);
	 if (r2 < r2Min)
	 {
	    r2Min = r2;
	    dest = jj;
	 }
      }
      if (dest != _cells[ii]._dest)
	 cout << "Fail "<<ii<<" fast " <<_cells[ii]._dest<<" brute "<<dest<<endl;
      
	 //assert(dest == _cells[ii]._dest);
      
   }
}

