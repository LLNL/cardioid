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

template<typename T>
inline void
allGather(vector<T>& data, size_t nDataPerTask, MPI_Comm comm)
{
   int myRank;
   MPI_Comm_rank(comm, &myRank);

   size_t localStart = myRank*nDataPerTask;
   typename vector<T>::const_iterator first = data.begin() + localStart;
   typename vector<T>::const_iterator last = first + nDataPerTask;
//   assert( data[localStart] == *first);

   vector<T> tmp(first, last);

   void* sendBuf = (void*) &tmp[0];
   void* recvBuf = (void*) &data[0];
   int nSend = nDataPerTask * sizeof(T);
   int nRecv = nDataPerTask * sizeof(T);
   
   MPI_Allgather(sendBuf, nSend, MPI_CHAR, recvBuf, nRecv, MPI_CHAR, comm);
}





KoradiTest::KoradiTest(AnatomyReader& anatomy,
		      int nCentersPerTask)
:_nCentersPerTask(nCentersPerTask),
 _indexToVector(anatomy._nx, anatomy._ny, anatomy._nz),
 _cells(anatomy._anatomy)
{
   MPI_Comm_size(MPI_COMM_WORLD, &_nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &_myRank);

   _alpha = 0.9;
   
   _localOffset = _myRank*_nCentersPerTask;

   _centers.resize(_nTasks*_nCentersPerTask);
   _radii.resize(_nTasks*_nCentersPerTask);
   _bias.resize(_nTasks*_nCentersPerTask, 0.0);
   _nbrDomains.resize(_nCentersPerTask);
   
   distributeCellsEvenly();
   pickInitialCenters();
   for (unsigned ii=0; ii<100; ++ii)
   {
      assignCells();
      moveCenters();
      _bias.assign(_bias.size(), 0);
   printStatistics();
   }
}

void KoradiTest::balanceStep()
{
   findNbrDomains();
   updateBias();
   assignCells();
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
   _cells.resize(nWant);
}


void KoradiTest::pickInitialCenters()
{
   assert(_cells.size() >= _nCentersPerTask);
   
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

   allGather(_centers, _nCentersPerTask, MPI_COMM_WORLD);
   if (_myRank ==0)
   {
      cout <<"Finished initial centers" <<endl;
   }
   
}

void KoradiTest::assignCells()
{
   calculateCellDestinations();
   exchangeCells();
}


void KoradiTest::calculateCellDestinations()
{
   // Without nbr domain info, bootstrap using a grid approach.
   if (_nbrDomains[0].size() == 0)
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
   }
   else // we know which domains might be close
   {
      for (unsigned ii=0; ii<_cells.size(); ++ii)
      {
	 int cellOwner = _cells[ii]._dest;
	 THREE_VECTOR rCell = _indexToVector(_cells[ii]._gid);
	 const vector<int> overLapList = _nbrDomains[cellOwner-_localOffset];
	 double r2Min = DIFFSQ(rCell, _centers[cellOwner]) - _bias[cellOwner];
	 _cells[ii]._dest = cellOwner;
	 for (unsigned jj=0; jj<overLapList.size(); ++jj)
	 {
	    double r2 = DIFFSQ(rCell, _centers[overLapList[jj]]) - _bias[overLapList[jj]];
	    if (r2 < r2Min)
	    {
	       r2Min = r2;
	       _cells[ii]._dest = overLapList[jj];
	    }
	 }
	 
      }
   }
}

void KoradiTest::exchangeCells()
{
   AnatomyCellDestSort sortByDest;
   sort(_cells.begin(), _cells.end(), sortByDest);
   vector<unsigned> dest(_cells.size());
   for (unsigned ii=0; ii<_cells.size(); ++ii)
      dest[ii] = _cells[ii]._dest/_nCentersPerTask;
   
   unsigned nLocal = _cells.size();
   unsigned capacity = max(10000ul, 3*_cells.size());
   _cells.resize(capacity);
   assignArray((unsigned char*) &(_cells[0]),
	       &nLocal,
	       capacity,
	       sizeof(AnatomyCell),
	       &(dest[0]),
	       0,
	       MPI_COMM_WORLD);
   _cells.resize(nLocal);
   sort(_cells.begin(), _cells.end(), sortByDest);
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
   allGather(_centers, _nCentersPerTask, MPI_COMM_WORLD);

   vector<double> oldRadii = _radii;
   computeRadii();
   for (unsigned ii=0; ii<_centers.size(); ++ii)
   {
      double rNew = _radii[ii];
      double rOld = oldRadii[ii];
      double update = 0.25*(rNew-rOld)*(rNew-rOld);
      _bias[ii] -= update;
//      cout <<"update " << ii<<": " << update<<endl;
   }
   


}

// We impose minimum radius of 1 units to ensure that if a domain
// happens to include no cells (or perhaps one cell) it still has a
// non-zero volume.
void KoradiTest::computeRadii()
{
   _radii.assign(_radii.size(), 0.0);
   
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

   allGather(_radii, _nCentersPerTask, MPI_COMM_WORLD);
}

void KoradiTest::printStatistics()
{
   computeRadii();

   vector<unsigned> nCells(_centers.size(), 0);
   for (unsigned ii=0; ii<_cells.size(); ++ii)
   {
      unsigned dest = _cells[ii]._dest;
      ++nCells[dest];
      assert( dest >= _localOffset && dest < _localOffset+_nCentersPerTask);
   }
   
   allGather(nCells, _nCentersPerTask, MPI_COMM_WORLD);

   int maxCells = *max_element(nCells.begin(), nCells.end());
   int minCells = *min_element(nCells.begin(), nCells.end());
   double maxRadius = *max_element(_radii.begin(), _radii.end());
   double minRadius = *min_element(_radii.begin(), _radii.end());
   double maxBias = *max_element(_bias.begin(), _bias.end());
   double minBias = *min_element(_bias.begin(), _bias.end());

   if (_myRank == 0)
   {
      cout << "min/max nCells, radius, bias  = "
	   << minCells << " " << maxCells << " "
	   << minRadius << " " << maxRadius << " "
	   << minBias << " " << maxBias << endl;
//       for (unsigned ii=0; ii<_centers.size(); ++ii)
//       {
// 	 cout << ii <<": " <<_bias[ii]<<" " << nCells[ii] <<" "<<_radii[ii]
// 	      <<" " << _nbrDomains[ii].size()<< " ";
	 
// 	 for (unsigned jj=0; jj<_nbrDomains[ii].size(); ++jj)
// 	    cout << _nbrDomains[ii][jj] <<",";
// 	 cout <<endl;
//       }

   }
}

void KoradiTest::updateBias()
{
   vector<double> load;
   computeLoad(load);
   allGather(load, _nCentersPerTask, MPI_COMM_WORLD);

   for (unsigned ii=0; ii<_nCentersPerTask; ++ii)
   {
      double localAverageLoad = 0;
      const vector<int>& overLapList = _nbrDomains[ii];
      for (unsigned jj=0; jj<overLapList.size(); ++jj)
	 localAverageLoad += load[overLapList[jj]];

      assert(overLapList.size() > 0);
      localAverageLoad /= (1.0*overLapList.size());
      double tmp = localAverageLoad/load[_localOffset+ii];
      tmp = pow(tmp, (2.0/3.0));
      double r = _radii[ii+_localOffset];
      _bias[ii+_localOffset] += _alpha*r*r*(tmp-1.0);
   }

   allGather(_bias, _nCentersPerTask, MPI_COMM_WORLD);

//    if (_myRank == 0)
//    {
//    }
   

}

void KoradiTest::findNbrDomains()
{
//    for (unsigned ii=0; ii<_nCentersPerTask; ++ii)
//       _nbrDomains[ii].clear();
//    for (unsigned ii=0; ii<_nCentersPerTask; ++ii)
//       for (unsigned jj=0; jj<_nCentersPerTask*_nTasks; ++jj)
//       {
// 	 if (ii==jj) continue;
// 	 _nbrDomains[ii].push_back(jj);
//       }
//    return;
   







   double deltaR = 2;
   for (unsigned ii=0; ii<_nCentersPerTask; ++ii)
      _nbrDomains[ii].clear();

   for (unsigned iCenter=0; iCenter<_nCentersPerTask; ++iCenter)
   {
      unsigned ii = iCenter+_localOffset;
      THREE_VECTOR& rC = _centers[ii];
      double rii = _radii[ii];
      for (unsigned jj=0; jj<_centers.size(); ++jj)
      {
	 if (jj == ii)
	    continue;
	 double r2 = DIFFSQ(rC, _centers[jj]);
	 double rjj = _radii[jj];
	 if (r2 < (rii+rjj+deltaR)*(rii+rjj+deltaR))
	     _nbrDomains[iCenter].push_back(jj);
      }
   }
}



// void KoradiTest::allGatherCenters()
// {
//    vector<THREE_VECTOR> tmp(_centers.size());
//    void* sendBuf = (void*)&(_centers[_localOffset]);
//    void* recvBuf = (void*)&(tmp[0]);
//    int nSend = _nCentersPerTask*sizeof(THREE_VECTOR);
//    int nRecv = _nCentersPerTask*sizeof(THREE_VECTOR);
   
//    MPI_Allgather(sendBuf, nSend, MPI_CHAR,
// 		 recvBuf, nRecv, MPI_CHAR, MPI_COMM_WORLD);
//    _centers = tmp;
// }

// void KoradiTest::allReduceRadii()
// {
//    vector<double> tmp(_radii.size());
//    void* sendBuf = (void*)&(_radii[0]);
//    void* recvBuf = (void*)&(tmp[0]);
//    int nSend = _radii.size();
   
//    MPI_Allreduce(sendBuf, recvBuf, nSend, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//    _radii = tmp;
// }

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

void KoradiTest::computeLoad(vector<double>& load)
{
   load.resize(_nTasks*_nCentersPerTask, 0);

   for (unsigned ii=0; ii<_cells.size(); ++ii)
   {
      load[_cells[ii]._dest] += 1;
   }
   allGather(load, _nCentersPerTask, MPI_COMM_WORLD);
}
