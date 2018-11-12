#include "Koradi.hh"

#include <mpi.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <sstream>
#include "ioUtils.h"
#include "GridAssignmentObject.h"
#include "mpiUtils.h"
#include "mpiTpl.hh"
#include "Drand48Object.hh"
#include "Simulate.hh"
#include "writeCells.hh"
using namespace std;


// ToDo
// 1.  The (test) results are not reproducible across different numbers
//     of tasks (when overprovisioned).  This is due to the way the
//     initial centers are chosen.
// 3.  Last vestige of THREE_VECTOR is in call to gao_nearestCenter.
//     Figure out how to get rid of it.
// 4.  Sort out whether it is better to balance to the local average or
//     the global average load.
// 5.  Need to implement an early bailout clause for the voronoi
//     balance.
// 6.  Need to handle the case where the tolerance is set to zero.
// 7.  There are some efficiency details to work out.  Make sure things
//     aren't calculated multiple times.
// 8.  Decide if there is any reason to support some of the feature that
//     I have played with during development such as the ability to
//     limit the size of the bias change, or allowing periodic
//     attempts to recondition the system by going to a pure voronoi
//     balancing.
// 9.  Think about whether there should be a minium radius that a domain
//     can have.  Any cell within that radius belongs unconditionally to
//     its nearest domain.



Koradi::Koradi(Anatomy& anatomy, const KoradiParms& parms)
:verbose_(parms.verbose),
 nCentersPerTask_(parms.nCentersPerTask),
 maxVoronoiSteps_(parms.maxVoronoiSteps),
 maxSteps_       (parms.maxSteps),
 outputRate_     (parms.outputRate),
 tolerance_      (parms.tolerance),
 nbrDeltaR_      (parms.nbrDeltaR),
 alphaStep_      (parms.alphaStep),
 indexToVector_  (anatomy.nx(), anatomy.ny(), anatomy.nz()),
 indexTo3Vector_ (anatomy.nx(), anatomy.ny(), anatomy.nz()),
 cells_          (anatomy.cellArray())
{
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks_);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);

   localOffset_ = myRank_*nCentersPerTask_;

   centers_.resize(nTasks_*nCentersPerTask_);
   radii_.resize(nTasks_*nCentersPerTask_);
   alpha_.resize(nTasks_*nCentersPerTask_, 1.0);
   load_.resize(nTasks_*nCentersPerTask_, 0.0);
   nbrDomains_.resize(nCentersPerTask_);

   if (centers_.size() == 1 ) maxVoronoiSteps_ = 1;

   distributeCellsEvenly();
   pickInitialCenters();
   voronoiBalance();

   for (unsigned ii=0; ii<maxSteps_; ++ii)
   {
      balanceStep();
      if (ii>0  && outputRate_>0 && ii%outputRate_ == 0)
      {
         stringstream name;
         name << "balance."<< setfill('0') <<setw(8) << ii;
         string fullname = name.str();
         if (myRank_ == 0)
            DirTestCreate(fullname.c_str());
         fullname += "/domains";
         writeCells(cells_, anatomy.nx(), anatomy.ny(), anatomy.nz(),
                    fullname.c_str());
      }
      
      if (verbose_)
         printStatistics();

      double maxLoad = *max_element(load_.begin(), load_.end());
      double imbalance = abs(maxLoad - targetLoad_)/targetLoad_;
      if (imbalance < tolerance_)
         break;
   }
   /*ewd DEBUG:  comment this out to save on I/O
   stringstream name;
   name << "balance.final";
   string fullname = name.str();
   if (myRank_ == 0)
   DirTestCreate(fullname.c_str());
   fullname += "/domains";
   writeCells(cells_, anatomy.nx(), anatomy.ny(), anatomy.nz(),
              fullname.c_str());
   */
}

void Koradi::balanceStep()
{
   findNbrDomains();
   biasAlpha();
   assignCells();
   moveCenters();
   computeRadii();
}

void Koradi::voronoiBalance()
{
   for (unsigned ii=0; ii<maxVoronoiSteps_; ++ii)
   {
      assignCells();
      moveCenters();
      computeRadii();
      if (verbose_) {
         computeLoad(load_);
         printStatistics();
      }

//       if (stats.imbalance < tolerance_)
//       break;
   }
}




void Koradi::distributeCellsEvenly()
{
   Long64 nLocal = cells_.size();
   Long64 nGlobal = 0;
   MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   if (myRank_ == 0) cout << "nGlobal = "<<nGlobal<<" nAve = " <<nGlobal/nTasks_/nCentersPerTask_<<endl;
   
   size_t nWant = nGlobal / nTasks_;
   unsigned nExtra = nGlobal % nTasks_;
   if (myRank_ < nExtra)
      ++nWant;

   cells_.resize(max(nWant, cells_.size()));
   distributeArray((unsigned char*)&(cells_[0]),
                   nLocal,
                   nWant,
                   sizeof(AnatomyCell),
                   MPI_COMM_WORLD);
   cells_.resize(nWant);
}


void Koradi::pickInitialCenters()
{
   assert(cells_.size() >= nCentersPerTask_);
   
   vector<int> indexArray(cells_.size());
   for (unsigned ii=0; ii< indexArray.size(); ++ii)
      indexArray[ii] = ii;
   Drand48Object rand(myRank_);
   random_shuffle(indexArray.begin(), indexArray.end(), rand);
   
   for (int ii=0; ii<nCentersPerTask_; ++ii)
   {
      Long64 gid = cells_[indexArray[ii]].gid_;
      centers_[ii+localOffset_] = indexToVector_(gid);
   }

   allGather(centers_, nCentersPerTask_, MPI_COMM_WORLD);
   if (myRank_ ==0)
   {
      cout <<"Finished initial centers" <<endl;
   }
   
}

void Koradi::assignCells()
{
   calculateCellDestinations();
   exchangeCells();
}


void Koradi::calculateCellDestinations()
{
   // Without nbr domain info, bootstrap using a grid approach.
   if (nbrDomains_[0].size() == 0)
   {
      GRID_ASSIGNMENT_OBJECT* gao = gao_init(centers_.size(),
                                             (const void*) &(centers_[0]),
                                             sizeof(Vector));
   
      for (unsigned ii=0; ii<cells_.size(); ++ii)
      {
         THREE_VECTOR r = indexTo3Vector_(cells_[ii].gid_);
         cells_[ii].dest_ = gao_nearestCenter(gao, r);
      }
      gao_destroy(gao);
   }
   else // we know which domains might be close
   {
      #pragma omp parallel for
      for (int ii=0; ii<cells_.size(); ++ii)
      {
         int cellOwner = cells_[ii].dest_;
         Vector rCell = indexToVector_(cells_[ii].gid_);
         assert(cellOwner >= localOffset_);
         assert(cellOwner < localOffset_+nCentersPerTask_);
         const vector<int>& overLapList = nbrDomains_[cellOwner-localOffset_];
         double r2Min = diffSq(rCell, centers_[cellOwner])/alpha_[cellOwner];
         cells_[ii].dest_ = cellOwner;
         for (unsigned jj=0; jj<overLapList.size(); ++jj)
         {
            double r2 = diffSq(rCell, centers_[overLapList[jj]])/alpha_[overLapList[jj]];
            if (r2 < r2Min)
            {
               r2Min = r2;
               cells_[ii].dest_ = overLapList[jj];
            }
         }
         
      }
   }
}

void Koradi::exchangeCells()
{
   AnatomyCellDestSort sortByDest;
   sort(cells_.begin(), cells_.end(), sortByDest);
   vector<unsigned> dest(cells_.size());
   for (unsigned ii=0; ii<cells_.size(); ++ii)
      dest[ii] = cells_[ii].dest_/nCentersPerTask_;
   
   unsigned nLocal = cells_.size();
   unsigned capacity = max(vector<AnatomyCell>::size_type(10000), 4*cells_.size());
   cells_.resize(capacity);
   assignArray((unsigned char*) &(cells_[0]),
               &nLocal,
               capacity,
               sizeof(AnatomyCell),
               &(dest[0]),
               0,
               MPI_COMM_WORLD);
   cells_.resize(nLocal);
   sort(cells_.begin(), cells_.end(), sortByDest);
}


void Koradi::moveCenters()
{
   // we include the current center location in the center of mass
   // calculation.  This means a domain with no cells won't suddenly
   // warp to the origin.
   vector<int> nCells(centers_.size(), 1);
   for (unsigned ii=0; ii<cells_.size(); ++ii)
   {
      int dest = cells_[ii].dest_;
      Vector r = indexToVector_(cells_[ii].gid_);
      ++nCells[dest];
      centers_[dest] += r;
   }

   for (unsigned ii=0; ii<nCentersPerTask_; ++ii)
      centers_[ii+localOffset_] /= double(nCells[ii+localOffset_]);
   allGather(centers_, nCentersPerTask_, MPI_COMM_WORLD);
}

// We impose minimum radius to ensure that if a domain happens to
// include no cells (or perhaps one cell) it still has a non-zero
// volume.
void Koradi::computeRadii()
{
   radii_.assign(radii_.size(), 0.);
   
   #pragma omp parallel for
   for (int ii=0; ii<cells_.size(); ++ii)
   {
      Vector v = indexToVector_(cells_[ii].gid_);
      int cellOwner = cells_[ii].dest_;
      assert(cellOwner/nCentersPerTask_ == myRank_);
      double r2 = diffSq(v,  centers_[cellOwner]);
      radii_[cellOwner] = max(radii_[cellOwner], r2);
   }

   for (unsigned ii=0; ii<radii_.size(); ++ii)
      radii_[ii] = sqrt(radii_[ii]);

   allGather(radii_, nCentersPerTask_, MPI_COMM_WORLD);
}

void Koradi::printStatistics()
{
   computeRadii();

   double maxLoad = *max_element(load_.begin(), load_.end());
   double minLoad = *min_element(load_.begin(), load_.end());
   double maxRadius = *max_element(radii_.begin(), radii_.end());
   double minRadius = *min_element(radii_.begin(), radii_.end());

   if (myRank_ == 0)
   {
      cout << "min/max load, radius = "
           << minLoad << " " << maxLoad << " "
           << minRadius << " " << maxRadius << endl;
   }
}

void Koradi::biasAlpha()
{
   computeLoad(load_);

   double globalAveLoad=0;
   for (unsigned ii=0; ii<load_.size(); ++ii)
      globalAveLoad += load_[ii];
   globalAveLoad /= load_.size();
   
   targetLoad_ = globalAveLoad;

   for (unsigned ii=0; ii<nCentersPerTask_; ++ii)
   {
      //Bigger Alpha = bigger radius
      if (load_[localOffset_+ii]/globalAveLoad < (1-tolerance_))
      {
         alpha_[localOffset_+ii] *= 1+alphaStep_;
      }
      else if (load_[localOffset_+ii]/globalAveLoad > (1+tolerance_))
      {
         alpha_[localOffset_+ii] /= 1+alphaStep_;
      }
   }

   allGather(alpha_, nCentersPerTask_, MPI_COMM_WORLD);
}

void Koradi::findNbrDomains()
{
   for (unsigned ii=0; ii<nCentersPerTask_; ++ii)
      nbrDomains_[ii].clear();

   for (unsigned iCenter=0; iCenter<nCentersPerTask_; ++iCenter)
   {
      unsigned ii = iCenter+localOffset_;
      const Vector& ci = centers_[ii];
      double rii = radii_[ii];
      for (unsigned jj=0; jj<centers_.size(); ++jj)
      {
         if (jj == ii)
            continue;
         Vector cij = ci - centers_[jj];
         double r2 = dot(cij, cij);
         double rjj = radii_[jj];
         if (r2 < (rii+rjj+nbrDeltaR_)*(rii+rjj+nbrDeltaR_))
             nbrDomains_[iCenter].push_back(jj);
      }
   }
}



void Koradi::bruteForceDistanceCheck()
{
   for (unsigned ii=0; ii<cells_.size(); ++ii)
   {
      double r2Min = 1e30;
      int dest = -1;
      Vector r = indexToVector_(cells_[ii].gid_);
      for (unsigned jj=0; jj<centers_.size(); ++jj)
      {
         Vector rij = r - centers_[jj];
         double r2 = dot(rij, rij);
         if (r2 < r2Min)
         {
            r2Min = r2;
            dest = jj;
         }
      }
      if (dest != cells_[ii].dest_)
         cout << "Fail "<<ii<<" fast " <<cells_[ii].dest_<<" brute "<<dest<<endl;
      
         //assert(dest == cells_[ii].dest_);
      
   }
}

void Koradi::computeLoad(vector<double>& load)
{
   load.assign(nTasks_*nCentersPerTask_, 0.);

   for (unsigned ii=0; ii<cells_.size(); ++ii)
   {
      load[cells_[ii].dest_] += 1;
   }
   allGather(load, nCentersPerTask_, MPI_COMM_WORLD);
}
