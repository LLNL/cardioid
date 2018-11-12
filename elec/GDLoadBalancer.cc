#include "GDLoadBalancer.hh"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <mpi.h>
#include "mpiUtils.h"
#include "GridPoint.hh"
#include "AnatomyCell.hh"
using namespace std;



////////////////////////////////////////////////////////////////////////////////
GDLoadBalancer::GDLoadBalancer(MPI_Comm comm, int npex, int npey, int npez):
    comm_(comm), npex_(npex), npey_(npey), npez_(npez), gapthresh_(5)
{
   MPI_Comm_size(comm_, &nTasks_);
   MPI_Comm_rank(comm_, &myRank_);

   npegrid_ = npex_*npey_*npez_;
   
   // reduced process grid dimensions
   tnx_ = tny_ = tnz_ = -1;
   
}
////////////////////////////////////////////////////////////////////////////////
GDLoadBalancer::~GDLoadBalancer()
{
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::initialDistByWorkFn(vector<AnatomyCell>& cells, int nx, int ny, int nz, double diffCost)
{
   // cells should contain only the cells with non-zero type indices, distributed
   // in x-y planes of grid points as follows:
   //    if (nTasks >= nz)     destRank = gid.z   (one plane per task for myRank < nz)
   //    else if (nTasks < nz) destRank = gid.z*nTasks/nz  (multiple planes per task)
   
   diffCost_ = diffCost;

   nx_ = nx;
   ny_ = ny;
   nz_ = nz;

   ownsZData_.resize(nz_);
   for (int iz=0; iz<nz_; iz++)
   {
      if (nTasks_ >= nz_) 
         ownsZData_[iz] = iz;
      else
         ownsZData_[iz] = iz*nTasks_/nz_;
   }

   // distribute grid points on process grid with uniform distribution in z-direction
   vector<int> pezind(nz_,-1);  // z-coordinate in process grid of data plane k
   peyind_.resize(npez_*ny_);
   peyind_.assign(npez_*ny_,-1);
   pexind_.resize(npez_*npey_*nx_);
   pexind_.assign(npez_*npey_*nx_,-1);

   vector<int> kpmin_z(npez_,-1);
   vector<int> kpmax_z(npez_,-1);
   vector<int> kpznum(npez_,0);
   vector<double> kpzwork(npez_,0);

   // zeroth order, calculate work of all groups of 4 z-planes
   vector<double> distWorkLoc(nz_,0.);
   vector<double> distWork(nz_);

   if (true)
   {
      const int zlen = 4;
      int zmin_loc = nz_;
      int zmax_loc = 0;
      vector<int> planecnt_loc(nz_,0);
      vector<int> planecnt(nz_,0);
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         Long64 gid = cells[ii].gid_;
         GridPoint gpt(gid,nx_,ny_,nz_);
         if (gpt.z < zmin_loc) zmin_loc = gpt.z;
         if (gpt.z > zmax_loc) zmax_loc = gpt.z;
         planecnt_loc[gpt.z]++;
      }
      int zmin, zmax;
      MPI_Allreduce(&zmin_loc, &zmin, 1, MPI_INT, MPI_MIN, comm_);
      MPI_Allreduce(&zmax_loc, &zmax, 1, MPI_INT, MPI_MAX, comm_);
      MPI_Allreduce(&planecnt_loc[0], &planecnt[0], nz_, MPI_INT, MPI_SUM, comm_);

      // store min/max values for distributeBox functions
      kpmin_z[0] = zmin;
      kpmax_z[npez_-1] = zmax;
      
      int deltaz = zmax-zmin+1;
      int npez_tmp = ( deltaz%zlen == 0 ? deltaz/zlen : deltaz/zlen + 1);
      
      for (int kp=0; kp<npez_tmp; kp++)
      {
         int tkpmin = zmin + kp*zlen;
         int tkpmax = tkpmin+zlen-1;
         if (tkpmax > zmax) tkpmax = zmax;

         // initial work (under)estimate for process plane kp
         int nkpcnt = 0;
         for (int iz=tkpmin; iz<=tkpmax; iz++)
            nkpcnt += planecnt[iz];         
         int tkpzwork = nkpcnt;

         if (myRank_ == 0)
            cout << "Pre-distributing plane " << kp << ", tkpmin = " << tkpmin << ", tkpmax = " << tkpmax << ", tkpzwork = " << tkpzwork << endl;

         if (myRank_ >= ownsZData_[tkpmin] && myRank_ <= ownsZData_[tkpmax])
         {
            double workFinal = distributePlaneByWorkFn(cells,tkpmin,tkpmax,tkpzwork,kp,false);
            if (myRank_ == ownsZData_[tkpmin])
               for (int kk=tkpmin; kk<=tkpmax; kk++)
                  distWorkLoc[kk] = workFinal;
         }
      }
   }
   MPI_Allreduce(&distWorkLoc[0], &distWork[0], nz_, MPI_DOUBLE, MPI_SUM, comm_);
   
   /* original first pass
   // distribute z planes across npez_ based on estimated work per plane
   distributeBoxByZWork(cells,kpmin_z,kpmax_z,kpznum,kpzwork,pezind);
   
   // now distribute planes over npey_ and npex_
   vector<double> distWork(nz_,0.);
   for (int kp=0; kp<npez_; kp++)
   {
      if (myRank_ == 0)
         cout << "Distributing plane " << kp << ", kpmin_z = " << kpmin_z[kp] << ", kpmax_z = " << kpmax_z[kp] << ", kpzwork = " << kpzwork[kp] << endl;
      double workFinal = distributePlaneByWorkFn(cells,kpmin_z[kp],kpmax_z[kp],kpzwork[kp],kp,true);
      for (int kk=kpmin_z[kp]; kk<=kpmax_z[kp]; kk++)
         distWork[kk] = workFinal;
   }
   */

   
   if (true) // try redistributing using work information from first pass
   {

      // distribute z planes across npez_ based on estimated work per plane
      distributeBoxByWorkDist(cells,kpmin_z,kpmax_z,kpznum,kpzwork,pezind,distWork);
      
      for (int kp=0; kp<npez_; kp++)
      {
         //double workFinal = distributePlaneByWorkFn(cells,kpmin_z[kp],kpmax_z[kp],kpzwork[kp],kp,true);

         if (myRank_ >= ownsZData_[kpmin_z[kp]] && myRank_ <= ownsZData_[kpmax_z[kp]])
         {
            double avgZwork = 0.0;
            for (int iz=kpmin_z[kp]; iz<=kpmax_z[kp]; iz++)
               avgZwork += distWork[iz];
            avgZwork /= (double)(kpmax_z[kp]-kpmin_z[kp]+1);
            if (myRank_ == ownsZData_[kpmin_z[kp]])
               cout << "Redistributing process plane " << kp << ", kpmin_z = " << kpmin_z[kp] << ", kpmax_z = " << kpmax_z[kp] << ", starting work = " << avgZwork/(npex_*npey_) << endl;
            double workFinal = distributePlaneByWorkFn(cells,kpmin_z[kp],kpmax_z[kp],avgZwork,kp,true);
         }
      }
   }

   //ewd DEBUG
   MPI_Barrier(MPI_COMM_WORLD);
   
   // we now have process coordinates for all local grid points, change dest_
   // index in cells array
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      Long64 gid = cells[ii].gid_;
      GridPoint gpt(gid,nx_,ny_,nz_);
      int kp = pezind[gpt.z];
      int jp = peyind_[ny_*kp + gpt.y];
      int ip = pexind_[nx_*npey_*kp + jp*nx_ + gpt.x];
      int peid = ip + jp*npex_ + kp*npex_*npey_;
      cells[ii].dest_ = peid;

      if (ip < 0 || jp < 0 || kp < 0)
         cout << "Grid point x,y,z = " << gpt.x << " " << gpt.y << " " << gpt.z << " has no peid:  ip = " << ip << ", jp = " << jp << ", kp = " << kp << endl;
      
      assert(ip >= 0 && ip < npex_);
      assert(jp >= 0 && jp < npey_);
      assert(kp >= 0 && kp < npez_);
      assert(peid >= 0);
   }

   // carry out communication to match computed distribution
   redistributeCells(cells);
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::distributeBoxByWorkDist(vector<AnatomyCell>& cells, vector<int>& kpmin_z,
                                             vector<int>& kpmax_z, vector<int>& kpznum, 
                                             vector<double>& kpzwork, vector<int>& pezind,
                                             vector<double>& workDist)
{
   // use calculated work information to redistribute planes and try again
   int minz = kpmin_z[0];
   int maxz = kpmax_z[npez_-1];
   double averageWork = 0.;
   for (int kk=minz; kk<=maxz; kk++)
      averageWork += workDist[kk];
   averageWork /= (double)npez_;
      
   int kpset = 0;
   int lastkk = minz;
   double lastwork = -1.0;
   double twork = 0.;
   double worksum = 0.;
   // add planes until work exceeds averageWork
   for (int kk=minz; kk<=maxz; kk++)
   {
      twork += workDist[kk];
      worksum += workDist[kk];
      if (worksum > averageWork*(kpset+1) && kpset < npez_)
      {
         kpmin_z[kpset] = lastkk;
         kpmax_z[kpset] = kk-1;
         if (lastwork > -1.)
            kpzwork[kpset] = lastwork;
         else
            kpzwork[kpset] = twork;
         twork = 0.;
         lastwork = -1.;
         lastkk = kk;
         kpset++;
      }
      else
      {
         lastwork = twork;
      }
   }
   kpmax_z[npez_-1] = maxz;
   if (kpset < npez_)
   {
      kpmin_z[npez_-1] = lastkk;
      kpzwork[npez_-1] = twork;
   }

   // save results in pezind
   for (int kp=0; kp<npez_; kp++)
      for (int iz=kpmin_z[kp]; iz<=kpmax_z[kp]; iz++)
         pezind[iz] = kp;

   return;
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::distributeBoxByZWork(vector<AnatomyCell>& cells, vector<int>& kpmin_z,
                                          vector<int>& kpmax_z, vector<int>& kpznum, 
                                          vector<double>& kpzwork, vector<int>& pezind)
{
   // calculate min/max boundaries of each plane
   vector<int> zmin_xloc(nz_,999999999);
   vector<int> zmin_yloc(nz_,999999999);
   vector<int> zmax_xloc(nz_,-999999999);
   vector<int> zmax_yloc(nz_,-999999999);
   vector<int> kpzcnt_loc(nz_,0);
   
   // determine which xy-planes this process owns data for
   int klocmax = -1;
   int klocmin = nz_+1;
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      if (gpt.z < klocmin) klocmin = gpt.z;
      if (gpt.z > klocmax) klocmax = gpt.z;
      if (gpt.x < zmin_xloc[gpt.z]) zmin_xloc[gpt.z] = gpt.x;
      if (gpt.x > zmax_xloc[gpt.z]) zmax_xloc[gpt.z] = gpt.x;
      if (gpt.y < zmin_yloc[gpt.z]) zmin_yloc[gpt.z] = gpt.y;
      if (gpt.y > zmax_yloc[gpt.z]) zmax_yloc[gpt.z] = gpt.y;
      kpzcnt_loc[gpt.z]++;
   }

   vector<int> zmin_x(nz_),zmin_y(nz_);
   vector<int> zmax_x(nz_),zmax_y(nz_);
   vector<int> kpzcnt(nz_);
   MPI_Allreduce(&zmin_xloc[0], &zmin_x[0], nz_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&zmin_yloc[0], &zmin_y[0], nz_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&zmax_xloc[0], &zmax_x[0], nz_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&zmax_yloc[0], &zmax_y[0], nz_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&kpzcnt_loc[0], &kpzcnt[0], nz_, MPI_INT, MPI_SUM, comm_);

   // divide planes across npez_ 
   vector<double> zwork(nz_);
   int minz = nz_+1;
   int maxz = -1;
   for (int kk=0; kk<nz_; kk++)
   {
      if (kpzcnt[kk] > 0)
      {
         int area = (zmax_x[kk]-zmin_x[kk]+1)*(zmax_y[kk]-zmin_y[kk]+1);
         zwork[kk] = costFunction(kpzcnt[kk], area, 1, diffCost_);
         if (kk < minz) minz = kk;
         if (kk > maxz) maxz = kk;
      }
      else
         zwork[kk] = 0.;
   }

   int kpavg = (maxz-minz+1)/npez_;
   if ((maxz-minz+1)%npez_ != 0) kpavg++;

   vector<int> kpmin_x(npez_,999999999);
   vector<int> kpmin_y(npez_,999999999);
   vector<int> kpmax_x(npez_,-999999999);
   vector<int> kpmax_y(npez_,-999999999);
   int kpset = 0;
   kpmin_z[kpset] = minz;
   for (int kk=minz; kk<=maxz; kk++)
   {
      kpznum[kpset] += kpzcnt[kk];
      if (zmin_x[kk] < kpmin_x[kpset]) kpmin_x[kpset] = zmin_x[kk];
      if (zmin_y[kk] < kpmin_y[kpset]) kpmin_y[kpset] = zmin_y[kk];
      if (zmax_x[kk] > kpmax_x[kpset]) kpmax_x[kpset] = zmax_x[kk];
      if (zmax_y[kk] > kpmax_y[kpset]) kpmax_y[kpset] = zmax_y[kk];
      
      if ((kk+1)%kpavg == 0)
      {
         kpmax_z[kpset] = kk;
         int area = (kpmax_x[kpset]-kpmin_x[kpset]+1)*(kpmax_y[kpset]-kpmin_y[kpset]+1);
         int height = (kk-kpmin_z[kpset]+1);
         kpzwork[kpset] = costFunction(kpznum[kpset], area, height, diffCost_);
         kpset++;
         if (kpset > npez_-1)
            kpset = npez_-1;
         else
            kpmin_z[kpset] = kk+1;
      }
   }

   {
      kpmax_z[npez_-1] = maxz;
      int area = (kpmax_x[npez_-1]-kpmin_x[npez_-1]+1)*(kpmax_y[npez_-1]-kpmin_y[npez_-1]+1);
      int height = (maxz-kpmin_z[npez_-1]+1);
      kpzwork[kpset] = costFunction(kpznum[npez_-1], area, height, diffCost_);
   }
   
   int maxwork = 0;
   for (int kp=0; kp<npez_; kp++)
      if (kpzwork[kp] > maxwork) maxwork = kpzwork[kp];

   vector<int> trialmin_z(npez_);
   vector<int> trialmax_z(npez_);
   vector<int> trialwork_z(npez_,-1);
   bool zworkConverged = false;
   double targetWork = maxwork;
   int viter = 0;
   while (!zworkConverged)
   {
      int kpset = 0;
      int xmin = 999999999;  int ymin = 999999999;
      int xmax = -999999999;  int ymax = -999999999;
      int lastkk = minz;
      int lastwork = -1;
      int knumsum = 0;
      // add planes until work exceeds targetWork
      for (int kk=minz; kk<=maxz; kk++)
      {
         knumsum += kpzcnt[kk];
            
         if (zmin_x[kk] < xmin) xmin = zmin_x[kk];
         if (zmin_y[kk] < ymin) ymin = zmin_y[kk];
         if (zmax_x[kk] > xmax) xmax = zmax_x[kk];
         if (zmax_y[kk] > ymax) ymax = zmax_y[kk];
         int area = (xmax-xmin+1)*(ymax-ymin+1);
         int height = (kk-lastkk+1);
         int twork = costFunction(knumsum, area, height, diffCost_);
         
         if (twork > targetWork && kpset < npez_)
         {
            trialmin_z[kpset] = lastkk;
            trialmax_z[kpset] = kk-1;
            if (lastwork > -1)
               trialwork_z[kpset] = lastwork;
            else
            {
               int tmpsum = 0;
               for (int tmpk=kk; tmpk<lastkk; tmpk++)
                  tmpsum += kpzcnt[tmpk];
               trialwork_z[kpset] = costFunction(tmpsum, area, kk-lastkk, diffCost_);
            }
            lastwork = -1;
            knumsum = 0;
            lastkk = kk;
            kpset++;
            xmin = zmin_x[kk];
            ymin = zmin_y[kk];
            xmax = zmax_x[kk];
            ymax = zmax_y[kk];
         }
         else
         {
            lastwork = twork;
         }
      }
      trialmax_z[npez_-1] = maxz;
      if (kpset < npez_)
         trialmin_z[npez_-1] = lastkk;
      
      int area = (xmax-xmin+1)*(ymax-ymin+1);
      int height = (maxz-trialmin_z[npez_-1]+1);
      trialwork_z[npez_-1] = costFunction(knumsum, area, height, diffCost_);

      int tworkmax = -1;
      int tworkmin = 999999999;
      for (int kp=0; kp<npez_; kp++)
      {
         if (trialwork_z[kp] > tworkmax) tworkmax = trialwork_z[kp];
         if (trialwork_z[kp] < tworkmin) tworkmin = trialwork_z[kp];
      }
      
      if (myRank_ == 0)
         cout << "Plane distribution: iteration " << viter++ << ", targetWork = " << targetWork << ", max work = " << tworkmax << ", tworkmin = " << tworkmin << endl;

      if ( (tworkmax < targetWork && tworkmax != tworkmin ) || tworkmin < 0 ) // keep going
      {
         for (int kp=0; kp<npez_; kp++)
         {
            kpmin_z[kp] = trialmin_z[kp];
            kpmax_z[kp] = trialmax_z[kp];
            kpzwork[kp] = trialwork_z[kp];
         }
      }
      else
      {
         if (myRank_ == 0)
         {
            cout << "Convergence reached." << endl;
            for (int kp=0; kp<npez_; kp++)
               cout << "  kp " << kp << ":  " << kpmin_z[kp] << " " << kpmax_z[kp] << ", avg work = " << kpzwork[kp]/(npey_*npex_) << endl;
         }
         
         // save results in pezind
         for (int kp=0; kp<npez_; kp++)
            for (int iz=kpmin_z[kp]; iz<=kpmax_z[kp]; iz++)
               pezind[iz] = kp;
         
         zworkConverged = true;
      }
      targetWork = 0.99*tworkmax;
   }
   return;
}
////////////////////////////////////////////////////////////////////////////////
double GDLoadBalancer::distributePlaneByWorkFn(vector<AnatomyCell>& cells, int zmin,
                                               int zmax, double zwork, int kp,
                                               bool assignpes)
{
   int kpdz = zmax - zmin + 1;
   vector<int> kpxyloc(nx_*ny_,0);
   int nloc = 0;
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      if (gpt.z >= zmin && gpt.z <= zmax)
      {
         kpxyloc[nx_*gpt.y+gpt.x]++;
         nloc++;
      }
   }
   //int size = nx_*ny_;
   //vector<int> kpxycnt(nx_*ny_);
   //MPI_Allreduce(&kpxyloc[0], &kpxycnt[0], size, MPI_INT, MPI_SUM, comm_);

   // sum kpxyloc over all tasks that own data for this process plane
   set<int> kpProcs;
   for (int iz=zmin; iz<=zmax; iz++)
   {
      int pe = ownsZData_[iz];
      kpProcs.insert(pe);
   }
   int firstPe = ownsZData_[zmin];
   
   vector<int> kpxycnt(nx_*ny_,0);
   int nkpProcs = kpProcs.size();
   for (set<int>::iterator zpeit = kpProcs.begin(); zpeit != kpProcs.end(); ++zpeit)
   {
      int size = nx_*ny_;
      vector<int> buf(size);
      if (myRank_ != *zpeit)
      {
         MPI_Status status;
         MPI_Recv(&buf[0],size,MPI_INT,*zpeit,myRank_,MPI_COMM_WORLD,&status);
         for (int ii=0; ii<size; ii++)
            kpxycnt[ii] += buf[ii];
      }
      else
      {
         for (set<int>::iterator destpeit = kpProcs.begin(); destpeit != kpProcs.end(); ++destpeit)
         {
            if (*destpeit != *zpeit)
               MPI_Send(&kpxyloc[0],size,MPI_INT,*destpeit,*destpeit,MPI_COMM_WORLD);
         }
         for (int ii=0; ii<size; ii++)
            kpxycnt[ii] += kpxyloc[ii];
      }
   }

   // kpxycnt should contain xy-distribution on all tasks that own data for this process plane

   // compute min/max values for all planes, strips
   vector<int> kpymin_x(ny_,999999999);
   vector<int> kpymax_x(ny_,-999999999);
   int xmin = 999999999;  int ymin = 999999999;
   int xmax = -999999999;  int ymax = -999999999;
   for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<nx_; ix++)
         if (kpxycnt[nx_*iy+ix] > 0)
         {
            if (ix < xmin) xmin = ix;
            if (ix > xmax) xmax = ix;
            if (iy < ymin) ymin = iy;
            if (iy > ymax) ymax = iy;
            if (ix < kpymin_x[iy]) kpymin_x[iy] = ix;
            if (ix > kpymax_x[iy]) kpymax_x[iy] = ix;
         }

   double targetWork = (double)zwork/(npex_*npey_);

   const int ninner = 64;
   vector<double> tWork(ninner);
   vector<double> maxWork(ninner);
   //ewd   double minW = 0.1*targetWork;
   double minW = 0.02*targetWork;
   double deltaW = targetWork;
   double lastWork = -1;
   const double threshold = 0.01;  // converge to within one percent
   for (int jj = 0; jj<40; jj++)
   {
      #pragma omp parallel for
      for (int ii=0; ii<ninner; ii++)
      {
         tWork[ii] = minW + ii*deltaW;
         double wsum;
         maxWork[ii] = calcMaxWork(tWork[ii],kpxycnt,ymin,ymax,kpdz,kp,wsum,false);
      }   
      
      double tmax = 1.E+19;
      int maxind = -1;
      for (int ii=0; ii<ninner; ii++)
      {
         if (maxWork[ii] < tmax)
         {
            tmax = maxWork[ii];
            maxind = ii;
         }
      }
      double change = abs(tmax - lastWork)/tmax;

      //ewd DEBUG
      //if (myRank_ == firstPe)
      //   cout << "GDLB.DEBUG: kp = " << kp << ", maxWork = " << tmax << ", lastWork = " << lastWork << ", change = " << change << ", threshold = " << threshold << endl;

      if (change < threshold) {
         targetWork = tmax;
         break;
      }
      lastWork = tmax;
      
      minW = (maxind > 0 ? tWork[maxind-1] : 0.1*tWork[0]);
      double maxW = (maxind < ninner-1 ? tWork[maxind+1] : 10.*tWork[ninner-1]);
      deltaW = (maxW - minW)/ninner;
   }

   
   double worksum;
   double wmax = calcMaxWork(targetWork,kpxycnt,ymin,ymax,kpdz,kp,worksum,assignpes);

   //ewd DEBUG
   //if (myRank_ == firstPe)
   //   cout << "    GDLB, mype = " << myRank_ << ", kp = " << kp << ", tworkmax = " << tworkmax << ", tworkmin = " << tworkmin << ", targetWork = " << targetWork << endl;
   //ewd DEBUG
      
   return worksum; // return total work
}

////////////////////////////////////////////////////////////////////////////////
double GDLoadBalancer::calcMaxWork(double targetWork, vector<int>& kpxycnt,int ymin,
                                   int ymax, int kpdz, int kp, double& twksum, bool assignpes)
{
   vector<int> trialmin_y(npey_,-1);
   vector<int> trialmax_y(npey_,-1);
   vector<int> trialmin_x(npex_*npey_,-1);
   vector<int> trialmax_x(npex_*npey_,-1);
   vector<int> tmin_x(npex_);
   vector<int> tmax_x(npex_);
   vector<double> tjpwork(npey_,-1.);
   int jpset = 0;
   trialmin_y[jpset] = ymin;
   trialmax_y[jpset] = ymin;
   twksum = 0.0;
   double laststripwork;
   for (int iy=ymin; iy<=ymax; iy++)
   {
      // try adding this iy strip to current jpset process strip
      vector<int> txcnt(nx_,0);
      int xsum = 0;
      for (int ty=trialmin_y[jpset]; ty<=iy; ty++)
      {
         for (int ix=0; ix<nx_; ix++)
         {
            txcnt[ix] += kpxycnt[nx_*ty+ix];
            xsum += kpxycnt[nx_*ty+ix];
         }
      }
      int dy = (iy-trialmin_y[jpset]+1);

      if (xsum == 0) // don't bother calculating empty distribution
         continue;
      
      xstripDistByWorkFn(txcnt,tmin_x,tmax_x,dy,kpdz,false);
      
      // calculate maximum work of (iy-trialmin_y) strips distributed over npex_ tasks
      double maxwork = -1.0;
      double stripwork = 0.0;
      for (int ip=0; ip<npex_; ip++)
      {
         if (!(tmin_x[ip] < 0 || tmax_x[ip] < 0))
         {
            int area = (tmax_x[ip]-tmin_x[ip]+1)*(iy-trialmin_y[jpset]+1);
            int tmpcnt = 0;
            for (int tx=tmin_x[ip]; tx<=tmax_x[ip]; tx++)
               tmpcnt += txcnt[tx];
            double twork = costFunction(tmpcnt, area, kpdz, diffCost_);
            if (twork > maxwork) maxwork = twork;
            stripwork += twork;
         }
      }

      if (maxwork > targetWork)
      {
         twksum += laststripwork;
         jpset++;
         if (jpset > npey_-1)
         {
            jpset = npey_-1;
            tjpwork[jpset] = maxwork;
         }
         else
         {
            trialmin_y[jpset] = iy;
         }
      }
      else {
         tjpwork[jpset] = maxwork;
         laststripwork = stripwork;
      }
      // save latest stats on this process grid strip (at jp,kp)
      trialmax_y[jpset] = iy;
      for (int ip=0; ip<npex_; ip++)
      {
         trialmin_x[jpset*npex_+ip] = tmin_x[ip];
         trialmax_x[jpset*npex_+ip] = tmax_x[ip];
      }            
   }
      
   // update npey_-1 term
   if (jpset == npey_-1)
   {
      trialmax_y[npey_-1] = ymax;
      twksum += laststripwork;
   }
   
   //recalculate work of npey_-1 term
   {
      double maxwork = -1.0;
      for (int ip=0; ip<npex_; ip++)
      {
         int jpind = (npey_-1)*npex_+ip;
         
         if (!(trialmin_x[jpind] < 0 || trialmax_x[jpind] < 0))
         {
            int tmpcnt = 0;
            for (int ty=trialmin_y[npey_-1]; ty<=ymax; ty++)
               for (int tx=trialmin_x[jpind]; tx<=trialmax_x[jpind]; tx++)
                  tmpcnt += kpxycnt[nx_*ty+tx];
            int area = (trialmax_x[jpind]-trialmin_x[jpind]+1)*(ymax-trialmin_y[npey_-1]+1);
            
            double twork = costFunction(tmpcnt, area, kpdz, diffCost_);
            if (twork > maxwork) maxwork = twork;
         }
      }
      if (maxwork > tjpwork[npey_-1]) tjpwork[npey_-1] = maxwork;
   }

      
   // all grid strips distributed, check work distribution to see if target work was exceeded
   double tworkmax = -1;
   double tworkmin = 999999999;
   for (int jp=0; jp<npey_; jp++)
   {
      if (tjpwork[jp] > tworkmax) tworkmax = tjpwork[jp];
      if (tjpwork[jp] < tworkmin) tworkmin = tjpwork[jp];
   }


   if (assignpes)
   {
      for (int jp=0; jp<npey_; jp++)
      {
         for (int iy=trialmin_y[jp]; iy<=trialmax_y[jp]; iy++)
            if (iy > -1)
               peyind_[ny_*kp+iy] = jp;
         
         for (int ip=0; ip<npex_; ip++)
            for (int ix=trialmin_x[jp*npex_+ip]; ix<=trialmax_x[jp*npex_+ip]; ix++)
               if (ix > -1)
                  pexind_[nx_*npey_*kp + nx_*jp + ix] = ip;
      }
   }

   return tworkmax;
}

////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::xstripDistByWorkFn(vector<int>& xcnt, vector<int>& pexmin, vector<int>& pexmax,
                                int dy, int dz, bool verbose)
{
   int xmin = 999999999;
   int xmax = -999999999;
   int xsum = 0;
   for (int ix=0; ix<nx_; ix++)
   {
      if (xcnt[ix] > 0 && xmin == 999999999)
         xmin = ix;
      if (xcnt[ix] > 0)
         xmax = ix;
      xsum += xcnt[ix];
   }
   for (int ip=0; ip<npex_; ip++)
   {
      pexmin[ip] = -1;
      pexmax[ip] = -2;
   }

   if (xsum == 0)
   {
      cout << "ERROR:  xstripDistByWorkFn called with no work, exiting." << endl;
      return;
   }
   
   int ngap = 0;
   int gapcnt = 0;
   int isumgap = -1;
   vector<double> bwork;
   vector<int> gapii(1,xmin);
   vector<int> bstart(1,xmin);
   vector<int> bend(1,xmax);
   vector<int> bcnt(1,0);
   for (int ii=xmin; ii<=xmax; ii++) 
   {
      bcnt[ngap] += xcnt[ii];
      if (isumgap == -1 && xcnt[ii] > 0)
         bstart[ngap] = ii;
      if (xcnt[ii] > 0)
         bend[ngap] = ii;
      isumgap += xcnt[ii];
      if (isumgap > -1 && xcnt[ii] == 0)
         gapcnt++;
      else
         gapcnt = 0;
      
      if (gapcnt > gapthresh_) {
         gapii.push_back(ii-gapcnt);
         int area = (bend[ngap]-bstart[ngap]+1)*dy;
         double tmpwork = costFunction(bcnt[ngap],area,dz,diffCost_);
         bwork.push_back(tmpwork);
         bstart.push_back(ii);
         bend.push_back(xmax);
         bcnt.push_back(0);
         ngap++;
         gapcnt = 0;
         isumgap = -1;
      }
   }
   
   int tmparea = (bend[ngap]-bstart[ngap]+1)*dy;
   double tmpwork = costFunction(bcnt[ngap],tmparea,dz,diffCost_);
   bwork.push_back(tmpwork);
   if (isumgap == -1 && ngap > 0)  // remove "gap" of empty space between last cells and box edge
      ngap--;
   gapii.push_back(xmax);
   
   // assign npex_ processors to each chunk, starting with the smallest chunks first
   vector<int> sortedWorkIndex(ngap+1);
   for (int ig=0; ig<=ngap; ig++)
      sortedWorkIndex[ig] = ig;
   sort(sortedWorkIndex.begin(),sortedWorkIndex.end(),GapSortDouble(bwork));
   vector<int> bnpes(ngap+1);
   int peleft = npex_;
   
   double worktot = 0.;
   for (int ig=0; ig<=ngap; ig++)
      worktot += bwork[ig];
   double avgwork = worktot/(double)npex_;
   for (int ig=0; ig<ngap; ig++)
   {
      int sind = sortedWorkIndex[ig];
      double gwork = bwork[sind];
      int gpes = gwork/avgwork + 1;
      bnpes[sind] = gpes;
      peleft -= gpes;
   }
   int sind = sortedWorkIndex[ngap];
   bnpes[sind] = peleft;

   //ewd DEBUG
   if (verbose)
   {
      cout << "XSTRIPDIST:  initial array:  "  << endl;
      for (int ix=0; ix<nx_; ix++)
         cout << "  " << ix << "   " << xcnt[ix] << endl;

      cout << "XSTRIPDIST:  ngap = " << ngap << endl;
      for (unsigned ig=0; ig<bstart.size(); ig++)
         cout << "   " << ig << "    " << bstart[ig] << "  " << bend[ig] << "   " << bwork[ig] << "   " << sortedWorkIndex[ig] << "   " << bnpes[ig] << endl;
   }
   //ewd DEBUG
   
      
   // loop through strip again and assign bnpes to each chunk of cells
   int ipset = 0;
   for (int ig=0; ig<=ngap; ig++)
   {
      pexmin[ipset] = bstart[ig];
      pexmax[ipset] = bstart[ig];
      double gavg = (double)bwork[ig]/(double)bnpes[ig] + 0.01;
      if (gavg < 1.01) gavg = 1.01;
      
      int gip = 0;
      int thisxcnt = 0;
      double worksum = 0.0;
      for (int ii=bstart[ig]; ii<=bend[ig]; ii++) 
      {
         thisxcnt += xcnt[ii];
         int area = (ii-pexmin[ipset])*dy;
         double thiswork = costFunction(thisxcnt,area,dz,diffCost_);

         if ((thiswork+worksum) > gavg*(gip+1)) {
            ipset++;
            worksum += thiswork;
            thisxcnt = 0;
            if (ipset > npex_-1)
               ipset = npex_-1;
            else
               pexmin[ipset] = ii;
            gip++;
         }
         if (xcnt[ii] > 0) 
            pexmax[ipset] = ii;
      }
      ipset++;
      if (ipset > npex_-1) ipset = npex_-1;
   }
   if (pexmin[npex_-1] > 0)
      pexmax[npex_-1] = xmax;

   //ewd DEBUG      
   if (verbose)
   {
      cout << "calculated distribution:" << endl;
      for (int ip=0; ip<npex_; ip++) {
         cout << "  " << ip << "    " << pexmin[ip] << "   " << pexmax[ip] << endl;
      }
   }
   //ewd DEBUG
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::initialDistByVol(vector<AnatomyCell>& cells, int nx, int ny, int nz)
{
   // cells should contain only the cells with non-zero type indices, distributed
   // in x-y planes of grid points as follows:
   //    if (nTasks >= nz)     destRank = gid.z   (one plane per task for myRank < nz)
   //    else if (nTasks < nz) destRank = gid.z*nTasks/nz  (multiple planes per task)
   
   int nlocal = cells.size();
   nx_ = nx;
   ny_ = ny;
   nz_ = nz;

   // distribute grid points on process grid with uniform distribution in z-direction
   vector<int> pezind(nz_,-1);  // z-coordinate in process grid of data plane k
   peyind_.resize(npez_*ny_);
   peyind_.assign(npez_*ny_,-1);
   pexind_.resize(npez_*npey_*nx_);
   pexind_.assign(npez_*npey_*nx_,-1);

   // calculate min/max boundaries of each plane
   vector<int> zmin_xloc(nz_,999999999);
   vector<int> zmin_yloc(nz_,999999999);
   vector<int> zmax_xloc(nz_,-999999999);
   vector<int> zmax_yloc(nz_,-999999999);
   vector<int> kpzcnt_loc(nz_,0);
   
   // determine which xy-planes this process owns data for
   int klocmax = -1;
   int klocmin = nz_+1;
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      if (gpt.z < klocmin) klocmin = gpt.z;
      if (gpt.z > klocmax) klocmax = gpt.z;
      if (gpt.x < zmin_xloc[gpt.z]) zmin_xloc[gpt.z] = gpt.x;
      if (gpt.x > zmax_xloc[gpt.z]) zmax_xloc[gpt.z] = gpt.x;
      if (gpt.y < zmin_yloc[gpt.z]) zmin_yloc[gpt.z] = gpt.y;
      if (gpt.y > zmax_yloc[gpt.z]) zmax_yloc[gpt.z] = gpt.y;
      kpzcnt_loc[gpt.z]++;
   }

   vector<int> zmin_x(nz_),zmin_y(nz_);
   vector<int> zmax_x(nz_),zmax_y(nz_);
   vector<int> kpzcnt(nz_);
   MPI_Allreduce(&zmin_xloc[0], &zmin_x[0], nz_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&zmin_yloc[0], &zmin_y[0], nz_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&zmax_xloc[0], &zmax_x[0], nz_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&zmax_yloc[0], &zmax_y[0], nz_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&kpzcnt_loc[0], &kpzcnt[0], nz_, MPI_INT, MPI_SUM, comm_);

   // divide planes across npez_ 
   vector<int> zvol(nz_);
   int minz = nz_+1;
   int maxz = -1;
   for (int kk=0; kk<nz_; kk++)
   {
      if (kpzcnt[kk] > 0)
      {
         zvol[kk] = (zmax_x[kk]-zmin_x[kk]+1)*(zmax_y[kk]-zmin_y[kk]+1);
         if (kk < minz) minz = kk;
         if (kk > maxz) maxz = kk;
      }
      else
         zvol[kk] = 0;
   }

   int kpavg = (maxz-minz+1)/npez_;
   if ((maxz-minz+1)%npez_ != 0) kpavg++;

   vector<int> kpmin_x(npez_,999999999);
   vector<int> kpmin_y(npez_,999999999);
   vector<int> kpmin_z(npez_,-1);
   vector<int> kpmax_x(npez_,-999999999);
   vector<int> kpmax_y(npez_,-999999999);
   vector<int> kpmax_z(npez_,-1);
   vector<int> kpzvol(npez_,0);
   int kpset = 0;
   kpmin_z[kpset] = minz;
   for (int kk=minz; kk<=maxz; kk++)
   {
      if (zmin_x[kk] < kpmin_x[kpset]) kpmin_x[kpset] = zmin_x[kk];
      if (zmin_y[kk] < kpmin_y[kpset]) kpmin_y[kpset] = zmin_y[kk];
      if (zmax_x[kk] > kpmax_x[kpset]) kpmax_x[kpset] = zmax_x[kk];
      if (zmax_y[kk] > kpmax_y[kpset]) kpmax_y[kpset] = zmax_y[kk];
      
      //if (kk > 0 && kk%kpavg == 0)
      if ((kk+1)%kpavg == 0)
      {
         kpmax_z[kpset] = kk;
         kpzvol[kpset] = (kpmax_x[kpset]-kpmin_x[kpset]+1)*(kpmax_y[kpset]-kpmin_y[kpset]+1)*(kk-kpmin_z[kpset]+1);
         kpset++;
         if (kpset > npez_-1)
            kpset = npez_-1;
         else
            kpmin_z[kpset] = kk+1;
      }
   }

   kpmax_z[npez_-1] = maxz;
   kpzvol[npez_-1] = (kpmax_x[npez_-1]-kpmin_x[npez_-1]+1)*(kpmax_y[npez_-1]-kpmin_y[npez_-1]+1)*(maxz-kpmin_z[npez_-1]+1);
   
   int maxvol = 0;
   for (int kp=0; kp<npez_; kp++)
      if (kpzvol[kp] > maxvol) maxvol = kpzvol[kp];

   // try recalculating planar distribution to achieve lower average and maximum volumes
   {
      vector<int> trialmin_z(npez_);
      vector<int> trialmax_z(npez_);
      vector<int> trialvol_z(npez_,-1);
      bool zvolConverged = false;
      double targetVol = maxvol;
      int viter = 0;
      while (!zvolConverged)
      {
         targetVol *= 0.99;
         int kpset = 0;
         int xmin = 999999999;  int ymin = 999999999;
         int xmax = -999999999;  int ymax = -999999999;
         int lastkk = minz;
         int lastvol = -1;
         // add planes until volume exceeds targetVol
         for (int kk=minz; kk<=maxz; kk++)
         {
            if (zmin_x[kk] < xmin) xmin = zmin_x[kk];
            if (zmin_y[kk] < ymin) ymin = zmin_y[kk];
            if (zmax_x[kk] > xmax) xmax = zmax_x[kk];
            if (zmax_y[kk] > ymax) ymax = zmax_y[kk];
            int tvol = (xmax-xmin+1)*(ymax-ymin+1)*(kk-lastkk+1);

            if (tvol > targetVol && kpset < npez_)
            {
               trialmin_z[kpset] = lastkk;
               trialmax_z[kpset] = kk-1;
               if (lastvol > -1)
                  trialvol_z[kpset] = lastvol;
               else
                  trialvol_z[kpset] = (xmax-xmin+1)*(ymax-ymin+1)*(kk-lastkk);                  
               lastvol = -1;
               lastkk = kk;
               kpset++;
               xmin = zmin_x[kk];
               ymin = zmin_y[kk];
               xmax = zmax_x[kk];
               ymax = zmax_y[kk];
            }
            else
            {
               lastvol = tvol;
            }
         }
         trialmax_z[npez_-1] = maxz;
         if (kpset < npez_)
            trialmin_z[npez_-1] = lastkk;
         trialvol_z[npez_-1] = (xmax-xmin+1)*(ymax-ymin+1)*(maxz-trialmin_z[npez_-1]+1);

         int tvolmax = -1;
         int tvolmin = 999999999;
         for (int kp=0; kp<npez_; kp++)
         {
            if (trialvol_z[kp] > tvolmax) tvolmax = trialvol_z[kp];
            if (trialvol_z[kp] < tvolmin) tvolmin = trialvol_z[kp];
         }
      
         if (myRank_ == 0)
            cout << "Plane distribution: iteration " << viter++ << ", targetVol = " << targetVol << ", max volume = " << tvolmax << ", tvolmin = " << tvolmin << endl;

         if ( (tvolmax < targetVol && tvolmax != tvolmin ) || tvolmin < 0 ) // keep going
         {
            for (int kp=0; kp<npez_; kp++)
            {
               kpmin_z[kp] = trialmin_z[kp];
               kpmax_z[kp] = trialmax_z[kp];
               kpzvol[kp] = trialvol_z[kp];
            }
         }
         else
         {
            if (myRank_ == 0)
            {
               cout << "Convergence reached." << endl;
               for (int kp=0; kp<npez_; kp++)
                  cout << "  kp " << kp << ":  " << kpmin_z[kp] << " " << kpmax_z[kp] << ", tot vol = " << kpzvol[kp] << ", avg vol = " << kpzvol[kp]/(npey_*npex_) << endl;
            }

            // save results in pezind
            for (int kp=0; kp<npez_; kp++)
               for (int iz=kpmin_z[kp]; iz<=kpmax_z[kp]; iz++)
                  pezind[iz] = kp;
            
            zvolConverged = true;
         }
      }
   }
   
   // now distribute planes over npey_ and npex_
   vector<int> distVol(nz_,0);
   for (int kp=0; kp<npez_; kp++)
   {
      if (myRank_ == 0)
         cout << "Distributing plane " << kp << ", kpzvol = " << kpzvol[kp] << ", kpmin_z = " << kpmin_z[kp] << ", kpmax_z = " << kpmax_z[kp] << endl;
      int volFinal;
      volFinal = distributePlaneByVol(cells,kpmin_z[kp],kpmax_z[kp],kpzvol[kp],kp);
      if (volFinal < 0)
      {
         if (myRank_ == 0)
            cout << "Increasing target volume and trying again...";
         int tmpvol = 2*kpzvol[kp];
         volFinal = distributePlaneByVol(cells,kpmin_z[kp],kpmax_z[kp],kpzvol[kp],kp);
      }
      for (int kk=kpmin_z[kp]; kk<=kpmax_z[kp]; kk++)
         distVol[kk] = volFinal;
   }


   // calculate volume of each plane
   if (false)
   {
      vector<int> cell_xmin(npegrid_,999999999);
      vector<int> cell_ymin(npegrid_,999999999);
      vector<int> cell_zmin(npegrid_,999999999);
      vector<int> cell_xmax(npegrid_,-999999999);
      vector<int> cell_ymax(npegrid_,-999999999);
      vector<int> cell_zmax(npegrid_,-999999999);
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         Long64 gid = cells[ii].gid_;
         GridPoint gpt(gid,nx_,ny_,nz_);
         int kp = pezind[gpt.z];
         int jp = peyind_[ny_*kp + gpt.y];
         int ip = pexind_[nx_*npey_*kp + jp*nx_ + gpt.x];
         int peid = ip + jp*npex_ + kp*npex_*npey_;

         if (gpt.x < cell_xmin[peid]) cell_xmin[peid] = gpt.x;
         if (gpt.y < cell_ymin[peid]) cell_ymin[peid] = gpt.y;
         if (gpt.z < cell_zmin[peid]) cell_zmin[peid] = gpt.z;
         if (gpt.x > cell_xmax[peid]) cell_xmax[peid] = gpt.x;
         if (gpt.y > cell_ymax[peid]) cell_ymax[peid] = gpt.y;
         if (gpt.z > cell_zmax[peid]) cell_zmax[peid] = gpt.z;
         
         if (ip < 0 || jp < 0 || kp < 0)
            cout << "Grid point x,y,z = " << gpt.x << " " << gpt.y << " " << gpt.z << " has no peid:  ip = " << ip << ", jp = " << jp << ", kp = " << kp << endl;
         
         assert(ip >= 0 && ip < npex_);
         assert(jp >= 0 && jp < npey_);
         assert(kp >= 0 && kp < npez_);
         assert(peid >= 0);
      }


      vector<int> cell_xmin_all(npegrid_);
      vector<int> cell_ymin_all(npegrid_);
      vector<int> cell_zmin_all(npegrid_);
      vector<int> cell_xmax_all(npegrid_);
      vector<int> cell_ymax_all(npegrid_);
      vector<int> cell_zmax_all(npegrid_);
      MPI_Allreduce(&cell_xmin[0], &cell_xmin_all[0], npegrid_, MPI_INT, MPI_MIN, comm_);
      MPI_Allreduce(&cell_ymin[0], &cell_ymin_all[0], npegrid_, MPI_INT, MPI_MIN, comm_);
      MPI_Allreduce(&cell_zmin[0], &cell_zmin_all[0], npegrid_, MPI_INT, MPI_MIN, comm_);
      MPI_Allreduce(&cell_xmax[0], &cell_xmax_all[0], npegrid_, MPI_INT, MPI_MAX, comm_);
      MPI_Allreduce(&cell_ymax[0], &cell_ymax_all[0], npegrid_, MPI_INT, MPI_MAX, comm_);
      MPI_Allreduce(&cell_zmax[0], &cell_zmax_all[0], npegrid_, MPI_INT, MPI_MAX, comm_);
      
      vector<int> slabvol(npez_,0);
      vector<int> slabvolsum(npez_,0);
      for (int ip=0; ip<npegrid_; ip++)
      {
         if (cell_xmin[ip] > 999999997) cell_xmin[ip] = 0;
         if (cell_ymin[ip] > 999999997) cell_ymin[ip] = 0;
         if (cell_zmin[ip] > 999999997) cell_zmin[ip] = 0;
         if (cell_xmax[ip] < 0) cell_xmax[ip] = -1;
         if (cell_ymax[ip] < 0) cell_ymax[ip] = -1;
         if (cell_zmax[ip] < 0) cell_zmax[ip] = -1;
         
         int procvol = (cell_xmax[ip]-cell_xmin[ip]+1)*(cell_ymax[ip]-cell_ymin[ip]+1)*(cell_zmax[ip]-cell_zmin[ip]+1);
         int kp = ip/(npex_*npey_);
         slabvol[kp] += procvol;
      }
      MPI_Allreduce(&slabvol[0], &slabvolsum[0], npez_, MPI_INT, MPI_SUM, comm_);

      for (int kp=0; kp<npez_; kp++)
         for (int kk=kpmin_z[kp]; kk<=kpmax_z[kp]; kk++)
            distVol[kk] = slabvolsum[kp];
   }
   
   int minBalVol = 999999999;
   int maxBalVol = -999999999;
   for (int kk=minz; kk<=maxz; kk++)
   {
      if (distVol[kk] < minBalVol) minBalVol = distVol[kk];
      if (distVol[kk] > maxBalVol) maxBalVol = distVol[kk];
   }
   const double imbalThreshold = 0.3;
   const double imbalCurrent = (double)(maxBalVol-minBalVol)/minBalVol;
   
   // redistribute planes with new information
   if (false && imbalCurrent > imbalThreshold)
   {
      if (myRank_ == 0)
         cout << "Planar imbalance = " << imbalCurrent << ", exceeding threshold of " << imbalThreshold << ", redistribute planes and repeat load balancing..." << endl;

      vector<int> trialmin_z(npez_);
      vector<int> trialmax_z(npez_);
      vector<int> trialvol_z(npez_,-1);

      double targetVol = 0.;
      for (int kk=minz; kk<=maxz; kk++)
         //targetVol += distVol[kk]*npex_*npey_;
         targetVol += distVol[kk];
      targetVol /= (double)npez_;
          
      int kpset = 0;
      int lastkk = minz;
      int lastvol = -1;
      int tvol = 0;
      // add planes until volume exceeds targetVol
      for (int kk=minz; kk<=maxz; kk++)
      {
         //tvol += distVol[kk]*npex_*npey_;
         tvol += distVol[kk];
         if (tvol > targetVol && kpset < npez_)
         {
            tvol = 0;
            trialmin_z[kpset] = lastkk;
            trialmax_z[kpset] = kk-1;
            if (lastvol > -1)
               trialvol_z[kpset] = lastvol;
            else
               trialvol_z[kpset] = tvol;
            lastvol = -1;
            lastkk = kk;
            kpset++;
         }
         else
         {
            lastvol = tvol;
         }
      }
      trialmax_z[npez_-1] = maxz;
      if (kpset < npez_)
         trialmin_z[npez_-1] = lastkk;
      if (trialvol_z[npez_-1] < 0)
         trialvol_z[npez_-1] = tvol;


      for (int kp=0; kp<npez_; kp++)
      {
         if (myRank_ == 0)
            cout << "Redistributing plane " << kp << ", trialvol = " << trialvol_z[kp] << ", zmin = " << trialmin_z[kp] << ", zmax = " << trialmax_z[kp] << endl;
         int volFinal;
         if (trialvol_z[kp] > 0)
         {
            volFinal = distributePlaneByVol(cells,trialmin_z[kp],trialmax_z[kp],trialvol_z[kp],kp);
            if (volFinal < 0)
            {
               if (myRank_ == 0)
                  cout << "Increasing target volume and trying again...";
               trialvol_z[kp] *= 2;
               volFinal = distributePlaneByVol(cells,trialmin_z[kp],trialmax_z[kp],trialvol_z[kp],kp);
            }
         }
         else
            if (myRank_ == 0)
               cout << "Redistributing:  plane " << kp << ", has zero target volume!" << endl;
      }

   }

   
   // we now have process coordinates for all local grid points, change dest_
   // index in cells array
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      Long64 gid = cells[ii].gid_;
      GridPoint gpt(gid,nx_,ny_,nz_);
      int kp = pezind[gpt.z];
      int jp = peyind_[ny_*kp + gpt.y];
      int ip = pexind_[nx_*npey_*kp + jp*nx_ + gpt.x];
      int peid = ip + jp*npex_ + kp*npex_*npey_;
      cells[ii].dest_ = peid;

      if (ip < 0 || jp < 0 || kp < 0)
         cout << "Grid point x,y,z = " << gpt.x << " " << gpt.y << " " << gpt.z << " has no peid:  ip = " << ip << ", jp = " << jp << ", kp = " << kp << endl;
      
      assert(ip >= 0 && ip < npex_);    //ewd DEBUG
      assert(jp >= 0 && jp < npey_);    //ewd DEBUG
      assert(kp >= 0 && kp < npez_);    //ewd DEBUG
      assert(peid >= 0);  //ewd DEBUG
   }

   // carry out communication to match computed distribution
   redistributeCells(cells);

}
////////////////////////////////////////////////////////////////////////////////
int GDLoadBalancer::distributePlaneByVol(vector<AnatomyCell>& cells, int zmin, int zmax, int zvol, int kp)
{
   vector<int> jpmin_y(npey_,-1);
   vector<int> jpmax_y(npey_,-1);      
   vector<int> jpmin_x(npex_*npey_,-1);
   vector<int> jpmax_x(npex_*npey_,-1);
   vector<int> jpvol(npey_,-1);      

   vector<int> oldjpmin_y(npey_,-1);
   vector<int> oldjpmax_y(npey_,-1);      
   vector<int> oldjpmin_x(npex_*npey_,-1);
   vector<int> oldjpmax_x(npex_*npey_,-1);
   vector<int> oldjpvol(npey_,-1);      
      
   int kpdz = zmax - zmin + 1;
   vector<int> kpxyloc(nx_*ny_,0);
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      if (gpt.z >= zmin && gpt.z <= zmax)
         kpxyloc[nx_*gpt.y+gpt.x]++;
   }
   int size = nx_*ny_;
   vector<int> kpxycnt(nx_*ny_);
   MPI_Allreduce(&kpxyloc[0], &kpxycnt[0], size, MPI_INT, MPI_SUM, comm_);

   // compute min/max values for all planes, strips
   vector<int> kpymin_x(ny_,999999999);
   vector<int> kpymax_x(ny_,-999999999);
   int xmin = 999999999;  int ymin = 999999999;
   int xmax = -999999999;  int ymax = -999999999;
   for (int iy=0; iy<ny_; iy++)
      for (int ix=0; ix<nx_; ix++)
         if (kpxycnt[nx_*iy+ix] > 0)
         {
            if (ix < xmin) xmin = ix;
            if (ix > xmax) xmax = ix;
            if (iy < ymin) ymin = iy;
            if (iy > ymax) ymax = iy;
            if (ix < kpymin_x[iy]) kpymin_x[iy] = ix;
            if (ix > kpymax_x[iy]) kpymax_x[iy] = ix;
         }

   //ewd DEBUG
   //if (myRank_ == 0)
   //   cout << "distributePlane DEBUG:  kp = " << kp << ", zmin = " << zmin << ", zmax = " << zmax << ", xmin = " << xmin << ", xmax = " << xmax << ", ymin = " << ymin << ", ymax = " << ymax << endl;
   
   double targetVol = (double)zvol/(npex_*npey_);
   bool xyvolConverged = false;
   int viter = 0;
   const int viterMax = 5000;
   while (!xyvolConverged && viter < viterMax)
   {
      vector<int> trialmin_y(npey_,-1);
      vector<int> trialmax_y(npey_,-1);
      vector<int> trialmin_x(npex_*npey_,-1);
      vector<int> trialmax_x(npex_*npey_,-1);
      vector<int> tmin_x(npex_);
      vector<int> tmax_x(npex_);
      vector<int> tjpvol(npey_,-1);
      int jpset = 0;
      trialmin_y[jpset] = ymin;
      trialmax_y[jpset] = ymin;
      for (int iy=ymin; iy<=ymax; iy++)
      {
         // try adding this iy strip to current jpset process strip
         vector<int> txcnt(nx_,0);
         for (int ty=trialmin_y[jpset]; ty<=iy; ty++)
            for (int ix=0; ix<nx_; ix++)
               txcnt[ix] += kpxycnt[nx_*ty+ix];

         xstripDistByVol(txcnt,tmin_x,tmax_x,false);
         
         // calculate maximum volume of (iy-trialmin_y) strips distributed over npex_ tasks
         int maxvol = -1;
         for (int ip=0; ip<npex_; ip++)
         {
            int tvol = (tmax_x[ip]-tmin_x[ip]+1)*(iy-trialmin_y[jpset]+1)*kpdz;
            if (tvol > maxvol) maxvol = tvol;
         }
         if (maxvol > targetVol)
         {
            jpset++;
            if (jpset > npey_-1)
            {
               jpset = npey_-1;
               tjpvol[jpset] = maxvol;
            }
            else
               trialmin_y[jpset] = iy;
         }
         else {
            tjpvol[jpset] = maxvol;
         }
         // save latest stats on this process grid strip (at jp,kp)
         trialmax_y[jpset] = iy;
         for (int ip=0; ip<npex_; ip++)
         {
            trialmin_x[jpset*npex_+ip] = tmin_x[ip];
            trialmax_x[jpset*npex_+ip] = tmax_x[ip];
         }            
      }
      
      // update npey_-1 term
      if (jpset == npey_-1)
         trialmax_y[npey_-1] = ymax;
         
      // all grid strips distributed, check volume distribution to see if target volume was exceeded
      int tvolmax = -1;
      int tvolmin = 999999999;
      for (int jp=0; jp<npey_; jp++)
      {
         if (tjpvol[jp] > tvolmax) tvolmax = tjpvol[jp];
         if (tjpvol[jp] < tvolmin) tvolmin = tjpvol[jp];
      }
      
      if ( (tvolmax < targetVol && tvolmax != tvolmin ) || tvolmin < 0 ) // keep going
      {
         // save current distribution, continue
         for (int jp=0; jp<npey_; jp++)
         {
            jpvol[jp] = tjpvol[jp];
            jpmin_y[jp] = trialmin_y[jp];
            jpmax_y[jp] = trialmax_y[jp];
            for (int ip=0; ip<npex_; ip++)
            {
               jpmin_x[jp*npex_+ip] = trialmin_x[jp*npex_+ip];
               jpmax_x[jp*npex_+ip] = trialmax_x[jp*npex_+ip];
            }
         }
      }
      else
      {
         if (viter == 0)
         {
            if (tvolmax != tvolmin)
            {
               if (myRank_ == 0)
                  cout << "Warning:  Plane " << kp << " could not reach target volume of " << targetVol << " on first iteration!  tvolmin = " << tvolmin << ", tvolmax = " << tvolmax << endl;

               // save results in peyind
               for (int jp=0; jp<npey_; jp++)
               {
                  for (int iy=trialmin_y[jp]; iy<=trialmax_y[jp]; iy++)
                     peyind_[ny_*kp+iy] = jp;
                  
                  for (int ip=0; ip<npex_; ip++)
                     for (int ix=trialmin_x[jp*npex_+ip]; ix<=trialmax_x[jp*npex_+ip]; ix++)
                        pexind_[nx_*npey_*kp + nx_*jp + ix] = ip;
               }

               //               return tvolmax;
               return -1;
            }
            //else
            //   if (myRank_ == 0)
            //      cout << "Plane " << kp << " already optimally converged, no further improvement possible:  tvolmin = " << tvolmin << ", tvolmax = " << tvolmax << endl;
         }
         
         // use previously stored distribution if available, otherwise save current one
         for (int jp=0; jp<npey_; jp++)
         {
            if (jpmin_y[jp] < 0) jpmin_y[jp] = trialmin_y[jp];
            if (jpmax_y[jp] < 0) jpmax_y[jp] = trialmax_y[jp];
            for (int ip=0; ip<npex_; ip++)
            {
               if (jpmin_x[jp*npex_+ip] < 0) jpmin_x[jp*npex_+ip] = trialmin_x[jp*npex_+ip];
               if (jpmax_x[jp*npex_+ip] < 0) jpmax_x[jp*npex_+ip] = trialmax_x[jp*npex_+ip];
            }
         }

         // compare maximum volume against max volume of last iteration
         int maxlast = -1;
         for (int jp=0; jp<npey_; jp++)
            if (jpvol[jp] > maxlast) maxlast = jpvol[jp];
         
         if (maxlast > 0 && maxlast < tvolmax)
         {
            tvolmax = maxlast;
            // restore last distribution
            for (int jp=0; jp<npey_; jp++)
            {
               jpmin_y[jp] = oldjpmin_y[jp];
               jpmax_y[jp] = oldjpmax_y[jp];
               for (int ip=0; ip<npex_; ip++)
               {
                  jpmin_x[jp*npex_+ip] = oldjpmin_x[jp*npex_+ip];
                  jpmax_x[jp*npex_+ip] = oldjpmax_x[jp*npex_+ip];
               }
            }
         }            
            
         if (myRank_ == 0)
            cout << "Strip convergence of plane " << kp << " reached at target volume per process = " << targetVol << ", tvolmax = " << tvolmax << endl;
         
         // save results in peyind
         for (int jp=0; jp<npey_; jp++)
         {
            for (int iy=jpmin_y[jp]; iy<=jpmax_y[jp]; iy++)
               peyind_[ny_*kp+iy] = jp;
            
            for (int ip=0; ip<npex_; ip++)
               for (int ix=jpmin_x[jp*npex_+ip]; ix<=jpmax_x[jp*npex_+ip]; ix++)
                  pexind_[nx_*npey_*kp + nx_*jp + ix] = ip;
         }               
         xyvolConverged = true;
      }      
      if (!xyvolConverged)
      {
         // store last distribution
         for (int jp=0; jp<npey_; jp++)
         {
            oldjpvol[jp] = jpvol[jp];
            oldjpmin_y[jp] = jpmin_y[jp];
            oldjpmax_y[jp] = jpmax_y[jp];
            for (int ip=0; ip<npex_; ip++)
            {
               oldjpmin_x[jp*npex_+ip] = jpmin_x[jp*npex_+ip];
               oldjpmax_x[jp*npex_+ip] = jpmax_x[jp*npex_+ip];
            }
         }

         viter++;
         targetVol *= 0.99;
      }
   }
   if (viter >= viterMax)
   {
      if (myRank_ == 0)
         cout << "ERROR:  distributePlane could not converge!" << endl;
      exit(1);
   }   
   
   return targetVol;
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::xstripDistByVol(vector<int>& xcnt, vector<int>& pexmin, vector<int>& pexmax, bool verbose)
{
   int xmin = 999999999;
   int xmax = -999999999;
   for (int ix=0; ix<nx_; ix++)
   {
      if (xcnt[ix] > 0 && xmin == 999999999)
         xmin = ix;
      if (xcnt[ix] > 0)
         xmax = ix;
   }
   for (int ip=0; ip<npex_; ip++)
   {
      pexmin[ip] = -1;
      pexmax[ip] = -1;
   }
   
   int ngap = 0;
   int gapcnt = 0;
   int isumgap = -1;
   vector<int> bvol;
   vector<int> gapii(1,xmin);
   vector<int> bstart(1,xmin);
   vector<int> bend(1,xmax);
   for (int ii=xmin; ii<=xmax; ii++) 
   {
      if (isumgap == -1 && xcnt[ii] > 0)
         bstart[ngap] = ii;
      if (xcnt[ii] > 0)
         bend[ngap] = ii;
      isumgap += xcnt[ii];
      if (isumgap > -1 && xcnt[ii] == 0)
         gapcnt++;
      else
         gapcnt = 0;
      
      if (gapcnt > gapthresh_) {
         gapii.push_back(ii-gapcnt);
         bvol.push_back(bend[ngap]-bstart[ngap]+1); // measuring "volume" w. x-dim only, y and z are constant
         bstart.push_back(ii);
         bend.push_back(xmax);
         ngap++;
         gapcnt = 0;
         isumgap = -1;
      }
   }
   bvol.push_back(bend[ngap]-bstart[ngap]+1);
   if (isumgap == -1 && ngap > 0)  // remove "gap" of empty space between last cells and box edge
      ngap--;
   gapii.push_back(xmax);

   // assign npex_ processors to each chunk, starting with the smallest chunks first
   vector<int> sortedVolIndex(ngap+1);
   for (int ig=0; ig<=ngap; ig++)
      sortedVolIndex[ig] = ig;
   sort(sortedVolIndex.begin(),sortedVolIndex.end(),GapSortInt(bvol));
   vector<int> bnpes(ngap+1);
   int peleft = npex_;
   
   int voltot = 0;
   for (int ig=0; ig<=ngap; ig++)
      voltot += bvol[ig];
   double avgvol = (double)voltot/(double)npex_;
   int avgvolint = (int)avgvol;
   for (int ig=0; ig<ngap; ig++)
   {
      int sind = sortedVolIndex[ig];
      int gvol = bvol[sind];
      int gpes = gvol/avgvol;
      if (abs(avgvol-(double)avgvolint) > 1.E-8 || gvol%avgvolint != 0) gpes++;
      bnpes[sind] = gpes;
      peleft -= gpes;
   }
   int sind = sortedVolIndex[ngap];
   bnpes[sind] = peleft;

   //ewd DEBUG
   if (verbose)
   {
      cout << "XSTRIPDIST:  initial array:  "  << endl;
      for (int ix=0; ix<nx_; ix++)
         cout << "  " << ix << "   " << xcnt[ix] << endl;

      cout << "XSTRIPDIST:  ngap = " << ngap << endl;
      for (int ig=0; ig<bstart.size(); ig++)
         cout << "   " << ig << "    " << bstart[ig] << "  " << bend[ig] << "   " << bvol[ig] << "   " << sortedVolIndex[ig] << "   " << bnpes[ig] << endl;
   }
   //ewd DEBUG
   
      
   // loop through strip again and assign bnpes to each chunk of cells
   int ipset = 0;
   for (int ig=0; ig<=ngap; ig++)
   {
      pexmin[ipset] = bstart[ig];
      pexmax[ipset] = bstart[ig];
      double gavg;

      //ewd DEBUG
      //cout << "DEBUG, ig = " << ig << ", bvol.size = " << bvol.size() << ", bnpes.size = " << bnpes.size() << ", bvol = " << bvol[ig] << ", bnpes = " << bnpes[ig] << endl;

      if (bnpes[ig] > 0)
      {
         if (bvol[ig]%bnpes[ig] == 0 && bvol[ig] > 0)
            gavg = (double)(bvol[ig]/bnpes[ig]);
         else
         {
            gavg = (double)bvol[ig]/(double)bnpes[ig] + 0.01;
            if (gavg < 1.01) gavg = 1.01;
         }
      
         int gip = 0;
         int volsum = 0;
         for (int ii=bstart[ig]; ii<=bend[ig]; ii++) 
         {
            volsum++;
            if (volsum > gavg*(gip+1)) {
               //if (volsum >= gavg*(gip+1)) {

               //ewd DEBUG
               //if (verbose && ngap == 3)
               //   cout << "XSTRIPDIST: gap " << ig << ", gavg = " << gavg << ", volsum = " << volsum << ", gavg*(gip+1) = " << gavg*(gip+1) << ", gip = " << gip << ", ipset = " << ipset << endl;
               
               ipset++;
               if (ipset > npex_-1)
                  ipset = npex_-1;
               else
                  pexmin[ipset] = ii;
               gip++;
            }
            if (xcnt[ii] > 0) 
               pexmax[ipset] = ii;
         }
         ipset++;
         if (ipset > npex_-1) ipset = npex_-1;
      }
   }
   if (pexmin[npex_-1] > 0)
      pexmax[npex_-1] = xmax;


   //ewd DEBUG      
   if (verbose)
   {
      cout << "calculated distribution:" << endl;
      for (int ip=0; ip<npex_; ip++) {
         cout << "  " << ip << "    " << pexmin[ip] << "   " << pexmax[ip] << endl;
      }
   }
   //ewd DEBUG



}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::redistributeCells(vector<AnatomyCell>& cells)
{
   bool testingOnly = (npegrid_ != nTasks_);
   if (!testingOnly)
   {
      sort(cells.begin(),cells.end(),AnatomyCell::destLessThan);
      Long64 nLocal = cells.size();
      vector<unsigned> dest(nLocal);
      for (unsigned ii=0; ii<cells.size(); ++ii)
         dest[ii] = cells[ii].dest_;

      // compute largest possible local size
      Long64 maxLocbuf = nLocal;
      Long64 nMax;
      MPI_Allreduce(&maxLocbuf, &nMax, 1, MPI_LONG_LONG, MPI_MAX, comm_);
      nMax *= 4;  //fudge factor
      //int nMax = 10*nx_*ny_*nz_/nTasks_;      // 10 is a fudge factor for safety
      cells.resize(nMax);
      unsigned nLocu = nLocal;
      assignArray((unsigned char*)&(cells[0]), &nLocu, cells.capacity(),
                  sizeof(AnatomyCell), &(dest[0]), 0, comm_);
      nLocal = nLocu;
      if (nLocal > nMax)
         cout << "GDLB::redistributeCells assertion failure on myRank = " << myRank_ << ", nLocal = " << nLocal << ", nMax = " << nMax << endl;
      assert(nLocal <= nMax);
      cells.resize(nLocal);
   }
   else  // testing case:  nTasks < target process grid size
   {
      Long64 nLocal = cells.size();
      vector<unsigned> dest(nLocal);

      // compress process grid onto nTasks
      if (tnx_ < 0 || tny_ < 0 || tnz_ < 0)
      {
         computeReducedProcGrid(nTasks_);
      }
         if (myRank_ == 0)
            cout << "GDLoadBalancer::redistributeCells:  reduced grid = " <<
                tnx_ << " x " << tny_ << " x " << tnz_ << " used for testing." << endl;

      // for testing, certain assumptions have to be fulfilled
      assert(npex_%2==0 && npey_%2==0 && npez_%2==0);
      assert(npex_%tnx_ == 0 && npey_%tny_ == 0 && npez_%tnz_ == 0);
        
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         int realdestpe = (int)cells[ii].dest_;
         GridPoint pept(realdestpe,npex_,npey_,npez_);
         int tix = pept.x*tnx_/npex_;
         int tiy = pept.y*tny_/npey_;
         int tiz = pept.z*tnz_/npez_;
         dest[ii] = tix + tiy*tnx_ + tiz*tnx_*tny_;
         cells[ii].sortind_ = dest[ii];
      }
      // need to sort dest, cells (by sort id, not dest)
      sort(dest.begin(),dest.end());
      sort(cells.begin(),cells.end(),AnatomyCell::indLessThan);

      Long64 nMax = 2*nx_*ny_*nz_/nTasks_;  // 2 is a fudge factor for safety
      cells.resize(nMax);
      unsigned nLocu = nLocal;
      assignArray((unsigned char*)&(cells[0]), &nLocu, cells.capacity(),
                  sizeof(AnatomyCell), &(dest[0]), 0, comm_);
      nLocal = nLocu;
      assert(nLocal <= nMax);
      cells.resize(nLocal);
   }
}

////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::setReducedProcGrid(int rnx, int rny, int rnz)
{
   tnx_ = rnx;
   tny_ = rny;
   tnz_ = rnz;
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::computeReducedProcGrid(int nTasks)
{
   int npes = nTasks;
   int npow = 0;
   while (npes > 1) {
      npes /= 2;
      npow++;
   }
   assert(ipow(2,npow) == nTasks);      // assume nTasks is a power of two
   tnx_ = tny_ = tnz_ = ipow(2,npow/3);
   int nrem = npow%3;
   if (nrem == 1)
      tnz_ *= 2;
   else if (nrem == 2)
   {
      tnz_ *= 2;
      tny_ *= 2;
   }
}

////////////////////////////////////////////////////////////////////////////////
int GDLoadBalancer::ipow(int a, int n)
{
   int b;
   if (n == 0) { b = 1; }
   else if (n < 0) { b = 0; }
   else
   {
      b = a;
      for (int ii=1; ii<n; ++ii)
         b *= a;
   }
   return b;
}    
////////////////////////////////////////////////////////////////////////////////
double GDLoadBalancer::costFunction(int nTissue, int area, int height, double a)
{
   int dz4 = (height+2);
   if (dz4 % 4  != 0) dz4 += (4-dz4%4);
   int bbVol = area*dz4;
   double cost = nTissue+a*bbVol;
   return cost;
}
