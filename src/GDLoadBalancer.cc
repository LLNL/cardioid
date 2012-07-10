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
#include <mpi.h>
#include "mpiUtils.h"
#include "GridPoint.hh"
#include "AnatomyCell.hh"
#include "ProcBox.hh"
using namespace std;



////////////////////////////////////////////////////////////////////////////////
GDLoadBalancer::GDLoadBalancer(MPI_Comm comm, int npex, int npey, int npez):
    comm_(comm), npex_(npex), npey_(npey), npez_(npez), gapthresh_(5)
{
   MPI_Comm_size(comm_, &nTasks_);
   MPI_Comm_rank(comm_, &myRank_);

   // only exchange load across faces
   nnbr_ = 6;

   // equivalent local neighbor index of neighboring point
   thatn_.resize(nnbr_);  
   thatn_[0] = 1; thatn_[1] = 0;
   thatn_[2] = 3; thatn_[3] = 2;
   thatn_[4] = 5; thatn_[5] = 4;

   // initialize arrays
   npegrid_ = npex_*npey_*npez_;
   togive_.resize(npegrid_,vector<int>(nnbr_,0));
   penbr_.resize(npegrid_,vector<int>(nnbr_,-1));

   peboxinfo_.resize(npegrid_);
   for (int ip=0; ip<npegrid_; ip++)
      peboxinfo_[ip] = new ProcBox(ip);

   // store process numbers of all neighbors
   for (int ip=0; ip<npegrid_; ip++)
   {
      // x,y,z coords of process ip
      GridPoint ipt(ip,npex_,npey_,npez_);

      // face to face exchanges
      int nbr;
      if (ipt.x > 0) {
         nbr = (ipt.x-1) + ipt.y*npex_ + ipt.z*npex_*npey_;
         penbr_[ip][0] = nbr;
      }
      if (ipt.x < npex_-1) {
         nbr = (ipt.x+1) + ipt.y*npex_ + ipt.z*npex_*npey_;
         penbr_[ip][1] = nbr;
      }
      if (ipt.y > 0) {
         nbr = ipt.x + (ipt.y-1)*npex_ + ipt.z*npex_*npey_;
         penbr_[ip][2] = nbr;
      }
      if (ipt.y < npey_-1) {
         nbr = ipt.x + (ipt.y+1)*npex_ + ipt.z*npex_*npey_;
         penbr_[ip][3] = nbr;
      }
      if (ipt.z > 0) {
         nbr = ipt.x + ipt.y*npex_ + (ipt.z-1)*npex_*npey_;
         penbr_[ip][4] = nbr;
      }
      if (ipt.z < npez_-1) {
         nbr = ipt.x + ipt.y*npex_ + (ipt.z+1)*npex_*npey_;
         penbr_[ip][5] = nbr;
      }
   }

   // reduced process grid dimensions
   tnx_ = tny_ = tnz_ = -1;
   
}
////////////////////////////////////////////////////////////////////////////////
GDLoadBalancer::~GDLoadBalancer()
{
   for (int ip=0; ip<npegrid_; ip++)
      delete peboxinfo_[ip];
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
      int volFinal = distributePlane(cells,kpmin_z[kp],kpmax_z[kp],kpzvol[kp],kp);
      for (int kk=kpmin_z[kp]; kk<=kpmax_z[kp]; kk++)
         distVol[kk] = volFinal;
   }


   // calculate volume of each plane
   {
      vector<int> cell_xmin(npegrid_,999999999);
      vector<int> cell_ymin(npegrid_,999999999);
      vector<int> cell_zmin(npegrid_,999999999);
      vector<int> cell_xmax(npegrid_,-999999999);
      vector<int> cell_ymax(npegrid_,-999999999);
      vector<int> cell_zmax(npegrid_,-999999999);
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         int gid = cells[ii].gid_;
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
         if (trialvol_z[kp] > 0)
            int volFinal = distributePlane(cells,trialmin_z[kp],trialmax_z[kp],trialvol_z[kp],kp);
         else
            if (myRank_ == 0)
               cout << "Redistributing:  plane " << kp << ", has zero target volume!" << endl;
      }

   }

   
   // we now have process coordinates for all local grid points, change dest_
   // index in cells array
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      int gid = cells[ii].gid_;
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

   if (myRank_ == 0)
      cout << "Load histogram after initial volume-weighted distribution:" << endl;
   nlocHistogram(cells);
}
////////////////////////////////////////////////////////////////////////////////
int GDLoadBalancer::distributePlane(vector<AnatomyCell>& cells, int zmin, int zmax, int zvol, int kp)
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

         xstripDist(txcnt,tmin_x,tmax_x,false);
         
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
                  cout << "ERROR:  Plane " << kp << " could not reach target volume of " << targetVol << " on first iteration!  tvolmin = " << tvolmin << ", tvolmax = " << tvolmax << endl;
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
         targetVol *= 0.995;
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
void GDLoadBalancer::xstripDist(vector<int>& xcnt, vector<int>& pexmin, vector<int>& pexmax, bool verbose)
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
   sort(sortedVolIndex.begin(),sortedVolIndex.end(),GapSort(bvol));
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
      unsigned nLocal = cells.size();
      vector<unsigned> dest(nLocal);
      for (unsigned ii=0; ii<cells.size(); ++ii)
         dest[ii] = cells[ii].dest_;

      // compute largest possible local size
      int maxLocbuf = nLocal;
      int nMax;
      MPI_Allreduce(&maxLocbuf, &nMax, 1, MPI_INT, MPI_MAX, comm_);
      nMax *= 4;  //fudge factor
      //int nMax = 10*nx_*ny_*nz_/nTasks_;      // 10 is a fudge factor for safety
      cells.resize(nMax);
      assignArray((unsigned char*)&(cells[0]), &nLocal, cells.capacity(),
                  sizeof(AnatomyCell), &(dest[0]), 0, comm_);
      if (nLocal > nMax)
         cout << "GDLB::redistributeCells assertion failure on myRank = " << myRank_ << ", nLocal = " << nLocal << ", nMax = " << nMax << endl;
      assert(nLocal <= nMax);
      cells.resize(nLocal);
   }
   else  // testing case:  nTasks < target process grid size
   {
      unsigned nLocal = cells.size();
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
         int realdestpe = cells[ii].dest_;
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

      int nMax = 2*nx_*ny_*nz_/nTasks_;  // 2 is a fudge factor for safety
      cells.resize(nMax);
      assignArray((unsigned char*)&(cells[0]), &nLocal, cells.capacity(),
                  sizeof(AnatomyCell), &(dest[0]), 0, comm_);
      assert(nLocal <= nMax);
      cells.resize(nLocal);
   }
}

////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::nlocHistogram(vector<AnatomyCell>& cells)
{
   // compute load histogram from data in cells
   histnloc_.resize(npegrid_,0);
   vector<int> mydata(npegrid_,0);
   for (unsigned ii=0; ii<cells.size(); ++ii)
      mydata[cells[ii].dest_]++;
   MPI_Allreduce(&mydata[0], &histnloc_[0], npegrid_, MPI_INT, MPI_SUM, comm_);
   if (myRank_ == 0)
      computeNlocHistogram();
    
}

////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::computeNlocHistogram()
{
   // compute histogram of current data distribution
   const int nhistmax = 100; // number of bins
   vector<int> phist(nhistmax,0);

   int maxnum = 0;
   int minnum = nx_*ny_*nz_;
   int maxpe = -1;
   for (int p=0; p<npegrid_; p++)
   {
      if (histnloc_[p] > maxnum) {
         maxnum = histnloc_[p];
         maxpe = p;
      }
      if (histnloc_[p] < minnum)
         minnum = histnloc_[p];
   }
   int nhist = maxnum - minnum + 1;
   if (nhist > nhistmax) nhist = nhistmax;
   int delta = (maxnum-minnum + 1)/nhist;
   if ((maxnum-minnum+1)%nhist !=0) delta++;
   for (int p=0; p<npegrid_; p++)
   {
      int bin = (histnloc_[p]-minnum)/delta;
      phist[bin]++;
      //ewd DEBUG: print process number of top bin
      //if (bin == nhist-1 && phist[bin] == 1)
      //  cout << "nlocHistogram: top bin pe " << p << ", nloc = " << histnloc_[p] << endl;
   }
   cout << "load balance histogram (ncells):  " << endl;
   for (int i=0; i<nhist; i++)
      cout << "  " << minnum+delta*i << " - " << minnum+delta*(i+1) << ":    " << phist[i] << endl;

   nloctot_ = 0;
   for (unsigned ii=0; ii<npegrid_; ++ii)
      nloctot_ += histnloc_[ii];
   double nlocavg_ = (double)nloctot_/(double)npegrid_; 
    
   cout << "total # of non-zero grid points = " << nloctot_ << ", avg. # per task = " << nlocavg_ << ", max pe = " << maxpe << " (" << maxnum << ")" << endl << endl;

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
      for (unsigned ii=1; ii<n; ++ii)
         b *= a;
   }
   return b;
}    
