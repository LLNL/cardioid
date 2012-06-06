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

      /*
      // face to face exchanges, allow for periodic connectivity
      int nbr;
      if (ipt.x > 0) 
         nbr = (ipt.x-1) + ipt.y*npex_ + ipt.z*npex_*npey_;
      else
         nbr = npex_-1 + ipt.y*npex_ + ipt.z*npex_*npey_;         
      penbr_[ip][0] = nbr;

      if (ipt.x < npex_-1) 
         nbr = (ipt.x+1) + ipt.y*npex_ + ipt.z*npex_*npey_;
      else
         nbr = 0 + ipt.y*npex_ + ipt.z*npex_*npey_;         
      penbr_[ip][1] = nbr;

      if (ipt.y > 0) 
         nbr = ipt.x + (ipt.y-1)*npex_ + ipt.z*npex_*npey_;
      else
         nbr = ipt.x + (npey_-1)*npex_ + ipt.z*npex_*npey_;
      penbr_[ip][2] = nbr;

      if (ipt.y < npey_-1) 
         nbr = ipt.x + (ipt.y+1)*npex_ + ipt.z*npex_*npey_;
      else
         nbr = ipt.x + ipt.z*npex_*npey_;         
      penbr_[ip][3] = nbr;

      if (ipt.z > 0) 
         nbr = ipt.x + ipt.y*npex_ + (ipt.z-1)*npex_*npey_;
      else
         nbr = ipt.x + ipt.y*npex_ + (npez_-1)*npex_*npey_;
      penbr_[ip][4] = nbr;

      if (ipt.z < npez_-1) 
         nbr = ipt.x + ipt.y*npex_ + (ipt.z+1)*npex_*npey_;
      else
         nbr = ipt.x + ipt.y*npex_;
      penbr_[ip][5] = nbr;

      */
   }

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
   vector<int> peyind(npez_*ny_,-1);
   vector<int> pexind(npez_*npey_*nx_,-1);

   // calculate min/max boundaries of each plane
   vector<int> zmin_xloc(nz_,99999999);
   vector<int> zmin_yloc(nz_,99999999);
   vector<int> zmax_xloc(nz_,-99999999);
   vector<int> zmax_yloc(nz_,-99999999);
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

   vector<int> kpmin_x(npez_,99999999);
   vector<int> kpmin_y(npez_,99999999);
   vector<int> kpmin_z(npez_,-1);
   vector<int> kpmax_x(npez_,-99999999);
   vector<int> kpmax_y(npez_,-99999999);
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
      
      if (kk > 0 && kk%kpavg == 0)
      {
         kpzvol[kpset] = (kpmax_x[kpset]-kpmin_x[kpset]+1)*(kpmax_y[kpset]-kpmin_y[kpset]+1)*(kk-kpmin_z[kpset]+1);
         kpmax_z[kpset] = kk;
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
         int xmin = 9999999;  int ymin = 9999999;
         int xmax = -9999999;  int ymax = -9999999;
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
         int tvolmin = 9999999;
         for (int kp=0; kp<npez_; kp++)
         {
            if (trialvol_z[kp] > tvolmax) tvolmax = trialvol_z[kp];
            if (trialvol_z[kp] < tvolmin) tvolmin = trialvol_z[kp];
         }
      
         if (myRank_ == 0)
            cout << "Plane distribution: iteration " << viter++ << ", targetVol = " << targetVol << ", max volume = " << tvolmax << ", tvolmin = " << tvolmin << endl;

         if (tvolmax < targetVol || tvolmin < 0)
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
   for (int kp=0; kp<npez_; kp++)
   {
      vector<int> jpmin_y(npey_,-1);
      vector<int> jpmax_y(npey_,-1);      
      vector<int> jpmin_x(npex_*npey_,-1);
      vector<int> jpmax_x(npex_*npey_,-1);
      
      int kpdz = kpmax_z[kp] - kpmin_z[kp] + 1;
      vector<int> kpxyloc(nx_*ny_,0);
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
         if (gpt.z >= kpmin_z[kp] && gpt.z <= kpmax_z[kp])
            kpxyloc[nx_*gpt.y+gpt.x]++;
      }
      int size = nx_*ny_;
      vector<int> kpxycnt(nx_*ny_);
      MPI_Allreduce(&kpxyloc[0], &kpxycnt[0], size, MPI_INT, MPI_SUM, comm_);

      // compute min/max values for all planes, strips
      vector<int> kpymin_x(ny_,9999999);
      vector<int> kpymax_x(ny_,-9999999);
      int xmin = 9999999;  int ymin = 9999999;
      int xmax = -9999999;  int ymax = -9999999;
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

      double targetVol = (double)kpzvol[kp]/(npex_*npey_);
      bool xyvolConverged = false;
      int viter = 0;
      while (!xyvolConverged)
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

            //if (false && myRank_ == 0 && kp == 20)
            if (false && myRank_ == 0 && kp == 0)
               xstripDist(txcnt,tmin_x,tmax_x,true);
            else
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
                  jpset = npey_-1;
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
         {
            int maxvol = -1;
            for (int ip=0; ip<npex_; ip++)
            {
               int tvol = (tmax_x[ip]-tmin_x[ip]+1)*(ymax-trialmin_y[npey_-1]+1)*kpdz;
               if (tvol > maxvol) maxvol = tvol;
            }
            trialmax_y[npey_-1] = ymax;
            tjpvol[npey_-1] = maxvol;
         }
         
         // all grid strips distributed, check volume distribution to see if target volume was exceeded
         int tvolmax = -1;
         int tvolmin = 9999999;
         for (int jp=0; jp<npey_; jp++)
         {
            if (tjpvol[jp] > tvolmax) tvolmax = tjpvol[jp];
            if (tjpvol[jp] < tvolmin) tvolmin = tjpvol[jp];
         }

         if (tvolmax < targetVol || tvolmin < 0)  // keep going
         {

            //ewd DEBUG
            if (false && myRank_ == 0)
            {
               cout << "Converging plane " << kp << ", iter " << viter << ", tvolmax = " << tvolmax << ", tvolmin = " << tvolmin << ", targetVol = " << targetVol << endl;
               for (int jp=0; jp<npey_; jp++)
                  cout << "   " << jp << "    " << tjpvol[jp] << endl;
            }
            //ewd DEBUG
 
           // save current distribution, continue
            for (int jp=0; jp<npey_; jp++)
            {
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

            //ewd DEBUG
            if (myRank_ == 0 && viter == 0)
               cout << "Plane " << kp << " could not reach target volume of " << targetVol << " on first iteration!  tvolmin = " << tvolmin << ", tvolmax = " << tvolmax << endl;
            //ewd DEBUG
            
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

            
            if (myRank_ == 0)
               cout << "Strip convergence of plane " << kp << " reached at target volume = " << targetVol << endl;

            // save results in peyind
            for (int jp=0; jp<npey_; jp++)
            {
               for (int iy=jpmin_y[jp]; iy<=jpmax_y[jp]; iy++)
                  peyind[ny_*kp+iy] = jp;

               for (int ip=0; ip<npex_; ip++)
                  for (int ix=jpmin_x[jp*npex_+ip]; ix<=jpmax_x[jp*npex_+ip]; ix++)
                     pexind[nx_*npey_*kp + nx_*jp + ix] = ip;
            }               
            xyvolConverged = true;
         }
         viter++;
         targetVol *= 0.99;
      }
   }

   // we now have process coordinates for all local grid points, change dest_
   // index in cells array
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      int gid = cells[ii].gid_;
      GridPoint gpt(gid,nx_,ny_,nz_);
      int kp = pezind[gpt.z];
      int jp = peyind[ny_*kp + gpt.y];
      int ip = pexind[nx_*npey_*kp + jp*nx_ + gpt.x];
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
   //volHistogram(cells);   
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::xstripDist(vector<int>& xcnt, vector<int>& pexmin, vector<int>& pexmax, bool verbose)
{
   int xmin = 9999999;
   int xmax = -9999999;
   for (int ix=0; ix<nx_; ix++)
   {
      if (xcnt[ix] > 0 && xmin == 9999999)
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
   for (int ii=xmin; ii<xmax; ii++) 
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
      double gavg = (double)bvol[ig]/(double)bnpes[ig] + 0.01;
      if (gavg < 1.01) gavg = 1.01;
      
      int gip = 0;
      int volsum = 0;
      for (int ii=bstart[ig]; ii<=bend[ig]; ii++) 
      {
         volsum++;
         if (volsum > gavg*(gip+1)) {

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
void GDLoadBalancer::initialDistribution(vector<AnatomyCell>& cells, int nx, int ny, int nz)
{
   // cells should contain only the cells with non-zero type indices, distributed
   // in x-y planes of grid points as follows:
   //    if (nTasks >= nz)     destRank = gid.z   (one plane per task for myRank < nz)
   //    else if (nTasks < nz) destRank = gid.z*nTasks/nz  (multiple planes per task)
   
   int nlocal = cells.size();
   nx_ = nx;
   ny_ = ny;
   nz_ = nz;

   // count total number of non-zero grid points
   MPI_Allreduce(&nlocal, &ntissue_, 1, MPI_INT, MPI_SUM, comm_);

   // determine which xy-planes this process owns data for
   int klocmax = -1;
   int klocmin = nz_+1;
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      if (gpt.z < klocmin) klocmin = gpt.z;
      if (gpt.z > klocmax) klocmax = gpt.z;
   }
   int nkloc = klocmax - klocmin + 1;
   if (nkloc < 0) nkloc = 0;

   //ewd DEBUG
   cout << "GDLBINIT1, myRank = " << myRank_ << ", nkloc = " << nkloc << ", klocmin = " << klocmin << ", klocmax = " << klocmax << endl;

   
   // calculate total number of non-zero points at each x-y plane, x row
   double xypeavg = (double)ntissue_/(double)npez_;   // avg. # of pts in each x-y plane of process grid
   vector<int> nxyplane_loc(nz_,0);
   vector<int> xrowcnt_loc(nz_*ny_,0);
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      nxyplane_loc[gpt.z]++;
      xrowcnt_loc[gpt.z*ny_+gpt.y]++;
   }

   // sum local nxyplane data (# of non-zero points/plane) across all tasks
   vector<int> nxyplane(nz);
   MPI_Allreduce(&nxyplane_loc[0], &nxyplane[0], nz, MPI_INT, MPI_SUM, comm_);

   // sum local xrowcnt data (# of non-zero points/grid plane row) across all tasks
   int nyz = ny_*nz_;
   vector<int> xrowcnt(nyz);  
   MPI_Allreduce(&xrowcnt_loc[0], &xrowcnt[0], nyz, MPI_INT, MPI_SUM, comm_);
  
   // distribute grid points on process grid with uniform distribution in z-direction
   vector<int> pezind(nz_);  // z-coordinate in process grid of data plane k
   vector<int> peyind(ny_);
   vector<int> pexind(npey_*nx_);
  
   vector<int> nxype(npez_,-1);
   vector<int> kpkmax(npez_,-1);
   vector<int> kpkmin(npez_);
   int kset = 0;
   int xysum = 0;
   int ksum = 0;
   for (int kk=0; kk<nz; kk++)
   {
      pezind[kk] = kset;
      xysum += nxyplane[kk];
      ksum += nxyplane[kk];
      if (xysum >= xypeavg*(kset+1)) {
         nxype[kset] = ksum;
         ksum = 0;
         kpkmax[kset] = kk+1;
         kpkmin[kset] = (kset > 0 ? kpkmax[kset-1] : 0);
         kset++;
      }
   }
   if (npez_ > 1)
      kpkmin[npez_-1] = kpkmax[npez_-2];
   else
      kpkmin[0] = 0;

   if (nxype[npez_-1] < 0)
      nxype[npez_-1] = ksum;
   if (kpkmax[npez_-1] < 0)
      kpkmax[npez_-1] = nz_;
   
   // loop over process grid planes owned by this process, compute distribution in x and y
   for (int kp = 0; kp < npez_; kp++)
   {
      vector<int> xdist_loc(npey_*nx_,0);
      vector<int> kpjsum(npey_,-1);
      vector<int> kpjmin(npey_,-1);
      vector<int> kpjmax(npey_,-1);
      if (kp >= pezind[klocmin] && kp <= pezind[klocmax])  // this task owns data for this slab
      {
         if (nxype[kp] > 0 )
         {
            // sum up all points in full process slab at each j
            vector<int> kprowsum(ny_,0);
            for (int jj=0; jj<ny_; jj++)
               for (int kk=kpkmin[kp]; kk<kpkmax[kp]; kk++)
                  kprowsum[jj] += xrowcnt[kk*ny_+jj];

            int jset = 0;
            int jsum = 0;
            int rowsum = 0;
            double kprowavg = (double)nxype[kp]/(double)npey_;
            for (int jj=0; jj<ny_; jj++)
            {
               peyind[jj] = jset;
               rowsum += kprowsum[jj];
               jsum += kprowsum[jj];
               if (rowsum >= kprowavg*(jset+1)) {
                  kpjsum[jset] = jsum;
                  jsum = 0;
                  kpjmax[jset] = jj+1;
                  kpjmin[jset] = (jset > 0 ? kpjmax[jset-1] : 0);
                  jset++;
               }
            }
            if (npey_ > 1)
               kpjmin[npey_-1] = kpjmax[npey_-2];
            else
               kpjmin[0] = 0;
            kpjmax[npey_-1] = ny_;
            if (kpjsum[npey_-1] < 0)
               kpjsum[npey_-1] = jsum;

            for (unsigned ii=0; ii<cells.size(); ++ii)
            {
               GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
               int kptmp = pezind[gpt.z];
               int jptmp = peyind[gpt.y];
               if (kptmp == kp)
                  xdist_loc[jptmp*nx_ + gpt.x]++;
            }
         }
         else {
            cout << "GDLoadBalancer::InitialDistribution WARNING:  nxype[kp] <= 0:  kp = " << kp << ", nxype = " << nxype[kp] << ", myRank = " << myRank_ << endl;
         }
      }
    
      bool sortXByVol = true;
      int xdsize = npey_*nx_;
      vector<int> xdistsum(xdsize);  
      MPI_Allreduce(&xdist_loc[0], &xdistsum[0], xdsize, MPI_INT, MPI_SUM, comm_);
      if (kp >= pezind[klocmin] && kp <= pezind[klocmax])  // this task owns data for this slab
      {
         for (int jp = 0; jp < npey_; jp++)
         {
            // for non-uniform geometries, e.g. hearts, previous versions of the load balancer assigned
            // cells which were highly separated in space, creating very large bounding boxes.  
            //
            // to avoid this, pre-calculate where all the gaps are in a given stripe and
            // assign pe boundaries to them
            //   a) pre-count number of gaps
            //   b) adjust average cells/pe to account for smaller domains caused by forced domain creation
            //      at gap boundaries
            //   c) if gap boundary falls at natural domain boundary, last processor won't get
            //      any work --> adjust xavg when this occurs
            
            int ngap = 0;
            int gapcnt = 0;
            int isumgap = -1;
            vector<int> bsize(1,0);
            vector<int> bvol;
            vector<int> gapii(1,0);
            int bstart = 0;
            int bend = 0;
            for (int ii=0; ii<nx_; ii++) 
            {
               if (isumgap == -1 && xdistsum[jp*nx_+ii] > 0)
                  bstart = ii;
               if (xdistsum[jp*nx_+ii] > 0)
                  bend = ii;
               isumgap += xdistsum[jp*nx_+ii];
               bsize[ngap] += xdistsum[jp*nx_+ii];
               if (isumgap > -1 && xdistsum[jp*nx_+ii] == 0)
                  gapcnt++;
               else
                  gapcnt = 0;

               if (gapcnt > gapthresh_) {
                  bsize.push_back(0);
                  gapii.push_back(ii);
                  bvol.push_back(bend-bstart+1); // measuring "volume" w. x-dim only, y and z are constant
                  bstart = ii;
                  ngap++;
                  gapcnt = 0;
                  isumgap = -1;
               }
            }
            bvol.push_back(bend-bstart+1);
            if (isumgap == -1)  // remove "gap" of empty space between last cells and box edge
               ngap--;
            gapii.push_back(nx_);

            // assign npex_ processors to each chunk, starting with the smallest chunks first
            double xavg = (double)kpjsum[jp]/(double)(npex_);
            int xavgint = (int)xavg;
            vector<int> sortedSizeIndex(ngap+1);
            for (int ig=0; ig<=ngap; ig++)
               sortedSizeIndex[ig] = ig;
            sort(sortedSizeIndex.begin(),sortedSizeIndex.end(),GapSort(bsize));
            vector<int> sortedVolIndex(ngap+1);
            for (int ig=0; ig<=ngap; ig++)
               sortedVolIndex[ig] = ig;
            sort(sortedVolIndex.begin(),sortedVolIndex.end(),GapSort(bvol));
            vector<int> bnpes(ngap+1);
            int peleft = npex_;

            if (sortXByVol)
            {
               int voltot = 0;
               for (int ig=0; ig<=ngap; ig++)
                  voltot += bvol[ig];

               //ewd DEBUG
               int sizetot = 0;
               for (int ig=0; ig<=ngap; ig++)
                  sizetot += bsize[ig];
               //ewd DEBUG

               
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
               bnpes[sortedVolIndex[ngap]] = peleft;
               
               //ewd DEBUG

               //manually compute dimensions of this process strip?
               //for (int jj=0; jj<ny_; jj++)


               for (int ig=0; ig<=ngap; ig++)
               {
                  int sind = sortedVolIndex[ig];
                  int gvol = bvol[sind];
                  int gpes = bnpes[sind];
                  int volmult = (klocmax-klocmin+1)*(kpjmax[jp]-kpjmin[jp]+1);
                  int stripvol = volmult*nx_;
                  int tvol = volmult*voltot;
                  cout << "GINIT2, pe = " << myRank_ << ", kp = " << kp << ", jp = " << jp << ", strip = 0:" << nx_ << " " << kpjmin[jp] << ":" << kpjmax[jp] << " " << klocmin << ":" << klocmax << ", strip vol = " << stripvol << ", sizetot = " << sizetot << ", bsize[ig] = " << bsize[ig] << ", voltot = " << tvol << ", gap " << ig << ", gvol = " << gvol*volmult << ", gpes = " << gpes << endl;

                  cout << "GCHECK, pe = " << myRank_ << ", kp = " << kp << ", jp = " << jp << ", dims for this pe strip: " << kpjmin[jp] << ":" << kpjmax[jp] << " " << klocmin << ":" << klocmax << ", jp_min = " << peyind[kpjmin[jp]] << ", jp_max = " << peyind[kpjmax[jp]] << ", kp_min = " << pezind[klocmin] << ", kp_max = " << pezind[klocmax] << endl;
               }
               //ewd DEBUG



               
               //ewd DEBUG
               /*
               ostringstream oss;
               oss << "XVOL kp = " << kp << ", jp = " << jp << ", ngaps = " << ngap << ": ";
               for (int ig=0; ig<=ngap; ig++)
                  oss << bvol[ig] << " ";
               oss << ", pes: ";
               for (int ig=0; ig<=ngap; ig++)
                  oss << bnpes[ig] << " ";
               oss << ", avgvol = " << avgvol;
               oss << endl;
               string outstr = oss.str();
               cout << outstr;
               */
               //ewd DEBUG
               
            }
            else
            {
               for (int ig=0; ig<ngap; ig++)
               {
                  int sind = sortedSizeIndex[ig];
                  int gsize = bsize[sind];
                  int gpes = gsize/xavg;
                  if (abs(xavg-(double)xavgint) > 1.E-8 || gsize%xavgint != 0) gpes++;
                  bnpes[sind] = gpes;
                  peleft -= gpes;
               }
               bnpes[sortedSizeIndex[ngap]] = peleft;
            }
            
            //ewd DEBUG
            //ostringstream oss;
            //oss << "XGAP  myRank = " << myRank_ << ", kp = " << kp << ", jp = " << jp << ", ngaps = " << ngap << ": ";
            //for (int ig=0; ig<ngap; ig++)
            //   oss << bsize[ig] << " ";
            //oss << ", sorted: ";
            //for (int ig=0; ig<ngap; ig++)
            //   oss << sortedSizeIndex[ig] << " ";
            //oss << ", pes: ";
            //for (int ig=0; ig<ngap; ig++)
            //   oss << bnpes[ig] << " ";
            //oss << ", xavg = " << xavg;
            //oss << endl;
            //string outstr = oss.str();
            //cout << outstr;   
            //ewd DEBUG

            // now loop through each strip again and assign bnpes to each chunk of cells
            if (sortXByVol)
            {
               int ipset = 0;
               for (int ig=0; ig<=ngap; ig++)
               {
                  double gavg = (double)bvol[ig]/(double)bnpes[ig];
                  int gip = 0;
                  int volsum = 0;
                  for (int ii=gapii[ig]; ii<gapii[ig+1]; ii++) 
                  {
                     if (xdistsum[jp*nx_+ii] > 0)
                        volsum++;
                     
                     if (volsum >= gavg*(gip+1)) {
                        ipset++;
                        if (ipset > npex_-1) ipset = npex_-1;
                        gip++;
                     }
                     pexind[jp*nx_+ii] = ipset;                     
                  }
                  ipset++;
                  if (ipset > npex_-1) ipset = npex_-1;
               }
            }
            else
            {
               int ipset = 0;
               for (int ig=0; ig<=ngap; ig++)
               {
                  int gsize = bsize[ig];
                  double gavg = (double)gsize/(double)bnpes[ig] + 0.5;
                  int gip = 0;
                  int blksum = 0;
                  for (int ii=gapii[ig]; ii<gapii[ig+1]; ii++) 
                  {
                     blksum += xdistsum[jp*nx_+ii];
                     if (blksum >= gavg*(gip+1)) {
                        ipset++;
                        if (ipset > npex_-1) ipset = npex_-1;
                        gip++;
                     }
                     pexind[jp*nx_+ii] = ipset;
                  }
               }
            }

            /*** EWD:  this was how we used to handle gaps, by just adjusting average at gap boundaries
            // it worked well overall, but couldn't handle difficult multi-gapped regions
            
            // use larger average values to compensate for smaller domains created by forced
            // processor transitions at gaps
            //ewd: comment this out for now
            //double xavg = (double)kpjsum[jp]/(double)(npex_-ngap);
            int ipset = 0;
            int isum = 0;
            int igap = 0;
            int blksum = 0;
            isumgap = 0;
            gapcnt = 0;
            for (int ii=0; ii<nx_; ii++)
            {
               // don't let ipset exceed process grid dimension
               if (ipset >= npex_)
                  ipset = npex_-1;
               
               pexind[jp*nx_+ii] = ipset;
               blksum += xdistsum[jp*nx_+ii];
               isum += xdistsum[jp*nx_+ii];

               // avoid creating domains with multiple disconnected regions
               if (isum > 0 && xdistsum[jp*nx_+ii] == 0)
                  igap++;
               else
                  igap = 0;

               // force a processor increment at every gap
               if (igap > gapthresh_) {
                  isum = 0;
                  ipset++;
                  igap = 0;
               }

               // keep running tally of gaps computed exactly as above
               isumgap += xdistsum[jp*nx_+ii];
               if (isumgap > 0 && xdistsum[jp*nx_+ii] == 0)
                  gapcnt++;
               else
                  gapcnt = 0;

               if (gapcnt > gapthresh_) {
                  gapcnt = 0;
                  isumgap = 0;

                  // gap falls right at natural processor boundary:  our current average will give the last
                  // processor no work, adjust xavg to distribute remaining work appropriately
                  if (blksum >= xavg*(ipset+1)) {
                     isum = 0;
                     ipset++;
                     int ipleft = npex_ - ipset;
                     if (ipleft > 0)
                        xavg *= (double)ipleft/(double)(ipleft+1.);
                  }
               }
               else if (blksum >= xavg*(ipset+1)) {
                  isum = 0;
                  ipset++;
               }
            }
            ***/
            
         }

         // we now have process coordinates for all local grid points, change dest_
         // index in cells array
         for (unsigned ii=0; ii<cells.size(); ++ii)
         {
            int gid = cells[ii].gid_;
            GridPoint gpt(gid,nx_,ny_,nz_);
            int tmpkp = pezind[gpt.z];
            if (tmpkp == kp)
            {
               int jp = peyind[gpt.y];
               int ip = pexind[jp*nx_+gpt.x];
               int peid = ip + jp*npex_ + kp*npex_*npey_;
               cells[ii].dest_ = peid;

               assert(ip >= 0 && ip < npex_);    //ewd DEBUG
               assert(jp >= 0 && jp < npey_);    //ewd DEBUG
               assert(kp >= 0 && kp < npez_);    //ewd DEBUG
               assert(peid >= 0);  //ewd DEBUG
            }
         }
      }
   }

   //ewd DEBUG
   //exit(1);
   
   // carry out communication to match computed distribution
   redistributeCells(cells);
   double t4 = MPI_Wtime();
   computeLoadInfo(cells);
   double t5 = MPI_Wtime();
   if (myRank_ == 0)
      cout << "computeLoadInfo timing:  " << t5-t4 << " sec" << endl;

   if (myRank_ == 0)
      cout << "Load histogram after initial distribution:" << endl;
   nlocHistogram(cells);
   volHistogram(cells);

   //ewd DEBUG
   exit(1);

   
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::gridVolMinLoop(vector<AnatomyCell>& cells)
{
   double t0 = MPI_Wtime();
   restrictMoves(cells,gapthresh_);
   nonzeroVolume(cells);
   double t1 = MPI_Wtime();
   if (myRank_ == 0)
      cout << "GDLB.volMinLoop:  init time = " << t1-t0 << endl;

   if (myRank_ == 0)
      cout << "Load histogram at start of volMinLoop:" << endl;
   nlocHistogram(cells);
   volHistogram(cells);

   bool balance = false;
   const int bprint = 1;
   const int maxiter = 100;
   const int maxStoredMoves = 1000000;
   double vmaxmult = 2.0;

   vector<MoveInfo> moves;

   int bcnt = 0;
   while (!balance && bcnt < maxiter)
   {
      double titer1 = MPI_Wtime();
      double pevolavg = 0.;
      for (unsigned ii=0; ii<npegrid_; ++ii)
         pevolavg += peboxinfo_[ii]->volume();
      pevolavg /= npegrid_;
      double vthresh = pevolavg*vmaxmult;

      int maxvol = -1;
      int maxpe = -1;
      for (unsigned ii=0; ii<npegrid_; ++ii)
         if (peboxinfo_[ii]->volume() > maxvol) {
            maxvol = peboxinfo_[ii]->volume();
            maxpe = ii;
         }
            
      double maxratio = maxvol / pevolavg;
      if (maxratio <= vmaxmult) balance = true;
      
      if (myRank_ == 0 && bcnt%bprint == 0)
         cout << "Volume minimizing loop iter " << bcnt << ":  avg. volume = " << pevolavg << ", max volume = " << maxvol << " (pe " << maxpe << "), vol_ratio = " << maxratio << endl;

      if (!balance)
      {
         int nMoves = 0;
         for (int ip=0; ip<npegrid_; ip++)
         {

            //ewd DEBUG
            /*
            if (myRank_ == 39)
            {
               GridPoint gpt(ip,npex_,npey_,npez_);
               cout << "TRIAL_HEAD:  pe " << ip << " (" << gpt.x << "," << gpt.y << "," << gpt.z << "), box = " <<
                   peboxinfo_[ip]->minPoint(0) << " " << peboxinfo_[ip]->maxPoint(0) << " " <<
                   peboxinfo_[ip]->minPoint(1) << " " << peboxinfo_[ip]->maxPoint(1) << " " <<
                   peboxinfo_[ip]->minPoint(2) << " " << peboxinfo_[ip]->maxPoint(2) << endl;
            }
            */
            //ewd DEBUG

            int pevol = peboxinfo_[ip]->volume();
            if (pevol > pevolavg)
            {
               for (int idim=0; idim<3; idim++)
               {
                  // try moving a plane in positive direction
                  {
                     int nbr = penbr_[ip][2*idim+1];
                     if (nbr > -1)
                     {
                        int val = peboxinfo_[ip]->maxAxisPoint(idim);
                        int acc = peboxinfo_[ip]->trialMove(idim,val,*peboxinfo_[nbr]); 
                        if (acc == 0)
                        {
                           nMoves++;
                           MoveInfo minfo(ip,nbr,idim,val);
                           moves.push_back(minfo);
                        }
                     }
                  }
                  // try moving a plane in negative direction
                  {
                     int nbr = penbr_[ip][2*idim];
                     if (nbr > -1)
                     {
                        int val = peboxinfo_[ip]->minAxisPoint(idim);
                        int acc = peboxinfo_[ip]->trialMove(idim,val,*peboxinfo_[nbr]); 
                        if (acc == 0)
                        {
                           nMoves++;
                           MoveInfo minfo(ip,nbr,idim,val);
                           moves.push_back(minfo);
                        }
                     }
                  }
               }
            }
         }

         if (myRank_ == 0)
            cout << "Loop iter " << bcnt << ", " << nMoves << " moves accepted." << endl;
         if (nMoves == 0) balance = true;

         // don't let moves array get too large
         //if (moves.size() > maxStoredMoves)
         {
            applyStoredMoves(moves,cells);
            moves.clear();
         }

         // carry out communication to match computed distribution
         double t6 = MPI_Wtime();
         redistributeCells(cells);
         computeLoadInfo(cells);
         double t7 = MPI_Wtime();
         if (myRank_ == 0)
            cout << "Redistribute and computeLoadInfo timing = " << t7-t6 << " sec" << endl;

      }
      bcnt++;
      double titer2 = MPI_Wtime();
      if (myRank_ == 0)
         cout << "volMinLoop iteration timing = " << titer2-titer1 << " sec" << endl;
   }
   double t2 = MPI_Wtime();
   if (myRank_ == 0)
      cout << "GDLB.volMinLoop:  timing = " << t2-t1 << endl;

   
   //ewd DEBUG
   /*
   {
      if (myRank_ == 0)
         for (unsigned ip=0; ip<npegrid_; ++ip)
         {
            GridPoint gpt(ip,npex_,npey_,npez_);
            cout << "PEBOX:  pe " << ip << " (" << gpt.x << "," << gpt.y << "," << gpt.z << "), box = " <<
                peboxinfo_[ip]->minPoint(0) << " " << peboxinfo_[ip]->maxPoint(0) << " " <<
                peboxinfo_[ip]->minPoint(1) << " " << peboxinfo_[ip]->maxPoint(1) << " " <<
                peboxinfo_[ip]->minPoint(2) << " " << peboxinfo_[ip]->maxPoint(2) << endl;
            
         }      
   }
   */
   //ewd DEBUG

   /*   
   //ewd DEBUG
   if (myRank_ == 0)
      cout << "Starting assignCells" << endl;
   assignCells(cells);
   */
   
   //ewd DEBUG
   if (myRank_ == 0)
      cout << "Starting applyStoredMoves, moves.size = " << moves.size() << endl;
   
   if (moves.size() > 0) {
      applyStoredMoves(moves,cells);
      moves.clear();
   }
   
   double t3 = MPI_Wtime();
   if (myRank_ == 0)
      cout << "GDLB.applyStoredMoves:  time = " << t3-t2 << endl;

   //ewd DEBUG
   /*
   {
      for (unsigned ii=0; ii<cells.size(); ++ii)
         if (cells[ii].dest_ == 19967) {
            GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_); 
            cout << "DEBUG_19967b, myRank = " << myRank_ << ", cell " << ii << ", gid = " << cells[ii].gid_ << ":  " << gpt.x << " " << gpt.y << " " << gpt.z << endl;
         }
   }
   */
   //ewd DEBUG

   
   if (myRank_ == 0)
      cout << "Load histogram after volMinLoop:" << endl;
   nlocHistogram(cells);
   volHistogram(cells);

   //ewd DEBUG
   /*
   for (unsigned ip=0; ip<npegrid_; ++ip)
   {
      int zcnt = 0;
      for (unsigned in=0; in<6; ++in)
         if (penbr_[ip][in] < 0)
            zcnt++;
      if (zcnt == 6)
      {
         for (unsigned ii=0; ii<cells.size(); ++ii)
         {
            int pe = cells[ii].dest_;
            if (pe == ip)
            {
               GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_); 
               cout << "DISCONNECTED CELLS:  myRank = " << myRank_ << ", dest pe = " << pe << ", gid = " << cells[ii].gid_ << ": " << gpt.x << " " << gpt.y << " " << gpt.z << endl;
            }
         }
      }
   }
   */
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::assignCells(vector<AnatomyCell>& cells)
{
   // carry out communication to match computed distribution
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_); 
      int nfound = 0;
      //bool found = false;
      for (int ip=0; ip<npegrid_; ip++)
      {
         //if (!found && peboxinfo_[ip]->containsPoint(gpt.x,gpt.y,gpt.z))
         if (peboxinfo_[ip]->containsPoint(gpt.x,gpt.y,gpt.z))
         {
            cells[ii].dest_ = ip;
            nfound++;
            //found = true;
         }
      }
      //if (!found)
      if (nfound == 0)
         cout << "ERROR:  cell " << ii << ": " << gpt.x << " " << gpt.y << " " << gpt.z << " could not be placed!" << endl;
      if (nfound > 1)
         cout << "ERROR:  cell " << ii << ": " << gpt.x << " " << gpt.y << " " << gpt.z << " fits in multiple processor domains!" << endl;
   }
   redistributeCells(cells);
   return;
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::applyStoredMoves(vector<MoveInfo> moves, vector<AnatomyCell>& cells)
{
   // to enable testing, we need to calculate which pes we own data for
   // (in production runs, this will just be myRank)
   vector<int> ownsDataLoc(npegrid_,0);
   for (unsigned ii=0; ii<cells.size(); ++ii)
      ownsDataLoc[cells[ii].dest_] = myRank_;  // who owns which processor's data
   vector<int> ownsData(npegrid_,0);
   MPI_Allreduce(&ownsDataLoc[0], &ownsData[0], npegrid_, MPI_INT, MPI_SUM, comm_);

   // carry out communication to match computed distribution
   for (int im=0; im<moves.size(); im++)
   {
      if (im%1000 == 0 && myRank_ == 0)
         cout << "Processing move " << im << "..." << endl;

      int pe = moves[im].srcpe;
      #pragma omp parallel for
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         if (cells[ii].dest_ == pe) {
            GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_); 
            if (moves[im].dim == 0 && gpt.x == moves[im].val)
               cells[ii].dest_ = moves[im].destpe;
            else if (moves[im].dim == 1 && gpt.y == moves[im].val)
               cells[ii].dest_ = moves[im].destpe;
            else if (moves[im].dim == 2 && gpt.z == moves[im].val)
               cells[ii].dest_ = moves[im].destpe;
         }
      }
   }
   redistributeCells(cells);
   return;
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::nonzeroVolume(vector<AnatomyCell>& cells)
{
   // do a single pass through pegrid, give any procs with zero volume an edge from
   // their neighbor with the largest volume
   for (unsigned ip=0; ip<npegrid_; ++ip)
   {
      if (peboxinfo_[ip]->volume() == 0)
      {
         int maxind = -1;
         int maxvol = -1;
         for (int in=0; in<6; in++) {
            int nbr = penbr_[ip][in];
            int nbrvol = (nbr > -1 ? peboxinfo_[nbr]->volume() : 0);
            if (nbrvol > maxvol)
            {
               maxvol = nbrvol;
               maxind = in;
            }
         }
         if (maxvol <= 0)
         {
            GridPoint gpt(ip,npex_,npey_,npez_); 
            if (myRank_ == 0)
               cout << "GDLB.gridVolMinLoop WARNING: proc " << ip << " (" << gpt.x << " " << gpt.y << " " << gpt.z << ") has no neighbors with volume > 0." << endl;
         }
         else
         {
            int nbr = penbr_[ip][maxind];
            int dim = maxind/2;
            int val;
            if (maxind%2 == 0)
               val = peboxinfo_[nbr]->maxAxisPoint(dim);
            else
               val = peboxinfo_[nbr]->minAxisPoint(dim);
            Plane* pptr = peboxinfo_[nbr]->getPlane(dim,val);

            //ewd DEBUG
            if (myRank_ == 0)
               cout << "NONZEROVOL1, proc " << ip << ", has zero volume, moving plane dim = " << dim << ", val = " << val << " from proc " << nbr << ", vol = " << peboxinfo_[nbr]->volume() << ", plane = " << pptr->min[0] << " " << pptr->max[0] << " " << pptr->min[1] << " " << pptr->max[1] << " " << pptr->min[2] << " " << pptr->max[2] << endl;
            //ewd DEBUG

            if (pptr != 0)
            {
               peboxinfo_[nbr]->removePlane(dim,val);
               peboxinfo_[ip]->addPlane(dim,pptr);
               peboxinfo_[nbr]->updateBoundaries();
               peboxinfo_[ip]->updateBoundaries();

               //ewd DEBUG
               if (myRank_ == 0)
               {
                  GridPoint gpt(ip,npex_,npey_,npez_);
                  cout << "NONZEROVOL2:  pe " << ip << " (" << gpt.x << "," << gpt.y << "," << gpt.z << "), box = " <<
                      peboxinfo_[ip]->minPoint(0) << " " << peboxinfo_[ip]->maxPoint(0) << " " <<
                      peboxinfo_[ip]->minPoint(1) << " " << peboxinfo_[ip]->maxPoint(1) << " " <<
                      peboxinfo_[ip]->minPoint(2) << " " << peboxinfo_[ip]->maxPoint(2) << endl;
               }
               //ewd DEBUG

               for (unsigned ii=0; ii<cells.size(); ++ii)
                  if (cells[ii].dest_ == nbr)
                  {
                     GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_); 
                     if (dim == 0 && gpt.x == val)
                        cells[ii].dest_ = ip;
                     else if (dim == 1 && gpt.y == val)
                        cells[ii].dest_ = ip;
                     else if (dim == 2 && gpt.z == val)
                        cells[ii].dest_ = ip;
                  }
            }
         }
      }
   }

   // carry out communication to match computed distribution
   redistributeCells(cells);
   computeLoadInfo(cells);
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::computeLoadInfo(vector<AnatomyCell>& cells)
{
   for (int ip=0; ip<npegrid_; ip++)
      delete peboxinfo_[ip];
   for (int ip=0; ip<npegrid_; ip++)
      peboxinfo_[ip] = new ProcBox(ip);

   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      int peid = cells[ii].dest_;
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      peboxinfo_[peid]->add3DCoord(gpt.x,gpt.y,gpt.z);
   }
   
   for (int ip=0; ip<npegrid_; ip++)
   {
      vector<int> ipcoords;
      if (peboxinfo_[ip]->nAdded() > 0)
         peboxinfo_[ip]->planesToCoords(ipcoords);

      vector<int> ipinfobuf(2,0);
      vector<int> ipinfo(2,0);
      if (peboxinfo_[ip]->nAdded() > 0)
      {
         ipinfobuf[0] = ipcoords.size();
         ipinfobuf[1] = myRank_;
      }
      MPI_Allreduce(&ipinfobuf[0], &ipinfo[0], 2, MPI_INT, MPI_SUM, comm_);
      int size = ipinfo[0];
      int owner = ipinfo[1];
      ipcoords.resize(size);
      MPI_Bcast(&ipcoords[0],size,MPI_INT,owner,comm_);

      for (int ii=0; ii<size; ii+=3)
      {
         int x = ipcoords[ii];
         int y = ipcoords[ii+1];
         int z = ipcoords[ii+2];
         peboxinfo_[ip]->add3DCoord(x,y,z);
      }
   }

   /* this is too slow, just communicate minimum data to describe each box
   // store processor box info on all tasks
   for (int ip=0; ip<npegrid_; ip++)
   {
      peboxinfo_[ip] = new ProcBox(ip);

      vector<int> ipcoords;
      if (ownsData[ip] == myRank_)
      {
         for (unsigned ii=0; ii<cells.size(); ++ii)
         {
            int peid = cells[ii].dest_;
            if (ip == peid)
            {
               GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
               ipcoords.push_back(gpt.x);
               ipcoords.push_back(gpt.y);
               ipcoords.push_back(gpt.z);
            }
         }
      }
      vector<int> ipinfobuf(2,0);
      vector<int> ipinfo(2,0);
      if (ownsData[ip] == myRank_)
      {
         ipinfobuf[0] = ipcoords.size();
         ipinfobuf[1] = myRank_;
      }
      MPI_Allreduce(&ipinfobuf[0], &ipinfo[0], 2, MPI_INT, MPI_SUM, comm_);
      int size = ipinfo[0];
      int owner = ipinfo[1];
      ipcoords.resize(size);
      MPI_Bcast(&ipcoords[0],size,MPI_INT,owner,comm_);

      for (int ii=0; ii<size; ii+=3)
      {
         int x = ipcoords[ii];
         int y = ipcoords[ii+1];
         int z = ipcoords[ii+2];
         peboxinfo_[ip]->add3DCoord(x,y,z);
      }
   }
   */


   /*
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      int peid = cells[ii].dest_;
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      peboxinfo_[peid]->add3DCoord(gpt.x,gpt.y,gpt.z);
   }

   // For now, let's manually copy this ProcBox object's data across all tasks.
   // It's hacky and involves a lot of communication, so we may have to rewrite it.
   for (unsigned ip=0; ip<npegrid_; ++ip)
   {
      vector<int> pesizebuf(4,0);
      vector<int> pesize(4,0);
      if (peboxinfo_[ip]->nAdded() > 0)
      {
         peboxinfo_[ip]->updateSize();
         pesizebuf[0] = peboxinfo_[ip]->xyzSize();
         pesizebuf[1] = peboxinfo_[ip]->xSize();
         pesizebuf[2] = peboxinfo_[ip]->ySize();
         pesizebuf[3] = peboxinfo_[ip]->zSize();
      }
      MPI_Allreduce(&pesizebuf[0], &pesize[0], 4, MPI_INT, MPI_SUM, comm_);

      int bufsize = pesize[0];
      vector<int> pedatabuf(bufsize,0);
      vector<int> pedata(bufsize,0);
      if (peboxinfo_[ip]->nAdded() > 0)
         peboxinfo_[ip]->getCoords(pedatabuf);
      MPI_Allreduce(&pedatabuf[0], &pedata[0], bufsize, MPI_INT, MPI_SUM, comm_);

      //ewd DEBUG
      if (pesize[0] != pesize[1]+pesize[2]+pesize[3])
         cout << "myRank_ = " << myRank_ << ", pesize mismatch:  " << pesize[0] << " " << pesize[1] << " " << pesize[2] << " " << pesize[3] << endl;
      
      // fill peboxinfo_[ip] with data
      int cnt = 0;
      for (int ii=0; ii<pesize[1]; ii++)
      {
         int x = pedata[cnt++];

         //ewd DEBUG
         assert(ip < npegrid_ && ip >= 0);
         assert(peboxinfo_.size() == npegrid_);
         assert(x >= 0 && x <= 99999);
         //ewd DEBUG
         
         peboxinfo_[ip]->addX(x);
      }
      for (int ii=0; ii<pesize[2]; ii++)
      {
         int y = pedata[cnt++];
         peboxinfo_[ip]->addY(y);
      }
      for (int ii=0; ii<pesize[3]; ii++)
      {
         int z = pedata[cnt++];
         peboxinfo_[ip]->addZ(z);
      }
      if (cnt > 0)
         peboxinfo_[ip]->updateSize();
   }
   */
   
   /*   
   // calculate volume, center of mass of each domain
   vector<int> pecnt_loc(npegrid_,0);
   vector<int> pevol_buf(npegrid_,0);
   vector<vector<double> > pecomloc_;
   vector<vector<int> > peminloc_;
   vector<vector<int> > pemaxloc_;
   pecomloc_.resize(3);
   peminloc_.resize(3);
   pemaxloc_.resize(3);
   pecom_.resize(3);
   pemin_.resize(3);
   pemax_.resize(3);
   for (unsigned ii=0; ii<3; ++ii)
   {
      pecomloc_[ii].assign(npegrid_,0.);
      peminloc_[ii].assign(npegrid_,100000000);
      pemaxloc_[ii].assign(npegrid_,-100000000);
      pecom_[ii].assign(npegrid_,0.);
      pemin_[ii].assign(npegrid_,0);
      pemax_[ii].assign(npegrid_,0);
   }
   
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      int peid = cells[ii].dest_;
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      pecomloc_[0][peid] += (double)gpt.x;
      pecomloc_[1][peid] += (double)gpt.y;
      pecomloc_[2][peid] += (double)gpt.z;
      pecnt_loc[peid]++;
      if (gpt.x > pemaxloc_[0][peid]) pemaxloc_[0][peid] = gpt.x;
      if (gpt.y > pemaxloc_[1][peid]) pemaxloc_[1][peid] = gpt.y;
      if (gpt.z > pemaxloc_[2][peid]) pemaxloc_[2][peid] = gpt.z;
      if (gpt.x < peminloc_[0][peid]) peminloc_[0][peid] = gpt.x;
      if (gpt.y < peminloc_[1][peid]) peminloc_[1][peid] = gpt.y;
      if (gpt.z < peminloc_[2][peid]) peminloc_[2][peid] = gpt.z;
   }

   //ewd DEBUG
   for (unsigned ii=0; ii<npegrid_; ++ii)
      if (pecnt_loc[ii] > 0)
         assert(ownsData[ii] == myRank_);
   //ewd DEBUG
   
   for (unsigned ii=0; ii<npegrid_; ++ii)
   {
      if (pecnt_loc[ii] > 0)
      {
         pecomloc_[0][ii] /= pecnt_loc[ii];
         pecomloc_[1][ii] /= pecnt_loc[ii];
         pecomloc_[2][ii] /= pecnt_loc[ii];
         //int minx = pbcDist(pemaxloc_[0][ii],peminloc_[0][ii],nx_);
         //int miny = pbcDist(pemaxloc_[1][ii],peminloc_[1][ii],ny_);
         //int minz = pbcDist(pemaxloc_[2][ii],peminloc_[2][ii],nz_);
         int minx = pemaxloc_[0][ii] - peminloc_[0][ii];
         int miny = pemaxloc_[1][ii] - peminloc_[1][ii];
         int minz = pemaxloc_[2][ii] - peminloc_[2][ii];
         pevol_buf[ii] = minx*miny*minz;
      }
   }
   
   // collect information from all tasks
   pevol_.assign(npegrid_,0);
   MPI_Allreduce(&pevol_buf[0], &pevol_[0], npegrid_, MPI_INT, MPI_SUM, comm_);
   MPI_Allreduce(&pecomloc_[0][0], &pecom_[0][0], npegrid_, MPI_DOUBLE, MPI_SUM, comm_);
   MPI_Allreduce(&pecomloc_[1][0], &pecom_[1][0], npegrid_, MPI_DOUBLE, MPI_SUM, comm_);
   MPI_Allreduce(&pecomloc_[2][0], &pecom_[2][0], npegrid_, MPI_DOUBLE, MPI_SUM, comm_);
   MPI_Allreduce(&pemaxloc_[0][0], &pemax_[0][0], npegrid_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&pemaxloc_[1][0], &pemax_[1][0], npegrid_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&pemaxloc_[2][0], &pemax_[2][0], npegrid_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&peminloc_[0][0], &pemin_[0][0], npegrid_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&peminloc_[1][0], &pemin_[1][0], npegrid_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&peminloc_[2][0], &pemin_[2][0], npegrid_, MPI_INT, MPI_MIN, comm_);

      */
}
////////////////////////////////////////////////////////////////////////////////
double GDLoadBalancer::pbcDist(double x1,double x2,int nx)
{
   double dist = abs(x1-x2);
   if (dist > 0.5*nx)
      dist = nx - dist;
   return dist;
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::balanceLoop(vector<AnatomyCell>& cells)
{
   // call balanceLoop w. default values
   const int bblock = 5;
   const int bthresh = 10;
   const int maxiter = 100000;
   balanceLoop(cells,bblock,bthresh,maxiter);
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::balanceLoop(vector<AnatomyCell>& cells, int bblock, int bthresh, int maxiter)
{
   bool balance = false;
   const int bprint = 100;
   
   // need to sync up information on each process' load on all tasks
   vector<int> nloc_buf_(npegrid_,0);
   vector<int> nloc_(npegrid_,0);

   // to allow for testing, don't assume that dest == myRank
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      int peid = cells[ii].dest_;
      nloc_buf_[peid]++;
   }
   MPI_Allreduce(&nloc_buf_[0], &nloc_[0], npegrid_, MPI_INT, MPI_SUM, comm_);

   //ewd:  prevent exchanges between neighbors separated by more than gapthresh cells
   //computeLoadInfo(cells);
   restrictMoves(cells,gapthresh_);

   // print out maximum bounding box volume
   vector<int>::iterator maxvol = max_element(pevol_.begin(),pevol_.end());
   int maxpe = maxvol - pevol_.begin();

   if (myRank_ == 0)
      cout << "GDLoadBalancer::balanceLoop start:  max bounding box volume = " << *maxvol << ", on pe " << maxpe << " (" << pemin_[0][maxpe] << ":" << pemax_[0][maxpe] << ", " << pemin_[1][maxpe] << ":" << pemax_[1][maxpe] << ", " << pemin_[2][maxpe] << ":" << pemax_[2][maxpe] << ")" << endl;



   
   nloctot_ = 0;
   for (unsigned ii=0; ii<npegrid_; ++ii)
      nloctot_ += nloc_[ii];
   nlocavg_ = (double)nloctot_/(double)npegrid_; 
  
   if (myRank_ == 0)
      cout << "Starting load balance loop: threshold = " << bthresh << ", inner loop size = " << bblock << ", max iterations = " << maxiter << ", nlocavg = " << nlocavg_ << endl;
  
   int bcnt = 0;
   while (!balance && bcnt < maxiter)
   {          
      for (int ip=0; ip<npegrid_; ip++)
      {
         for (int b=0; b<bblock; b++)
         {
            if (nloc_[ip] < nlocavg_)
            {
               // loop over ip's neighbors, check if any are above average
               int maxval = 0;
               int nbrpmax = -1;
               int thisn = -1;
               for (int n=0; n<nnbr_; n++)
               {
                  int nbr = penbr_[ip][n];
                  if (nbr > -1 && nloc_[nbr] > nlocavg_ && nloc_[nbr] > maxval)
                  {
                     maxval= nloc_[nbr];
                     nbrpmax = nbr;
                     thisn = n;
                  }
               }
               if (nbrpmax == -1)  // move to any neighbor w. more points than us
               {
                  for (int n=0; n<nnbr_; n++)
                  {
                     int nbr = penbr_[ip][n];
                     if (nbr > -1 && nloc_[nbr] > nloc_[ip] && nloc_[nbr] > maxval)
                     {
                        maxval= nloc_[nbr];
                        nbrpmax = nbr;
                        thisn = n;
                     }
                  }
               }

               // add grid point to ip, subtract from nbrpmax
               if (nbrpmax > -1)
               {
                  nloc_[ip]++;
                  togive_[ip][thisn]--;  // we are owed one grid point from neighbor thisn
                  nloc_[nbrpmax]--;
                  togive_[nbrpmax][thatn_[thisn]]++;  // we owe one grid point to neighbor thatn_
               }
            }
            else if (nloc_[ip] > nlocavg_)
            {
               // loop over ip's neighbors, check if any are below average
               int minval = nloctot_;
               int nbrpmin = -1;
               int thisn = -1;
               for (int n=0; n<nnbr_; n++)
               {
                  int nbr = penbr_[ip][n];
                  if (nbr > -1 && nloc_[nbr] < nlocavg_ && nloc_[nbr] < minval)
                  {
                     minval= nloc_[nbr];
                     nbrpmin = nbr;
                     thisn = n;
                  }
               }
               if (nbrpmin == -1)  // move to any neighbor w. fewer points than us
               {
                  for (int n=0; n<nnbr_; n++)
                  {
                     int nbr = penbr_[ip][n];
                     if (nbr > -1 && nloc_[nbr] < nloc_[ip] && nloc_[nbr] < minval)
                     {
                        minval= nloc_[nbr];
                        nbrpmin = nbr;
                        thisn = n;
                     }
                  }
               }
               // subtract grid point from ip, add to nbrpmin
               if (nbrpmin > -1)
               {
                  nloc_[ip]--;
                  togive_[ip][thisn]++;  // we owe one grid point to neighbor thisn
                  nloc_[nbrpmin]++;
                  togive_[nbrpmin][thatn_[thisn]]--;  // we are owed one grid point from neighbor thatn_
               }
            }
         }
      }
      int maxnum = 0;
      int minnum = nloctot_;
      int maxp = -1;
      for (int p=0; p<npegrid_; p++)
      {
         if (nloc_[p] > maxnum) {
            maxnum = nloc_[p];
            maxp = p;
         }
         if (nloc_[p] < minnum)
            minnum = nloc_[p];
      }
      if ((maxnum - minnum) < bthresh)
      {
         balance = true;
         if (myRank_ == 0)
         {
            cout << "load balance iteration " << bcnt << ":  " << maxnum << " - " << minnum << " = " << maxnum-minnum << ", threshold = " << bthresh << endl;
            nlocHistogram(nloc_);
            cout << endl << " *** Load balance achieved in " << bcnt << " iterations. ***" << endl << endl;
         }
      }

      bcnt++;
      if (bcnt%bprint == 0 && myRank_ == 0)
      {
         cout << "load balance iteration " << bcnt << ":  " << maxnum << " - " << minnum << " = " << maxnum-minnum << ", threshold = " << bthresh << endl;
         nlocHistogram(nloc_);
      }

   }
   
   // to enable testing, we need to calculate which pes we own data for
   // (in production runs, this will just be myRank)
   vector<int> ownsDataLoc(npegrid_,0);
   for (unsigned ii=0; ii<cells.size(); ++ii)
      ownsDataLoc[cells[ii].dest_] = myRank_;  // who owns which processor's data
   vector<int> ownsData(npegrid_);
   MPI_Allreduce(&ownsDataLoc[0], &ownsData[0], npegrid_, MPI_INT, MPI_SUM, comm_);

   // test that no processor's data is shared
   haveData_.resize(npegrid_,0);
   for (unsigned ii=0; ii<cells.size(); ++ii)
      haveData_[cells[ii].dest_] = 1;  
   vector<int> testOwnership(npegrid_);
   MPI_Allreduce(&haveData_[0], &testOwnership[0], npegrid_, MPI_INT, MPI_SUM, comm_);
   for (unsigned ip=0; ip<npegrid_; ++ip)
   {
      if (testOwnership[ip] > 1)
      {
         if (myRank_ == 0)
            cout << "ERROR in GDLoadBalancer::balanceLoop:  incorrect data distribution!" << endl;
         MPI_Abort(comm_,2);
      }
   }
   myRankList_.clear();
   for (unsigned ip=0; ip<npegrid_; ++ip)
      if (haveData_[ip] == 1)
         myRankList_.push_back(ip);


   // for non-testing case, each task owns its own data (even if it doesn't currently
   // have any)
   if (myRankList_.size() < 1)
      myRankList_.push_back(myRank_);
   
   // use togive_ array to redistribute data
   bool allDataSent = false;
   int ccnt = 0;
   while (!allDataSent && ccnt < 1000)
   {
      //ewd DEBUG
      if (myRank_ == 0)
         cout << "Calling diffuseDest ccnt = " << ccnt << endl;

      allDataSent = diffuseDest(cells);
      ccnt++;
   }
   // calculate load histogram directly from cells
   if (myRank_ == 0)
      cout << "Load histogram after " << ccnt << " comm iterations:" << endl;
   nlocHistogram(cells);
   volHistogram(cells);
  
   if (myRank_ == 0)
      cout << "GDLoadBalancer::balanceLoop stop:  max bounding box volume = " << *maxvol << ", on pe " << maxpe << " (" << pemin_[0][maxpe] << ":" << pemax_[0][maxpe] << ", " << pemin_[1][maxpe] << ":" << pemax_[1][maxpe] << ", " << pemin_[2][maxpe] << ":" << pemax_[2][maxpe] << ")" << endl;

   if (myRank_ == 0)
      if (!balance)
         cout << "balanceLoop did not achieve balance after " << maxiter << " iterations." << endl;
  
}
////////////////////////////////////////////////////////////////////////////////
bool GDLoadBalancer::diffuseDest(vector<AnatomyCell>& cells)
{
   bool allDataSent_ = false;
   int allSentFlag = 1;
   for (int ir=0; ir<myRankList_.size(); ir++)
   {
      int ip = myRankList_[ir];
      //ewd DEBUG
      //if (myRank_ == 0)
      //   cout << "DIFFDEST, myRank = " << myRank_ << ", ir = " << ir << ", size = " << myRankList_.size() << ", ip = " << ip << endl;

      vector<int> locgid_;
      vector<int> locind_;
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         if (cells[ii].dest_ == ip)
         {
            locgid_.push_back(cells[ii].gid_);
            locind_.push_back(ii);
         }
      }
      int ngidloc = locgid_.size();

      vector<int> keep(ngidloc);
      for (int i=0; i<ngidloc; i++)
         keep[i] = 1;

      vector<int> xgid(ngidloc);
      vector<int> ygid(ngidloc);
      vector<int> zgid(ngidloc);
      vector<int> xgid_sort(ngidloc);
      vector<int> ygid_sort(ngidloc);
      vector<int> zgid_sort(ngidloc);
      for (int i=0; i<ngidloc; i++)
      {
         int gid = locgid_[i];
         GridPoint gpt(gid,nx_,ny_,nz_);
         xgid[i] = gpt.z + gpt.y*nz_ + gpt.x*nz_*ny_;
         ygid[i] = gpt.x + gpt.z*nx_ + gpt.y*nx_*nz_;
         zgid[i] = gpt.y + gpt.x*ny_ + gpt.z*ny_*nx_;
         xgid_sort[i] = xgid[i];
         ygid_sort[i] = ygid[i];
         zgid_sort[i] = zgid[i];
      }
      sort(xgid_sort.begin(),xgid_sort.end());
      sort(ygid_sort.begin(),ygid_sort.end());
      sort(zgid_sort.begin(),zgid_sort.end());
      
      int nleft = ngidloc;
      for (int n=0; n<nnbr_; n++)
      {
         if (togive_[ip][n] > 0)      // have data to transfer to this neighbor
         {
            int offset = 0;
            int s = 0;
            while (s < togive_[ip][n] && nleft > 0)
            {
               int tid = -1;
               if (n==0) {
                  for (int j=0; j<ngidloc; j++)
                     if (xgid[j] == xgid_sort[s+offset]) // send this id
                        tid = j;
               }
               else if (n==1) {
                  for (int j=0; j<ngidloc; j++)
                     if (xgid[j] == xgid_sort[ngidloc-1-s-offset]) // send this id
                        tid = j;
               }
               else if (n==2) {
                  for (int j=0; j<ngidloc; j++)
                     if (ygid[j] == ygid_sort[s+offset]) // send this id
                        tid = j;
               }
               else if (n==3) {
                  for (int j=0; j<ngidloc; j++)
                     if (ygid[j] == ygid_sort[ngidloc-1-s-offset]) // send this id
                        tid = j;
               }
               else if (n==4) {
                  for (int j=0; j<ngidloc; j++)
                     if (zgid[j] == zgid_sort[s+offset]) // send this id
                        tid = j;
               }
               else if (n==5) {
                  for (int j=0; j<ngidloc; j++)
                     if (zgid[j] == zgid_sort[ngidloc-1-s-offset]) // send this id
                        tid = j;
               }
            
               if (tid > -1 && keep[tid] == 1) // hasn't been moved yet
               {
                  int thisgid = locgid_[tid];
                  int ipn = penbr_[ip][n];
                  assert(ipn > -1);
                  int ii = locind_[tid];
                  cells[ii].dest_ = ipn;
                  keep[tid] = 0;
                  nleft--;
                  s++;
               }
               else
               {
                  offset++;
               }
            }
            togive_[ip][n] -= s;
            if (togive_[ip][n] > 0)
               allSentFlag = 0;

            //ewd DEBUG
            //{
            //   allSentFlag = 0;
            //   int givesum = togive_[ip][0]+togive_[ip][1]+togive_[ip][2]+togive_[ip][3]+togive_[ip][4]+togive_[ip][5];
            //   cout << "GDLB, ip = " << ip << ", n = " << n << ", togive_[ip][n] = " << togive_[ip][n] << ", ngidloc = " << ngidloc << ", nleft = " << nleft << ", penbr = " << penbr_[ip][n] << ", givesum = " << givesum << endl;
            //}               
            //ewd DEBUG
         }
      }
   }

   // carry out assigned communication 
   redistributeCells(cells);

   // see if anyone still has unsent data
   int allSent;
   MPI_Allreduce(&allSentFlag, &allSent, 1, MPI_INT, MPI_SUM, comm_);
   if (allSent == nTasks_)
      allDataSent_ = true;
   return allDataSent_;
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
      computeReducedProcGrid(nTasks_);
      //if (myRank_ == 0)
      //   cout << "GDLoadBalancer::redistributeCells:  reduced grid = " <<
      //       tnx_ << " x " << tny_ << " x " << tnz_ << " used for testing." << endl;

      // for testing, certain assumptions have to be fulfilled
      assert(npex_%2==0 && npey_%2==0 && npez_%2==0);
      assert(npex_%tnx_ == 0 && npey_%tny_ == 0 && npez_%tnz_ == 0);
        
      for (unsigned ii=0; ii<cells.size(); ++ii)
      {
         int gid = cells[ii].dest_;
         GridPoint gpt(gid,npex_,npey_,npez_);
         int tix = gpt.x*tnx_/npex_;
         int tiy = gpt.y*tny_/npey_;
         int tiz = gpt.z*tnz_/npez_;
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
void GDLoadBalancer::restrictMoves(vector<AnatomyCell>& cells, int ngap)
{
   for (int ip=0; ip<npegrid_; ip++)
   {
      for (int idim = 0; idim < 3; idim++)
      {
         int nbrl = penbr_[ip][2*idim];
         int nbrr = penbr_[ip][2*idim+1];
         if (peboxinfo_[ip]->volume() > 0)
         {
            if (nbrl > -1 && peboxinfo_[nbrl]->volume() > 0)
            {
               int edgegapl = peboxinfo_[ip]->minPoint(idim) - peboxinfo_[nbrl]->maxPoint(idim);
               if (edgegapl > ngap && edgegapl < 10000000)
               {
                  penbr_[ip][2*idim] = -2;
                  penbr_[nbrl][thatn_[2*idim]] = -2;
               }
            }
            if (nbrr > -1 && peboxinfo_[nbrr]->volume() > 0)
            {
               int edgegapr = peboxinfo_[nbrr]->minPoint(idim) - peboxinfo_[ip]->maxPoint(idim);
               if (edgegapr > ngap && edgegapr < 10000000)
               {
                  penbr_[ip][2*idim+1] = -2;
                  penbr_[nbrr][thatn_[2*idim+1]] = -2;
               }
            }
         }
      }
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
void GDLoadBalancer::nlocHistogram(vector<int>& nloc)
{
   // compute load histogram from data in nloc
   assert(nloc.size() == npegrid_);
   histnloc_.resize(npegrid_,0);
   for (int ip=0; ip<npegrid_; ip++)
      histnloc_[ip] = nloc[ip];
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
void GDLoadBalancer::volHistogram(vector<AnatomyCell>& cells)
{
   //ewd: this uses procbox array that was already computed and distributed
   //if (myRank_ == 0)
   //   computeVolHistogram();    

   //ewd:  calculate histogram directly from current cells array
   computeVolHistogram(cells);    

}

////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::computeVolHistogram(vector<AnatomyCell>& cells)
{
   // compute volumes from cells
   vector<int> peminx(npegrid_,99999999);
   vector<int> peminy(npegrid_,99999999);
   vector<int> peminz(npegrid_,99999999);
   vector<int> pemaxx(npegrid_,-99999999);
   vector<int> pemaxy(npegrid_,-99999999);
   vector<int> pemaxz(npegrid_,-99999999);
   vector<int> nloc(npegrid_,0);
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      int peid = cells[ii].dest_;
      GridPoint gpt(cells[ii].gid_,nx_,ny_,nz_);
      if (gpt.x < peminx[peid]) peminx[peid] = gpt.x;
      if (gpt.y < peminy[peid]) peminy[peid] = gpt.y;
      if (gpt.z < peminz[peid]) peminz[peid] = gpt.z;
      if (gpt.x > pemaxx[peid]) pemaxx[peid] = gpt.x;
      if (gpt.y > pemaxy[peid]) pemaxy[peid] = gpt.y;
      if (gpt.z > pemaxz[peid]) pemaxz[peid] = gpt.z;
      nloc[peid]++;

      //ewd DEBUG
      if (peid == 8128)
         cout << "VOLDEBUG, pe 8128, gpt " << ii << ", gid = " << cells[ii].gid_ << ", x = " << gpt.x << ", y = " << gpt.y << ", z = " << gpt.z << endl;

   }
   vector<int> peminx_all(npegrid_);
   vector<int> peminy_all(npegrid_);
   vector<int> peminz_all(npegrid_);
   vector<int> pemaxx_all(npegrid_);
   vector<int> pemaxy_all(npegrid_);
   vector<int> pemaxz_all(npegrid_);
   vector<int> nall(npegrid_);
   MPI_Allreduce(&peminx[0], &peminx_all[0], npegrid_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&peminy[0], &peminy_all[0], npegrid_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&peminz[0], &peminz_all[0], npegrid_, MPI_INT, MPI_MIN, comm_);
   MPI_Allreduce(&pemaxx[0], &pemaxx_all[0], npegrid_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&pemaxy[0], &pemaxy_all[0], npegrid_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&pemaxz[0], &pemaxz_all[0], npegrid_, MPI_INT, MPI_MAX, comm_);
   MPI_Allreduce(&nloc[0], &nall[0], npegrid_, MPI_INT, MPI_SUM, comm_);

   if (myRank_ == 0)
   {
      const int nhistmax = 100; // number of bins
      vector<int> phist(nhistmax,0);
      int maxvol = 0;
      int minvol = nx_*ny_*nz_;
      int maxvolip;
      
      vector<int> pevol_all(npegrid_);
      for (int ip=0; ip<npegrid_; ip++)
      {
         int tvol = 0;
         if (nall[ip] > 0)
            tvol = (pemaxx_all[ip]-peminx_all[ip]+1)*(pemaxy_all[ip]-peminy_all[ip]+1)*(pemaxz_all[ip]-peminz_all[ip]+1);
         if (tvol > maxvol)
         {
            maxvol = tvol;
            maxvolip = ip;
         }
         if (tvol < minvol) minvol = tvol;
         pevol_all[ip] = tvol;
      }
      
      int nhist = maxvol - minvol + 1;
      if (nhist > nhistmax) nhist = nhistmax;
      int delta = (maxvol-minvol + 1)/nhist;
      if ((maxvol-minvol+1)%nhist !=0) delta++;
      for (int ip=0; ip<npegrid_; ip++)
      {
         int pvol = pevol_all[ip];
         int bin = (pvol-minvol)/delta;
         phist[bin]++;
      }
      cout << "load balance histogram (volume):  " << endl;
      for (int i=0; i<nhist; i++)
         cout << "  " << minvol+delta*i << " - " << minvol+delta*(i+1) << ":    " << phist[i] << endl;

      int voltot_ = 0;
      for (unsigned ip=0; ip<npegrid_; ++ip)
         voltot_ += pevol_all[ip];
      double volavg_ = (double)voltot_/(double)npegrid_; 
      cout << "total assigned volume = " << voltot_ << ", avg. volume = " << volavg_ << ", max volume = " << maxvol << " (pe " << maxvolip << ")" << endl << endl;
   }
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::computeVolHistogram()
{
   // compute histogram of current data distribution
   const int nhistmax = 100; // number of bins
   vector<int> phist(nhistmax,0);
   int maxvol = 0;
   int minvol = nx_*ny_*nz_;
   for (int p=0; p<npegrid_; p++)
   {
      int pvol = peboxinfo_[p]->volume();
      if (pvol > maxvol)
         maxvol = pvol;
      if (pvol < minvol)
         minvol = pvol;
   }
   int nhist = maxvol - minvol + 1;
   if (nhist > nhistmax) nhist = nhistmax;
   int delta = (maxvol-minvol + 1)/nhist;
   if ((maxvol-minvol+1)%nhist !=0) delta++;
   for (int ip=0; ip<npegrid_; ip++)
   {
      int pvol = peboxinfo_[ip]->volume();
      int bin = (pvol-minvol)/delta;
      phist[bin]++;
      //ewd DEBUG: print process number of top bin
      if (bin == nhist-1 && phist[bin] == 1)
        cout << "nlocHistogram: top bin pe " << ip << ", maxvol = " << peboxinfo_[ip]->volume() << endl;
   }
   cout << "load balance histogram (volume):  " << endl;
   for (int i=0; i<nhist; i++)
      cout << "  " << minvol+delta*i << " - " << minvol+delta*(i+1) << ":    " << phist[i] << endl;

   int voltot_ = 0;
   for (unsigned ii=0; ii<npegrid_; ++ii)
      voltot_ += peboxinfo_[ii]->volume();
   double volavg_ = (double)voltot_/(double)npegrid_; 

      //ewd DEBUG
   /*
   for (int ip=0; ip<npegrid_; ip++)
   {
      int pvol = peboxinfo_[ip]->volume();
      if (pvol > 1.5*volavg_) {
         GridPoint gpt(ip,npex_,npey_,npez_);
         cout << "LARGEVOL:  pvol = " << pvol << ", pe " << ip << " (" << gpt.x << "," << gpt.y << "," << gpt.z << "), box = " <<
             peboxinfo_[ip]->minPoint(0) << " " << peboxinfo_[ip]->maxPoint(0) << " " <<
             peboxinfo_[ip]->minPoint(1) << " " << peboxinfo_[ip]->maxPoint(1) << " " <<
             peboxinfo_[ip]->minPoint(2) << " " << peboxinfo_[ip]->maxPoint(2) << endl;

         for (int in=0; in<6; in++)
         {
            int nbr = penbr_[ip][in];
            if (nbr > -1)
            {
               int nbrvol = peboxinfo_[nbr]->volume();
               GridPoint ngpt(nbr,npex_,npey_,npez_);
               cout << "LARGEVOL_NBR:  nbrvol = " << nbrvol << ", nbr " << nbr << " (" << ngpt.x << "," << ngpt.y << "," << ngpt.z << "), box = " <<
                   peboxinfo_[nbr]->minPoint(0) << " " << peboxinfo_[nbr]->maxPoint(0) << " " <<
                   peboxinfo_[nbr]->minPoint(1) << " " << peboxinfo_[nbr]->maxPoint(1) << " " <<
                   peboxinfo_[nbr]->minPoint(2) << " " << peboxinfo_[nbr]->maxPoint(2) << endl;
            }
         }
      }
   }
   */
   //ewd DEBUG


   
   cout << "total assigned volume = " << voltot_ << ", avg. volume = " << volavg_ << ", max volume = " << maxvol << endl << endl;

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
