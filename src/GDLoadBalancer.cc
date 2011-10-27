#include "GDLoadBalancer.hh"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <map>
#include <mpi.h>
#include "mpiUtils.h"
#include "GridPoint.hh"
#include "AnatomyCell.hh"
#include "Timer.hh"
using namespace std;



////////////////////////////////////////////////////////////////////////////////
GDLoadBalancer::GDLoadBalancer(MPI_Comm comm, int npex, int npey, int npez):
    comm_(comm), npex_(npex), npey_(npey), npez_(npez)
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
  //locgid_.resize(npegrid_,vector<int>(0));
  //locind_.resize(npegrid_,vector<int>(0));
  penbr_.resize(npegrid_,vector<int>(nnbr_,-1));

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

}
////////////////////////////////////////////////////////////////////////////////
GDLoadBalancer::~GDLoadBalancer()
{

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
  
  vector<int> nxype(npez_);
  vector<int> kpkmax(npez_);
  vector<int> kpkmin(npez_);
  int kset = 0;
  int xysum = 0;
  int ksum = 0;
  for (int kk=0; kk<nz; kk++)
  {
    pezind[kk] = kset;
    xysum += nxyplane[kk];
    ksum += nxyplane[kk];
    if (xysum > xypeavg*(kset+1)) {
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
  nxype[npez_-1] = ksum;
  kpkmax[npez_-1] = nz_;

  // loop over process grid planes owned by this process, compute distribution in x and y
  //for (int kploc = pezind[klocmin]; kploc <= pezind[klocmax]; kploc++)
  for (int kp = 0; kp < npez_; kp++)
  {
    vector<int> xdist_loc(npey_*nx_,0);
    vector<int> kpjsum(npey_,-1);
    if (kp >= pezind[klocmin] && kp <= pezind[klocmax])  // this task owns data for this slab
    {
      if (nxype[kp] > 0 )
      {
        // sum up all points in full process slab at each j
        vector<int> kprowsum(ny_,0);
        for (int jj=0; jj<ny_; jj++)
          for (int kk=kpkmin[kp]; kk<kpkmax[kp]; kk++)
            kprowsum[jj] += xrowcnt[kk*ny_+jj];
    
        vector<int> kpjmin(npey_);
        vector<int> kpjmax(npey_);
        int jset = 0;
        int jsum = 0;
        int rowsum = 0;
        double kprowavg = (double)nxype[kp]/(double)npey_;
        for (int jj=0; jj<ny_; jj++)
        {
          peyind[jj] = jset;
          rowsum += kprowsum[jj];
          jsum += kprowsum[jj];
          if (rowsum > kprowavg*(jset+1)) {
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
    
    int xdsize = npey_*nx_;
    vector<int> xdistsum(xdsize);  
    MPI_Allreduce(&xdist_loc[0], &xdistsum[0], xdsize, MPI_INT, MPI_SUM, comm_);
    
    if (kp >= pezind[klocmin] && kp <= pezind[klocmax])  // this task owns data for this slab
    {
      for (int jp = 0; jp < npey_; jp++)
      {
        double xavg = (double)kpjsum[jp]/(double)npex_;
        int iset = 0;
        int isum = 0;
        int blksum = 0;
        for (int ii=0; ii<nx_; ii++)
        {
          pexind[jp*nx_+ii] = iset;
          blksum += xdistsum[jp*nx_+ii];
          isum += xdistsum[jp*nx_+ii];
          if (blksum > xavg*(iset+1)) {
            isum = 0;
            iset++;
          }
        }
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

  // carry out communication to match computed distribution
  redistributeCells(cells);

  if (myRank_ == 0)
    cout << "Load histogram after initial distribution:" << endl;
  loadHistogram(cells);

}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::balanceLoop(vector<AnatomyCell>& cells)
{
  // call balanceLoop w. default values
  const int bblock = 5;
  const int bthresh = 10;
  const int maxiter = 100000000;
  balanceLoop(cells,bblock,bthresh,maxiter);
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::balanceLoop(vector<AnatomyCell>& cells, int bblock, int bthresh, int maxiter)
{
  bool balance = false;
  const int bprint = 10000;
  
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
      loadHistogram(nloc_);
        cout << endl << " *** Load balance achieved in " << bcnt << " iterations. ***" << endl << endl;
      }
    }

    bcnt++;
    if (bcnt%bprint == 0 && myRank_ == 0)
    {
      cout << "load balance iteration " << bcnt << ":  " << maxnum << " - " << minnum << " = " << maxnum-minnum << ", threshold = " << bthresh << endl;
      loadHistogram(nloc_);
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
  
  // use togive_ array to redistribute data
  bool allDataSent = false;
  int ccnt = 0;
  while (!allDataSent && ccnt < 1000)
  {
    allDataSent = diffuseDest(cells);
    ccnt++;
  }
  // calculate load histogram directly from cells
  if (myRank_ == 0)
    cout << "Load histogram after " << ccnt << " comm iterations:" << endl;
  loadHistogram(cells);
  
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
      xgid[i] = gpt.x + gpt.y*nx_ + gpt.z*nx_*ny_;
      ygid[i] = gpt.y + gpt.z*ny_ + gpt.x*ny_*nz_;
      zgid[i] = gpt.z + gpt.x*nz_ + gpt.y*nz_*nx_;
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

    int nMax = 2*nx_*ny_*nz_/nTasks_;      // 2 is a fudge factor for safety
    cells.resize(nMax);
    assignArray((unsigned char*)&(cells[0]), &nLocal, cells.capacity(),
                sizeof(AnatomyCell), &(dest[0]), 0, comm_);
    assert(nLocal <= nMax);
    cells.resize(nLocal);
  }
  else  // testing case:  nTasks < target process grid size
  {
    unsigned nLocal = cells.size();
    vector<unsigned> dest(nLocal);

    // compress process grid onto nTasks
    computeReducedProcGrid(nTasks_);
    if (myRank_ == 0)
      cout << "assignCellsToTask:  nTasks < process grid size, reduced grid = " <<
          tnx_ << " x " << tny_ << " x " << tnz_ << " used for testing." << endl;

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
void GDLoadBalancer::loadHistogram(vector<AnatomyCell>& cells)
{
  // compute load histogram from data in cells
  histloc_.resize(npegrid_,0);
  vector<int> mydata(npegrid_,0);
  for (unsigned ii=0; ii<cells.size(); ++ii)
    mydata[cells[ii].dest_]++;
  MPI_Allreduce(&mydata[0], &histloc_[0], npegrid_, MPI_INT, MPI_SUM, comm_);
  if (myRank_ == 0)
    computeHistogram();
    
}

////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::loadHistogram(vector<int>& nloc)
{
  // compute load histogram from data in nloc
  assert(nloc.size() == npegrid_);
  histloc_.resize(npegrid_,0);
  for (int ip=0; ip<npegrid_; ip++)
    histloc_[ip] = nloc[ip];
  if (myRank_ == 0)
    computeHistogram();
}

////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::computeHistogram()
{
  // compute histogram of current data distribution
  const int nhistmax = 100; // number of bins
  vector<int> phist(nhistmax,0);

  int maxnum = 0;
  int minnum = nx_*ny_*nz_;
  for (int p=0; p<npegrid_; p++)
  {
    if (histloc_[p] > maxnum)
      maxnum = histloc_[p];
    if (histloc_[p] < minnum)
      minnum = histloc_[p];
  }
  int nhist = maxnum - minnum + 1;
  if (nhist > nhistmax) nhist = nhistmax;
  int delta = (maxnum-minnum + 1)/nhist;
  if ((maxnum-minnum+1)%nhist !=0) delta++;
  for (int p=0; p<npegrid_; p++)
  {
    int bin = (histloc_[p]-minnum)/delta;
    phist[bin]++;
    //ewd DEBUG: print process number of top bin
    //if (bin == nhist-1 && phist[bin] == 1)
    //  cout << "loadHistogram: top bin pe " << p << ", nloc = " << histloc_[p] << endl;
  }
  cout << "load balance histogram:  " << endl;
  for (int i=0; i<nhist; i++)
    cout << "  " << minnum+delta*i << " - " << minnum+delta*(i+1) << ":    " << phist[i] << endl;

  nloctot_ = 0;
  for (unsigned ii=0; ii<npegrid_; ++ii)
    nloctot_ += histloc_[ii];
  double nlocavg_ = (double)nloctot_/(double)npegrid_; 
    
  cout << "total # of non-zero grid points = " << nloctot_ << ", avg. # per task = " << nlocavg_ << endl << endl;

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
