#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <map>
#if USE_MPI
#include <mpi.h>
#endif
#include "GDLoadBalancer.h"
#include "Control.h"
#include "Timer.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
GDLoadBalancer::GDLoadBalancer(int npex, int npey, int npez):
    npex_(npex), npey_(npey), npez_(npez), nnbr_(6)
{
  // initialize arrays
  npegrid_ = npex_*npey_*npez_;
  nloc_.resize(npegrid_);
  for (int ip=0; ip<npegrid_; ip++)
    nloc_[ip] = 0;

  tosend_.resize(npegrid_);
  for (int ip=0; ip<npegrid_; ip++)
    tosend_[ip].resize(nnbr_);
  for (int ip=0; ip<npegrid_; ip++)
    for (int n=0; n<nnbr_; n++)
      tosend_[ip][n] = 0;
  
  loctype_.resize(npegrid_);
  locgid_.resize(npegrid_);
  for (int ip=0; ip<npegrid_; ip++)
  {
    loctype_[ip].resize(0);
    locgid_[ip].resize(0);
  }
    
  // store process numbers of all neighbors
  penbr_.resize(npegrid_);
  for (int ip=0; ip<npegrid_; ip++)
    penbr_[ip].resize(nnbr_);

  for (int ip=0; ip<npegrid_; ip++)
  {
    for (int j=0; j<nnbr_; j++)
      penbr_[ip][j] = -1;


    // want x,y,z coords of process ip
    const int npexy = npex_*npey_;
    int pz = ip/npexy;
    int xy0 = ip;
    while (xy0 >= npexy) xy0 -= npexy;
    int py = xy0/npex_;
    int px = xy0 - py*npex_;  
    
    int nbr;
    if (px > 0) {
      nbr = (px-1) + py*npex_ + pz*npex_*npey_;
      penbr_[ip][0] = nbr;
    }
    if (px < npex_-1) {
      nbr = (px+1) + py*npex_ + pz*npex_*npey_;
      penbr_[ip][1] = nbr;
    }
    if (py > 0) {
      nbr = px + (py-1)*npex_ + pz*npex_*npey_;
      penbr_[ip][2] = nbr;
    }
    if (py < npey_-1) {
      nbr = px + (py+1)*npex_ + pz*npex_*npey_;
      penbr_[ip][3] = nbr;
    }
    if (pz > 0) {
      nbr = px + py*npex_ + (pz-1)*npex_*npey_;
      penbr_[ip][4] = nbr;
    }
    if (pz < npez_-1) {
      nbr = px + py*npex_ + (pz+1)*npex_*npey_;
      penbr_[ip][5] = nbr;
    }
  }

  // equivalent neighbor index of neighboring point
  thatn_.resize(nnbr_);  
  thatn_[0] = 1; thatn_[1] = 0;
  thatn_[2] = 3; thatn_[3] = 2;
  thatn_[4] = 5; thatn_[5] = 4;
}
////////////////////////////////////////////////////////////////////////////////
GDLoadBalancer::~GDLoadBalancer()
{

}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::initialDistribution(vector<int>& types, int nx, int ny, int nz)
{
  int npts = types.size();
  if (npts != nx*ny*nz)
    cout << "WARNING:  GDLoadBalancer::initialDistribution called with npts != nx*ny*nz!" << endl;
  
  // count non-zero grid points
  ntissue_ = 0;
  for (int i=0; i<npts; i++)
    if (types[i] > 0)
      ntissue_++;

  // store list of which process owns each grid points
  //ewd: this is probably gratuitous, but good for debugging
  gpe_.resize(npts);
  for (int i=0; i<npts; i++)
    gpe_[i] = -1;

  // calculate total number of non-zero points at each x-y plane, x row
  int xypeavg = ntissue_/npez_;   // avg. # of pts in each x-y plane of process grid
  vector<int> nxyplane(nz);
  for (int k=0; k<nz; k++)
    nxyplane[k] = 0;
  vector<vector<int> > nxyrow;
  nxyrow.resize(nz);
  for (int k=0; k<nz; k++)
    nxyrow[k].resize(ny);
  for (int k=0; k<nz; k++)
    for (int j=0; j<ny; j++)
      nxyrow[k][j] = 0;

  for (int k=0; k<nz; k++)
    for (int j=0; j<ny; j++)
      for (int i=0; i<nx; i++)
      {
        int gid = i + j*nx + k*nx*ny;
        if (types[gid] > 0) {
          nxyplane[k]++;
          nxyrow[k][j]++;
        }
      }
  
    
  // distribute grid points on process grid with uniform distribution in z-direction
  vector<int> nxype(npez_);
  vector<int> kmax(npez_);
  int kset = 0;
  int xysum = 0;
  int ksum = 0;
  for (int k=0; k<nz; k++)
  {
    xysum += nxyplane[k];
    ksum += nxyplane[k];
    if (xysum > xypeavg*(kset+1)) {
      nxype[kset] = ksum;
      ksum = 0;
      kmax[kset] = k+1;
      kset++;
    }
  }
  kmax[npez_-1] = nz;

  int k0 = 0;
  vector<vector<int> > jmax;
  vector<vector<int> > nrowxpe;
  jmax.resize(npez_);
  nrowxpe.resize(npez_);
  for (int kp=0; kp<npez_; kp++)
  {
    jmax[kp].resize(npey_);
    nrowxpe[kp].resize(npey_);
  }
  
  for (int kp=0; kp<npez_; kp++)
  {
    int rowavg = nxype[kp]/npey_;  // avg. # of non-zero pts in each row of this plane
    vector<int> nrowk(ny);
    for (int j=0; j<ny; j++)
      nrowk[j] = 0;
    for (int j=0; j<ny; j++)
      for (int k=k0; k<kmax[kp]; k++)
        nrowk[j] += nxyrow[k][j];

    int jset = 0;
    int jsum = 0;
    int rowsum = 0;
    for (int j=0; j<ny; j++)
    {
      rowsum += nrowk[j];
      jsum += nrowk[j];
      if (rowsum > rowavg*(jset+1)) {
        nrowxpe[kp][jset] = jsum;
        jsum = 0;
        jmax[kp][jset] = j;
        jset++;
      }
    }
    jmax[kp][npey_-1] = ny-1;
    k0 = kmax[kp];
  }
    
  k0 = 0;    
  for (int kp=0; kp<npez_; kp++)
  {
    int j0 = 0;
    for (int jp=0; jp<npey_; jp++)
    {
      int bavg = nrowxpe[kp][jp]/npex_;
      int iset = 0;
      int isum = 0;
      for (int k=k0; k<kmax[kp]; k++)
      {
        for (int j=j0; j<jmax[kp][jp]; j++)
        {
          for (int i=0; i<nx; i++)
          {
            int gid = i + j*nx + k*nx*ny;
            int pex = iset;
            int pey = jp;
            int pez = kp;
            int peid = pex + pey*npex_ + pez*npex_*npey_;
            if (peid >= npegrid_ || peid < 0) {
              cout << "ERROR:  npegrid_ = " << npegrid_ << ", peid = " << peid << ", pex = " << pex << ", pey = " << pey << ", pez = " << pez << ", gid = " << gid << ", kset = " << kset << endl;
              exit(-1);
            }
            else {
              int t = types[gid];
              if (t > 0) {
                nloc_[peid]++;
                loctype_[peid].push_back(t);
                locgid_[peid].push_back(gid);
                
                if (gpe_[gid] != -1)
                {
                  cout << "ERROR:  grid point " << gid << " assigned more than once!" << ", peid = " << peid << endl;
                  exit(-1);
                }
                gpe_[gid] = peid;
                
                isum++;
                if (isum > bavg*(iset+1) && iset < (npex_-1))
                  iset++;
              }
            }
          }
        }
      }
      j0 = jmax[kp][jp];   // j0 = previous jmax
    }
    k0 = kmax[kp];   // k0 = previous kmax
  }

  nloctot_ = 0;
  for (int i=0; i<npegrid_; i++)
    nloctot_ += nloc_[i];
  nlocavg_ = (double)nloctot_/(double)npegrid_;

  loadHistogram();

  //ewd DEBUG
  // print out number of non-zero grid points at each xyplane
  cout << endl;
  cout << "X-Y data grid distribution (avg points per plane = " << ntissue_/nz << ", avg points per pe plane = " << ntissue_/npez_ << ")" << endl;
  for (int k=0; k<nz; k++)
  {
    cout << " x-y plane " << k << ", non-zero grid pts = " << nxyplane[k] << endl;
  }
  cout << endl;
  //ewd DEBUG
  
  //ewd DEBUG
  // calculate total number of non-zero points at each x-y plane of pe grid
  vector<int> xype_hist(npez_);
  for (int kp=0; kp<npez_; kp++)
    xype_hist[kp] = 0;
  for (int kp=0; kp<npez_; kp++)
    for (int jp=0; jp<npey_; jp++)
      for (int ip=0; ip<npex_; ip++)
      {
        int pe = ip + jp*npex_ + kp*npex_*npey_;
        xype_hist[kp] += nloc_[pe];
      }
  cout << endl;
  cout << "X-Y process grid distribution (avg points per plane = " << ntissue_/npez_ << ")" << endl;
  for (int kp=0; kp<npez_; kp++)
  {
    cout << " x-y pe plane " << kp << ", non-zero grid pts = " << xype_hist[kp] << endl;
  }
  cout << endl;
  //ewd DEBUG

}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::balanceLoop()
{
  // call balanceLoop w. default values
  const int bblock = 5;
  const int bthresh = 10;
  const int maxiter = 100000000;
  balanceLoop(bblock,bthresh,maxiter);
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::balanceLoop(int bblock, int bthresh, int maxiter)
{
  bool balance = false;
  const int bprint = 10000;
  bool rattle_meta = false;   // allow moves to procs w. nlocavg_ to break metastable pinning
  
  cout << "Starting load balance loop: threshold = " << bthresh << ", inner loop size = " << bblock << ", max iterations = " << maxiter << endl;
  
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
          bool allavg = true;
          for (int n=0; n<nnbr_; n++)
          {
            int nbr = penbr_[ip][n];
            if (nbr > -1 && nloc_[nbr] > nlocavg_ && nloc_[nbr] > maxval)
            {
              maxval= nloc_[nbr];
              nbrpmax = nbr;
              thisn = n;
            }
            if (nbr > -1 && nloc_[nbr] != nlocavg_)
              allavg = false;
          }
          // add grid point to ip, subtract from nbrpmax
          if (nbrpmax > -1)
          {
            nloc_[ip]++;
            tosend_[ip][thisn]--;  // we are owed one grid point from neighbor thisn
            nloc_[nbrpmax]--;
            tosend_[nbrpmax][thatn_[thisn]]++;  // we owe one grid point to neighbor thatn_
          }
          else if (rattle_meta && allavg)  // avoid metastable vacancy pinning
          {
            int tnbr = ((int)(nlocavg_-nloc_[ip]) >= nnbr_ ? nnbr_ : (int)(nlocavg_-nloc_[ip]));
            for (int n=0; n<tnbr; n++)
            {
              // take points from tnbr neighbors
              int tip = penbr_[ip][n];
              if (tip > -1) 
              {
                nloc_[ip]++;
                tosend_[ip][n]--;  // we are owed one grid point from neighbor n
                nloc_[tip]--;
                tosend_[tip][thatn_[n]]++;  // we owe one grid point to neighbor thatn_
              }
            }
          }
        }
        else if (nloc_[ip] > nlocavg_)
        {
          // loop over ip's neighbors, check if any are below average
          int minval = nloctot_;
          int nbrpmin = -1;
          int thisn = -1;
          bool allavg = true;
          for (int n=0; n<nnbr_; n++)
          {
            int nbr = penbr_[ip][n];
            if (nbr > -1 && nloc_[nbr] > nlocavg_ && nloc_[nbr] < minval)
            {
              minval= nloc_[nbr];
              nbrpmin = nbr;
              thisn = n;
            }
            if (nbr > -1 && nloc_[nbr] != nlocavg_)
              allavg = false;
          }
          // subtract grid point from ip, add to nbrpmin
          if (nbrpmin > -1)
          {
            nloc_[ip]--;
            tosend_[ip][thisn]++;  // we owe one grid point to neighbor thisn
            nloc_[nbrpmin]++;
            tosend_[nbrpmin][thatn_[thisn]]--;  // we are owed one grid point from neighbor thatn_
          }
          else if (rattle_meta && allavg)  // avoid metastable vacancy pinning
          {
            int tnbr = ((int)(nloc_[ip]-nlocavg_) >= nnbr_ ? nnbr_ : (int)(nloc_[ip]-nlocavg_));
            for (int n=0; n<tnbr; n++)
            {
              // give points to tnbr neighbors
              int tip = penbr_[ip][n];
              if (tip > -1) 
              {
                nloc_[ip]--;
                tosend_[ip][n]++;  // we owe one grid point to neighbor n
                nloc_[tip]++;
                tosend_[tip][thatn_[n]]--;  // we are owed one grid point from neighbor thatn_
              }
            }
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
      cout << endl << " *** Load balance achieved in " << bcnt << " iterations. ***" << endl << endl;
    }
    bcnt++;

    if (bcnt%bprint == 0)
    {
      cout << "load balance iteration " << bcnt << ":  " << maxnum << " - " << minnum << " = " << maxnum-minnum << ", threshold = " << bthresh << endl;
      loadHistogram();
    }
        

  }

  /*
  //ewd:  where should this take place?  ProcessGrid?
  // use tosend_ array to achieve calculated load balance, move gids between neighbors
  for (int ip=0; ip<npegrid_; ip++)
  {
    for (int n=0; n<nnbr_; n++)
    {
      if (tosend_[ip][n] > 0)      // send data
      {
        

      }
      else if (tosend_[ip][n] > 0) // receive data
      {
        
      }
    }
  }
  */
  
  if (!balance)
    cout << "balanceLoop did not converge after " << maxiter << " iterations." << endl;
  loadHistogram();
}
////////////////////////////////////////////////////////////////////////////////
void GDLoadBalancer::loadHistogram()
{
  // compute histogram of current data distribution
  const int nhist = 100; // number of bins
  vector<int> phist(nhist);
  {
    for (int i=0; i<nhist; i++)
      phist[i] = 0;
    int maxnum = 0;
    for (int p=0; p<npegrid_; p++)
      if (nloc_[p] > maxnum)
        maxnum = nloc_[p];
      
    int delta = maxnum/nhist;
    if (maxnum%nhist !=0) delta++;
    for (int p=0; p<npegrid_; p++)
    {
      int bin = nloc_[p]/delta;
      phist[bin]++;

      //ewd DEBUG: print process number of top bin
      if (bin == nhist-1 && phist[bin] == 1)
        cout << "loadHistogram: top bin pe " << p << ", nloc = " << nloc_[p] << endl;
      

    }
    cout << "load balance histogram:  " << endl;
    for (int i=0; i<nhist; i++)
      cout << "  " << delta*i << " - " << delta*(i+1) << ":    " << phist[i] << endl;

    cout << "total # of non-zero grid points = " << nloctot_ << ", avg. # per task = " << nlocavg_ << endl << endl;
  }
}
