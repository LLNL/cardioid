#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "ProcessGrid3D.h"
using namespace std;

#ifdef USE_MPI
#include <mpi.h> 
#else
typedef int MPI_Comm;
#endif

////////////////////////////////////////////////////////////////////////////////
ProcessGrid3D::ProcessGrid3D(MPI_Comm comm, int np0, int np1, int np2):
  np0_(np0), np1_(np1), np2_(np2)
{
  int npes, mype;
#if USE_MPI
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);
#else
  npes = 1;
  mype = 0;
#endif
  npes_ = np0_*np1_*np2_;
  if (npes == npes_)
  {
    mype_ = mype;
#if USE_MPI
    int ierr = MPI_Comm_dup(comm,&comm_);
    assert(ierr == 0);
#else
    comm_ = 0;
#endif
  }
  else if (npes_ < npes)
  {
    // comm_ = subcommunicator of comm
    vector<int> pmap_;
    pmap_.resize(npes_);
    for (int i=0; i<npes_; i++)
      pmap_[i] = i;
#if USE_MPI
    MPI_Group c_group, subgroup;
    MPI_Comm_group(comm,&c_group);
    MPI_Group_incl(c_group,npes_,&pmap_[0],&subgroup);
    MPI_Comm_create(comm,subgroup,&comm_);
    MPI_Group_free(&c_group);
    MPI_Group_free(&subgroup);
    if (mype < npes_)
      MPI_Comm_rank(comm_, &mype_);
    else
      mype_ = mype;
#else
    comm_ = 0;
#endif
  }
  else
  {
    if (mype == 0)
      cout << "ERROR in ProcessGrid3D:  process grid too large for communicator!" << endl;
    exit(1);
  }
  active_ = false;
  if (mype_ < npes_)
    active_ = true;
  
  // calculate coordinates of local process (x fastest)
  //   e.g. for an 8x4x2 process grid, pe 11 = (3,1,0), pe 52 = (4,2,1)
  if (active_)
  {
    const int np01 = np0_*np1_;
    ip2_ = mype_/np01;
    int xy0 = mype_;
    while (xy0 >= np01) xy0 -= np01;
    ip1_ = xy0/np0_;
    ip0_ = xy0 - ip1_*np0_;  
  }
  else
  {
    ip0_ = -1;
    ip1_ = -1;
    ip2_ = -1;
  }
}
////////////////////////////////////////////////////////////////////////////////
ProcessGrid3D::~ProcessGrid3D()
{
  //if (comm_ != MPI_COMM_NULL)
  //  MPI_Comm_free(&comm_);
}
////////////////////////////////////////////////////////////////////////////////
MPI_Comm ProcessGrid3D::subcomm(int nsubx, int nsuby, int nsubz, int istart, int jstart, int kstart)
{
  // carve out a subcommunicator containing a subset of the process grid
  int nprocs = nsubx*nsuby*nsubz;
  vector<int> pmap;
  pmap.resize(nprocs);
  // build pmap
  int p = 0;
  for (int k=kstart; k<nsubz+kstart; k++) {
    for (int j=jstart; j<nsuby+jstart; j++) {
      for (int i=istart; i<nsubx+istart; i++) {
        pmap[p] = i + j*np0_ + k*np0_*np1_;
        p++;
      }
    }
  }

  MPI_Comm sub;
#if USE_MPI
  MPI_Group c_group, subgroup;
  MPI_Comm_group(comm_,&c_group);
  MPI_Group_incl(c_group,nprocs,&pmap[0],&subgroup);
  MPI_Comm_create(comm_,subgroup,&sub);
  MPI_Group_free(&c_group);
  MPI_Group_free(&subgroup);
#else
  sub = 0;
#endif
  return sub;
}
////////////////////////////////////////////////////////////////////////////////
int ProcessGrid3D::dsend(int destpe, int dind, vector<double>& vec, vector<int>& id)
{
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
int ProcessGrid3D::drecv(int srcpe, int dind, vector<double>& vec, vector<int>& id)
{
  return 0;
}
