#include "GridRouter.hh"
#include <iostream>
#include <cassert>
#include <cmath>
#include "CommTable.hh"
#include "Grid3DStencil.hh"
#include "DomainInfo.hh"

using namespace std;


GridRouter::GridRouter(vector<Long64>& gid, int nx, int ny, int nz, MPI_Comm comm)
: comm_(comm)
{
  // input:
  //   gid: vector of grid ids (gid = i+j*nx+k*nx*ny) owned by this task
  //   pectr:  cartesian coordinate of center of sphere enclosing all grid points
  //   perad:  radius of enclosing sphere
  // output:
  //   nSend:  number of tasks to send data to
  //   sendRank: vector of task ranks to send to
  //   sendOffset: vector of indices of local data for each send
  //       note: assumes packing such that sendRank[i] gets the elements
  //       sendOffset[i]:sendOffset[i+1]-1
    
  int npes, mype;
  MPI_Comm_size(comm_, &npes);
  MPI_Comm_rank(comm_, &mype);  

  double deltaR = 2.0;
  

  DomainInfo dInfo(gid, nx, ny, nz);
  vector<double> center = dInfo.center();
  double radius = dInfo.radius();
  
  // distribute process centers and radii to all tasks
  int ctrsize = 4*npes;
  vector<double> sumbuf(ctrsize);  // center coord + radius for each pe
  for (int i=0; i<ctrsize; i++)
    sumbuf[i] = 0.;
  sumbuf[4*mype+0] = center[0];
  sumbuf[4*mype+1] = center[1];
  sumbuf[4*mype+2] = center[2];
  sumbuf[4*mype+3] = radius;
  vector<double> pectr(ctrsize, 0.0);
  MPI_Allreduce(&sumbuf[0],&pectr[0],ctrsize,MPI_DOUBLE,MPI_SUM,comm_);

  // compute all tasks which overlap with this process
  vector<int> penbrs(0);
  for (int ip=0; ip<npes; ip++)
  {
     if (ip == mype)
	continue;
     double distsq = 0.0;
     for (int i=0; i<3; i++)
	distsq += (pectr[4*ip+i]-center[i])*(pectr[4*ip+i]-center[i]);
     assert(distsq > 0.0);
     double dist = sqrt(distsq);
     // this process is a potential neighbor, add it to the list
     if (dist <= radius + pectr[4*ip+3] + deltaR)
	penbrs.push_back(ip);
  }
  
  // send mpi rank, local data size to all neighbors
  int nnbrs = penbrs.size();
  int locsize = gid.size();
  MPI_Request sendReq[nnbrs];
  MPI_Request recvReq[nnbrs];
  vector<int> recvinfobuf(2*nnbrs);
  vector<int> sendinfobuf(2*nnbrs);
  for (int in=0; in<nnbrs; in++)
  {
    int recvpe = penbrs[in];
    MPI_Irecv(&recvinfobuf[2*in],2,MPI_INT,recvpe,1,comm_,recvReq+in);
  }  
  for (int in=0; in<nnbrs; in++)
  {
    int destpe = penbrs[in];
    sendinfobuf[2*in+0] = mype;
    sendinfobuf[2*in+1] = locsize;
    MPI_Isend(&sendinfobuf[2*in],2,MPI_INT,destpe,1,comm_,sendReq+in);
  }  
  MPI_Waitall(nnbrs, sendReq, MPI_STATUSES_IGNORE);
  MPI_Waitall(nnbrs, recvReq, MPI_STATUSES_IGNORE);

  // send local data to neighbors
  int nrecvtot = 0;
  for (int in=0; in<nnbrs; in++)
    nrecvtot += recvinfobuf[2*in+1];

  vector<Long64> recvdatabuf(nrecvtot);
  vector<int> recvoff(nnbrs);
  int offset = 0;
  for (int in=0; in<nnbrs; in++)
  {
    int recvpe = penbrs[in];
    int recvsize = recvinfobuf[2*in+1];
    MPI_Irecv(&recvdatabuf[offset],recvsize,MPI_LONG_LONG,recvpe,1,comm_,recvReq+in);
    recvoff[in] = offset;
    offset += recvsize;
  }  
  for (int in=0; in<nnbrs; in++)
  {
    int destpe = penbrs[in];
    MPI_Isend(&gid[0],locsize,MPI_LONG_LONG,destpe,1,comm_,sendReq+in);
  }  
  MPI_Waitall(nnbrs, sendReq, MPI_STATUSES_IGNORE);
  MPI_Waitall(nnbrs, recvReq, MPI_STATUSES_IGNORE);

  // recvdatabuf contains all gids of all possible neighbors, make accompanying
  // list of which pes own them
  vector<int> recvdatanbr(nrecvtot);
  for (int in=0; in<nnbrs; in++)
  {
    int insize = recvinfobuf[2*in+1];
    int inoff = recvoff[in];
    for (int i=0; i<insize; i++)
      recvdatanbr[inoff+i] = in;
  }      

  // use stencil to determine which gids to send/recv
  instencil_.resize(nnbrs);
  vector<int> sendthis(nnbrs);
  for (int ig=0; ig<locsize; ig++)
  {
    // neighbor process may own multiple grid points in stencil of ig
    // we only want to send it once, so just set send flag
    for (int in=0; in<nnbrs; in++)
      sendthis[in]=0;
    
    Grid3DStencil stencil(gid[ig], nx, ny, nz);
    for (int is=0; is<stencil.nStencil(); is++)
      for (int ib=0; ib<nrecvtot; ib++)
        if (stencil.nbrGid(is) == recvdatabuf[ib])
          sendthis[recvdatanbr[ib]] = 1;

    for (int in=0; in<nnbrs; in++)
      if (sendthis[in] == 1)
        instencil_[in].push_back(ig);
  }

  nSend_ = 0;
  int datasize = 0;
  sendRank_.clear();
  sendOffset_.clear();
  for (int in=0; in<instencil_.size(); in++)
  {
    if (instencil_[in].size() > 0)
    {
      sendRank_.push_back(recvinfobuf[2*in]);
      sendOffset_.push_back(datasize);
      datasize += instencil_[in].size();
      sendIndex_.resize(datasize);
      for (int is=0; is<instencil_[in].size(); is++)
        sendIndex_[sendOffset_[nSend_]+is] = instencil_[in][is];
      nSend_++;
    }
  }
}

CommTable GridRouter::commTable() const
{
   return CommTable(sendRank_, sendOffset_, comm_);
}


