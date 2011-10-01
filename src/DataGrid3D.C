#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "ProcessGrid3D.h"
#include "DataGrid3D.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
DataGrid3D::DataGrid3D(ProcessGrid3D& pgrid, int ng0, int ng1, int ng2):
    pgrid_(pgrid), ng0_(ng0), ng1_(ng1), ng2_(ng2)
{
  nloc_ = 0;
  nstate_ = 0;
  const int np0 = pgrid_.np0();
  const int np1 = pgrid_.np1();
  const int np2 = pgrid_.np2();
  assert(np0*np1*np2 != 0);
  nb0_ = (ng0_%np0 == 0 ? ng0_/np0 : ng0_/np0 + 1);
  nb1_ = (ng1_%np1 == 0 ? ng1_/np1 : ng1_/np1 + 1);
  nb2_ = (ng2_%np2 == 0 ? ng2_/np2 : ng2_/np2 + 1);
  assert(nb0_*nb1_*nb2_ != 0);

  //ewd DEBUG
  if (pgrid_.mype() == 0)
    cout << "DataGrid constructor:  grid = " << ng0 << " x " << ng1 << " x " << ng2 << ", nb0 = " << nb0_ << ", nb1 = " << nb1_ << ", nb2 = " << nb2_ << endl;


}
////////////////////////////////////////////////////////////////////////////////
DataGrid3D::~DataGrid3D()
{
}
////////////////////////////////////////////////////////////////////////////////
int DataGrid3D::gridCoords(int gid, int &ig0, int &ig1, int &ig2)
{
  const int ng01 = ng0_*ng1_;
  ig2 = gid/ng01;
  int xy0 = gid;
  while (xy0 >= ng01) xy0 -= ng01;
  ig1 = xy0/ng0_;
  ig0 = xy0 - ig1*ng0_;  
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
void DataGrid3D::addState(string name)
{
  statename_.push_back(name);
  nstate_++;
}
////////////////////////////////////////////////////////////////////////////////
int DataGrid3D::initSendLocalData(string state, vector<double>& locdata, int g0)
{
  // initLocalData distributes blocks of grid data using a simple parallel
  // distribution without regard to load balance
  //
  // Note: this routine is only called by the I/O tasks
  
  int dind = stateInd(state);
  if (dind == -1)
  {
    if (pgrid_.mype() == 0)
      cout << "ERROR in DataGrid3D::initSendLocalData, state " << state << " not found!" << endl;
    return 1;
  }

  // uniform grid distribution:  compute destination processor from global grid id
  vector<vector<double> > chunk;
  vector<vector<int> > chunkid;
  chunk.resize(pgrid_.npes());
  chunkid.resize(pgrid_.npes());
  int ip0,ip1,ip2;
  for (int ig=0; ig<locdata.size(); ig++)
  {
    int gid = g0+ig;
    int ierr = initProcLoc(gid,ip0,ip1,ip2);
    if (ierr > 0) return ierr;
    int destpe = pgrid_.pe(ip0,ip1,ip2);
    chunk[destpe].push_back(locdata[ig]);
    chunkid[destpe].push_back(gid);
  }

  // send data
  for (int ip=0; ip<pgrid_.npes(); ip++) 
    pgrid_.dsend(ip,dind,chunk[ip],chunkid[ip]);
  
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
int DataGrid3D::initRecvLocalData(string state, int srcpe)
{
  // call this from every task:  need to know which pes are sending data
  // (i.e. which tasks are reading input/checkpoint data from disk)

  int dind = stateInd(state);
  if (dind == -1)
  {
    if (pgrid_.mype() == 0)
      cout << "ERROR in DataGrid3D::initRecvLocalData, state " << state << " not found!" << endl;
    return 1;
  }

  vector<double> recvchunk;
  vector<int> recvid;
  pgrid_.drecv(srcpe,dind,recvchunk,recvid);
  for (int i=0; i<recvchunk.size(); i++)
    locval_[dind].push_back(recvchunk[i]);

  if (dind == 0) // assume all grid data is given in same order, only store gid once
    for (int i=0; i<recvid.size(); i++)
      locgid_.push_back(recvid[i]);

  //ewd DEBUG: this should always be true by definition, but let's check just in case
  assert(recvid.size() == recvchunk.size());

  return 0;
}
////////////////////////////////////////////////////////////////////////////////
int DataGrid3D::initProcLoc(int gid, int &p0, int &p1, int &p2)
{
  int ig0, ig1, ig2;
  int ierr = gridCoords(gid,ig0,ig1,ig2);
  p0 = ig0/nb0_;
  p1 = ig1/nb1_;
  p2 = ig2/nb2_;
  return ierr;
}
////////////////////////////////////////////////////////////////////////////////
int DataGrid3D::stateInd(string state)
{
  int dind = -1;
  for (int i=0; i<statename_.size(); i++)
    if (statename_[i] == state)
      dind = i;
  return dind;
}
////////////////////////////////////////////////////////////////////////////////
