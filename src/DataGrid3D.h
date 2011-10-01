#ifndef DATAGRID3D_H
#define DATAGRID3D_H

#include <vector>
#include "ProcessGrid3D.h"
using namespace std;

class DataGrid3D {
  private:

  ProcessGrid3D& pgrid_;
  int ng0_,ng1_,ng2_;
  int nb0_,nb1_,nb2_;
  int nloc_, nstate_;
  vector<vector<double> > locval_;   // ndata_ values at each local grid point
  vector<int> locgid_;               // global grid id of each local grid point
  vector<string> statename_;         // name of each nstate_ value
  bool active_;

  int initProcLoc(int gid, int &p0, int &p1, int &p2);
  int stateInd(string state);
  
  public:

  DataGrid3D(ProcessGrid3D& pgrid, int ng0, int ng1, int ng2);
  ~DataGrid3D();
  void addState(string name);
  int gridCoords(int gid, int &ig0, int &ig1, int &ig2);
  int initSendLocalData(string state, vector<double>& locdata, int g0);
  int initRecvLocalData(string state, int srcpe);
};
#endif
