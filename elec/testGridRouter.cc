#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#if USE_MPI
#include <mpi.h>
#endif
#include "Control.h"
#include "Heart.h"
#include "TT04Model.h"
#include "SimpleInputParser.h"
#include "ProcessGrid3D.h"
#include "DataGrid3D.h"
#include "GDLoadBalancer.h"
#include "Timer.h"
#include "GridPoint.h"
#include "GridRouter.h"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char** argv)
{

  const bool realdata_ = true;

  Control ctrl;
  Heart heart;
  map<std::string,Timer> tmap;
  
  int npes, mype;
#if USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);  
#else
  npes = 1;
  mype = 0;
#endif
  
  // get input file name from command line argument
  if (argc != 2 && argc != 1)
  {
    if (mype == 0)
      cout << "Usage:  bigHeart [input file]" << endl;
    exit(1);
  }

  // parse input file
  string inputfile("input");
  if (argc > 1)
  {
    string argfile(argv[1]);
    inputfile = argfile;
    cout << "argc = " << argc << endl;
  }
  SimpleInputParser* sip = new SimpleInputParser();
  if (sip->readInput((char*)inputfile.c_str(),&ctrl))
  {
    if (mype == 0)
      cout << "Input parsing error!" << endl;
    return 1;
  }
  delete sip;

  // compute process grid info 
  int npex = ctrl.npegrid_ep[0];
  int npey = ctrl.npegrid_ep[1];
  int npez = ctrl.npegrid_ep[2];

  assert(npex*npey*npez == npes);
  
  int nx,ny,nz,npts;
  nx = ctrl.ngrid_ep[0];
  ny = ctrl.ngrid_ep[1];
  nz = ctrl.ngrid_ep[2];
  npts = nx*ny*nz;
    
  vector<int> types(npts);
  GDLoadBalancer* loadbal;

  if (mype == 0)
  {
    types.resize(npts);
    tmap["randomload"].start();
    for (int i=0; i<npts; i++) {
      if (drand48() < 0.2)
        types[i] = 10;
      else
        types[i] = 0;
    }
    tmap["randomload"].stop();

    loadbal = new GDLoadBalancer(npex,npey,npez);
    loadbal->initialDistribution(types,nx,ny,nz);
    loadbal->balanceLoop();
  }

  // send grid data to pes
  vector<int> locgid;
  if (mype == 0)
    locgid = loadbal->locgid(0);
  for (int ip=1; ip<npes; ip++)
  {
    if (mype == 0)
    {
      vector<int> gid = loadbal->locgid(ip);
      int sendsize = gid.size();
      MPI_Send(&sendsize,1,MPI_INT,ip,ip,MPI_COMM_WORLD);
      MPI_Send(&gid[0],sendsize,MPI_INT,ip,ip,MPI_COMM_WORLD);
    }
    if (mype == ip)
    {
      // receive local data
      int locsize;
      MPI_Recv(&locsize,1,MPI_INT,0,mype,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      locgid.resize(locsize);
      MPI_Recv(&locgid[0],locsize,MPI_INT,0,mype,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
  }

  // calculate center, radius of process domain
  GridPoint ip0(locgid[0],nx,ny,nz);
  int xmin = ip0.x;
  int xmax = ip0.x;
  int ymin = ip0.y;
  int ymax = ip0.y;
  int zmin = ip0.z;
  int zmax = ip0.z;
  for (int i=0; i<locgid.size(); i++)
  {
    GridPoint ipt(locgid[i],nx,ny,nz);
    if (ipt.x < xmin) xmin = ipt.x;
    if (ipt.x > xmax) xmax = ipt.x;
    if (ipt.y < ymin) ymin = ipt.y;
    if (ipt.y > ymax) ymax = ipt.y;
    if (ipt.z < zmin) zmin = ipt.z;
    if (ipt.z > zmax) zmax = ipt.z;
  }

  vector<double> center(3);
  center[0] = 0.5*(double)(xmin+xmax+1);
  center[1] = 0.5*(double)(ymin+ymax+1);
  center[2] = 0.5*(double)(zmin+zmax+1);

  double radsq = 0.0;
  for (int i=0; i<locgid.size(); i++)
  {
    GridPoint ipt(locgid[i],nx,ny,nz);
    double dsq = (center[0]-(double)ipt.x)*(center[0]-(double)ipt.x)
        + (center[1]-(double)ipt.y)*(center[1]-(double)ipt.y)
        + (center[2]-(double)ipt.z)*(center[2]-(double)ipt.z);
    if (dsq > radsq) radsq = dsq;
    dsq = (center[0]-(double)ipt.x-1.)*(center[0]-(double)ipt.x-1.)
        + (center[1]-(double)ipt.y)*(center[1]-(double)ipt.y)
        + (center[2]-(double)ipt.z)*(center[2]-(double)ipt.z);
    if (dsq > radsq) radsq = dsq;
    dsq = (center[0]-(double)ipt.x)*(center[0]-(double)ipt.x)
        + (center[1]-(double)ipt.y-1.)*(center[1]-(double)ipt.y-1.)
        + (center[2]-(double)ipt.z)*(center[2]-(double)ipt.z);
    if (dsq > radsq) radsq = dsq;
    dsq = (center[0]-(double)ipt.x)*(center[0]-(double)ipt.x)
        + (center[1]-(double)ipt.y)*(center[1]-(double)ipt.y)
        + (center[2]-(double)ipt.z-1.)*(center[2]-(double)ipt.z-1.);
    if (dsq > radsq) radsq = dsq;
    dsq = (center[0]-(double)ipt.x-1.)*(center[0]-(double)ipt.x-1.)
        + (center[1]-(double)ipt.y-1.)*(center[1]-(double)ipt.y-1.)
        + (center[2]-(double)ipt.z)*(center[2]-(double)ipt.z);
    if (dsq > radsq) radsq = dsq;
    dsq = (center[0]-(double)ipt.x-1.)*(center[0]-(double)ipt.x-1.)
        + (center[1]-(double)ipt.y)*(center[1]-(double)ipt.y)
        + (center[2]-(double)ipt.z-1.)*(center[2]-(double)ipt.z-1.);
    if (dsq > radsq) radsq = dsq;
    dsq = (center[0]-(double)ipt.x)*(center[0]-(double)ipt.x)
        + (center[1]-(double)ipt.y-1.)*(center[1]-(double)ipt.y-1.)
        + (center[2]-(double)ipt.z-1.)*(center[2]-(double)ipt.z-1.);
    if (dsq > radsq) radsq = dsq;
    dsq = (center[0]-(double)ipt.x-1.)*(center[0]-(double)ipt.x-1.)
        + (center[1]-(double)ipt.y-1.)*(center[1]-(double)ipt.y-1.)
        + (center[2]-(double)ipt.z-1.)*(center[2]-(double)ipt.z-1.);
    if (dsq > radsq) radsq = dsq;
  }
  double radius = sqrt(radsq);
  
  if (mype == 0)
    tmap["grid_router"].start();
  GridRouter(locgid,nx,ny,nz,center,radius,MPI_COMM_WORLD);
  if (mype == 0)
    tmap["grid_router"].stop();


  
  if (mype == 0)
  {
    cout.setf(ios::fixed,ios::floatfield);
    for ( map<std::string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
      double time = (*i).second.real();
      cout << "timing name=" << setw(15) << (*i).first << ":   time=" << setprecision(6) << setw(12) << time << " sec" << endl;
    }
  }

#if USE_MPI
  MPI_Finalize();
#endif



  return 0;
}
