#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <map>
#if WITH_MPI
#include <mpi.h>
#endif

using std::cout;
using std::endl;
using std::string;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char** argv)
{

  cout << "Hello world." << endl;
  return 0;
}


//ewd:  testLoadBalancer is now obsolete -- I've replaced it with a Hello world test code to use
// as a model for adding independent test codes to the make system.  If we need it, I can go through and
// update the objects with the appropriate interfaces.

/*
#include "Control.h"
#include "Heart.h"
#include "TT04Model.h"
#include "SimpleInputParser.h"
#include "ProcessGrid3D.h"
#include "DataGrid3D.h"
#include "GDLoadBalancer.h"
#include "Timer.h"

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

  GDLoadBalancer* loadbal = new GDLoadBalancer(npex,npey,npez);

  /////
  //ewd: for now, just hack in a simple load
  /////
  if (mype == 0)
  {


    // use real heart data
    int nx,ny,nz,npts;
    ifstream isanat,isthet,isphi;
    if (realdata_)
    {
      tmap["simpleload"].start();
      string anatfile("heart_data.anatomy.dat");
      string thetafile("heart_data.theta.dat");
      string phifile("heart_data.phi.dat");
      isanat.open(anatfile.c_str(),ifstream::in);
      isthet.open(thetafile.c_str(),ifstream::in);
      isphi.open(phifile.c_str(),ifstream::in);
  
      // read header
      string tmp;
      isanat >> nx >> ny >> nz;
      getline(isanat,tmp);  // skip blank line
  
      npts = nx*ny*nz;
      assert(npts > 0);
      cout << anatfile << ":  header grid = " << nx << " x " << ny << " x " << nz << " (" << npts*sizeof(int) << " bytes)" << endl;
    }
    else
    {
      nx = ctrl.ngrid_ep[0];
      ny = ctrl.ngrid_ep[1];
      nz = ctrl.ngrid_ep[2];
      npts = nx*ny*nz;
    }
    
    vector<int> types(npts);
    const int readprint = 1000000;
    if (realdata_)
    {
      int ii = 0;
      int x,y,z;
      while (isanat.good()) {
        isanat >> x >> y >> z >> types[ii++];
        if (ii%readprint == 0)
          cout << "reading anatomy line " << ii << "..." << endl;
      }    
      tmap["simpleload"].stop();
    }
    else
    {
      tmap["randomload"].start();
      for (int i=0; i<npts; i++) {
        if (drand48() < 0.2)
          types[i] = 10;
        else
          types[i] = 0;
      }
      tmap["randomload"].stop();
    }

    tmap["assign_init"].start();
    loadbal->initialDistribution(types,nx,ny,nz);
    tmap["assign_init"].stop();

    tmap["balance"].start();
    loadbal->balanceLoop();
    tmap["balance"].stop();

    // check that every point was assigned to a pe and no point was assigned twice
    vector<int> mypt(npts);
    for (int i=0; i<npts; i++)
      mypt[i] = -1;

    vector<vector<int> > locgid = loadbal->locgid();
    for (int ip=0; ip<npex*npey*npez; ip++)
    {
      int ngidloc = locgid[ip].size();
      for (int i=0; i<ngidloc; i++)
      {
        int id = locgid[ip][i];
        if (mypt[id] != -1)
          cout << "ERROR: point " << id << " defined more than once, ip = " << ip << endl;
        mypt[id] = ip;
      }
    }
    for (int i=0; i<npts; i++)
      if (mypt[i] == -1 && types[i] > 0)
        cout << "ERROR: point " << i << " not assigned to any tasks!" << endl;
    
    cout.setf(ios::fixed,ios::floatfield);
    for ( map<std::string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
      double time = (*i).second.real();
      cout << "timing name=" << setw(15) << (*i).first << ":   time=" << setprecision(6) << setw(12) << time << " sec" << endl;
    }
    
  }
#if USE_MPI
  MPI_Finalize();
#endif
  delete loadbal;

  return 0;
}

*/
                                                                                                                
