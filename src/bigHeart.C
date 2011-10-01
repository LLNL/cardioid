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
#include "Control.h"
#include "Heart.h"
#include "TT04Model.h"
#include "SimpleInputParser.h"
#include "ProcessGrid3D.h"
#include "DataGrid3D.h"
#include "Timer.h"
#include "GDLoadBalancer.h"

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
  int MPI_COMM_WORLD = 0;
#endif
  
  // get input file name from command line argument
  if (argc != 2 && argc != 1)
  {
    if (mype == 0)
      cout << "Usage:  bigHeart [input file]" << endl;
    exit(-1);
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

    cout.setf(ios::fixed,ios::floatfield);
    for ( map<std::string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
      double time = (*i).second.real();
      cout << "timing name=" << setw(15) << (*i).first << ":   time=" << setprecision(6) << setw(12) << time << " sec" << endl;
    }
    
    // compute data decomposition, redistribute data until good load balance is achieved


    // move data to tasks:  placeholder grid objects constructed
  
    assert(ctrl.npegrid_ep.size() == 3);
    ProcessGrid3D* pgridep = new ProcessGrid3D(MPI_COMM_WORLD,ctrl.npegrid_ep[0],ctrl.npegrid_ep[1],ctrl.npegrid_ep[2]);

    cout << "mype = " << mype << ", grid.pe = " << pgridep->mype()
         << ", grid coord = " << pgridep->ip0() << " " << pgridep->ip1() << " "
         << pgridep->ip2() << ", active = " << pgridep->active() << endl;

    DataGrid3D* dgridep = new DataGrid3D(*pgridep,ctrl.ngrid_ep[0],ctrl.ngrid_ep[1],ctrl.ngrid_ep[2]);
  
    // compute comm tables


    
    // create integrator object, run time steps



    
  }
#if USE_MPI  
  MPI_Finalize();
#endif
  
  return 0;
}
