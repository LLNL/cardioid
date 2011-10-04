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
#include "AnatomyReader.h"
#include "mpiUtils.h"

Control parseCommandLineAndReadInputFile(int argc, char** argv, int rank);

using std::cout;
using std::endl;
using std::string;

// Sorry about this.
MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char** argv)
{

   
  const bool realdata_ = true;

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

  Control ctrl = parseCommandLineAndReadInputFile(argc, argv, mype);

//   AnatomyReader reader("anatomy#", MPI_COMM_WORLD);
//   Long64 nCells = reader._nx*reader._ny*reader._nz;
  
//   for (int ii=0; ii<reader._anatomy.size(); ++ii)
//      reader._anatomy[ii]._dest = 0;

//   unsigned nLocal = reader._anatomy.size();
  
//   vector<unsigned> dest(nLocal, 0);
  
//   reader._anatomy.resize(nCells);
//   assignArray((unsigned char*)&(reader._anatomy[0]),
// 	      &nLocal,
// 	      reader._anatomy.capacity(),
// 	      sizeof(AnatomyCell),
// 	      &(dest[0]),
// 	      0,
// 	      MPI_COMM_WORLD);
//   reader._anatomy.resize(nLocal);
//   cout << "task " << mype << " size " << reader._anatomy.size() <<endl;
//   for (unsigned ii=0; ii<reader._anatomy.size(); ++ii)
//      cout << reader._anatomy[ii]._cellType << " "
// 	  << reader._anatomy[ii]._theta << " "
// 	  << reader._anatomy[ii]._phi << " "
// 	  << endl;
  
//   exit(0);

  
  
  // compute process grid info 
  assert(ctrl.npegrid_ep.size() == 3);
  int npex = ctrl.npegrid_ep[0];
  int npey = ctrl.npegrid_ep[1];
  int npez = ctrl.npegrid_ep[2];

  ProcessGrid3D* pgridep = new ProcessGrid3D(MPI_COMM_WORLD,npex,npey,npez);
  DataGrid3D* dgridep = new DataGrid3D(*pgridep,ctrl.ngrid_ep[0],ctrl.ngrid_ep[1],ctrl.ngrid_ep[2]);
  
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

    // compute data decomposition, redistribute data until load balance is achieved
    tmap["assign_init"].start();
    loadbal->initialDistribution(types,nx,ny,nz);
    tmap["assign_init"].stop();

    tmap["balance"].start();
    loadbal->balanceLoop();
    tmap["balance"].stop();

    // move data to tasks
    
    
  
    // compute comm tables


    
    // create integrator object, run time steps



    

    cout.setf(ios::fixed,ios::floatfield);
    for ( map<std::string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
      double time = (*i).second.real();
      cout << "timing name=" << setw(15) << (*i).first << ":   time=" << setprecision(6) << setw(12) << time << " sec" << endl;
    }
    
  }
#if USE_MPI  
  MPI_Finalize();
#endif

  delete pgridep;
  delete dgridep;
  delete loadbal;
  
  return 0;
}


Control parseCommandLineAndReadInputFile(int argc, char** argv, int rank)
{
   // get input file name from command line argument
   if (argc != 2 && argc != 1)
   {
      if (rank == 0)
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
   SimpleInputParser sip;
   Control ctrl;
   if (sip.readInput(inputfile.c_str(), &ctrl))
   {
      if (rank == 0)
	 cout << "Input parsing error!" << endl;
      exit(1);
   }
   return ctrl;
}
