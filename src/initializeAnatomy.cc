#include "initializeAnatomy.hh"

#include <iostream>
#include <cassert>

#include "Simulate.hh"
#include "AnatomyReader.hh"
#include "object_cc.hh"
#include "Anatomy.hh"

using namespace std;

namespace
{
   void readUsingPio(Simulate& sim, OBJECT* obj, MPI_Comm comm);
}

   

void initializeAnatomy(Simulate& sim, const string& name, MPI_Comm comm)
{
   OBJECT* obj = object_find(name.c_str(), "ANATOMY");

   string method;
   objectGet(obj, "method", method, "pio");

   if (method == "pio")
      readUsingPio(sim, obj, comm);
   else if (method == "simple")
      // We can wire in the simple load code that Erik originally wrote
      // here if we still need it.
      assert(1==0);
   else
      assert(1==0);
   double dx, dy, dz;
   objectGet(obj, "dx", dx, "0.2");
   objectGet(obj, "dy", dy, "0.2");
   objectGet(obj, "dz", dz, "0.2");
   sim.anatomy_.dx() = dx;
   sim.anatomy_.dy() = dy;
   sim.anatomy_.dz() = dz;
}


namespace {
   void readUsingPio(Simulate& sim, OBJECT* obj, MPI_Comm comm)
   {
      int myRank;
      MPI_Comm_rank(comm, &myRank);
      
      string fileName;
      objectGet(obj, "fileName", fileName, "snapshot.initial/atatomy#");
      
      if (myRank==0) cout << "Starting read" <<endl;

      AnatomyReader reader(fileName, comm, sim);
      if (myRank==0) cout << "Finished read" <<endl;
   }
}



// Here is the simple load code jut cut and pasted from main.  It would
// take about half an hour or so to get it working.  It needs to be
// modified to get its parameters from object.data and populate Simulate
// the Simulate::cells_ instead of types.
#if 0
void simpleLoad()
{
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
}
#endif
