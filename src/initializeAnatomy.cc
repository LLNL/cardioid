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
