#ifndef SIMULATE_HH
#define SIMULATE_HH

#include <map>
#include <string>
#include <vector>

#include "Timer.hh"
#include "AnatomyCell.hh"
#include "LocalGrid.hh"
#include "Array3d.hh"

class Diffusion;

// Kitchen sink class for heart simulation.  This is probably a great
// example of how *not* to do design, but we need to get something
// running before we can understand what good design would be.

// This class intentionally exposes all the members as public items.
// This is poor encapsulation, but it seems silly at this point to make
// get/set calls for all of the externals while we are still figuring
// out what we are doing.

class Simulate
{
 public:

   std::map<std::string, Timer> tmap_;

   int loop_;
   int maxLoop_;
   
   int nx_, ny_, nz_; // global size of grid


   unsigned nCellLocal_;
   unsigned nCellRemote_;
   std::vector<AnatomyCell> cell_;
   std::vector<double> VmArray_;

   Diffusion* diffusion_;
   
};

#endif
