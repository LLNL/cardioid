#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <map>
#include <mpi.h>
#include <cstdlib>
#include <sstream>

#include "Simulate.hh"
#include "PerformanceTimers.hh"
#include "mpiUtils.h"
#include "pio.h"

#include "initializeSimulate.hh"
#include "simulationLoop.hh"
#include "heap.h"
#include "object_cc.hh"

using namespace std;

namespace
{
   void parseCommandLineAndReadInputFile(int argc, char** argv, int rank);
   void printBanner();
}



// Sorry about this.
MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char** argv)
{
   int npes, mype;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);  

   profileStart("Total");
   profileSetPrintOrder("Total");
   profileSetPrintOrder("Assignment");
   profileSetPrintOrder("Loop");
   profileSetPrintOrder("");
   profileSetRefTimer("Loop");

   heap_start(100);

   if (mype == 0)
     printBanner();
   parseCommandLineAndReadInputFile(argc, argv, mype);
   
   Simulate sim;
   initializeSimulate("simulate", sim);

   //ewd:  turn on mpiP
   MPI_Pcontrol(1);

   profileStart("Loop");
   simulationLoop(sim);  
   profileStop("Loop");
   
   //ewd:  turn off mpiP
   MPI_Pcontrol(0);

   profileStop("Total");
   if (mype == 0)
   {
      profileDumpTimes(cout);
      cout << "\n" << endl;
   }
   profileDumpStats(cout);
   
   stringstream dirname;
   dirname << "snapshot."<<setfill('0')<<setw(8)<<sim.loop_;
   profileDumpAll(dirname.str());
   MPI_Finalize();
   
   return 0;
}


namespace
{
void parseCommandLineAndReadInputFile(int argc, char** argv, int rank)
{
   // get input file name from command line argument
   if (argc != 2 && argc != 1)
   {
      if (rank == 0)
         cout << "Usage:  bigHeart [input file]" << endl;
      exit(-1);
   }

   // parse input file
   string inputfile("object.data");
   if (argc > 1)
   {
      string argfile(argv[1]);
      inputfile = argfile;
      //if (rank == 0)
      //cout << "argc = " << argc << endl;
   }

   object_compilefile(inputfile.c_str());
}
}


namespace 
{
void printBanner()
{
   cout <<
      "Cardioid v1.0\n"
      "Unclassified/Code in development Distribution\n"
      "LLNL-CODE-508771\n\n";
}
}
