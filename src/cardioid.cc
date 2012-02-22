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
#include <omp.h>

#include "Simulate.hh"
#include "PerformanceTimers.hh"
#include "mpiUtils.h"
#include "pio.h"
#include "units.h"

#include "initializeSimulate.hh"
#include "simulationLoop.hh"
#include "heap.h"
#include "object_cc.hh"
#include "Version.hh"

using namespace std;

namespace
{
   void parseCommandLineAndReadInputFile(int argc, char** argv, MPI_Comm comm);
   void printBanner();
}

/**
 * Units:
 *
 *  For Cardioid, the internal and external units are the same.  The
 *  seven fundamental units are:
 *  - length:             millimeter
 *  - mass:               microgram
 *  - time:               millisecond
 *  - current:            milliamp
 *  - temperature:        Kelvin
 *  - amount:             mol
 *  - luminous intensity: candella
 *
 *  This seven imply derived units as follows:
 *  - voltage:      millivolts
 *  - conductivity: Siemens
 *  - capacitance:  millifarad
 *  - charge:       microcoulomb
 *  
 */

// Sorry about this.
MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char** argv)
{
   int npes, mype;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);  

   // See units above.
   units_internal(1e-3, 1e-9, 1e-3, 1e-3, 1, 1, 1);
   units_external(1e-3, 1e-9, 1e-3, 1e-3, 1, 1, 1);
   
   profileStart("Total");
   profileSetPrintOrder("Total");
   profileSetPrintOrder("Assignment");
   profileSetPrintOrder("Loop");
   profileSetPrintOrder("");
   profileSetRefTimer("Loop");
   heap_start(100);

   if (mype == 0)
     printBanner();
   parseCommandLineAndReadInputFile(argc, argv, MPI_COMM_WORLD);
   
   Simulate sim;
   initializeSimulate("simulate", sim);

   //ewd:  turn on mpiP
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Pcontrol(1);

   profileStart("Loop");
   if ( !sim.parallelDiffusionReaction_) simulationLoop(sim);  
   else  simulationLoopParallelDiffusionReaction(sim);
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
void parseCommandLineAndReadInputFile(int argc, char** argv, MPI_Comm comm)
{
   int myRank;
   MPI_Comm_rank(comm, &myRank);      

   // get input file name from command line argument
   if (argc != 2 && argc != 1)
   {
      if (myRank == 0)
         cout << "Usage:  cardioid [input file]" << endl;
      exit(-1);
   }

   // parse input file
   string inputfile("object.data");
   if (argc > 1)
   {
      string argfile(argv[1]);
      inputfile = argfile;
   }

   if (myRank == 0)
   {
      object_compilefile(inputfile.c_str());
      printf("\nContents of object database:\n"
              "----------------------------------------------------------------------\n");
      object_print_all(stdout);
      printf("----------------------------------------------------------------------\n"
             "End of object database\n\n");
   }
   object_Bcast(0, MPI_COMM_WORLD);
}
}


namespace 
{
void printBanner()
{
   cout <<
      "Cardioid v1.0\n"
      "Unclassified/Code in development Distribution\n"
      "LLNL-CODE-508771\n"
      "Do not redistribute without permission\n\n";

   Version::getInstance().versionPrint(cout);
}
}
