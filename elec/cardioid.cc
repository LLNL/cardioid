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
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include "Simulate.hh"
#include "PerformanceTimers.hh"
#include "mpiUtils.h"
#include "pio.h"
#include "units.h"
#include "ioUtils.h"

#include "initializeSimulate.hh"
#include "simulationLoop.hh"
#include "heap.h"
#include "object_cc.hh"
#include "Version.hh"

#ifdef HPM
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif

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
 *  - amount:             nanomol
 *  - luminous intensity: candella
 *
 *  These seven imply derived units as follows:
 *  - voltage:       millivolts
 *  - conductivity:  Siemens
 *  - capacitance:   millifarad
 *  - charge:        microcoulomb
 *  - concentration: millimolar
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

#ifdef USE_CUDA
   //FIXME!!! Have to figure out a better way to do binding!
   cudaSetDevice(mype % 4);
   // A ugly way to trigger the default CUDA context
   int *d_i;
   cudaMalloc(&d_i, sizeof(int));
   cudaFree(d_i);
#endif
   // See units above.
   units_internal(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);
   units_external(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);
   
   profileInit();
   profileStart("Total");
   // heap_start moved to initializeSimulate so that the size can be set
   // in the input deck.
//   heap_start(500);

   if (mype == 0)
     printBanner();
   parseCommandLineAndReadInputFile(argc, argv, MPI_COMM_WORLD);

   
   timestampBarrier("Starting initializeSimulate", MPI_COMM_WORLD);
   Simulate sim;
   initializeSimulate("simulate", sim);
   timestampBarrier("Finished initializeSimulate", MPI_COMM_WORLD);

   //ewd:  turn on mpiP
   //MPI_Barrier(MPI_COMM_WORLD);
   //MPI_Pcontrol(1);

#ifdef HPM
  HPM_Start("Loop"); 
#endif 
   timestampBarrier("Starting Simulation Loop", MPI_COMM_WORLD);
   profileStart_HW("Loop");
   switch (sim.loopType_)
   {
     case Simulate::omp:
      simulationLoop(sim);
      break;
     case Simulate::pdr:
     // printf("Cardioid pdr ptr=%p %p %p\n",sim.diffusion_->blockIndex(),sim.diffusion_->dVmBlock(),sim.diffusion_->VmBlock());  fflush(stdout); 
      simulationLoopParallelDiffusionReaction(sim);
      break;
     default:
      assert(false);
   }
   profileStop_HW("Loop");
   timestampBarrier("Finished Simulation Loop", MPI_COMM_WORLD);
#ifdef HPM
  HPM_Stop("Loop"); 
#endif 
   
   //ewd:  turn off mpiP
   //MPI_Pcontrol(0);

   profileStop("Total");
   //profileSetRefTimer("00:Loop");
   //profileSetPrintOrder("Total");
   //profileSetPrintOrder("Assignment");
   //profileSetPrintOrder("Loop");
   //profileSetPrintOrder("parallelDiffReac");
   //profileSetPrintOrder("DiffusionLoop");
   //profileSetPrintOrder("ReactionLoop");
   //profileSetPrintOrder("Dummy");
   //profileSetPrintOrder("HaloExchange");
   //profileSetPrintOrder("HaloExchMove2Buf");
   //profileSetPrintOrder("Integrator");
   //profileSetPrintOrder("ReactionWait");
   //profileSetPrintOrder("reactionL2Arrive");
   //profileSetPrintOrder("reactionL2Rest");
   //profileSetPrintOrder("Reaction");
   //profileSetPrintOrder("Reaction_nonGate");
   //profileSetPrintOrder("GateNonGateBarrier");
   //profileSetPrintOrder("Reaction_Gate");

   //profileSetPrintOrder("");
   if (mype == 0)
   {
      //profileDumpTimes(cout);
      //cout << "\n" << endl;
   }
   //profileDumpStats(cout);
   stringstream dirname;
   dirname << "snapshot."<<setfill('0')<<setw(12)<<sim.loop_;
   //profileDumpAll(dirname.str());
   heap_deallocate();
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
   for (int argvCursor=1; argvCursor<argc; argvCursor++)
   {
      if (argv[argvCursor] == "-h" ||
          argv[argvCursor] == "--help")
      {
         if (myRank == 0)
            cout << "Usage:  cardioid [input files]" << endl;
         exit(1);
      }
   }

   // parse input file
   vector<string> objectFiles;
   string restartFile("restart");
   if (argc == 1)
   {
      objectFiles.push_back("object.data");
   }
   else
   {
      for (int argvCursor=1; argvCursor<argc; argvCursor++)
      {
         objectFiles.push_back(argv[argvCursor]);
      }
   }

   if (myRank == 0)
   {
      bool inputDeckExists=true;
      for (int ifile=0; ifile<objectFiles.size(); ifile++)
      {
         if (filetest(objectFiles[ifile].c_str(), S_IFREG) != 0)
         {
            printf("objectfile=%s does not exist or wrong type\n", objectFiles[ifile].c_str());
            inputDeckExists=false;
         }
      }
      assert(inputDeckExists);

      string inputFiles = objectFiles[0];
      for (int ifile=1; ifile<objectFiles.size(); ifile++)
      {
         inputFiles += " " + objectFiles[ifile];
      }
      if (filetest(restartFile.c_str(), S_IFREG) == 0)
         inputFiles += " " + restartFile;

      object_set("files", inputFiles.c_str());
      
      object_compile();
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
      "Cardioid v2.0.0\n"
      "LLNL-CODE-764041\n"
      "\n";

   Version::getInstance().versionPrint(cout);
}
}
