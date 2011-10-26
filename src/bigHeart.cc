#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <map>
#include <mpi.h>

#include "Simulate.hh"
#include "Timer.hh"
#include "mpiUtils.h"

#include "initializeAnatomy.hh"
#include "assignCellsToTasks.hh"
#include "simulationLoop.hh"
#include "heap.h"
#include "object_cc.hh"
#include "diffusionFactory.hh"
#include "reactionFactory.hh"
#include "getRemoteCells.hh"
#include "Anatomy.hh"

using namespace std;

void parseCommandLineAndReadInputFile(int argc, char** argv, int rank);



// Sorry about this.
MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

int main(int argc, char** argv)
{
   int npes, mype;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);  
   heap_start(100);
   
   parseCommandLineAndReadInputFile(argc, argv, mype);
   
   Simulate sim;
   OBJECT* simObj = object_find("simulate", "SIMULATE");
   string nameTmp;
   
   objectGet(simObj, "anatomy", nameTmp, "anatomy");
   initializeAnatomy(sim, nameTmp, MPI_COMM_WORLD);
   
   objectGet(simObj, "decomposition", nameTmp, "decomposition");
   assignCellsToTasks(sim, nameTmp, MPI_COMM_WORLD);
   
   getRemoteCells(sim, MPI_COMM_WORLD);
   
   objectGet(simObj, "diffusion", nameTmp, "diffusion");
   sim.diffusion_ = diffusionFactory(nameTmp, sim.anatomy_);
   
   objectGet(simObj, "reaction", nameTmp, "reaction");
   sim.reaction_ = reactionFactory(nameTmp, sim.anatomy_);
   
   simulationLoop(sim);  
   
   if (mype == 0)
   {
      
      cout.setf(ios::fixed,ios::floatfield);
      for ( map<std::string,Timer>::iterator i = sim.tmap_.begin(); i != sim.tmap_.end(); i++ )
      {
	 double time = (*i).second.real();
	 cout << "timing name=" << setw(15) << (*i).first << ":   time=" << setprecision(6) << setw(12) << time << " sec" << endl;
      }
      
   }
   MPI_Finalize();
   
   return 0;
}


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
      cout << "argc = " << argc << endl;
   }

   object_compilefile(inputfile.c_str());
   

}
