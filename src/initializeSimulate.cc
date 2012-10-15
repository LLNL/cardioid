#include "initializeSimulate.hh"

#include <iostream>
#include <sstream>

#include "object_cc.hh"
#include "Simulate.hh"

#include "initializeAnatomy.hh"
#include "assignCellsToTasks.hh"
#include "diffusionFactory.hh"
#include "reactionFactory.hh"
#include "stimulusFactory.hh"
#include "sensorFactory.hh"
#include "getRemoteCells.hh"
#include "Anatomy.hh"
#include "mpiUtils.h"
#include "PerformanceTimers.hh"
#include "readSnapshotCellList.hh"
#include "hardwareInfo.h"
#include "pio.h"

using namespace std;


namespace
{
   /** Makes up a core list with one thread per core, except on BGQ
    *  where we have four threads per core. */
   void buildCoreList(unsigned& nCores, vector<unsigned>& cores)
   {
      int threadsPerCore = 1;
      #ifdef BGQ
      threadsPerCore = 4;
      #endif
      assert(cores.size() == 0);
      int nThreads = nCores*threadsPerCore;
      cores.resize(nThreads);
      for (unsigned ii=0; ii<nThreads; ++ii)
         cores[ii] = ii/threadsPerCore;
   }
}

Long64 findGlobalMinGid(const Anatomy& anatomy)
{
   Long64 minGid=anatomy.nx()*anatomy.ny()*anatomy.nz();
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
      if (anatomy.gid(ii) < minGid)
         minGid = anatomy.gid(ii); 
   Long64 globalMin;
   MPI_Allreduce(&minGid, &globalMin, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
   return globalMin;
}  

namespace
{
   void writeTorusMap(MPI_Comm comm, const string& filename)
   {
      if ( hi_hasTorus() == 0)
         return;

      int myRank;
      MPI_Comm_rank(comm, &myRank);

      int nDims = hi_nTorusDim();
      int coord[nDims];
      hi_torusCoords(coord);
      
      PFILE* pfile = Popen(filename.c_str(), "w", comm);
      PioSet(pfile, "ngroup", 1);
      if (myRank == 0)
      {
         int size[nDims];
         hi_torusSize(size);
         Pprintf(pfile, "# mpiRank coords[%d]\n", nDims);
         Pprintf(pfile, "# torus size:");
         for (unsigned ii=0; ii<nDims; ++ii)
            Pprintf(pfile, " %4d", size[ii]);
         Pprintf(pfile, "\n");
      }
      
      Pprintf(pfile, "%7d", myRank);
      for (unsigned ii=0; ii<nDims; ++ii)
         Pprintf(pfile, " %4d", coord[ii]);
      Pprintf(pfile, "\n");

      Pclose(pfile);      
   }
}


void initializeSimulate(const string& name, Simulate& sim)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   sim.name_ = name;

   OBJECT* obj = objectFind(name, "SIMULATE");

   objectGet(obj, "checkRanges", sim.checkRange_.on, "1");
   objectGet(obj, "VmMin", sim.checkRange_.vMin, "-110");
   objectGet(obj, "VmMax", sim.checkRange_.vMax, " 60");
   objectGet(obj, "loop", (int&)sim.loop_, "0"); // cast away volatile
   objectGet(obj, "maxLoop", sim.maxLoop_, "1000");
   objectGet(obj, "dt", sim.dt_, "0.01", "t");
   objectGet(obj, "time", sim.time_, "0", "t");
   objectGet(obj, "printRate", sim.printRate_, "20");
   objectGet(obj, "globalSyncRate", sim.globalSyncRate_, "-1");
   objectGet(obj, "printGid", sim.printGid_, "-1");
   objectGet(obj, "snapshotRate", sim.snapshotRate_, "100");
   objectGet(obj, "checkpointRate", sim.checkpointRate_, "-1");
   objectGet(obj, "nFiles", sim.nFiles_, "0");
   {
      int tmp; objectGet(obj, "profileAllCounters", tmp, "0");
      if (tmp == 1)
         profileSetVerbosity(true);
      else
         profileSetVerbosity(false);
   }
   {
      int tmp;         objectGet(obj, "writeTorusMap", tmp, "1");
      string filename; objectGet(obj, "torusMapFile", filename, "torusMap");
      if (tmp == 1)
         writeTorusMap(MPI_COMM_WORLD, filename);
   }
   
   timestampBarrier("initializing anatomy", MPI_COMM_WORLD);
   string nameTmp;
   objectGet(obj, "anatomy", nameTmp, "anatomy");
   initializeAnatomy(sim.anatomy_, nameTmp, MPI_COMM_WORLD);
   sim.nx_ = sim.anatomy_.nx();
   sim.ny_ = sim.anatomy_.ny();
   sim.nz_ = sim.anatomy_.nz();

   objectGet(obj, "stateFile", sim.stateFilename_);
   for (unsigned ii=0; ii<sim.stateFilename_.size(); ++ii)
   {
      string& name = sim.stateFilename_[ii];
      if (name[name.size()-1] != '#')
         name += "#";
   }

   // read either the loopType or parallelDiffusionReaction
   string tmp; objectGet(obj, "loopType", tmp, "omp");
   if (tmp == "omp")
       sim.loopType_ = Simulate::omp;
   else if (tmp == "pdr")
      sim.loopType_ = Simulate::pdr;
   else if (tmp == "allSkate")
      sim.loopType_ = Simulate::allSkate;
   else if (tmp == "lag")
      sim.loopType_ = Simulate::lag;
   else
      assert(false);
   if (object_testforkeyword(obj, "parallelDiffusionReaction"))
   {
      // backward compatibility
      int tmp; 
      objectGet(obj, "parallelDiffusionReaction", tmp, "0");
      if (tmp == 1)
         sim.loopType_ = Simulate::pdr;
      else
         sim.loopType_ = Simulate::omp;
   }
   timestampBarrier("assigning cells to tasks", MPI_COMM_WORLD);
   string decompositionName;
   objectGet(obj, "decomposition", decompositionName, "decomposition");
   LoadLevel loadLevel = assignCellsToTasks(sim, decompositionName, MPI_COMM_WORLD);

   stringstream stream;
   stream << loadLevel.nDiffusionCores;
   string defaultNDiffusionCores  = stream.str();

   unsigned nDiffusionCores;
   objectGet(obj, "nDiffusionCores", nDiffusionCores, defaultNDiffusionCores);

   vector<unsigned> diffusionCores;
   objectGet(obj, "diffusionThreads", diffusionCores);

   if (sim.loopType_ == Simulate::pdr || sim.loopType_ == Simulate::lag)
   {
      // diffusionThreads overrides nDiffusionCores, but when no thread
      // list is specified, we use 1 core.
      if (diffusionCores.size() == 0)
         buildCoreList(nDiffusionCores, diffusionCores);
      ThreadServer& threadServer = ThreadServer::getInstance();
      sim.diffusionThreads_ = threadServer.getThreadTeam(diffusionCores);
      if (getRank(0) == 0)
         cout << "Diffusion Threads: " << sim.diffusionThreads_ << endl;
      sim.reactionThreads_ = threadServer.getThreadTeam(vector<unsigned>());
      if (getRank(0) == 0)
         cout << "Reaction Threads: " << sim.reactionThreads_ << endl;
   }
   if (sim.loopType_ == Simulate::allSkate)
   {
      ThreadServer& threadServer = ThreadServer::getInstance();
      sim.reactionThreads_ = threadServer.getThreadTeam(vector<unsigned>());
      sim.diffusionThreads_ = sim.reactionThreads_;
      if (getRank(0) == 0)
      {
         cout << "Reaction Threads: " << sim.reactionThreads_ << endl;
         cout << "Diffusion Threads: " << sim.diffusionThreads_ << endl;
      }
   }
   
   
   timestampBarrier("building reaction object", MPI_COMM_WORLD);
   objectGet(obj, "reaction", nameTmp, "reaction");
   sim.reaction_ = reactionFactory(nameTmp, sim.anatomy_, sim.reactionThreads_);

   sim.printIndex_ = -1;
   // -2 -> print index 0 rank 0
   if (sim.printGid_ == -2 && myRank == 0)
   {
      sim.printGid_ = sim.anatomy_.gid(0);
   }
   // -1 -> global min gid
   if ( sim.printGid_ == -1 )    
      sim.printGid_ = findGlobalMinGid(sim.anatomy_);
   
   for (unsigned ii=0; ii<sim.anatomy_.nLocal(); ii++)
      if (sim.anatomy_.gid(ii) == sim.printGid_)
      {
         sim.printIndex_ = ii;
         break;
      }
   
   sim.printFile_=NULL; 
   if (sim.printIndex_ >=0)
   {
      sim.printFile_=fopen("data","a"); 
      printf(                "#   Loop     Time         gid            Vm(t)              dVm_r(t-h)             dVm_d(t-h)\n");
      fprintf(sim.printFile_,"#   Loop     Time         gid            Vm(t)              dVm_r(t-h)             dVm_d(t-h)\n");
   }
   
   
   getRemoteCells(sim, decompositionName, MPI_COMM_WORLD);

   timestampBarrier("building diffusion object", MPI_COMM_WORLD);
   objectGet(obj, "diffusion", nameTmp, "diffusion");
   sim.diffusion_ = diffusionFactory(nameTmp, sim.anatomy_, sim.diffusionThreads_,
                                     sim.reactionThreads_,
                                     sim.loopType_,loadLevel.stencil);
   
   timestampBarrier("building stimulus object", MPI_COMM_WORLD);
   vector<string> names;
   objectGet(obj, "stimulus", names);
   for (unsigned ii=0; ii<names.size(); ++ii)
      sim.stimulus_.push_back(stimulusFactory(names[ii], sim.anatomy_));

   timestampBarrier("building sensor object", MPI_COMM_WORLD);
   names.clear();
   objectGet(obj, "sensor", names);
   for (unsigned ii=0; ii<names.size(); ++ii)
      sim.sensor_.push_back(sensorFactory(names[ii], sim));

   // let user specify a filename containing the list of cell gids they want in the snapshots
   // (need to do this after data is distributed, so we can store just the local subset of points
   // on each task)
   string snapshotCLFile;
   objectGet(obj, "snapshotCellList", snapshotCLFile, "");
   if (snapshotCLFile != "")
      readSnapshotCellList(snapshotCLFile, sim,obj);
}
