#include "hardwareInfo.h"

#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#include "ddcMalloc.h"

#ifdef BGL
#include <bglpersonality.h>
#include <rts.h>
#endif

#ifdef BGP
#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#endif

typedef struct ioData_st
{
   MPI_Comm comm;
   int      nIoTasks;
   int*     ioTaskList;
} IO_DATA;

typedef struct hardwareInfo_st
{
   int pset;
   int rankInPset;
   int cpuId;
   int mpiRank;
} HARDWARE_INFO;



static IO_DATA* _ioData = NULL;
static int      _ioDataSize = 0;
static int      _ioDataCapacity = 0;
static int      _capacityIncrement = 10;

static int findComm(MPI_Comm comm);
static int mkIoData(MPI_Comm comm);
static int mkIoDataGeneric(MPI_Comm comm);
#ifdef BGL
static int mkIoDataBgl(MPI_Comm comm);
static void torusCoordsBgl(int* x, int* y, int* z, int* t);
static void torusSizeBgl(int* x, int* y, int* z);
#endif
#ifdef BGP
static int mkIoDataBgp(MPI_Comm comm);
static void torusCoordsBgp(int* x, int* y, int* z, int* t);
static void torusSizeBgp(int* x, int* y, int* z);
#endif
int hInfoSortFunction(const void* av, const void* bv); // changed from static
int intSortFunction(const void* av, const void* bv); // changed from static


int hi_nIoTasks(MPI_Comm comm)
{
   int ii = findComm(comm);
   return _ioData[ii].nIoTasks;
}

const int* hi_ioTaskList(MPI_Comm comm)
{
   int ii = findComm(comm);
   return _ioData[ii].ioTaskList;
}

int hi_hasTorus(void)
{
   #ifdef BGL
     return 1;
   #endif
   #ifdef BGP
     return 1;
   #endif
   
   return 0;
}

void hi_torusCoords(int* x, int* y, int* z, int* t)
{
   *x=0; *y=0; *z=0; *t=0;

   #ifdef BGL
   torusCoordsBgl(x, y, z, t);
   #endif

   #ifdef BGP
   torusCoordsBgp(x, y, z, t);
   #endif
   
   return;
}

void hi_torusSize(int* x, int* y, int* z)
{
   *x=0; *y=0; *z=0;

   #ifdef BGL
   torusSizeBgl(x, y, z);
   #endif

   #ifdef BGP
   torusSizeBgp(x, y, z);
   #endif
   
   return;
}


// "Private" functions

/** Returns the index in _ioData for the given comm.  This routine
 *  always succeeds since it creates an entry if none is already
 *  present. */
int findComm(MPI_Comm comm)
{
   for (int ii=0; ii<_ioDataSize; ++ii)
   {
      int result;
      MPI_Comm_compare(_ioData[ii].comm, comm, &result);
      if (result == MPI_IDENT)
	 return ii;
   }
   return mkIoData(comm);
}

/** Computes which tasks on the given comm are available as I/O tasks.
 * Returns the index in _ioData where the computed information is
 * stored.  Other than making sure there is enough storage, all of the
 * real work is delegated to a platform specific function that knows the
 * details of which tasks can optimally perform I/O. */
int mkIoData(MPI_Comm comm)
{
   if (_ioDataSize == _ioDataCapacity)
   {
      _ioDataCapacity += _capacityIncrement;
      _ioData = ddcRealloc(_ioData, _ioDataCapacity*sizeof(_ioData));
   }

#ifdef BGL
   return mkIoDataBgl(comm);
#endif
#ifdef BGP
   return mkIoDataBgp(comm);
#endif
   return mkIoDataGeneric(comm);
}

/** In the generic case all tasks are eligible to do I/O. */
int mkIoDataGeneric(MPI_Comm comm)
{
   assert(_ioDataCapacity > _ioDataSize+1);
   int commSize;
   MPI_Comm_size(comm, &commSize);
   _ioData[_ioDataSize].comm = comm;
   _ioData[_ioDataSize].nIoTasks = commSize;
   _ioData[_ioDataSize].ioTaskList = ddcMalloc(commSize*sizeof(int));
   for (int ii=0; ii<commSize; ++ii)
      _ioData[_ioDataSize].ioTaskList[ii] = ii;
   ++_ioDataSize;
   return _ioDataSize-1;
}

/** For BG/L we want only one task per pset doing I/O.  (All nodes in a
 * pset are connected to the same I/O node.)  For an arbitrary comm we
 * don't know how many psets will be represented or which ranks from a
 * given pset will be present.  So for each pset represented in comm, we
 * choose the task with the lowest rank in the pset and the lowest CPU
 * id. */
#ifdef BGL
int mkIoDataBgl(MPI_Comm comm)
{
   assert(_ioDataCapacity > _ioDataSize+1);
   int commSize;
   MPI_Comm_size(comm, &commSize);
   int myRank;
   MPI_Comm_rank(comm, &myRank);

   MPI_Datatype hInfoType;
   MPI_Type_contiguous(4, MPI_INT, &hInfoType);
   MPI_Type_commit(&hInfoType);
   
   BGLPersonality p;
   rts_get_personality(&p, sizeof(p));
   HARDWARE_INFO* hInfo = ddcMalloc(commSize * sizeof(HARDWARE_INFO));
   hInfo[myRank].pset = BGLPersonality_psetNum(&p);
   hInfo[myRank].rankInPset = BGLPersonality_rankInPset(&p);
   hInfo[myRank].cpuId = rts_get_processor_id();
   hInfo[myRank].mpiRank = myRank;

   HARDWARE_INFO myInfo = hInfo[myRank];
   MPI_Allgather(&myInfo, 1, hInfoType, hInfo, 1, hInfoType, comm);

   qsort(hInfo, commSize, sizeof(HARDWARE_INFO), hInfoSortFunction);
   unsigned nIoTasks = 1;
   for (unsigned ii=1; ii<commSize; ++ii)
      if (hInfo[ii].pset != hInfo[ii-1].pset)
	 ++nIoTasks;

   _ioData[_ioDataSize].comm = comm;
   _ioData[_ioDataSize].nIoTasks = nIoTasks;
   _ioData[_ioDataSize].ioTaskList = ddcMalloc(nIoTasks*sizeof(int));

   _ioData[_ioDataSize].ioTaskList[0] = hInfo[0].mpiRank;
   unsigned cnt = 1;
   for (unsigned ii=1; ii<commSize; ++ii)
      if (hInfo[ii].pset != hInfo[ii-1].pset)
	 _ioData[_ioDataSize].ioTaskList[cnt++] = hInfo[ii].mpiRank;
   assert(cnt == nIoTasks);
   qsort(_ioData[_ioDataSize].ioTaskList, nIoTasks, sizeof(int), intSortFunction);

   ++_ioDataSize;
   
   ddcFree(hInfo);
   MPI_Type_free(&hInfoType);
   return _ioDataSize-1;
}
#endif

#ifdef BGL
void torusCoordsBgl(int* x, int* y, int* z, int* t)
{
   BGLPersonality p;
   rts_get_personality(&p, sizeof(p));
   *x = p.xCoord;
   *y = p.yCoord;
   *z = p.zCoord;
   *t = rts_get_processor_id();
}
#endif

#ifdef BGL
void torusSizeBgl(int* x, int* y, int* z)
{
   BGLPersonality p;
   rts_get_personality(&p, sizeof(p));
   *x = p.xSize;
   *y = p.ySize;
   *z = p.zSize;
}
#endif

/** For BG/P it takes more than one task per I/O node to saturate the
 *  outbound connection.  Also, limiting ourselves to one task/file per
 *  I/O node too severely limits the number of files we can use to hit
 *  all OST on a large lustre files system.
 *
 *  All nodes in a pset are connected to the same I/O node.  For an
 *  arbitrary comm we don't know how many psets will be represented or
 *  which ranks from a given pset will be present.  So for each pset
 *  represented in comm, we choose a selection of the tasks on the pset
 *  attempting to choose all the tasks on different nodes in order to
 *  attempt to maximize available bandwidth on the internal network.  */
#ifdef BGP
int mkIoDataBgp(MPI_Comm comm)
{
   const unsigned ioTasksPerPset = 4;
   
   assert(_ioDataCapacity > _ioDataSize+1);
   int commSize;
   MPI_Comm_size(comm, &commSize);
   int myRank;
   MPI_Comm_rank(comm, &myRank);

   MPI_Datatype hInfoType;
   MPI_Type_contiguous(4, MPI_INT, &hInfoType);
   MPI_Type_commit(&hInfoType);
   
   _BGP_Personality_t p;
   Kernel_GetPersonality(&p, sizeof(p));
   HARDWARE_INFO* hInfo = ddcMalloc(commSize * sizeof(HARDWARE_INFO));
   hInfo[myRank].pset = p.Network_Config.PSetNum;
   hInfo[myRank].rankInPset = p.Network_Config.RankInPSet;
   hInfo[myRank].cpuId = Kernel_PhysicalProcessorID();
   hInfo[myRank].mpiRank = myRank;

   HARDWARE_INFO myInfo = hInfo[myRank];
   MPI_Allgather(&myInfo, 1, hInfoType, hInfo, 1, hInfoType, comm);

   qsort(hInfo, commSize, sizeof(HARDWARE_INFO), hInfoSortFunction);
   unsigned nPsets = 1;
   for (unsigned ii=1; ii<commSize; ++ii)
      if (hInfo[ii].pset != hInfo[ii-1].pset)
	 ++nPsets;

   _ioData[_ioDataSize].comm = comm;
   _ioData[_ioDataSize].ioTaskList = ddcMalloc(nPsets * ioTasksPerPset * sizeof(int));

   unsigned nIoTasks = 0;
   unsigned iBegin = 0;
   for (unsigned ii=0; ii<commSize; ++ii)
   {
      if ( ii==commSize-1 || hInfo[ii+1].pset != hInfo[iBegin].pset )
      {
	 unsigned psetSize = ii - iBegin + 1;
	 unsigned stride = psetSize/ioTasksPerPset;
	 if (psetSize%ioTasksPerPset != 0)
	    ++stride;
	 for (unsigned jj=iBegin; jj<ii+1; jj+=stride)
	 {
	    _ioData[_ioDataSize].ioTaskList[nIoTasks] = hInfo[jj].mpiRank;
	    ++nIoTasks;
	 }
	 iBegin = ii+1;
      }
   }
   _ioData[_ioDataSize].nIoTasks = nIoTasks;

   //ddt
/*    if (myRank == 0) */
/*    { */
/*       printf("I/O task info\n" */
/* 	     "-------------------------------------------\n" */
/* 	     "comm size = %d\n" */
/* 	     "nIoTasks  = %d\n" */
/* 	     "Task list:\n" */
/* //	      123456 12345678 123456 1234567890 123456 */
/* 	     " Index  mpiRank   pset rankInPset  cpuId\n" */
/* 	     , commSize, nIoTasks); */
/*       for (unsigned ii=0; ii<nIoTasks; ++ii) */
/* 	 for (unsigned jj=0; jj<commSize; ++jj) */
/* 	    if ( _ioData[_ioDataSize].ioTaskList[ii] == hInfo[jj].mpiRank ) */
/* 	       printf("%6d %8d %6d %10d %6d\n", */
/* 		      ii, hInfo[jj].mpiRank, hInfo[jj].pset, hInfo[jj].rankInPset, */
/* 		      hInfo[jj].cpuId); */
/*    } */
   //ddt end
   
   ++_ioDataSize;
   
   ddcFree(hInfo);
   MPI_Type_free(&hInfoType);
   return _ioDataSize-1;
}
#endif

#ifdef BGP
void torusCoordsBgp(int* x, int* y, int* z, int* t)
{
   _BGP_Personality_t p;
   Kernel_GetPersonality(&p, sizeof(p));
   *x = p.Network_Config.Xcoord;
   *y = p.Network_Config.Ycoord;
   *z = p.Network_Config.Zcoord;
   *t = Kernel_PhysicalProcessorID();
}
#endif

#ifdef BGP
void torusSizeBgp(int* x, int* y, int* z)
{
   _BGP_Personality_t p;
   Kernel_GetPersonality(&p, sizeof(p));
   *x = p.Network_Config.Xnodes;
   *y = p.Network_Config.Ynodes;
   *z = p.Network_Config.Znodes;
}
#endif

/** Sorts HARDWARE_INFO by pset then by rankInPset then by cpuId. */
int hInfoSortFunction(const void* av, const void* bv)
{
	const HARDWARE_INFO* a = (const HARDWARE_INFO*)av;
	const HARDWARE_INFO* b = (const HARDWARE_INFO*)bv;
	if (a->pset > b->pset) return 1;
	if (a->pset < b->pset) return -1;

	if (a->rankInPset > b->rankInPset) return 1;
	if (a->rankInPset < b->rankInPset) return -1;

 	if (a->cpuId > b->cpuId) return 1;
	if (a->cpuId < b->cpuId) return -1;

	return 0;
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
