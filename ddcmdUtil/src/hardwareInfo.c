#include "hardwareInfo.h"

#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <inttypes.h>

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

#ifdef BGQ
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>
#endif

typedef struct ioData_st
{
   MPI_Comm comm;
   int      nIoTasks;
   int*     ioTaskList;
} IO_DATA;

// This structure has a corresponding MPI_Datatype.
typedef struct hardwareInfo_st
{
   uint64_t pset;
   uint64_t rankInPset;
   int cpuId;
   int mpiRank;
} HARDWARE_INFO;

MPI_Datatype hi_MPIType(void)
{
   static int initialized = 0;
   static MPI_Datatype hiType;
   if (initialized == 0)
   {
      HARDWARE_INFO hi;
      int n=2;
      int blkcnt[n];
      MPI_Aint disp[n];
      MPI_Datatype types[n];
      blkcnt[0] = 2;
      blkcnt[1] = 2;
      MPI_Address(&hi.pset, &disp[0]);
      MPI_Address(&hi.cpuId, &disp[1]);
      types[0] = MPI_LONG_LONG;
      types[1] = MPI_INT;
      for (int i = n-1; i >= 0; i--)
         disp[i] -= disp[0];
      MPI_Type_struct(n, blkcnt, disp, types, &hiType);
      MPI_Type_commit(&hiType);
      initialized = 1;
   }
   return hiType;
}



static IO_DATA* _ioData = NULL;
static int      _ioDataSize = 0;
static int      _ioDataCapacity = 0;
static int      _capacityIncrement = 10;

static int findComm(MPI_Comm comm);
static int mkIoData(MPI_Comm comm);
static int mkIoDataGeneric(MPI_Comm comm);
#ifdef BGL
static int mkIoDataBgl(MPI_Comm comm);
static void torusCoordsBgl(int* coord);
static void torusSizeBgl(int* size);
static int coreIdBgl(void);
#endif
#ifdef BGP
static int mkIoDataBgp(MPI_Comm comm);
static void torusCoordsBgp(int* coord);
static void torusSizeBgp(int* size);
static int coreIdBgp(void);
#endif
#ifdef BGQ
static int mkIoDataBgq(MPI_Comm comm);
static void torusCoordsBgq(int* coord);
static void torusSizeBgq(int* size);
static int coreIdBgq(void);
static int hw_info_to_io_task_list(int pid,int np,MPI_Comm comm,
											  HARDWARE_INFO self,int ioTasksPerPset,
											  int **ioTaskList_ptr);
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
#ifdef BGQ
   return 1;
#endif
   
   return 0;
}

int hi_nTorusDim(void)
{
#ifdef BGL
   return 3;
#endif
#ifdef BGP
   return 3;
#endif
#ifdef BGQ
   return 5;
#endif
   
   return 0;
}

void hi_torusCoords(int* coord)
{
#ifdef BGL
   torusCoordsBgl(coord);
#endif

#ifdef BGP
   torusCoordsBgp(coord);
#endif
   
#ifdef BGQ
   torusCoordsBgq(coord);
#endif
   
   return;
}

void hi_torusSize(int* size)
{
#ifdef BGL
   torusSizeBgl(size);
#endif

#ifdef BGP
   torusSizeBgp(size);
#endif

#ifdef BGQ
   torusSizeBgq(size);
#endif

   return;
}

int hi_coreId(void)
{
#ifdef BGL
   return coreIdBgl();
#endif
#ifdef BGP
   return coreIdBgp();
#endif
#ifdef BGQ
   return coreIdBgq();
#endif
   return 0;
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
#ifdef BGQ
   return mkIoDataBgq(comm);
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

   MPI_Datatype hInfoType = hi_MPIType();
   
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
void torusCoordsBgl(int* coord)
{
   BGLPersonality p;
   rts_get_personality(&p, sizeof(p));
   coord[0] = p.xCoord;
   coord[1] = p.yCoord;
   coord[2] = p.zCoord;
}
#endif

#ifdef BGL
void torusSizeBgl(int* size)
{
   BGLPersonality p;
   rts_get_personality(&p, sizeof(p));
   size[0] = p.xSize;
   size[1] = p.ySize;
   size[2] = p.zSize;
}
#endif

#ifdef BGL
int coreIdBgl(void)
{
   return rts_get_processor_id();
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

   MPI_Datatype hInfoType = hi_MPIType();
   
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
void torusCoordsBgp(int* coord)
{
   _BGP_Personality_t p;
   Kernel_GetPersonality(&p, sizeof(p));
   coord[0] = p.Network_Config.Xcoord;
   coord[1] = p.Network_Config.Ycoord;
   coord[2] = p.Network_Config.Zcoord;
}
#endif

#ifdef BGP
void torusSizeBgp(int* size)
{
   _BGP_Personality_t p;
   Kernel_GetPersonality(&p, sizeof(p));
   size[0] = p.Network_Config.Xnodes;
   size[1] = p.Network_Config.Ynodes;
   size[2] = p.Network_Config.Znodes;
}
#endif

#ifdef BGP
int coreIdBgp(void)
{
   return Kernel_PhysicalProcessorID();
}
#endif

#ifdef BGQ
uint64_t getIdFromCoord(int coord[5])
{
   uint64_t tmp = 0;
   for (unsigned ii=0; ii<5; ++ii)
      tmp = (tmp<<8) + coord[ii];
   tmp = tmp<<1;
   tmp += 1;
   return tmp;
}
#endif
#ifdef BGQ
uint64_t getMyBridgeNodeId(Personality_t p)
{
   int coord[5];
   coord[0] = p.Network_Config.cnBridge_A;
   coord[1] = p.Network_Config.cnBridge_B;
   coord[2] = p.Network_Config.cnBridge_C;
   coord[3] = p.Network_Config.cnBridge_D;
   coord[4] = p.Network_Config.cnBridge_E;
   return getIdFromCoord(coord);
}
#endif
#ifdef BGQ
uint64_t getCoordId(Personality_t p)
{
   int coord[5];
   torusCoordsBgq(coord);
   return getIdFromCoord(coord);
}
#endif

/** BG/Q doesn't actually have psets, but to make this code look as much
 *  as possible like the code for L and P, we can use a similar
 *  concept.  Each BG/Q compute node is connected to a bridge node
 *  (which is itself a compute node) that sends data to the actual I/O
 *  node.  The mapping between compute nodes and bridge nodes and
 *  between bridge nodes and I/O nodes is set up at boot time and is not
 *  dynamic during the run.  Hence we can treat the set of all nodes
 *  with the same bridge node as a pset and the id of the processor as
 *  the rank in the pset.
 *
 *  We will use one little trick however:  Instead of copying IBM's
 *  method of turning coords into an id directly, we're going to add one
 *  to the id to make the ids start at 1 instead of zero.  Then we will
 *  make the rank of any bridge node 0. That will make the bridge node
 *  always have the lowest rank in the pset.
 */
#ifdef BGQ
int mkIoDataBgq(MPI_Comm comm)
{
   const unsigned ioTasksPerPset = 2;
   
   assert(_ioDataCapacity > _ioDataSize+1);
   int commSize;
   MPI_Comm_size(comm, &commSize);
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   
   Personality_t p;
   int rc = Kernel_GetPersonality(&p, sizeof(p));
   assert(rc == 0);

   HARDWARE_INFO myInfo;
   myInfo.pset = getMyBridgeNodeId(p);
   myInfo.rankInPset = getCoordId(p);
   if (myInfo.pset == myInfo.rankInPset)
	  myInfo.rankInPset = 0;
   myInfo.cpuId = Kernel_ProcessorID(); // 0-63
   myInfo.mpiRank = myRank;

	{
	  int n_io_tasks,*io_task_list;
	  n_io_tasks =
		 hw_info_to_io_task_list(myRank,commSize,comm,
										 myInfo,ioTasksPerPset,&io_task_list);
	  _ioData[_ioDataSize].comm = comm;
	  _ioData[_ioDataSize].nIoTasks = n_io_tasks;
	  _ioData[_ioDataSize].ioTaskList = io_task_list;
	}

   ++_ioDataSize;   
   return _ioDataSize-1;
}
#endif

#ifdef BGQ
void torusCoordsBgq(int* coord)
{
   Personality_t pers;
   int rc = Kernel_GetPersonality(&pers, sizeof(pers));
   assert(rc == 0);
   coord[0] = pers.Network_Config.Acoord;
   coord[1] = pers.Network_Config.Bcoord;
   coord[2] = pers.Network_Config.Ccoord;
   coord[3] = pers.Network_Config.Dcoord;
   coord[4] = pers.Network_Config.Ecoord;
}
#endif

#ifdef BGQ
void torusSizeBgq(int* size)
{
   Personality_t pers;
   int rc = Kernel_GetPersonality(&pers, sizeof(pers));
   assert(rc == 0);
   size[0] = pers.Network_Config.Anodes;
   size[1] = pers.Network_Config.Bnodes;
   size[2] = pers.Network_Config.Cnodes;
   size[3] = pers.Network_Config.Dnodes;
   size[4] = pers.Network_Config.Enodes;
}
#endif

#ifdef BGQ
int coreIdBgq(void)
{
   int coreId = Kernel_ProcessorCoreID();   // 0-15   
   int procID = Kernel_ProcessorID();       // 0-63
   return coreId;
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

#ifdef BGQ

/*
  Sort the HARDWARE_INFO data (self argument) among processors
  in comm, in parallel, and figure out which tasks belong to
  which psets, and then determine which tasks should be IO
  tasks. The number of IO tasks is returned. The last argument
  is set to an array of the ranks of the IO tasks on return.
 */
static int hw_info_to_io_task_list(int pid,int np,MPI_Comm comm,
											  HARDWARE_INFO self,int ioTasksPerPset,
											  int **ioTaskList_ptr)
{
  struct
  {
    uint64_t bitbucket[(sizeof(HARDWARE_INFO)+sizeof(uint64_t)-1)/sizeof(uint64_t)];
    HARDWARE_INFO info;
  } item,tmpitem;
  const int nfields = sizeof(item.bitbucket)/sizeof(*item.bitbucket);
  
  item.info = self;
  
  /* Stable multi-radix sort of infos... */
  {
    uint64_t
      *bits = item.bitbucket,mybits[3],orall[3],andall[3];
    const int mybits_len = sizeof(mybits)/sizeof(*mybits);
    const int wordbits = sizeof(*mybits)*8;
    int ibit = 0;
    
    /* Figure out which bits differ between ranks, and
       create a packed bit representation of the sorting keys */ {
      int i;
		
      mybits[0] = item.info.pset;        /* Most significant bits to sort on */
      mybits[1] = item.info.rankInPset;
      mybits[2] = item.info.cpuId;       /* Least significant bits to sort on */
      
      MPI_Allreduce(mybits,orall,3,MPI_UNSIGNED_LONG_LONG,MPI_BOR,comm);
      MPI_Allreduce(mybits,andall,3,MPI_UNSIGNED_LONG_LONG,MPI_BAND,comm);
      
      for(i = 0; i<nfields; i++)
		  bits[i] = 0;
      for(i = mybits_len-1; i>=0; i--)
		  {
			 assert(wordbits == 64);
			 uint64_t mask = orall[i] ^ andall[i];
			 uint64_t data = mybits[i];
			 while(mask)
				{
				  if(mask & 1)
					 {
						bits[ibit/wordbits] |= (data & 1)<<(ibit%wordbits);
						ibit++;
					 }
				  data >>= 1;
				  mask >>= 1;
				}
		  }
    }
	 
    /* ibit now holds the number of bits in the packed representation */
    
    /* This is the actual radix sort, treating nbits at a time */
	 {
      const int nbits = 4,nbuckets = (1<<nbits);
      const uint64_t bucketmask = nbuckets-1;
      int buckets[nbuckets],groupidx[nbuckets],groupsum[nbuckets];
      
      for(int i = 0; i<nbuckets; i++)
		  buckets[i] = 0;
		
      while(ibit > 0)
		  {
			 const int ibucket = bits[0] & bucketmask;
			 
			 /* Figure out sorting index on last nbits bits */
			 {
				buckets[ibucket] = 1;
				MPI_Allreduce(buckets,groupsum,nbuckets,MPI_INT,MPI_SUM,comm);
				MPI_Scan(buckets,groupidx,nbuckets,MPI_INT,MPI_SUM,comm);
				buckets[ibucket] = 0;
			 }
			 
			 /* Shift out the last we just figured out the sorting index on */
			 {
				uint64_t carry = 0;
				int i;
				for(i = mybits_len-1; i>=0; i--)
				  {
					 const uint64_t tmp = bits[i] & bucketmask;
					 bits[i] = (bits[i]>>nbits) | (carry<<(wordbits-nbits));
					 carry = tmp;
				  }
			 }
			 
			 /* Send my info data to correct rank */
			 {
				MPI_Request req;
				int i,destrank = groupidx[ibucket] - 1;
				for(i = 0; i<ibucket; i++)
				  destrank += groupsum[i];	  
				
				MPI_Irecv(&tmpitem,sizeof(tmpitem),MPI_BYTE,MPI_ANY_SOURCE,19,comm,&req);
				MPI_Send(&item,sizeof(item),MPI_BYTE,destrank,19,comm);
				MPI_Wait(&req,MPI_STATUS_IGNORE);
				item = tmpitem;
			 }
		  
			 ibit -= nbits;
		  }
    }
  }
  
  /* Data is sorted... */
  
  int group_pid,group_size;
  /* Figure out group boundaries and sizes */
  {
    int igroup,group_leader;
	 
    MPI_Request req;
    uint64_t prev_pset = -1;
    int tag,tagpid;
    if(pid > 0)
      MPI_Irecv(&prev_pset,1,MPI_UNSIGNED_LONG_LONG,pid-1,17,comm,&req);
    if(pid < np-1)
      MPI_Send(&item.info.pset,1,MPI_UNSIGNED_LONG_LONG,pid+1,17,comm);
    if(pid > 0)
      MPI_Wait(&req,MPI_STATUS_IGNORE);
    tag = (pid == 0 || item.info.pset != prev_pset);
    
    /* Find out which group I am in */
	 {
      igroup = 0;
      MPI_Exscan(&tag,&igroup,1,MPI_INT,MPI_SUM,comm);
    }
	 
    tagpid = tag*pid;
    group_size = 0;
    MPI_Exscan(&tagpid,&group_leader,1,MPI_INT,MPI_MAX,comm);
    if(tag == 1)
		{
		  MPI_Irecv(&group_size,1,MPI_INT,MPI_ANY_SOURCE,15,comm,&req);
		  if(pid > 0)
			 {
				int prev_group_size = pid - group_leader;
				MPI_Send(&prev_group_size,1,MPI_INT,group_leader,15,comm);
			 }
		  group_leader = pid;
		}
    if(pid == np-1)
		{
		  int prev_group_size = np - group_leader;
		  MPI_Send(&prev_group_size,1,MPI_INT,group_leader,15,comm);
		}
    if(tag == 1)
      MPI_Wait(&req,MPI_STATUS_IGNORE);
    group_pid = pid - group_leader;
    
    /* Send group sizes to all */
	 {
      const uint64_t one = 1;
      uint64_t igroup_isize = 0,tmp;
      if(tag == 1) igroup_isize = (((uint64_t) igroup) << 32) | group_size;
      MPI_Scan(&igroup_isize,&tmp,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm);
      group_size = tmp & ((one<<32)-1);
    }
  }
  
  /* Set up ioTaskList array */
  {
    int i_am_io_task,io_task_list_idx;
    int nIoTasks,*ioTaskList;
	 
    /* Very weird formula, but let's be compatible for now.
       I believe that the formula does the right thing
       only wgeb ioTasksPerPet == 2.
    */
	 {
      int stride = group_size / ioTasksPerPset;
      if(group_size % ioTasksPerPset != 0) stride += 1;
      i_am_io_task = ((group_pid % stride) == 0);
    }
    
    io_task_list_idx = 0;
    MPI_Exscan(&i_am_io_task,&io_task_list_idx,1,MPI_INT,MPI_SUM,comm);
    MPI_Allreduce(&i_am_io_task,&nIoTasks,1,MPI_INT,MPI_SUM,comm);
    ioTaskList = (int *) ddcMalloc(sizeof(*ioTaskList) * nIoTasks);
    {
      int *tmpList = (int *) ddcMalloc(sizeof(*tmpList) * nIoTasks);
      int i;
      for(i = 0; i<nIoTasks; i++) tmpList[i] = 0;
      if(i_am_io_task) tmpList[io_task_list_idx] = item.info.mpiRank;
      MPI_Allreduce(tmpList,ioTaskList,nIoTasks,MPI_INT,MPI_SUM,comm);
      ddcFree(tmpList);
    }
    *ioTaskList_ptr = ioTaskList;
    return nIoTasks;
  }
}
#endif

/* Local Variables: */
/* tab-width: 3 */
/* End: */
