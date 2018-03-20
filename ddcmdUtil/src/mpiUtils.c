// $Id$

#include "mpiUtils.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <stdio.h>

#ifdef WITH_MPI
#include "ddcMalloc.h"
#include "utilities.h"
#include "lessThan.h"
#include "external.h"
#endif

typedef unsigned long long long64;

#ifdef WITH_MPI
void timestampBarrier(const char* msg, MPI_Comm comm)
{
   int myId;
   MPI_Comm_rank(comm, &myId);
   MPI_Barrier(comm);
   if (myId == 0) timestamp_anyTask(msg);
   MPI_Barrier(comm);
}
#endif // #ifdef WITH_MPI


/** The caller is responsible to ensure that the data array has
 *  sufficient capacity to hold nWant elements. */
#ifdef WITH_MPI
void distributeArray(unsigned char* data,
		     const unsigned nHave,
		     const unsigned nWant,
		     const unsigned width,
		     MPI_Comm comm)
{
   int nTasks;   MPI_Comm_size(comm, &nTasks);
   int myId;     MPI_Comm_rank(comm, &myId);
   
   unsigned* recvBuf = ddcMalloc(nTasks*2*sizeof(unsigned));
   unsigned sendBuf[2];
   sendBuf[0] = nHave;
   sendBuf[1] = nWant;
   
   MPI_Allgather(sendBuf, 2, MPI_INT, recvBuf, 2, MPI_INT, comm);

   long64* endWant = ddcMalloc(nTasks*sizeof(long64));
   long64  nGlobal = recvBuf[0];
   endWant[0] = recvBuf[1];
   for (int ii=1; ii<nTasks; ++ii)
   {
      nGlobal += recvBuf[(2*ii)];
      endWant[ii] = endWant[ii-1] + recvBuf[(2*ii)+1];
   }

   assert(nGlobal == endWant[nTasks-1]);
   
   long64 myBase=0;
   for (int ii=0; ii<myId; ++ii)
      myBase += recvBuf[(2*ii)];
   
   unsigned* dest = ddcMalloc(nHave*sizeof(unsigned));
   int iTask = 0;
   for (unsigned ii=0; ii<nHave; ++ii)
   {
      while (myBase+ii >= endWant[iTask])
	 ++iTask;
      dest[ii] = iTask;
   }

   assert(iTask < nTasks);
   
   ddcFree(recvBuf);
   ddcFree(endWant);   

   unsigned nLocal = nHave;
   assignArray(data, &nLocal, nWant, width, dest, 0, comm);
   assert(nLocal == nWant);

   ddcFree(dest);
}
#endif // #ifdef WITH_MPI


/** Returns 1 if a is sorted in ascending order, 0 otherwise */
static int sorted(const unsigned* a, unsigned n)
{
   for (unsigned ii=1; ii<n; ++ii)
      if (a[ii] < a[ii-1])
	 return 0;
   return 1;
}


#ifdef WITH_MPI
void assignArray(unsigned char* data,
		 unsigned* nLocal,
		 const unsigned capacity,
		 const unsigned width,
		 const unsigned* dest,
		 int verbose,
		 MPI_Comm comm)
{
   const int msgInfoTag=1;
   const int msgTag=2;

   int nTasks;   MPI_Comm_size(comm, &nTasks);
   int myId;     MPI_Comm_rank(comm, &myId);

   if ( sorted(dest, *nLocal) == 0)
      assert(1==0);
   
   if (verbose)
      timestampBarrier("Negotiating routing", comm);
   int* isTarget = ddcMalloc(nTasks*sizeof(int));
   int* nSenders = ddcMalloc(nTasks*sizeof(int));
   for (int ii=0; ii<nTasks; ++ii)
      isTarget[ii] = 0;
   for (unsigned ii=0; ii<*nLocal; ++ii)
      isTarget[dest[ii]] = 1;
   unsigned nMsgSend = 0;
   for (int ii=0; ii<nTasks; ++ii)
      nMsgSend += isTarget[ii];
   int* offset = ddcMalloc( (nMsgSend+1)*sizeof(int) );
   offset[0] = 0;
   offset[nMsgSend] = *nLocal;
   unsigned cnt=0;
   for (unsigned ii=1; ii<*nLocal; ++ii)
      if (dest[ii] != dest[ii-1])
      {
	 offset[++cnt] = ii;
      }
   assert(cnt == nMsgSend-1 || *nLocal == 0); 
   
   MPI_Allreduce(isTarget, nSenders, nTasks, MPI_INT, MPI_SUM, comm);
   unsigned nMsgRecv = nSenders[myId];
   ddcFree(isTarget);
   ddcFree(nSenders);

   MPI_Request* recvRequest = ddcMalloc(nMsgRecv*sizeof(MPI_Request));
   MPI_Request* sendRequest = ddcMalloc(nMsgSend*sizeof(MPI_Request));
   int* msgInfo = ddcMalloc(2*nMsgRecv*sizeof(int));   
   int* infoBuf = ddcMalloc(2*nMsgSend*sizeof(int));
   for (unsigned ii=0; ii<nMsgRecv; ++ii)
      MPI_Irecv(msgInfo+2*ii, 2, MPI_INT, MPI_ANY_SOURCE, msgInfoTag, comm, recvRequest+ii);
   for (unsigned ii=0; ii<nMsgSend; ++ii)
   {
      cnt = offset[ii + 1] - offset[ii];
      unsigned target = dest[offset[ii]];
      infoBuf[2*ii] = myId;
      infoBuf[2*ii+1] = cnt;
      MPI_Isend(infoBuf+2*ii, 2, MPI_INT, target, msgInfoTag, comm, sendRequest+ii);
   }
   MPI_Waitall(nMsgSend, sendRequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(nMsgRecv, recvRequest, MPI_STATUSES_IGNORE);

   qsort(msgInfo, nMsgRecv, 2*sizeof(int), intLess);
   
   if (verbose)
      timestampBarrier("Communicating data", comm);
   unsigned nRecv = 0;
   for (unsigned ii=0; ii<nMsgRecv; ++ii)
      nRecv += msgInfo[2*ii+1];

   if (nRecv > capacity)
   {
      printf("Error (task %d): Insufficient capacity to receive messages in assignment.\n"
	     "                 capacity = %u\n"
	     "                 nRecv = %u\n", myId, capacity, nRecv);
      MPI_Abort(comm, -1);
   }

   unsigned char* sendBuf = ddcMalloc(*nLocal*width);
   memcpy(sendBuf, data, *nLocal*width);
      
   unsigned index = 0;
   for (unsigned ii=0; ii<nMsgRecv; ++ii)
   {
      unsigned source = msgInfo[2*ii];
      cnt = msgInfo[2*ii+1];
      int sendSize = cnt * width;
      MPI_Irecv(data + index, sendSize, MPI_BYTE, source, msgTag, comm, recvRequest+ii);
      index += sendSize;
   }

   index = 0;
   for (unsigned ii=0; ii<nMsgSend; ++ii)
   {
      int sendSize  = (offset[ii + 1] - offset[ii])*width;
      unsigned target = dest[offset[ii]];
      MPI_Isend(sendBuf + index, sendSize, MPI_BYTE, target, msgTag, comm, sendRequest+ii);
      index += sendSize;
   }
   MPI_Waitall(nMsgSend, sendRequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(nMsgRecv, recvRequest, MPI_STATUSES_IGNORE);
   
   if (verbose)
      timestampBarrier("End receiving data", comm);
   ddcFree(sendBuf);
   *nLocal = nRecv;
   ddcFree(msgInfo);
   ddcFree(infoBuf);
   ddcFree(offset);
   ddcFree(recvRequest);
   ddcFree(sendRequest);
}
#endif // #ifdef WITH_MPI


int  getRank(int flag)
{
#ifdef WITH_MPI
   int rank;
   MPI_Comm comm = MPI_COMM_NULL;
   switch (flag)
   {
     case -1:
      comm = MPI_COMM_WORLD;
      break;
     case 0:
      comm = COMM_LOCAL;
      break;
     default:
      assert(1==0);
   }
   MPI_Comm_rank(comm, &rank);
   return rank;
#else
   return 0;
#endif // #ifdef WITH_MPI
}


int getSize(int flag)
{
#ifdef WITH_MPI
   int size;
   MPI_Comm comm = MPI_COMM_NULL;
   switch (flag)
   {
     case -1:
      comm = MPI_COMM_WORLD;
      break;
     case 0:
      comm = COMM_LOCAL;
      break;
     default:
      assert(1==0);
   }
   MPI_Comm_size(comm, &size);
   return size;
#else
   return 1;
#endif // #ifdef WITH_MPI
}

void abortAll_core(int rc,char *file,int line)
{
#ifdef WITH_MPI
  int pid = getRank(0);
  printf("abortAll(%d) called by rank %d from %s:%d\n",
			rc,pid,file,line);
   MPI_Abort(MPI_COMM_WORLD, rc);
#endif
   exit(rc); 
}


void WAIT(int flag)
{
#ifdef WITH_MPI
   MPI_Comm comm = MPI_COMM_NULL;
   switch (flag)
   {
     case -1:
      comm = MPI_COMM_WORLD;
      break;
     case 0:
      comm = COMM_LOCAL;
      break;
     default:
      assert(1==0);
   }
   MPI_Barrier(comm);
#endif // ifdef WITH_MPI
   return;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
