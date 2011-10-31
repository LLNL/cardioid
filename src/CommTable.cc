// $Id$

#include "CommTable.hh"

#include <cassert>

using std::vector;

/** Given the send-side routing information, the constructor figures out
 *  the corresponding recv side information, nRecv, sourceTask, and
 *  recvOffset.
 */
CommTable::CommTable(
   const vector<int>& sendTask, const vector<int>& sendOffset, MPI_Comm comm)
: _sendTask(sendTask), _sendOffset(sendOffset), _comm(comm)
{
   const int tag = 141414;
   int nTasks;
   MPI_Comm_size(_comm, &nTasks);
   int myId;
   MPI_Comm_rank(_comm, &myId);
   
   // How many tasks will we receive data from?
   int* msgCount = new int[nTasks];
   int* recvCnt =   new int[nTasks];
   for (int ii=0; ii<nTasks; ++ii)
   {
      msgCount[ii] = 0;
      recvCnt[ii] = 1;
   }
   for (unsigned ii=0; ii<_sendTask.size(); ++ii)
   {
      msgCount[_sendTask[ii]] += 1;
      assert(msgCount[_sendTask[ii]] == 1);
   }
   int nRecv;
   MPI_Reduce_scatter(msgCount, &nRecv, recvCnt, MPI_INT, MPI_SUM, _comm);
   
   delete [] recvCnt;
   delete [] msgCount;
   _recvTask.resize(nRecv);
   _recvOffset.resize(nRecv+1);

   // Find out which tasks we will be exchanging data with.
   int* recvBuf = new int[2*_recvTask.size()];
   int* sendBuf = new int[2*_sendTask.size()];
   MPI_Request* recvReq = new MPI_Request[_recvTask.size()];
   MPI_Request* sendReq = new MPI_Request[_sendTask.size()];
	
   for (unsigned ii=0; ii<_recvTask.size(); ++ii)
   {
      MPI_Irecv(recvBuf+2*ii, 2, MPI_INT, MPI_ANY_SOURCE, tag, _comm, recvReq+ii);
   }
   for (unsigned ii=0; ii<_sendTask.size(); ++ii)
   {
      sendBuf[2*ii] = myId;
      sendBuf[2*ii+1] = _sendOffset[ii+1] - _sendOffset[ii];
      int dest = _sendTask[ii];
      MPI_Isend(sendBuf+2*ii, 2, MPI_INT, dest, tag, _comm, sendReq+ii);
   }
   MPI_Waitall(_sendTask.size(), sendReq, MPI_STATUSES_IGNORE);
   MPI_Waitall(_recvTask.size(), recvReq, MPI_STATUSES_IGNORE);
	
   // set recvTask and recvOffset arrays.
   _recvOffset[0] = 0;
   for (unsigned ii=0; ii<_recvTask.size(); ++ii)
   {
      _recvTask[ii] = recvBuf[2*ii];
      _recvOffset[ii+1] = _recvOffset[ii] + recvBuf[2*ii+1];
   }

   delete[] recvBuf;
   delete[] sendBuf;
   delete[] recvReq;
   delete[] sendReq;
}


CommTable::~CommTable()
{
}

void CommTable::execute(void* sendBufV, void* recvBufV, unsigned width)
{
   char* sendBuf = (char*)sendBufV;
   char* recvBuf = (char*)recvBufV;
   
   const int tag = 151515;
   MPI_Request recvReq[_recvTask.size()];
   MPI_Request sendReq[_sendTask.size()];

   for (unsigned ii=0; ii<_recvTask.size(); ++ii)
   {
      assert(recvBuf);
      unsigned sender = _recvTask[ii];
      unsigned nItems = _recvOffset[ii+1] - _recvOffset[ii];
      unsigned len = nItems * width;
      char* recvPtr = recvBuf + _recvOffset[ii]*width;
      MPI_Irecv(recvPtr, len, MPI_CHAR, sender, tag, _comm, recvReq+ii);
   }

   for (unsigned ii=0; ii<_sendTask.size(); ++ii)
   {
      assert(sendBuf);
      unsigned target = _sendTask[ii];
      unsigned nItems = _sendOffset[ii+1] - _sendOffset[ii];
      unsigned len = nItems * width;
      char* sendPtr = sendBuf + _sendOffset[ii]*width;
      MPI_Isend(sendPtr, len, MPI_CHAR, target, tag, _comm, sendReq+ii);
   }

   MPI_Waitall(_sendTask.size(), sendReq, MPI_STATUS_IGNORE);
   MPI_Waitall(_recvTask.size(), recvReq, MPI_STATUS_IGNORE);
}

unsigned CommTable::nRemote()
{
   return _recvOffset.back();
}

