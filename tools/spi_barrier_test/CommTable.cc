// $Id: CommTable.cc 1028 2012-09-18 16:44:28Z kim70 $

#include "CommTable.hh"

#include <cassert>
#include <iostream>
using namespace std;

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
   int myId;
   MPI_Comm_size(_comm, &nTasks);
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
   _recvIdx.resize(nRecv);

   // Find out which tasks we will be exchanging data with.
   int* recvBuf = new int[4*_recvTask.size()];
   int* sendBuf = new int[4*_sendTask.size()];
   MPI_Request* recvReq = new MPI_Request[_recvTask.size()];
   MPI_Request* sendReq = new MPI_Request[_sendTask.size()];
        
   for (unsigned ii=0; ii<_recvTask.size(); ++ii)
   {
      MPI_Irecv(recvBuf+3*ii, 3, MPI_INT, MPI_ANY_SOURCE, tag, _comm, recvReq+ii);
   }

   for (unsigned ii=0; ii<_sendTask.size(); ++ii)
   {
      sendBuf[3*ii] = myId;
      sendBuf[3*ii+1] = _sendOffset[ii+1] - _sendOffset[ii];
      sendBuf[3*ii+2] = ii;
      int dest = _sendTask[ii];
      MPI_Isend(sendBuf+3*ii, 3, MPI_INT, dest, tag, _comm, sendReq+ii);
   }
   MPI_Waitall(_sendTask.size(), sendReq, MPI_STATUSES_IGNORE);
   MPI_Waitall(_recvTask.size(), recvReq, MPI_STATUSES_IGNORE);
        
   // hfwen: Need to add a Barrier here
   MPI_Barrier(_comm);

 
   // set recvTask and recvOffset arrays.
   _recvOffset[0] = 0;
   for (unsigned ii=0; ii<_recvTask.size(); ++ii)
   {
      _recvTask[ii] = recvBuf[3*ii];
   #ifdef debug
       cout << _recvTask[ii] << " ";
   #endif
      _recvOffset[ii+1] = _recvOffset[ii] + recvBuf[3*ii+1];
      _recvIdx[ii] = recvBuf[3*ii+2];
   }

   #ifdef SPI
   //send out the recv offset so that each node knows where to put data 
   for (unsigned ii=0; ii<_sendTask.size(); ++ii)
   {
      MPI_Irecv(sendBuf+ii*4, 4, MPI_INT, MPI_ANY_SOURCE, tag, _comm, sendReq+ii);
   }
   #ifdef debug
   cout << endl << "dest=" ;
   #endif
   for (unsigned ii=0; ii<_recvTask.size(); ++ii)
   {
      recvBuf[4*ii] = myId;
      recvBuf[4*ii+1] = _recvOffset[ii];
      recvBuf[4*ii+2] = ii;
      recvBuf[4*ii+3] = _recvIdx[ii];
      int dest = _recvTask[ii];
   #ifdef debug
      cout << dest << ":" << _recvOffset[ii]  << " " ;
   #endif
      MPI_Isend(recvBuf+ii*4, 4, MPI_INT, dest, tag, _comm, recvReq+ii);
   }
   #ifdef debug
   cout << endl;
   #endif
//   MPI_Waitall(_recvTask.size(), recvReq, MPI_STATUSES_IGNORE);
   MPI_Waitall(_sendTask.size(), sendReq, MPI_STATUSES_IGNORE);
   MPI_Waitall(_recvTask.size(), recvReq, MPI_STATUSES_IGNORE);
 
   // set recvTask and recvOffset arrays.
   _putTask.resize(_sendTask.size());
   _putOffset.resize(_sendTask.size());
   _putCntOffset.resize(_sendTask.size());
   _putIdx.resize(_sendTask.size());
   for (unsigned ii=0; ii<_sendTask.size(); ++ii)
   {
      _putTask[ii] = sendBuf[4*ii];
      _putOffset[ii] = sendBuf[4*ii+1];
      _putCntOffset[ii] = sendBuf[4*ii+2];
      _putIdx[ii] = sendBuf[4*ii+3];
   }
   #endif

   MPI_Barrier(_comm);

   _offsets = 0;
   #ifdef SPI
   _offsets = new int*[5];  //need to pass these offsets to spi_implementation
   _offsets[0]=&(_sendOffset[0]);
   _offsets[1]=&(_putOffset[0]);
   _offsets[2]=&(_putCntOffset[0]);
   _offsets[3]=&(_recvOffset[0]);
   _offsets[4]=&(_recvTask[0]);
   #endif
   
   delete[] recvBuf;
   delete[] sendBuf;
   delete[] recvReq;
   delete[] sendReq;

}

void CommTable::dump_put()
{
  for(int ii=0 ; ii < _putTask.size() ; ++ii )
    cout << _putTask[ii] <<"=" << _sendTask[_putIdx[ii]] << ":" << _putOffset[ii] << ":" <<  _putCntOffset[ii] << " ";
  cout<<endl;
}

CommTable::~CommTable()
{
  delete [] _offsets;
}

unsigned CommTable::nRemote()
{
   return _recvOffset.back();
}

unsigned CommTable::nNbrs()
{
   return _sendTask.size();
}

vector<int> CommTable::msgSize()
{
   vector<int> tmp(_sendOffset.size()-1);
   for (unsigned ii=1; ii<_sendOffset.size(); ++ii)
      tmp[ii-1] = _sendOffset[ii] - _sendOffset[ii-1];
   return tmp;
}
