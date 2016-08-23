// $Id:
#ifndef COMM_TABLE_HH
#define COMM_TABLE_HH

#include <mpi.h>
#include <vector>
#include <iostream>
#include <inttypes.h>

/** Implements a communications table
 *
 *  Note that this class handles only the communication, it does not
 *  handle the packing and unpacking of the send and receive buffers.
 *
 *  LIMITATIONS:
 *
 *  This implementation strives for simplicity of implementation at the
 *  expense of optimization.  One consequence is that this method uses
 *  more memory for intermediate storage that is necessary for at least
 *  some cases.
 *
 *  The send and receive buffers are assumed to composed of
 *  non-overlapping sections.  This means whem multiple copies of the
 *  same data are sent to multiple receivers there will be several
 *  copies in the send buffer.
 *
 *  All sends and receives are processed together.  This class does
 *  not exploit potential memory savings from doing only one send or
 *  recv at a time to reduce intermediate storage space.
 */

class CommTable
{
  public:

   CommTable(const std::vector<int>& sendTask,
             const std::vector<int>& sendOffset,
             MPI_Comm comm);
   ~CommTable();

   void dump_put();
   inline uint32_t sendSize() const { return _sendOffset[_sendTask.size()]; };
   inline uint32_t recvSize() const { return _recvOffset[_recvTask.size()]; };

   unsigned nRemote();
   unsigned nNbrs();
   std::vector<int> msgSize(); // nItems in each msg that will be sent


  MPI_Comm  _comm;

  std::vector<int> _sendTask;   // ranks to send to
  std::vector<int> _recvTask;   // ranks to recv from
  std::vector<int> _putTask; // offsets to sendBuf. offsets are size n+1

  std::vector<int> _sendOffset; // offsets to sendBuf. offsets are size n+1
  std::vector<int> _recvOffset; // offsets to recvBuf
  std::vector<int> _putOffset; 
  std::vector<int> _putCntOffset; 
  std::vector<int> _putIdx; 
  std::vector<int> _recvIdx; 
  int** _offsets;

};


#endif // #ifndef COMM_TABLE_HH


/* Local Variables: */
/* tab-width: 3 */
/* End: */
