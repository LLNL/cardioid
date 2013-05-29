#ifndef MPITPL
#define MPITPL

#include <mpi.h>
#include <vector>

/** Templatized wrappers around MPI functions.  Intended to provide
 *  simplified calling sequences in C++, especially for stl objects. */

template<typename T>
inline void
allGather(std::vector<T>& data, MPI_Comm comm);

template<typename T>
inline void
allGather(std::vector<T>& data, size_t nDataPerTask, MPI_Comm comm);

template<typename T>
inline void
allReduce(std::vector<T>& data, MPI_Datatype type, MPI_Op op, MPI_Comm comm);



template<typename T>
inline void
allGather(std::vector<T>& data, MPI_Comm comm)
{
   allGather(data, 1, comm);
}

template<typename T>
inline void
allGather(std::vector<T>& data, size_t nDataPerTask, MPI_Comm comm)
{
   int myRank;
   MPI_Comm_rank(comm, &myRank);

   size_t localStart = myRank*nDataPerTask;
   typename std::vector<T>::const_iterator first = data.begin() + localStart;
   typename std::vector<T>::const_iterator last = first + nDataPerTask;
//   assert( data[localStart] == *first);

   std::vector<T> tmp(first, last);

   void* sendBuf = (void*) &tmp[0];
   void* recvBuf = (void*) &data[0];
   int nSend = nDataPerTask * sizeof(T);
   int nRecv = nDataPerTask * sizeof(T);
   
   MPI_Allgather(sendBuf, nSend, MPI_CHAR, recvBuf, nRecv, MPI_CHAR, comm);
}

template<typename T>
inline void
allReduce(std::vector<T>& data, MPI_Datatype type, MPI_Op op, MPI_Comm comm)
{
   std::vector<T> tmp(data);
   void* sendBuf = (void*) &tmp[0];
   void* recvBuf = (void*) &data[0];
   int nSend = data.size();
   MPI_Allreduce(sendBuf, recvBuf, nSend, type, op, comm);
}


#endif
