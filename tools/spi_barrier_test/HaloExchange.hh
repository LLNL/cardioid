#ifndef HALO_EXCHANGE_HH
#define HALO_EXCHANGE_HH

#include <vector>
#include <cassert>
#include "CommTable.hh"
#include <iostream>
#include "PerformanceTimers.hh"
using namespace std;

#ifdef SPI
#include "spi_impl.h"
#endif

/**
 *  There are two ways of using this class.  The simplest is to call
 *  execute().  This method performs all aspects of the HaloExchange
 *  including tranferring data to the send buffer and leaving the remote
 *  data at the end of the data vector that is passed in.  If necessary
 *  the vector will be expanded to make room for the remote elements.
 *  If you call execute() you should not count on any aspect of the
 *  internal state of the HaloExchange
 *
 *  The other way of calling this class is to do practically everything
 *  yourself.  First you need to fill the send buffer.  You can either
 *  call fillSendBuffer, or you can get the sendMap and sendBuffer and
 *  do it yourself.  Then call startComm() to get the communication going.
 *  Next call wait().  Finally you'll need to getRecvBuf() and put the
 *  data where you actually want it to be.  When you're working in this
 *  "advanced" mode the class receives data into an internal buffer
 *  which must be copied out manually before starting the next halo
 *  exchange. 
 */

template <class T, class Allocator = std::allocator<T> >
class HaloExchangeBase
{
 public:
   HaloExchangeBase(const std::vector<int>& sendMap, const CommTable* comm)
     : width_(sizeof(T)), commTable_(comm), sendMap_(sendMap)
   {
     sendBuf_ = new T[commTable_->sendSize()];
     recvBuf_ = new T[commTable_->recvSize()*2]; // why 2*?
   };

   ~HaloExchangeBase()
   {
     delete [] recvBuf_;
     delete [] sendBuf_;
   };

   T* getSendBuf() {return sendBuf_;}
   const vector<int>& getSendMap() const {return sendMap_;}
   
   void fillSendBuffer(std::vector<T, Allocator>& data)
   {
      startTimer(PerformanceTimers::haloMove2BufTimer);
      // fill send buffer
      assert(sendMap_.size() == commTable_->sendSize()); 
      for (unsigned ii=0; ii<sendMap_.size(); ++ii) { sendBuf_[ii]=data[sendMap_[ii]]; }
      stopTimer(PerformanceTimers::haloMove2BufTimer);
   };

//   void set_recv_buf(int bw) { for (unsigned ii=0; ii<commTable_->recvSize();++ii) {recv_buf_[ii+bw*commTable_->recvSize()]=-1;} }

   void dumpRecvBuf()
   {
      cout << " recv is : ";
      for (unsigned ii=0; ii<(commTable_->recvSize()) * 2;++ii) {cout << recvBuf_[ii] << " ";}
      cout << endl;
   }
     
   void dumpSsendBuf()
   {
      cout << " send is : ";
      for (unsigned ii=0; ii<commTable_->sendSize();++ii) {cout << sendBuf_[ii] << " ";}
      cout << endl;
   }

   //virtual void execute(std::vector<T>& data, unsigned nLocal);
   //virtual void complete(std::vector<T>& data, unsigned nLocal);
 protected:
   
   unsigned width_;
   const CommTable* commTable_;
   std::vector<int> sendMap_;
   T* sendBuf_;
   T* recvBuf_;
};

#ifdef SPI
// spi version
template <class T, class Allocator = std::allocator<T> >
class HaloExchange : public HaloExchangeBase<T, Allocator>
{
  public:
   HaloExchange(const std::vector<int>& sendMap, const CommTable* comm)
   : HaloExchangeBase<T, Allocator>(sendMap,comm), bw_(1)
  {
    //create mapping table
//    mapping_table(&spiHdl_);
    myID=mapping_table_new(&spiHdl_);
 
    //setup base address table
    //search and allocate
    setup_bat(&spiHdl_,(void*)recvBuf_,commTable_->recvSize()*width_);
 
    //setup injection memory fifo
    //search,allocate,initialize,activate
    setup_inj_fifo(&spiHdl_);

    //setup descriptor
    setup_descriptors(commTable_->_offsets, &(commTable_->_putTask[0]),
                      &(commTable_->_putIdx[0]), commTable_->_putTask.size(),
                      (void*)sendBuf_, commTable_->sendSize(), &spiHdl_,width_);

    barrier();

  };

   ~HaloExchange()
   {
     free_spi(&spiHdl_);
   };

   T* getRecvBuf() { return &(recvBuf_[bw_*commTable_->recvSize()]); };

   void startComm()
   {
     bw_=1-bw_;
#pragma omp critical 
     execute_spi_alter(&spiHdl_,commTable_->_putTask.size(),bw_);
   }
   
   
   void execute(vector<T, Allocator>& data, int nLocal)
   {
      fillSendBuffer(data);
      data.resize(nLocal + commTable_->recvSize());
      startComm();
      wait();
      for (int ii=0; ii<commTable_->recvSize(); ++ii)
         data[nLocal+ii] = recvBuf_[ii];
   };

//    void execute()
//    {
//      bw_=1-bw_;
//      execute_spi_alter(&spiHdl_,commTable_->_putTask.size(),bw_);
//    };
   void wait()
   {
#pragma omp critical 
      //complete_spi_alter(&spiHdl_, commTable_->_recvTask.size(), commTable_->_offsets[3], bw_, width_ );
      complete_spi_alter_monitor(&spiHdl_, commTable_->_recvTask.size(), commTable_->_offsets[3], commTable_->_offsets[4], bw_, width_, myID );
   };

//    void execute3() {execute_spi(&spiHdl_,commTable_->_putTask.size());};
//    void execute2() {execute_spi_2(&spiHdl_,commTable_->_putTask.size());};
//    void complete2() {complete_spi(&spiHdl_, commTable_->_recvTask.size());};
//    void dump_mapping_table() { spi_dump_mapping( &spiHdl_); };
   void barrier() {global_sync(&spiHdl_);};
   void barrierWithTimeout(uint64_t timeout) {global_sync_2(&spiHdl_,timeout);};


 private:
   
   int bw_;
   spi_hdl_t spiHdl_;
   uint32_t myID;
};

#else // not SPI

// MPI version
template <class T, class Allocator = std::allocator<T> >
class HaloExchange : public HaloExchangeBase<T, Allocator>
{
  public:
   using HaloExchangeBase<T, Allocator>::commTable_;
   using HaloExchangeBase<T, Allocator>::sendBuf_;
   using HaloExchangeBase<T, Allocator>::recvBuf_;
   using HaloExchangeBase<T, Allocator>::width_;

   HaloExchange(const std::vector<int>& sendMap, const CommTable* comm) 
   : HaloExchangeBase<T, Allocator>(sendMap,comm),
     recvReq_(comm->recvSize()),
     sendReq_(comm->sendSize())
   {};

   T* getRecvBuf() {return recvBuf_;}
   void execute(vector<T>& data, int nLocal)
   {
      fillSendBuffer(data);
      T* tmp = recvBuf_;
      data.resize(nLocal + commTable_->recvSize());
      recvBuf_ = (&data[nLocal]);
      startComm();
      wait();
      recvBuf_ = tmp;
   }
   
   void startComm()
   {
#pragma omp critical 
      {
         
      char* sendBuf = (char*)sendBuf_;
      char* recvBuf = (char*)recvBuf_;

      MPI_Request* recvReq = &recvReq_[0];
      const int tag = 151515;
      for (unsigned ii=0; ii< commTable_->_recvTask.size(); ++ii)
      {
         assert(recvBuf);
         unsigned sender = commTable_->_recvTask[ii];
         unsigned nItems = commTable_->_recvOffset[ii+1] - commTable_->_recvOffset[ii];
         unsigned len = nItems * width_;
         char* recvPtr = recvBuf + commTable_->_recvOffset[ii]*width_;
         MPI_Irecv(recvPtr, len, MPI_CHAR, sender, tag, commTable_->_comm, recvReq+ii);
      }
 
      MPI_Request* sendReq = &sendReq_[0];
      for (unsigned ii=0; ii<commTable_->_sendTask.size(); ++ii)
      {
         assert(sendBuf);
         unsigned target = commTable_->_sendTask[ii];
         unsigned nItems = commTable_->_sendOffset[ii+1] - commTable_->_sendOffset[ii];
         unsigned len = nItems * width_;
         char* sendPtr = sendBuf + commTable_->_sendOffset[ii]*width_;
         MPI_Isend(sendPtr, len, MPI_CHAR, target, tag, commTable_->_comm, sendReq+ii);
      }
      }
      
   };

   void wait()
   {
#pragma omp critical              
   {                                            
   
      MPI_Waitall(commTable_->_sendTask.size(), &sendReq_[0], MPI_STATUS_IGNORE);
      MPI_Waitall(commTable_->_recvTask.size(), &recvReq_[0], MPI_STATUS_IGNORE);
   }
   };

   void barrier()
   {
      MPI_Barrier(commTable_->_comm);
   }
   

 private:
   std::vector<MPI_Request> recvReq_;
   std::vector<MPI_Request> sendReq_;
   
};

#endif // ifdef SPI

#endif
