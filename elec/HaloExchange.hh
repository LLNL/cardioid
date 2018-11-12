#ifndef HALO_EXCHANGE_HH
#define HALO_EXCHANGE_HH

#include <vector>
#include <cassert>
#include "CommTable.hh"
#include <iostream>
#include "PerformanceTimers.hh"
#include "DeviceFor.hh"
#include "cudaNVTX.h"

#ifdef SPI
#include "spi_impl.h"
#endif

/**
 *  There are two ways of using this class.  The simplest is to call
 *  execute().  This method performs all aspects of the HaloExchange
 *  including tranferring data to the send buffer and leaving the remote
 *  data in the internal recieve buffer.  Users are required to copy the
 *  data out themselves.  If you call execute() you should not count on
 *  any aspect of the internal state of the HaloExchange.
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

template <class T>
class HaloExchangeBase
{
 public:
   HaloExchangeBase(const std::vector<int>& sendMap, const CommTable* comm)
     : width_(sizeof(T)), commTable_(comm), sendMap_(sendMap)
   {
      sendBuf_.resize(commTable_->sendSize());
   };

   virtual ~HaloExchangeBase() {};

   rw_mgarray_ptr<T>   getSendBuf() {return sendBuf_;}
   ro_mgarray_ptr<int> getSendMap() const {return sendMap_;}

   int recvSize() { return commTable_->recvSize(); }


   void execute(rw_mgarray_ptr<T> _data, int nLocal)
   {
      fillSendBuffer(_data);
      startComm();
      wait();
      {
         rw_array_ptr<T> data = _data.useOn(CPU);
         ro_array_ptr<T> recvBuf = getRecvBuf().useOn(CPU);
         for (int ii=0; ii<recvSize(); ii++)
         {
            data[nLocal+ii] = recvBuf[ii];
         }
      }
   }
   virtual void fillSendBuffer(ro_mgarray_ptr<T> _externalData)
   {
      startTimer(PerformanceTimers::haloMove2BufTimer);

      ro_array_ptr<int> sendMap = sendMap_.readonly(CPU);
      wo_array_ptr<T> sendBuf = sendBuf_.writeonly(CPU);
      ro_array_ptr<T> externalData = _externalData.useOn(CPU);
      HOST_PARALLEL_FORALL(sendMap.size(), ii,
                           sendBuf[ii] = externalData[sendMap[ii]]);
      stopTimer(PerformanceTimers::haloMove2BufTimer);
         
   }
   virtual ro_mgarray_ptr<T> getRecvBuf() = 0;
   virtual void startComm() = 0;
   virtual void wait() = 0;
   virtual void barrier() = 0;

 protected:
   
   unsigned width_;
   const CommTable* commTable_;
   lazy_array<int> sendMap_;
   lazy_array<T> sendBuf_;
};



#ifdef SPI
// spi version
template <class T>
class HaloExchange : public HaloExchangeBase<T>
{
  public:
   using HaloExchangeBase<T>::commTable_;
   using HaloExchangeBase<T>::sendBuf_;
   using HaloExchangeBase<T>::sendMap_;
   using HaloExchangeBase<T>::width_;
   HaloExchange(const std::vector<int>& sendMap, const CommTable* comm)
   : HaloExchangeBase<T>(sendMap,comm), bw_(1)
  {
     recvBuf_.resize(commTable_->recvSize()*2);

     rw_array_ptr<T> sendBuf = sendBuf_.useOn(CPU);
     rw_array_ptr<T> recvBuf = recvBuf_.useOn(CPU);
     
     //create mapping table
     //mapping_table(&spiHdl_);
     myID=mapping_table_new(&spiHdl_);
 
     //setup base address table
     //search and allocate
     setup_bat(&spiHdl_,(void*)recvBuf.raw(),commTable_->recvSize()*width_);
 
     //setup injection memory fifo
     //search,allocate,initialize,activate
     setup_inj_fifo(&spiHdl_);

     
     //setup descriptor
     setup_descriptors(commTable_->_offsets, &(commTable_->_putTask[0]),
                       &(commTable_->_putIdx[0]), commTable_->_putTask.size(),
                       (void*)sendBuf.raw(), commTable_->sendSize(), &spiHdl_,width_);

     barrier();
  };

   ~HaloExchange()
   {
     free_spi(&spiHdl_);
   };

   virtual void startComm()
   {
     bw_=1-bw_;
#pragma omp critical
     {
        sendBuf_.readonly(CPU);
        recvBuf_.writeonly(CPU);
        execute_spi_alter(&spiHdl_,commTable_->_putTask.size(),bw_);
     }
     
   }

   virtual void wait()
   {
#pragma omp critical
      {
         //complete_spi_alter(&spiHdl_, commTable_->_recvTask.size(), commTable_->_offsets[3], bw_, width_ );
         complete_spi_alter_monitor(&spiHdl_, commTable_->_recvTask.size(), commTable_->_offsets[3], commTable_->_offsets[4], bw_, width_, myID );
      }
   };

//    void execute3() {execute_spi(&spiHdl_,commTable_->_putTask.size());};
//    void execute2() {execute_spi_2(&spiHdl_,commTable_->_putTask.size());};
//    void complete2() {complete_spi(&spiHdl_, commTable_->_recvTask.size());};
//    void dump_mapping_table() { spi_dump_mapping( &spiHdl_); };
   virtual void barrier() {global_sync(&spiHdl_);};
   void barrierWithTimeout(uint64_t timeout) {global_sync_2(&spiHdl_,timeout);};
   
   virtual ro_mgarray_ptr<T> getRecvBuf() { return recvBuf_.readonly().slice(recvSize()*bw_,recvSize()); }



 private:
   lazy_array<T> recvBuf_; //double buffer
   int bw_;
   spi_hdl_t spiHdl_;
   uint32_t myID;
};

#else // not SPI

// MPI version
template <class T>
class HaloExchange : public HaloExchangeBase<T>
{
  public:
   using HaloExchangeBase<T>::commTable_;
   using HaloExchangeBase<T>::sendBuf_;
   using HaloExchangeBase<T>::sendMap_;
   using HaloExchangeBase<T>::width_;

   HaloExchange(const std::vector<int>& sendMap, const CommTable* comm) 
   : HaloExchangeBase<T>(sendMap,comm),
     recvReq_(comm->recvSize()),
     sendReq_(comm->sendSize())
   {
      recvBuf_.resize(comm->recvSize());
   }
   
   virtual void startComm()
   {
#pragma omp critical 
      {
         const char* sendBuf = (const char*)sendBuf_.readonly(CPU).raw();
         char* recvBuf = (char*)recvBuf_.writeonly(CPU).raw();

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
            const char* sendPtr = sendBuf + commTable_->_sendOffset[ii]*width_;
            MPI_Isend(sendPtr, len, MPI_CHAR, target, tag, commTable_->_comm, sendReq+ii);
         }
      }
   };
   
   virtual void wait()
   {
#pragma omp critical              
      {                                            
         MPI_Waitall(commTable_->_sendTask.size(), &sendReq_[0], MPI_STATUS_IGNORE);
         MPI_Waitall(commTable_->_recvTask.size(), &recvReq_[0], MPI_STATUS_IGNORE);
      }
   };

   virtual void barrier()
   {
      MPI_Barrier(commTable_->_comm);
   }
   
   virtual ro_mgarray_ptr<T> getRecvBuf()
   {
      return recvBuf_;
   }
   
 private:
   lazy_array<T> recvBuf_;
   std::vector<MPI_Request> recvReq_;
   std::vector<MPI_Request> sendReq_;
};

#endif // ifdef SPI

template <class T>
class HaloExchangeDevice : public HaloExchange<T>
{
   
 public:
   using HaloExchangeBase<T>::sendBuf_;
   using HaloExchangeBase<T>::sendMap_;
   
   HaloExchangeDevice(const std::vector<int>& sendMap, const CommTable* comm)
     : HaloExchange<T>(sendMap, comm)
   {};

   ~HaloExchangeDevice()
   {
   };

   void fillSendBuffer(ro_mgarray_ptr<T> _data)
   {
      startTimer(PerformanceTimers::haloMove2BufTimer);

      ro_array_ptr<int> sendMap = sendMap_.readonly(DEFAULT_COMPUTE_SPACE);
      wo_array_ptr<T> sendBuf = sendBuf_.writeonly(DEFAULT_COMPUTE_SPACE);
      ro_array_ptr<T> data = _data.useOn(DEFAULT_COMPUTE_SPACE);
      DEVICE_PARALLEL_FORALL(sendMap.size(), ii,
                             sendBuf[ii] = data[sendMap[ii]]);

      stopTimer(PerformanceTimers::haloMove2BufTimer);
   };

   virtual void startComm()
   {
      PUSH_RANGE("Halo_Exchange", 1);
      HaloExchange<T>::startComm();
   }
   
   virtual void wait()
   {
      HaloExchange<T>::wait();
      POP_RANGE();
   }
};


//#endif

#endif
