#ifndef HALO_EXCHANGE_HH
#define HALO_EXCHANGE_HH

#include <vector>
#include <cassert>
#include "CommTable.hh"
#include <iostream>

using namespace std;

#ifdef SPI
#include "spi_impl.h"
extern "C"{
uint32_t mapping_table(spi_hdl_t* spi_hdl);
uint32_t mapping_table_new(spi_hdl_t* spi_hdl);
void spi_dump_mapping(spi_hdl_t* spi_hdl);
void setup_bat(spi_hdl_t* spi_hdl,void* recv_buf,uint32_t recv_buf_size);
uint32_t setup_inj_fifo(spi_hdl_t* spi_hdl);
uint32_t setup_descriptors(int** offsets, int* dest, int* putIdx, int put_size , void* send_buf, uint32_t send_buf_size, spi_hdl_t* spi_hdl, int width);
void execute_spi(spi_hdl_t* spi_hdl, uint32_t put_size);
void execute_spi_2(spi_hdl_t* spi_hdl, uint32_t put_size);
void execute_spi_alter(spi_hdl_t* spi_hdl, uint32_t put_size,int bw);
void complete_spi(spi_hdl_t* spi_hdl, uint32_t recv_size);
void complete_spi_alter(spi_hdl_t* spi_hdl, uint32_t recv_size, int* recv_offset, int bw, int width);
void free_spi(spi_hdl_t* spi_hdl);
void global_sync(spi_hdl_t* spi_hdl);
}
#endif

template <class T>
class HaloExchange
{
 public:
   HaloExchange(const std::vector<int>& sendMap, const CommTable* comm)
     : width_(sizeof(T)), commTable_(comm), sendMap_(sendMap)
   {
     //
     send_buf_ = new T[commTable_->sendSize()];
   };

   ~HaloExchange()
   {
     delete [] send_buf_;
   };

   void move2Buf(std::vector<T>& data)
   {
     // fill send buffer
     assert(sendMap_.size() == commTable_->sendSize()); 
     for (unsigned ii=0; ii<sendMap_.size(); ++ii) { send_buf_[ii]=data[sendMap_[ii]]; }
   };

   void set_recv_buf(int bw) { for (unsigned ii=0; ii<commTable_->recvSize();++ii) {recv_buf_[ii+bw*commTable_->recvSize()]=-1;} }

   void dump_recv_buf()
   {
     cout << " recv is : ";
     for (unsigned ii=0; ii<(commTable_->recvSize()) * 2;++ii) {cout << recv_buf_[ii] << " ";}
     cout << endl;
   }
     
   void dump_send_buf()
   {
     cout << " send is : ";
     for (unsigned ii=0; ii<commTable_->sendSize();++ii) {cout << send_buf_[ii] << " ";}
     cout << endl;
   }

   //virtual void execute(std::vector<T>& data, unsigned nLocal);
   //virtual void complete(std::vector<T>& data, unsigned nLocal);
 protected:
   
   unsigned width_;
   const CommTable* commTable_;
   std::vector<int> sendMap_;
   T* send_buf_;
   T* recv_buf_;
};

#ifdef SPI
template <class T>
class spi_HaloExchange : public HaloExchange<T>
{
  public:

  spi_hdl_t spi_hdl;
  int bw;

  public:
   spi_HaloExchange(const std::vector<int>& sendMap, const CommTable* comm)
     : HaloExchange<T>(sendMap,comm), bw(1)
  {
    recv_buf_ = new T[commTable_->recvSize()*2];

    //create mapping table
//    mapping_table(&spi_hdl);
    mapping_table_new(&spi_hdl);
 
    //setup base address table
    //search and allocate
    setup_bat(&spi_hdl,(void*)recv_buf_,commTable_->recvSize()*width_);
 
    //setup injection memory fifo
    //search,allocate,initialize,activate
    setup_inj_fifo(&spi_hdl);

    //setup descriptor
    setup_descriptors(commTable_->_offsets, const_cast<int*>(&(commTable_->_putTask[0])), const_cast<int*>(&(commTable_->_putIdx[0])), commTable_->_putTask.size(), (void*)send_buf_, commTable_->sendSize(), &spi_hdl,width_);

    barrier();

  };
  ~spi_HaloExchange()
   {
     free_spi(&spi_hdl);
     delete [] recv_buf_;
   };

   T* get_recv_buf_() { return &(recv_buf_[bw*commTable_->recvSize()]); };

   void execute(vector<T>& Data,int nLocal) {
     move2Buf(Data);
     bw=1-bw;
     execute_spi_alter(&spi_hdl,commTable_->_putTask.size(),bw);
   };

   void execute()
   {
     bw=1-bw;
     execute_spi_alter(&spi_hdl,commTable_->_putTask.size(),bw);
   };
   void complete() {complete_spi_alter(&spi_hdl, commTable_->_recvTask.size(), commTable_->_offsets[3], bw, width_ );};

   void execute3() {execute_spi(&spi_hdl,commTable_->_putTask.size());};
   void execute2() {execute_spi_2(&spi_hdl,commTable_->_putTask.size());};
   void complete2() {complete_spi(&spi_hdl, commTable_->_recvTask.size());};
   void dump_mapping_table() { spi_dump_mapping( &spi_hdl); };
   void barrier() {global_sync(&spi_hdl);};

};
#endif

template <class T>
class mpi_HaloExchange : public HaloExchange<T>
{
  public:
   using HaloExchange<T>::commTable_;
   using HaloExchange<T>::send_buf_;
   using HaloExchange<T>::recv_buf_;
   using HaloExchange<T>::width_;

   int t1,t2,t3,t4,t5,t6,t7;

   mpi_HaloExchange(const std::vector<int>& sendMap, const CommTable* comm) 
    : HaloExchange<T>(sendMap,comm)
   {t5=t6=t7=0;};

   T* get_recv_buf_() { return this->recv_buf_; }

   void execute(vector<T>& Data,int nLocal)
   {
      move2Buf(Data);

      Data.resize(nLocal + commTable_->recvSize());
      recv_buf_ = (&Data[nLocal]);
      char* sendBuf = (char*)send_buf_;
      char* recvBuf = (char*)recv_buf_;

   
      const int tag = 151515;
      MPI_Request recvReq[commTable_->recvSize()];
      MPI_Request sendReq[commTable_->sendSize()];

      for (unsigned ii=0; ii< commTable_->_recvTask.size(); ++ii)
      {
         assert(recvBuf);
         unsigned sender = commTable_->_recvTask[ii];
         unsigned nItems = commTable_->_recvOffset[ii+1] - commTable_->_recvOffset[ii];
         unsigned len = nItems * width_;
         char* recvPtr = recvBuf + commTable_->_recvOffset[ii]*width_;
         MPI_Irecv(recvPtr, len, MPI_CHAR, sender, tag, commTable_->_comm, recvReq+ii);
      }
 
      for (unsigned ii=0; ii<commTable_->_sendTask.size(); ++ii)
      {
         assert(sendBuf);
         unsigned target = commTable_->_sendTask[ii];
         unsigned nItems = commTable_->_sendOffset[ii+1] - commTable_->_sendOffset[ii];
         unsigned len = nItems * width_;
         char* sendPtr = sendBuf + commTable_->_sendOffset[ii]*width_;
         MPI_Isend(sendPtr, len, MPI_CHAR, target, tag, commTable_->_comm, sendReq+ii);
      }
 
      MPI_Waitall(commTable_->_sendTask.size(), sendReq, MPI_STATUS_IGNORE);
      MPI_Waitall(commTable_->_recvTask.size(), recvReq, MPI_STATUS_IGNORE);
   };
   void complete() {};

   void set_recv_buf(vector<T>& Data,int nLocal)
   {  
      Data.resize(nLocal + commTable_->recvSize());
      recv_buf_ = (&Data[nLocal]);
   }

   void execute2()
   {
      char* sendBuf = (char*)send_buf_;
      char* recvBuf = (char*)recv_buf_;

      const int tag = 151515;
      MPI_Request recvReq[commTable_->recvSize()];
      MPI_Request sendReq[commTable_->sendSize()];

      timebase(t1);
      for (unsigned ii=0; ii< commTable_->_recvTask.size(); ++ii)
      {
         assert(recvBuf);
         unsigned sender = commTable_->_recvTask[ii];
         unsigned nItems = commTable_->_recvOffset[ii+1] - commTable_->_recvOffset[ii];
         unsigned len = nItems * width_;
         char* recvPtr = recvBuf + commTable_->_recvOffset[ii]*width_;
         MPI_Irecv(recvPtr, len, MPI_CHAR, sender, tag, commTable_->_comm, recvReq+ii);
      }
 
      timebase(t2);
      for (unsigned ii=0; ii<commTable_->_sendTask.size(); ++ii)
      {
         assert(sendBuf);
         unsigned target = commTable_->_sendTask[ii];
         unsigned nItems = commTable_->_sendOffset[ii+1] - commTable_->_sendOffset[ii];
         unsigned len = nItems * width_;
         char* sendPtr = sendBuf + commTable_->_sendOffset[ii]*width_;
         MPI_Isend(sendPtr, len, MPI_CHAR, target, tag, commTable_->_comm, sendReq+ii);
      }
 
      timebase(t3);
      MPI_Waitall(commTable_->_sendTask.size(), sendReq, MPI_STATUS_IGNORE);
      MPI_Waitall(commTable_->_recvTask.size(), recvReq, MPI_STATUS_IGNORE);
      timebase(t4);

      t5+=t2-t1;
      t6+=t3-t2;
      t7+=t4-t3;
   };
};

#endif
