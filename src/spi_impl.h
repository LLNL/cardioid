#ifndef SPI_IMPL_H
#define SPI_IMPL_H

#include <inttypes.h>

typedef struct spi_hdl_s
{
    void* inj_fifo_hdl; //fifo_pointer
    void* inj_fifo_subgrp; // fifo_subgroup
    void* inj_fifo ;  //fifo_addr
    void* recv_cnt;
    void* bat_hdl;
    void* mapping_hdl;
    void* cr;
    uint32_t giID;
    uint32_t nNodes;
    uint32_t n_bat_ids;
    uint32_t bat_subgrp_id;
    uint32_t n_fifo_ids;
    uint32_t free_bat_id[256];
    uint32_t free_fifo_id[256];
    uint64_t send_buf_base_pa;
    void* mem_region_hdl[4];  //recv_buf, recv_cnt, inj_fifo, send_buf
    void* barrier_hdl;
    void* barrier_hdl_B;
  double* recvBuf0;
  double* recvBuf1;
  unsigned recvBufSize;
} spi_hdl_t;


#ifdef __cplusplus
extern "C" {
#endif
   uint32_t mapping_table(spi_hdl_t* spi_hdl);
   uint32_t mapping_table_new(spi_hdl_t* spi_hdl);
   void spi_dump_mapping(spi_hdl_t* spi_hdl);
   void setup_bat(spi_hdl_t* spi_hdl,void* recv_buf,uint32_t recv_buf_size);
   void setup_GI(spi_hdl_t* spi_hdl);
   uint32_t setup_inj_fifo(spi_hdl_t* spi_hdl);
   uint32_t setup_descriptors(int** offsets, const int* dest, const int* putIdx, int put_size , void* send_buf, uint32_t send_buf_size, spi_hdl_t* spi_hdl, int width);
   void execute_spi(spi_hdl_t* spi_hdl, uint32_t put_size);
   void execute_spi_2(spi_hdl_t* spi_hdl, uint32_t put_size);
   void execute_spi_alter(spi_hdl_t* spi_hdl, uint32_t put_size,int bw);
   void complete_spi(spi_hdl_t* spi_hdl, uint32_t recv_size);
   void complete_spi_alter(spi_hdl_t* spi_hdl, uint32_t recv_size, int32_t* recv_offset, int bw, int width);
   void complete_spi_alter_monitor(spi_hdl_t* spi_hdl, uint32_t recv_size, int32_t* recv_offset, int32_t* recv_task,int bw, int width, uint32_t myID);
   void free_spi(spi_hdl_t* spi_hdl);
   void global_sync(spi_hdl_t* spi_hdl);
   void global_sync_A(spi_hdl_t* spi_hdl);
   void global_sync_2(spi_hdl_t* spi_hdl,uint64_t timeout);

#ifdef __cplusplus
}
#endif


#endif
