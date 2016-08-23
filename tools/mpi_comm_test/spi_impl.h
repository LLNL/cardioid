typedef struct spi_hdl_s {
    void* inj_fifo_hdl; //fifo_pointer
    void* inj_fifo_subgrp; // fifo_subgroup
    void* inj_fifo ;  //fifo_addr
    void* recv_cnt;
    void* bat_hdl;
    void* mapping_hdl;
    uint32_t n_bat_ids;
    uint32_t bat_subgrp_id;
    uint32_t n_fifo_ids;
    uint32_t free_bat_id[256];
    uint32_t free_fifo_id[256];
    uint64_t send_buf_base_pa;
    void* mem_region_hdl[4];  //recv_buf, recv_cnt, inj_fifo, send_buf
    void* barrier_hdl;
  } spi_hdl_t;
