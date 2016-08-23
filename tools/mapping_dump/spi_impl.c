
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef SPI
#include <spi/include/kernel/gi.h>
#include <spi/include/kernel/location.h>
#include <spi/include/kernel/memory.h>
#include <spi/include/kernel/MU.h>
#include <spi/include/kernel/collective.h>
#include <spi/include/mu/GIBarrier.h>
#include <spi/include/mu/RecFifo.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/Descriptor_inlines.h>
#include <spi/include/mu/Descriptor.h>
#include <spi/include/mu/Addressing_inlines.h>
#include <spi/include/mu/Classroute.h>
#include <hwi/include/bqc/classroute.h>
#include <spi/include/l2/lock.h>
#include <hwi/include/bqc/dcr_support.h>
#include <hwi/include/bqc/nd_500_dcr.h>
#include "node_personality.h"
#include "spi_impl.h"

/* pt2pt constants */
#define PT2PT_MISC1                                     \
  ( MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE   |            \
    MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |                    \
    MUHWI_PACKET_DO_NOT_DEPOSIT )                                               

#define PT2PT_MISC2 MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC

/* adaptive routing may affect algorithms that check for the put
 *    completion by checking a tag in the final part of the message */
/* #define PT2PT_MISC1                                  \ */
/*   ( MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE   |         \ */
/*     MUHWI_PACKET_USE_DYNAMIC_ROUTING |                       \ */
/*     MUHWI_PACKET_DO_NOT_DEPOSIT )                                            */

/* #define PT2PT_MISC2 MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC */


#define PT2PT_FIFO_MAP                          \
( MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM |          \
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP )

//#define spi_debug_fifo
//#define spi_debug_bat
//#define spi_debug_desc
//#define spi_debug_exec
//#define spi_debug_compl
#define MAX_WAIT 1600000000
#define MAX_PUT_NUM 128
#define SizeInjMemFifo 128


uint32_t compute_node_id( Personality_t* personality );
uint32_t is_power_of_two( uint32_t value );
void* aligned_malloc( size_t size, size_t align_size ); 
void aligned_free( void *ptr ); 
uint32_t int_log2( uint32_t val );
uint32_t mapping_table(spi_hdl_t* spi_hdl);
uint32_t mapping_table_new(spi_hdl_t* spi_hdl);


uint32_t compute_node_id( Personality_t* personality )
{
  uint32_t id;

  Personality_Networks_t *net = &(personality->Network_Config);
  id = ((((((( (net->Acoord * net->Bnodes) + net->Bcoord)
	     * net->Cnodes) + net->Ccoord)
	   * net->Dnodes) + net->Dcoord)
	 * net->Enodes) + net->Ecoord);

  return id;
}

void spi_dump_mapping(spi_hdl_t* spi_hdl)
{
  MUHWI_Destination_t* mapping = (MUHWI_Destination_t*)(spi_hdl->mapping_hdl);
  for(int id=0;id<32;id++)
    printf("id=%d, %d,%d,%d,%d,%d\n",id, (mapping)[id].Destination.A_Destination, (mapping)[id].Destination.B_Destination, (mapping)[id].Destination.C_Destination, (mapping)[id].Destination.D_Destination, (mapping)[id].Destination.E_Destination);
};

void setup_bat(spi_hdl_t* spi_hdl,void* recv_buf,uint32_t recv_buf_size)
{
  uint32_t grp;
  uint32_t nids;
  void* recv_cnt;
  uint32_t* free_ids=spi_hdl->free_bat_id;
  //uint32_t free_ids[256];
  uint32_t ii,size_in_bytes,rc;

  spi_hdl->barrier_hdl = malloc(sizeof(MUSPI_GIBarrier_t));
  MUSPI_GIBarrierInit( (MUSPI_GIBarrier_t*)(spi_hdl->barrier_hdl), 0 ); //ready global barrier
  Kernel_MemoryRegion_t* recv_buf_mem=(Kernel_MemoryRegion_t*)malloc(sizeof(Kernel_MemoryRegion_t));
  Kernel_MemoryRegion_t* recv_cnt_mem=(Kernel_MemoryRegion_t*)malloc(sizeof(Kernel_MemoryRegion_t));
  MUSPI_BaseAddressTableSubGroup_t* bat=(MUSPI_BaseAddressTableSubGroup_t*)malloc(sizeof(MUSPI_BaseAddressTableSubGroup_t));

  
  nids=0;
  for(grp=0;grp<BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE && nids==0;grp++) 
  {
    rc = Kernel_QueryBaseAddressTable (grp, &nids, free_ids);
    assert(rc == 0);
  }
  grp--;
  assert(nids > 3 );
  {
  #ifdef spi_debug_bat
    printf("Slots in BAT are found in group %d.\n Those are ",grp);
    for(ii=0;ii<nids;ii++) printf("%d ",free_ids[ii]);
    printf("\n");
  #endif

    rc = Kernel_AllocateBaseAddressTable( grp, bat, 3, free_ids, 0 );
    assert(rc==0);

    //in order to setup the base address table, we need allocate memory first
    size_in_bytes = recv_buf_size*2;
    rc = Kernel_CreateMemoryRegion( recv_buf_mem, (void*) recv_buf, size_in_bytes );
    assert( rc == 0 );

  #ifdef spi_debug_bat
    printf("memory region info:\n");
    printf("physical base address %p :\n",recv_buf_mem->BasePa);
    printf("Virtural base address %p :\n",recv_buf_mem->BaseVa);
    printf("Pointer address %p :\n",recv_buf);
    printf("size : %x \n",recv_buf_mem->Bytes);
  #endif
 
    rc = MUSPI_SetBaseAddress( bat, free_ids[0], (uint64_t)recv_buf_mem->BasePa + (uint64_t)recv_buf - (uint64_t)recv_buf_mem->BaseVa);
    assert(rc==0);
    rc = MUSPI_SetBaseAddress( bat, free_ids[1], (uint64_t)recv_buf_mem->BasePa + (uint64_t)recv_buf - (uint64_t)recv_buf_mem->BaseVa + recv_buf_size);
    assert(rc==0);
 
 
    size_in_bytes = sizeof(uint64_t)*MAX_PUT_NUM*2;
    recv_cnt = malloc(size_in_bytes);
    memset(recv_cnt,0x0,size_in_bytes);
    rc = Kernel_CreateMemoryRegion( recv_cnt_mem, recv_cnt, size_in_bytes );
    assert(rc==0);

  #ifdef spi_debug_bat
    printf("memory region info:\n");
    printf("cnt physical base address %p :\n",recv_cnt_mem->BasePa);
    printf("cnt Virtural base address %p :\n",recv_cnt_mem->BaseVa);
    printf("cnt Pointer address %p :\n",recv_cnt);
    printf("cnt size : %x \n",recv_cnt_mem->Bytes);
  #endif
 
    //rc = MUSPI_SetBaseAddress( bat, free_ids[1], (uint64_t)recv_cnt_mem->BasePa + (uint64_t)recv_cnt - (uint64_t)recv_cnt_mem->BaseVa );
    rc = MUSPI_SetBaseAddress( bat, free_ids[2], MUSPI_GetAtomicAddress((uint64_t)recv_cnt_mem->BasePa + (uint64_t)recv_cnt - (uint64_t)recv_cnt_mem->BaseVa ,MUHWI_ATOMIC_OPCODE_STORE_ADD));
    assert(rc==0);

 
  #ifdef spi_debug_bat
    printf("BAT physical physical base address of %d = %x \n", free_ids[0], MUSPI_GetBaseAddress( bat, free_ids[0] ));
    printf("BAT physical base address of %d= %x \n",free_ids[1], MUSPI_GetBaseAddress( bat, free_ids[1] ));
  #endif
  }
  spi_hdl -> n_bat_ids = nids;
  spi_hdl -> bat_subgrp_id = grp;
  spi_hdl -> recv_cnt = recv_cnt;
  spi_hdl->bat_hdl=(void*)bat;
  spi_hdl->mem_region_hdl[0]=(void*) recv_buf_mem;
  spi_hdl->mem_region_hdl[1]=(void*) recv_cnt_mem;
  _bgq_msync();
}

uint32_t setup_inj_fifo(spi_hdl_t* spi_hdl)
{
  uint32_t grp;
  uint32_t nids;
  uint32_t* free_ids=spi_hdl->free_fifo_id;
  uint32_t ii,rc;

  MUSPI_InjFifo_t*    IF;
  MUSPI_InjFifoSubGroup_t* IF_subgroup=(MUSPI_InjFifoSubGroup_t*)malloc(sizeof(MUSPI_InjFifoSubGroup_t));
  Kernel_MemoryRegion_t* inj_fifo_mem=(Kernel_MemoryRegion_t*)malloc(sizeof(Kernel_MemoryRegion_t));


  grp=-1;
  nids=0;
  #ifdef spi_debug_fifo
  printf("searching free inj hdw fifos ...");
  #endif
  for(grp=0;((grp<BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE) && (nids==0));grp++) 
  {
    //printf("%d ",grp);
    rc = Kernel_QueryInjFifos( grp , &nids, free_ids );
    assert(rc==0);
  }
  #ifdef spi_debug_fifo
  printf("\n");
  #endif
  grp--;

  if ( nids == 0 ) { printf("No free injFifo is found. Please use \"export MUSPI_NUMINJFIFOS=2\".  \n"); }
  else
  {

  #ifdef spi_debug_fifo
    printf("Slots in injFifo are found.\n Those are ");
    for(ii=0;ii<nids;ii++) printf("%d ",free_ids[ii]);
    printf("\n");
  #endif

    // reset the fifo attributes 
    Kernel_InjFifoAttributes_t inj_fifo_attrs[BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP];
    memset( inj_fifo_attrs, 0x0, sizeof( inj_fifo_attrs ) );

  #ifdef spi_debug_fifo
    printf("Need only one. Allocating id=%d subgroup=%d\n",free_ids[0],grp);
  #endif

    rc = Kernel_AllocateInjFifos( grp, IF_subgroup, 1, free_ids, inj_fifo_attrs );
    assert( rc == 0 );

    //in order to setup the base address table, we need allocate memory first
    //size_in_bytes = put_size;
    //rc = Kernel_CreateMemoryRegion( send_buf_mem, (void*) send_buf, size_in_bytes );
    //trace_assert( rc == 0, "recv_buf_mem: failed Kernel_CreateMemoryRegion\n" );
    //printf("memory region info:\n");
    
    void* inj_fifo_addr = aligned_malloc(64*(SizeInjMemFifo+2),64);
    rc = Kernel_CreateMemoryRegion( inj_fifo_mem, (void*) inj_fifo_addr, 64*SizeInjMemFifo );
    assert(rc==0);

  #ifdef spi_debug_fifo
    printf("fifo:the base address is %p\n", inj_fifo_mem->BasePa);
    printf("fifo:the base address is %p\n", inj_fifo_mem->BaseVa);
    printf("fifo:the base address is %p\n", inj_fifo_addr);
  #endif

    rc = Kernel_InjFifoInit( IF_subgroup, free_ids[0], inj_fifo_mem, (uint64_t)inj_fifo_addr - (uint64_t)inj_fifo_mem->BaseVa, 64*SizeInjMemFifo-1);
    assert( rc == 0 );
 
  #ifdef spi_debug_fifo
    printf("Injifo Init done\n");
  #endif

    rc = Kernel_InjFifoActivate( IF_subgroup,1, free_ids, KERNEL_INJ_FIFO_DEACTIVATE );
    assert( rc == 0);

  #ifdef spi_debug_fifo
    printf("Fifo DeActivated. Will be activated by setup_descriptors call\n");
  #endif

    IF = MUSPI_IdToInjFifo( free_ids[0], IF_subgroup );
    IF->hw_injfifo->descCount=0;

  #ifdef spi_debug_fifo
    printf("fifo:the tail value of the fifo is %p\n", MUSPI_getTailVa(&(IF->_fifo)));
    printf("fifo:free space %d\n", IF -> freeSpace);
    printf("fifo:Injected  %d\n", IF -> numInjected);
    printf("fifo:descCount :%d \n", MUSPI_getHwDescCount(IF));
    printf("fifo:HwTail,HwHead: %p,%p\n",MUSPI_getHwTail(&(IF->_fifo)),MUSPI_getHwHead(&(IF->_fifo)));
  #endif

    spi_hdl->n_fifo_ids = nids;
    spi_hdl->inj_fifo_subgrp = (void*)IF_subgroup;
    spi_hdl->inj_fifo_hdl = (void*)IF;
    spi_hdl->mem_region_hdl[2] = (void*)inj_fifo_mem;
    spi_hdl->inj_fifo = inj_fifo_addr;
  }
  _bgq_msync();
}

uint32_t setup_descriptors(int** offsets, int* dest, int* putIdx, int put_size , void* send_buf, uint32_t send_buf_size, spi_hdl_t* spi_hdl, int width)
{

  MUHWI_Descriptor_t* mu_desc;  
  MUSPI_BaseAddressTableSubGroup_t* bat = (MUSPI_BaseAddressTableSubGroup_t*) spi_hdl->bat_hdl;
  MUHWI_Destination_t* mapping = spi_hdl->mapping_hdl;


  uint32_t* free_bat_ids = spi_hdl->free_bat_id;
  uint64_t send_addr;
  uint64_t send_size;
  void* inj_fifo = spi_hdl -> inj_fifo;
  int ii,rc,kk;

  int* send_offset = offsets[0];
  int* put_offset = offsets[1];
  int* put_cnt_offset = offsets[2];

  Kernel_MemoryRegion_t* send_buf_mem=(Kernel_MemoryRegion_t*)malloc(sizeof(Kernel_MemoryRegion_t));
  rc = Kernel_CreateMemoryRegion( send_buf_mem, send_buf, send_buf_size );
  assert(rc==0);
  spi_hdl->mem_region_hdl[3] = (void*)send_buf_mem;
  

  for(int jj=0;jj<2;jj++)
  for(ii=0;ii<put_size;ii++)
  {
    send_addr = (uint64_t)send_buf_mem->BasePa + (uint64_t)send_buf - (uint64_t)send_buf_mem->BaseVa + send_offset[putIdx[ii]]*width;
    send_size = width*(send_offset[putIdx[ii]+1]-send_offset[putIdx[ii]]);
    assert(send_size>0);
    #ifdef spi_debug_desc
    printf("%dth descriptor setting\n",ii);
    printf("send_addr=%x size=%d\n",send_addr,send_size);
    //printf("sending... ");
    //for(int ii=0;ii<send_size/width;ii++) printf("%f ",((double*)send_addr)[ii]);
    //printf("\n");
    #endif

    mu_desc =  (MUHWI_Descriptor_t*)( (char*)inj_fifo + 64*(ii+jj*put_size) );
    //printf("mu_desc address %p\n",mu_desc);
    MUSPI_DescriptorZeroOut ( mu_desc );
    MUSPI_SetPrefetchOnly( mu_desc, MUHWI_DESCRIPTOR_PRE_FETCH_ONLY_NO );
    MUSPI_SetPayload( mu_desc, send_addr  , send_size );
    MUSPI_SetTorusInjectionFIFOMap( mu_desc, PT2PT_FIFO_MAP );
    MUSPI_SetDestination( mu_desc, mapping[dest[ii]] );

    #ifdef spi_debug_desc
    printf("dest=%d\n",dest[ii]);
    #endif
    
    /* set the network packet */
    MUSPI_SetDataPacketType( mu_desc, MUHWI_PT2PT_DATA_PACKET_TYPE );
    MUSPI_SetMessageUnitPacketType( mu_desc, MUHWI_PACKET_TYPE_PUT );
    MUSPI_SetHintsABCD( mu_desc, 0 ); // 0 -> machine determine
    MUSPI_SetPt2PtMisc1( mu_desc, PT2PT_MISC1 );
    MUSPI_SetPt2PtMisc2( mu_desc, PT2PT_MISC2 );
    MUSPI_SetInterrupt( mu_desc, MUHWI_DESCRIPTOR_DO_NOT_INTERRUPT_ON_PACKET_ARRIVAL );
    MUSPI_SetSize( mu_desc );
 
    /* set the MU info */
    MUSPI_SetRecPayloadBaseAddressInfo( mu_desc, 8*(bat->subgrpid) + free_bat_ids[jj], put_offset[ii]*width); //(uint64_t)recv_buf_mem.BasePa + (uint64_t)recv_buf - (uint64_t)recv_buf_mem.BaseVa  );
    //MUSPI_SetRecCounterBaseAddressInfo( mu_desc, 1, MUSPI_GetAtomicOffsetFromBaseAddress ( &bat, 1, (uint64_t)recv_cnt - (uint64_t)recv_cnt_mem.BaseVa + (uint64_t)recv_cnt_mem.BasePa, MUHWI_ATOMIC_OPCODE_STORE_ADD ) );
    MUSPI_SetRecCounterBaseAddressInfo( mu_desc, 8*(bat->subgrpid) + free_bat_ids[2], (put_cnt_offset[ii]+ MAX_PUT_NUM *jj)*8);
    MUSPI_SetPacing( mu_desc, MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED );
  _bgq_msync();
    
  #ifdef spi_debug_desc
    printf("destination subgroup id is %d,%d and the offset is %d,%d\n",8*(bat->subgrpid) + free_bat_ids[jj],8*(bat->subgrpid) + free_bat_ids[2],put_offset[ii],put_cnt_offset[ii]);
    printf("%dth descriptor is ready\n",ii);
    char tmp_str[1024];
    MUSPI_DescriptorDumpHex( tmp_str, mu_desc );
  #endif
  }

  MUSPI_InjFifo_t* IF = (MUSPI_InjFifo_t*)(spi_hdl -> inj_fifo_hdl);
  MUSPI_InjFifoSubGroup_t* IF_subgroup=(MUSPI_InjFifoSubGroup_t*)(spi_hdl -> inj_fifo_subgrp);
//  rc = Kernel_InjFifoActivate( IF_subgroup,1, spi_hdl->free_fifo_id, KERNEL_INJ_FIFO_DEACTIVATE );
//  MUSPI_setHwHead(&(IF->_fifo),put_size*64*2); 
//  MUSPI_setHwTail(&(IF->_fifo),put_size*64*2); 
//  rc = Kernel_InjFifoActivate( IF_subgroup,1, spi_hdl->free_fifo_id, KERNEL_INJ_FIFO_ACTIVATE );
  assert(rc==0);
  #ifdef spi_debug_desc
    printf("All descriptors are ready\n",ii);
  #endif

  _bgq_msync();
}

void global_sync(spi_hdl_t *spi_hdl)
{
  _bgq_msync();
//  MUSPI_GIBarrier_t barrier;
//  MUSPI_GIBarrierInit( &barrier, 0 );
  MUSPI_GIBarrierEnterAndWait((MUSPI_GIBarrier_t*)spi_hdl->barrier_hdl );
}

void execute_spi(spi_hdl_t* spi_hdl, uint32_t put_size)
{
    printf("initiate sending %d messages\n",put_size);
    MUSPI_InjFifo_t* IF = (MUSPI_InjFifo_t*)(spi_hdl -> inj_fifo_hdl);
    MUSPI_InjFifoAdvanceDescMultiple ( IF, put_size );
    printf("initiated :%d \n", IF->numInjected);
    printf("initiated:descCount :%d \n", MUSPI_getHwDescCount(IF));
    printf("initialed:HwTail,HwHead,Start: %p,%p,%p\n",MUSPI_getHwTail(&(IF->_fifo)),MUSPI_getHwHead(&(IF->_fifo)),MUSPI_getHwStart(&(IF->_fifo)));
    //rewind Multiple?
    //Advance Multiple?
}

void execute_spi_alter(spi_hdl_t* spi_hdl, uint32_t put_size,int bw)
{
    int rc;
#ifdef spi_debug_exec
    printf("initiate sending %d messages with bw=%d\n",put_size,bw);
#endif
    MUSPI_InjFifo_t* IF = (MUSPI_InjFifo_t*)(spi_hdl -> inj_fifo_hdl);
    MUSPI_InjFifoSubGroup_t* IF_subgroup=(MUSPI_InjFifoSubGroup_t*)(spi_hdl -> inj_fifo_subgrp);
    rc = Kernel_InjFifoActivate( IF_subgroup,1, spi_hdl->free_fifo_id, KERNEL_INJ_FIFO_DEACTIVATE );
    assert(rc==0);
    _bgq_msync();
#ifdef spi_debug_exec
    printf("InjFifo Deactivated \n");
#endif
    IF->numInjected = put_size;
    IF->hw_injfifo->descCount=0;
#ifdef spi_debug_exec
    printf("before:HwTail,HwHead,Start: %p,%p,%p\n",MUSPI_getHwTail(&(IF->_fifo)),MUSPI_getHwHead(&(IF->_fifo)),MUSPI_getHwStart(&(IF->_fifo)));
#endif
    if(bw == 0 )
    {
      IF->_fifo.hwfifo->pa_head = IF->_fifo.hwfifo->pa_start;
      IF->_fifo.hwfifo->pa_tail = IF->_fifo.hwfifo->pa_start + 64*put_size;
    }
    else
    {
      IF->_fifo.hwfifo->pa_head = IF->_fifo.hwfifo->pa_start + 64*put_size;
      IF->_fifo.hwfifo->pa_tail = IF->_fifo.hwfifo->pa_start + 2*64*put_size;
    }
#ifdef spi_debug_exec
    printf("after:HwTail,HwHead,Start: %p,%p,%p\n",MUSPI_getHwTail(&(IF->_fifo)),MUSPI_getHwHead(&(IF->_fifo)),MUSPI_getHwStart(&(IF->_fifo)));
    printf("InjFifo Activating...");
#endif
    _bgq_msync();
    rc = Kernel_InjFifoActivate( IF_subgroup,1, spi_hdl->free_fifo_id, KERNEL_INJ_FIFO_ACTIVATE );
    assert(rc==0);
#ifdef spi_debug_exec
    printf("Done \n");
#endif

}
void complete_spi_alter(spi_hdl_t* spi_hdl, uint32_t recv_size, uint32_t* recv_offset, int bw, int width)
{
  volatile uint64_t* recv_cnt = spi_hdl -> recv_cnt;
  //checking if fifo is done.
  #ifdef spi_debug_compl
  printf("completion check...");
  #endif
  MUSPI_InjFifo_t*    IF = (MUSPI_InjFifo_t*)spi_hdl->inj_fifo_hdl;
  while(MUSPI_getHwTail(&(IF->_fifo)) != MUSPI_getHwHead(&(IF->_fifo)));
  #ifdef spi_debug_compl
  printf("head=tail done...");
  #endif

//  printf("recv size ");
//  for(int ii=0;ii<recv_size;ii++)
//    printf("%d ",(recv_offset[ii+1]-recv_offset[ii])*width);
//  printf("\n");

  uint64_t knt=0,sum=1;
  for(knt=0;(knt<MAX_WAIT && sum != 0);knt++)
  {
    sum=0;
    for(int ii=0;ii<recv_size;ii++)
    {
      int recved=recv_cnt[ii+bw*MAX_PUT_NUM]+(recv_offset[ii+1]-recv_offset[ii])*width;
//      if(knt%10==0) printf("%d ",recved);
      sum+= (recved == 0 ? 0 : 1);
    }
//    printf("\n");
//    if(bw==0) for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii] == 0 ? 1 : 0);
//    else for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii+MAX_PUT_NUM] == 0 ? 1 : 0);
  }
  assert(knt <  MAX_WAIT); 
//  //for(int ii=0;ii<recv_size;ii++) printf("%x ",recv_cnt[ii+MAX_PUT_NUM*bw]);
  //printf("\n");
  #ifdef spi_debug_compl
  printf("counter touched\n");
  #endif


  #ifdef spi_debug_compl
  //for(int ii=0;ii<recv_size;ii++)
  //{
  //  if(bw==0)printf("%dth packet is at %x\n",ii,recv_cnt[ii]);
  //  else printf("%dth packet is at %x\n",ii,recv_cnt[ii+MAX_PUT_NUM]);
  //}
  #endif

  for(int ii=0;ii<recv_size;ii++)
    recv_cnt[ii+MAX_PUT_NUM*bw]=0;
  _bgq_msync();

}


void execute_spi_2(spi_hdl_t* spi_hdl, uint32_t put_size)
{
    printf("initiate sending %d messages\n",put_size);
    MUSPI_InjFifo_t* IF = (MUSPI_InjFifo_t*)(spi_hdl -> inj_fifo_hdl);
    //MUSPI_InjFifoAdvanceDescMultiple ( IF, put_size );
    IF->numInjected = put_size;
    IF->hw_injfifo->descCount=0;
    IF->_fifo.hwfifo->pa_head = IF->_fifo.hwfifo->pa_start;
     uint64_t tmp=IF->_fifo.hwfifo->pa_head;
    _bgq_msync();

    printf("initiated :%d \n", IF->numInjected);
    printf("initiated:descCount :%d \n", MUSPI_getHwDescCount(IF));
    printf("initialed:HwTail,HwHead,Start: %p,%p,%p\n",MUSPI_getHwTail(&(IF->_fifo)),tmp,MUSPI_getHwStart(&(IF->_fifo)));
    //rewind Multiple?
    //Advance Multiple?
}

void complete_spi(spi_hdl_t* spi_hdl, uint32_t recv_size)
{
  volatile uint64_t* recv_cnt = spi_hdl -> recv_cnt;
  //checking if fifo is done.
  printf("completion check\n");
  MUSPI_InjFifo_t*    IF = (MUSPI_InjFifo_t*)spi_hdl->inj_fifo_hdl;
  //for(int ii=0;(ii<1000000 && MUSPI_getHwTail(&(IF->_fifo)) != MUSPI_getHwHead(&(IF->_fifo)));ii++ ) ;
  while(MUSPI_getHwTail(&(IF->_fifo)) != MUSPI_getHwHead(&(IF->_fifo)));
  printf("complete:Injected :%d \n", IF->numInjected);
  printf("complete:descCount :%d \n", MUSPI_getHwDescCount(IF));

  uint64_t knt=0,sum=1;
  for(knt=0;(knt<MAX_WAIT && sum != 0);knt++)
  {
    sum=0;
    for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii] == 0 ? 1 : 0);
    for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii+MAX_PUT_NUM] == 0 ? 1 : 0);
  }
  assert(knt <  MAX_WAIT); 
  _bgq_msync();


  for(int ii=0;ii<recv_size;ii++)
  {
    printf("%dth packet is at %x\n",ii,recv_cnt[ii]);
    printf("%dth packet is at %x\n",ii,recv_cnt[ii+MAX_PUT_NUM]);
  }

  for(int ii=0;ii<recv_size;ii++)
    recv_cnt[ii]=0;
}

void free_spi(spi_hdl_t* spi_hdl)
{
  uint32_t rc;
  rc = Kernel_InjFifoActivate( (MUSPI_InjFifoSubGroup_t *)(spi_hdl->inj_fifo_subgrp), 1, spi_hdl->free_fifo_id, KERNEL_INJ_FIFO_DEACTIVATE );
  assert(rc==0);
  rc = Kernel_DestroyMemoryRegion((Kernel_MemoryRegion_t*)spi_hdl->mem_region_hdl[0]);
  assert(rc==0);
  rc = Kernel_DestroyMemoryRegion((Kernel_MemoryRegion_t*)spi_hdl->mem_region_hdl[1]);
  assert(rc==0);
  rc = Kernel_DestroyMemoryRegion((Kernel_MemoryRegion_t*)spi_hdl->mem_region_hdl[2]);
  assert(rc==0);
  rc = Kernel_DestroyMemoryRegion((Kernel_MemoryRegion_t*)spi_hdl->mem_region_hdl[3]);
  assert(rc==0);
  rc = Kernel_DeallocateInjFifos ((MUSPI_InjFifoSubGroup_t *)(spi_hdl->inj_fifo_subgrp), 1, spi_hdl->free_fifo_id);
  assert(rc==0);
  rc = Kernel_DeallocateBaseAddressTable( (MUSPI_BaseAddressTableSubGroup_t *)(spi_hdl->bat_hdl), 3, spi_hdl->free_bat_id);
  assert(rc==0);

  free(spi_hdl->recv_cnt);
  aligned_free(spi_hdl->inj_fifo);
  free(spi_hdl->bat_hdl);
  free(spi_hdl->mapping_hdl);
  free(spi_hdl->inj_fifo_subgrp);

}

uint32_t is_power_of_two( uint32_t value )
{
  return ( value & ( ~value + 1 ) ) == value;
}


/* align_size has to be a power of two */
void* aligned_malloc( size_t size, size_t align_size ) 
{
  char *ptr, *ptr2, *aligned_ptr;
  uint32_t align_mask = align_size - 1;

  ptr = (char*) malloc( size + align_size + sizeof( uint32_t ) );
  if ( ptr == NULL ) 
    return ( NULL );

  ptr2 = ptr + sizeof( uint32_t );
  aligned_ptr = ptr2 + (align_size - ( (size_t) ptr2 & align_mask ) );


  ptr2 = aligned_ptr - sizeof( uint32_t );
  *( (uint32_t *) ptr2 ) = ( uint32_t )( aligned_ptr - ptr );

  return aligned_ptr;
}


uint32_t int_log2( uint32_t val )
{
  uint32_t r = 0;
  
  while ( val >>= 1 )
    r++;

  return r;
} 


void aligned_free( void *ptr ) 
{
  uint32_t *ptr2= (uint32_t*) ptr - 1;
  ptr -= *ptr2;
  free(ptr);
}

uint32_t mapping_table(spi_hdl_t* spi_hdl)
{
  uint32_t a, b, c, d, e, id;
  uint32_t rc;
  node_personality node; 

  rc = Kernel_GetPersonality( &node.personality, sizeof( node.personality ) );
  assert(rc == 0);

  node.node_id = compute_node_id( &node.personality );
  
  node.coord[CR_AXIS_A] = node.a = node.personality.Network_Config.Acoord;
  node.coord[CR_AXIS_B] = node.b = node.personality.Network_Config.Bcoord;
  node.coord[CR_AXIS_C] = node.c = node.personality.Network_Config.Ccoord;
  node.coord[CR_AXIS_D] = node.d = node.personality.Network_Config.Dcoord;
  node.coord[CR_AXIS_E] = node.e = node.personality.Network_Config.Ecoord;

  node.dim[CR_AXIS_A] = node.a_nodes = node.personality.Network_Config.Anodes;
  node.dim[CR_AXIS_B] = node.b_nodes = node.personality.Network_Config.Bnodes;
  node.dim[CR_AXIS_C] = node.c_nodes = node.personality.Network_Config.Cnodes;
  node.dim[CR_AXIS_D] = node.d_nodes = node.personality.Network_Config.Dnodes;
  node.dim[CR_AXIS_E] = node.e_nodes = node.personality.Network_Config.Enodes;

  node.nodes =  node.a_nodes * node.b_nodes * node.c_nodes *
    node.d_nodes * node.e_nodes;
 
  /* allocate the mapping table */
  MUHWI_Destination_t* mapping = malloc( node.nodes * sizeof(  MUHWI_Destination_t ) );
  assert(mapping);
  spi_hdl->mapping_hdl = (void*) mapping; 


  for( a = 0; a < node.a_nodes; a++ )
    for( b = 0; b < node.b_nodes; b++ )
      for( c = 0; c < node.c_nodes; c++ )
        for( d = 0; d < node.d_nodes; d++ )
          for( e = 0; e < node.e_nodes; e++ )
            {
              Personality_Networks_t *net = &(node.personality.Network_Config);
              id = ((((((( (a * net->Bnodes) + b)
                         * net->Cnodes) + c )
                       * net->Dnodes) + d )
                     * net->Enodes) + e );

              (mapping)[id].Destination.A_Destination = a;
              (mapping)[id].Destination.B_Destination = b;
              (mapping)[id].Destination.C_Destination = c;
              (mapping)[id].Destination.D_Destination = d;
              (mapping)[id].Destination.E_Destination = e;
            }

  #ifdef spi_debug
  printf( "anodes %u bnodes %u cnodes %u dnodes %u enodes %u\n",
             node.a_nodes, node.b_nodes, node.c_nodes, node.d_nodes,
             node.e_nodes );
  printf("completing init_mapping_table\n" );
  #endif

  return node.node_id;

};

#if 1
uint32_t mapping_table_new(spi_hdl_t* spi_hdl)
{
  uint32_t a, b, c, d, e, id;
  uint32_t rc,ii;
  uint32_t node_id = Kernel_GetRank();
  BG_JobCoords_t sblock;
  Kernel_JobCoords(&sblock); 
  assert(sblock.corner.core == 0);
  
//  node.dim[CR_AXIS_A] = node.a_nodes = sblock.shape.a;
//  node.dim[CR_AXIS_B] = node.b_nodes = sblock.shape.b;
//  node.dim[CR_AXIS_C] = node.c_nodes = sblock.shape.c;
//  node.dim[CR_AXIS_D] = node.d_nodes = sblock.shape.d;
//  node.dim[CR_AXIS_E] = node.e_nodes = sblock.shape.e;

  uint64_t nodes = sblock.shape.a*sblock.shape.b*sblock.shape.c*sblock.shape.d*sblock.shape.e;
  /* allocate the mapping table */
  MUHWI_Destination_t* mapping = malloc( nodes * sizeof(  MUHWI_Destination_t ) );
  assert(mapping);
  spi_hdl->mapping_hdl = (void*) mapping; 

  /* allocate BG_CoordinateMapping */
  BG_CoordinateMapping_t* bg_map = malloc( nodes * sizeof(BG_CoordinateMapping_t));
  assert(bg_map);
  uint64_t totN;
  Kernel_RanksToCoords(nodes * sizeof(BG_CoordinateMapping_t),bg_map,&totN);
  printf("my node id=%d\n",node_id);
  if(node_id==0) printf("total nodes in this block :%d\n",nodes);
  if(node_id==0) printf("number of node being used :%d\n",totN);
  if(node_id==0) printf("the corner node is %d,%d,%d,%d,%d\n", sblock.corner.a,sblock.corner.b,sblock.corner.c,sblock.corner.d,sblock.corner.e);
//  assert(totN == nodes);

  for( ii =0 ; ii < nodes ; ii++)
  {
    (mapping)[ii].Destination.A_Destination = bg_map[ii].a + sblock.corner.a;
    (mapping)[ii].Destination.B_Destination = bg_map[ii].b + sblock.corner.b;
    (mapping)[ii].Destination.C_Destination = bg_map[ii].c + sblock.corner.c;
    (mapping)[ii].Destination.D_Destination = bg_map[ii].d + sblock.corner.d;
    (mapping)[ii].Destination.E_Destination = bg_map[ii].e + sblock.corner.e;
  }


  #ifdef spi_debug
  printf("my id is %d\n",node_id);
  printf("geom is %d x %d x %d x %d x %d\n", sblock.shape.a,sblock.shape.b,sblock.shape.c,sblock.shape.d,sblock.shape.e);
  printf("my rel coord is %d x %d x %d x %d x %d\n", bg_map[node_id].a,bg_map[node_id].b,bg_map[node_id].c,bg_map[node_id].d,bg_map[node_id].e);

  printf("completing init_mapping_table\n" );
  #endif

  free(bg_map);
  return node_id;

};


uint32_t coord_list(int *myid, int *mytotN, int* mapping)
{
  uint32_t a, b, c, d, e, id;
  uint32_t rc,ii;
  uint32_t node_id = Kernel_GetRank();
  *myid = node_id;
  BG_JobCoords_t sblock;
  Kernel_JobCoords(&sblock); 
  assert(sblock.corner.core == 0);
  
  uint64_t nodes = sblock.shape.a*sblock.shape.b*sblock.shape.c*sblock.shape.d*sblock.shape.e;
  uint64_t totN;
  /* allocate the mapping table */
  assert(mapping);

  /* allocate BG_CoordinateMapping */
  BG_CoordinateMapping_t* bg_map = malloc( nodes * sizeof(BG_CoordinateMapping_t));
  assert(bg_map);
  Kernel_RanksToCoords(nodes * sizeof(BG_CoordinateMapping_t),bg_map,&totN);
  *mytotN = totN;
  if(node_id==0) printf("total nodes in this block :%d\n",nodes);
  if(node_id==0) printf("number of node being used :%d\n",totN);
  if(node_id==0) printf("the corner node is %d,%d,%d,%d,%d\n", sblock.corner.a,sblock.corner.b,sblock.corner.c,sblock.corner.d,sblock.corner.e);
//  assert(totN == nodes);

  for( ii =0 ; ii < totN ; ii++)
  {
    mapping[ii*5+0] = bg_map[ii].a + sblock.corner.a;
    mapping[ii*5+1] = bg_map[ii].b + sblock.corner.b;
    mapping[ii*5+2] = bg_map[ii].c + sblock.corner.c;
    mapping[ii*5+3] = bg_map[ii].d + sblock.corner.d;
    mapping[ii*5+4] = bg_map[ii].e + sblock.corner.e;
  }


  #ifdef spi_debug
  printf("my id is %d\n",node_id);
  printf("geom is %d x %d x %d x %d x %d\n", sblock.shape.a,sblock.shape.b,sblock.shape.c,sblock.shape.d,sblock.shape.e);
  printf("my rel coord is %d x %d x %d x %d x %d\n", bg_map[node_id].a,bg_map[node_id].b,bg_map[node_id].c,bg_map[node_id].d,bg_map[node_id].e);

  printf("completing init_mapping_table\n" );
  #endif

  free(bg_map);
  return node_id;

};
#endif
//MUSPI_GIBarrier_t barrier;
  //MUSPI_GIBarrierInit( &barrier, 0 );
  //MUSPI_GIBarrierEnterAndWait( &barrier );
  //printf("GI Barrier done\n");
#endif
