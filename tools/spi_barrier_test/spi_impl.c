#include "spi_impl.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef SPI
#define NEW_GI

// for some reason (maybe some include file calculates
// something in an assert?!), compiling these include files with -DNDEBUG
// leads to a crash. So let's undef it here and redfine it later
// if it was defined 
// jlf 8/2/2012
#ifdef NDEBUG
#undef NDEBUG
#define NNDEBUG
#endif

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
#include <hwi/include/bqc/classroute.h>

#ifdef NNDEBUG
#define NDEBUG
#endif

#include "node_personality.h"

static unsigned nCompleteCalls = 0;
static unsigned nExecuteCalls = 0;
static int staticId;

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
#define MAX_WAIT   1600000000
#define MAX_WAIT_SE 200000000
#define MAX_RECV_NUM 128
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

void setup_GI(spi_hdl_t* spi_hdl)
{
  uint32_t nClassRoutes;
  uint32_t classRouteIds[256];
  size_t   sizeofClassRouteIds;
  uint32_t ii,rc;

  if (spi_hdl->nNodes==1) 
  {
    printf("Single node run\n");
    return;
  }

  rc=Kernel_QueryGlobalInterruptClassRoutes ( &nClassRoutes, classRouteIds, 256*4 );
  //printf("# of free Id:%d\n",nClassRoutes);
  assert(nClassRoutes>0);
  assert(classRouteIds[0]==1);
  //printf("free Ids are ");
//  for(ii=0;ii<nClassRoutes;ii++) printf("%d ",classRouteIds[ii]);
  //printf("\n");

  Kernel_AllocateGlobalInterruptClassRoute( classRouteIds[0], NULL);
  spi_hdl->giID = classRouteIds[0];

  Kernel_SetGlobalInterruptClassRoute (spi_hdl->giID, (spi_hdl->cr));
  Kernel_CheckGlobalInterruptClassRoute (spi_hdl->giID, (spi_hdl->cr));

  //copy form MUSPI_GIBarrierInitMU
  global_sync_A(spi_hdl);

  /* Step 2: Reset the control register to the initial state */
  rc = MUSPI_GIBarrierInitMU1(spi_hdl->giID);
  assert(rc==0);

  global_sync_A(spi_hdl);

  rc = MUSPI_GIBarrierInitMU2(spi_hdl->giID,MAX_WAIT);

  assert(rc==0);

  global_sync_A(spi_hdl);

  spi_hdl->barrier_hdl_B = malloc(sizeof(MUSPI_GIBarrier_t));

  MUSPI_GIBarrierInit( (MUSPI_GIBarrier_t*)(spi_hdl->barrier_hdl_B), spi_hdl->giID ); 
}

void setup_bat(spi_hdl_t* spi_hdl,void* recv_buf,uint32_t recv_buf_size)
{
   nCompleteCalls = 0;
   nExecuteCalls = 0;

  uint32_t grp;
  uint32_t nids;
  void* recv_cnt;
  uint32_t* free_ids=spi_hdl->free_bat_id;
  //uint32_t free_ids[256];
  uint32_t ii,size_in_bytes,rc;
  uint32_t classRouteId;

  spi_hdl->recvBuf0 = recv_buf;
  spi_hdl->recvBuf1 = (double*)(((char*) recv_buf) + recv_buf_size);

  spi_hdl->recvBufSize = recv_buf_size;

  Kernel_GetGlobalBarrierUserClassRouteId ( &classRouteId );
  spi_hdl->barrier_hdl = malloc(sizeof(MUSPI_GIBarrier_t));
  MUSPI_GIBarrierInit( (MUSPI_GIBarrier_t*)(spi_hdl->barrier_hdl), classRouteId ); 

  //printf("Barrier Init status:%x\n",((MUSPI_GIBarrier_t*)(spi_hdl->barrier_hdl))->state);
  //printf("Barrier ID :%d\n",classRouteId);

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
 
 
    size_in_bytes = sizeof(uint64_t)*MAX_RECV_NUM*2;
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

  #ifdef NEW_GI
  setup_GI(spi_hdl);
  #endif
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

uint32_t setup_descriptors(int** offsets, const int* dest, const int* putIdx, int put_size , void* send_buf, uint32_t send_buf_size, spi_hdl_t* spi_hdl, int width)
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
    MUSPI_SetRecCounterBaseAddressInfo( mu_desc, 8*(bat->subgrpid) + free_bat_ids[2], (put_cnt_offset[ii]+ MAX_RECV_NUM *jj)*8);
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

void global_sync_A(spi_hdl_t *spi_hdl)
{
  _bgq_msync();
//  MUSPI_GIBarrier_t barrier;
//  MUSPI_GIBarrierInit( &barrier, 0 );
  MUSPI_GIBarrierEnterAndWait((MUSPI_GIBarrier_t*)spi_hdl->barrier_hdl );
}

void global_sync(spi_hdl_t *spi_hdl)
{
  _bgq_msync();
//  MUSPI_GIBarrier_t barrier;
//  MUSPI_GIBarrierInit( &barrier, 0 );
  if(spi_hdl->nNodes == 1)
  {
    printf("global sync:single node run\n");
    return;
  }
  #ifdef NEW_GI
  MUSPI_GIBarrierEnterAndWait((MUSPI_GIBarrier_t*)spi_hdl->barrier_hdl_B );
  #else
  MUSPI_GIBarrierEnterAndWait((MUSPI_GIBarrier_t*)spi_hdl->barrier_hdl );
  #endif
}

void global_sync_2(spi_hdl_t *spi_hdl,uint64_t timeout)
{
  _bgq_msync();
  //MUSPI_GIBarrier_t barrier;
  //  MUSPI_GIBarrierInit( &barrier, 0 );
    int32_t rc= MUSPI_GIBarrierEnterAndWaitWithTimeout((MUSPI_GIBarrier_t*)spi_hdl->barrier_hdl,timeout);
      if( rc!=0 )
        {
            printf("barrier failed. It could be timeout\n");
                assert(rc==0);
        }
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
  assert(bw == nExecuteCalls%2);
  ++nExecuteCalls;
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
void complete_spi_alter(spi_hdl_t* spi_hdl, uint32_t recv_size, int32_t* recv_offset, int bw, int width)
{
  ++nCompleteCalls;
  volatile uint64_t* recv_cnt = spi_hdl -> recv_cnt;
  assert(recv_size < MAX_RECV_NUM ); 
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

  volatile uint64_t sum=1;
  uint64_t knt=0;
  uint64_t maxWait = MAX_WAIT/recv_size;
  for(knt=0;(knt<maxWait && sum != 0);knt++)
  {
    sum=0;
    for(int ii=0;ii<recv_size;ii++)
    {
      volatile int recved=recv_cnt[ii+bw*MAX_RECV_NUM]+(recv_offset[ii+1]-recv_offset[ii])*width;
//      if(knt%10==0) printf("%d ",recved);
      sum+= (recved == 0 ? 0 : 1);
    }
//    printf("\n");
//    if(bw==0) for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii] == 0 ? 1 : 0);
//    else for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii+MAX_RECV_NUM] == 0 ? 1 : 0);
  }
  assert(knt <  maxWait); 
//  //for(int ii=0;ii<recv_size;ii++) printf("%x ",recv_cnt[ii+MAX_RECV_NUM*bw]);
  //printf("\n");
  #ifdef spi_debug_compl
  printf("counter touched\n");
  #endif


  #ifdef spi_debug_compl
  //for(int ii=0;ii<recv_size;ii++)
  //{
  //  if(bw==0)printf("%dth packet is at %x\n",ii,recv_cnt[ii]);
  //  else printf("%dth packet is at %x\n",ii,recv_cnt[ii+MAX_RECV_NUM]);
  //}
  #endif

  for(int ii=0;ii<recv_size;ii++)
    recv_cnt[ii+MAX_RECV_NUM*bw]=0;
  _bgq_msync();

}

void complete_spi_alter_monitor(spi_hdl_t* spi_hdl, uint32_t recv_size, int32_t* recv_offset, int32_t* recv_task,int bw, int width, uint32_t myID)
{
  staticId = myID;
  assert(bw == nCompleteCalls%2);
  ++nCompleteCalls;

  volatile uint64_t* recv_cnt = spi_hdl -> recv_cnt;
  uint64_t knt=0;
  //checking if fifo is done.
  #ifdef spi_debug_compl
  printf("completion check...");
  #endif
  MUSPI_InjFifo_t*    IF = (MUSPI_InjFifo_t*)spi_hdl->inj_fifo_hdl;
  for(knt=0; knt<MAX_WAIT_SE &&  (MUSPI_getHwTail(&(IF->_fifo)) != MUSPI_getHwHead(&(IF->_fifo))); knt++);
  assert(knt<MAX_WAIT_SE);
  
  #ifdef spi_debug_compl
  printf("head=tail done...");
  #endif

//  printf("recv size ");
//  for(int ii=0;ii<recv_size;ii++)
//    printf("%d ",(recv_offset[ii+1]-recv_offset[ii])*width);
//  printf("\n");

  volatile uint64_t sum=1;
  uint64_t maxWait = recv_size==0 ?  100 : MAX_WAIT/recv_size;
  for(knt=0;(knt<maxWait && sum != 0);knt++)
  {
    sum=0;
    for(int ii=0;ii<recv_size;ii++)
    {
      volatile int recved=recv_cnt[ii+bw*MAX_RECV_NUM]+(recv_offset[ii+1]-recv_offset[ii])*width;
//      if(knt%10==0) printf("%d ",recved);
      sum+= (recved == 0 ? 0 : 1);
    }
//    printf("\n");
//    if(bw==0) for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii] == 0 ? 1 : 0);
//    else for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii+MAX_RECV_NUM] == 0 ? 1 : 0);
  }

  if ( knt == maxWait )
  {
    printf("node:%d have reached MAX_WAIT. startCalls = %u waitCalls = %u\n",myID, nExecuteCalls, nCompleteCalls);
    printf("node:%d expecting %d messages\n", myID, recv_size);
    for(int ii=0;ii<recv_size;ii++)
    {
      volatile int recved=recv_cnt[ii+bw*MAX_RECV_NUM]+(recv_offset[ii+1]-recv_offset[ii])*width;
//      if(knt%10==0) printf("%d ",recved);
      if(recved != 0)
      {
	printf("node:%d did NOT receive %d bytes from node:%d on step %u while expecting %d\n",myID,recved,recv_task[ii],nCompleteCalls, (recv_offset[ii+1]-recv_offset[ii])*width); 
/* 	{ */
/* 	  double* rBuf = spi_hdl->recvBuf0; */
/* 	  if (bw = 1) rBuf = spi_hdl->recvBuf1; */
/* 	  printf("node:%d recvBuf[%d] = %f\n", myID, recv_offset[ii], rBuf[recv_offset[ii]]); */
/* 	  for (int jj = recv_offset[ii]; jj<recv_offset[ii+1]; ++jj) */
/* 	    if (rBuf[jj] != rBuf[recv_offset[ii]]) */
/* 		printf("node:%d Inconsistent recv buffer at %d\n", myID, jj); */
/* 	} */
/* 	{ */
/* 	  double* rBuf = spi_hdl->recvBuf0; */
/* 	  if (bw = 0) rBuf = spi_hdl->recvBuf1; */
/* 	  printf("node:%d recvBufOther[%d] = %f\n", myID, recv_offset[ii], rBuf[recv_offset[ii]]); */
/* 	  for (int jj = recv_offset[ii]; jj<recv_offset[ii+1]; ++jj) */
/* 	    if (rBuf[jj] != rBuf[recv_offset[ii]]) */
/* 		printf("node:%d Inconsistent other recv buffer at %d\n", myID, jj); */
/* 	} */
	
         MUHWI_Destination_t* mapping = spi_hdl->mapping_hdl;
         printf("node:%d=(%d,%d,%d,%d,%d)\n",myID,(mapping)[myID].Destination.A_Destination,(mapping)[myID].Destination.B_Destination,(mapping)[myID].Destination.C_Destination,(mapping)[myID].Destination.D_Destination,(mapping)[myID].Destination.E_Destination);
         printf("node:%d=(%d,%d,%d,%d,%d)\n",recv_task[ii],(mapping)[recv_task[ii]].Destination.A_Destination,(mapping)[recv_task[ii]].Destination.B_Destination,(mapping)[recv_task[ii]].Destination.C_Destination,(mapping)[recv_task[ii]].Destination.D_Destination,(mapping)[recv_task[ii]].Destination.E_Destination);
     } 
     else printf("node:%d did receive all the packets from node:%d\n",myID,recv_task[ii]);
    }
    //wait more for others to quit
    for(knt=0;knt<maxWait ;knt++) { }
  }
  if (knt >= maxWait)
     //while(1){};
     sleep(1800);  // die after 30 minutes
  assert(knt <  maxWait); 


//  //for(int ii=0;ii<recv_size;ii++) printf("%x ",recv_cnt[ii+MAX_RECV_NUM*bw]);
  //printf("\n");
  #ifdef spi_debug_compl
  printf("counter touched\n");
  #endif


  #ifdef spi_debug_compl
  //for(int ii=0;ii<recv_size;ii++)
  //{
  //  if(bw==0)printf("%dth packet is at %x\n",ii,recv_cnt[ii]);
  //  else printf("%dth packet is at %x\n",ii,recv_cnt[ii+MAX_RECV_NUM]);
  //}
  #endif

  for(int ii=0;ii<recv_size;ii++)
    recv_cnt[ii+MAX_RECV_NUM*bw]=0;
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

  volatile uint64_t sum=1;
  uint64_t knt=0;
  for(knt=0;(knt<MAX_WAIT && sum != 0);knt++)
  {
    sum=0;
    for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii] == 0 ? 1 : 0);
    for(int ii=0;ii<recv_size;ii++) sum+= (recv_cnt[ii+MAX_RECV_NUM] == 0 ? 1 : 0);
  }
  assert(knt <  MAX_WAIT); 
  _bgq_msync();


  for(int ii=0;ii<recv_size;ii++)
  {
    printf("%dth packet is at %x\n",ii,recv_cnt[ii]);
    printf("%dth packet is at %x\n",ii,recv_cnt[ii+MAX_RECV_NUM]);
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
  if (spi_hdl->nNodes!=1)
  {
    rc = Kernel_DeallocateGlobalInterruptClassRoute ( spi_hdl->giID );
    assert(rc==0);
  }

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

uint32_t mapping_table_new(spi_hdl_t* spi_hdl)
{
  uint32_t a, b, c, d, e, id;
  uint32_t rc,ii;
  uint32_t node_id = Kernel_GetRank();
  BG_JobCoords_t sblock;
  Kernel_JobCoords(&sblock); 
  assert(sblock.corner.core == 0);
  //Personality_t Mypersonality;
  
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
  spi_hdl->nNodes=totN;
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

  //determine the root of GI
  {
    uint32_t Ac = (int)(sblock.shape.a/2) + sblock.corner.a;
    uint32_t Bc = (int)(sblock.shape.b/2) + sblock.corner.b;
    uint32_t Cc = (int)(sblock.shape.c/2) + sblock.corner.c;
    uint32_t Dc = (int)(sblock.shape.d/2) + sblock.corner.d;
    uint32_t Ec = (int)(sblock.shape.e/2) + sblock.corner.e;

    uint32_t myCoord[5]={bg_map[node_id].a + sblock.corner.a,
                         bg_map[node_id].b + sblock.corner.b,
                         bg_map[node_id].c + sblock.corner.c,
                         bg_map[node_id].d + sblock.corner.d,
                         bg_map[node_id].e + sblock.corner.e};

    spi_hdl->cr = malloc(sizeof(ClassRoute_t));
    uint16_t input,output;
    //rc = Kernel_GetPersonality( &Mypersonality, sizeof( Mypersonality ) );
    //assert(rc == 0);
    //((ClassRoute_t*)(spi_hdl->cr))->input=Mypersonality.Network_Config.PrimordialClassRoute.GlobIntUpPortInputs;
    //((ClassRoute_t*)(spi_hdl->cr))->output=Mypersonality.Network_Config.PrimordialClassRoute.GlobIntUpPortOutputs;

    //A direction
    input=BGQ_CLASS_INPUT_LINK_LOCAL;
    #define nestedif(xxx,yyy,ppp,mmm,aaa) \
    if ( myCoord[xxx] < yyy ) \
    { \
      output=BGQ_CLASS_LINK_##ppp; \
      if (bg_map[node_id].##aaa>0) input=(BGQ_CLASS_LINK_##mmm | input);  \
    } \
    else if ( myCoord[xxx] > yyy) \
    { \
      output=BGQ_CLASS_LINK_##mmm; \
      if ( bg_map[node_id].##aaa<sblock.shape.##aaa-1) input=(BGQ_CLASS_LINK_##ppp | input);  \
    } \
    else if ( myCoord[xxx] == yyy ) \
    { \
      if (bg_map[node_id].##aaa>0) input=(BGQ_CLASS_LINK_##mmm | input);  \
      if (bg_map[node_id].##aaa<sblock.shape.##aaa-1) input=(BGQ_CLASS_LINK_##ppp | input);  

    #define close_brace }

    nestedif(0,Ac,AP,AM,a)
      nestedif(1,Bc,BP,BM,b)
        nestedif(2,Cc,CP,CM,c)
          nestedif(3,Dc,DP,DM,d)
            nestedif(4,Ec,EP,EM,e)

    close_brace
    close_brace
    close_brace
    close_brace
    close_brace

    #undef close_brace
    #undef nestedif
  }
      

//  node.dim[CR_AXIS_A] = node.a_nodes = sblock.shape.a;
//  node.dim[CR_AXIS_B] = node.b_nodes = sblock.shape.b;
//  node.dim[CR_AXIS_C] = node.c_nodes = sblock.shape.c;
//  node.dim[CR_AXIS_D] = node.d_nodes = sblock.shape.d;
//  node.dim[CR_AXIS_E] = node.e_nodes = sblock.shape.e;


  #ifdef spi_debug
  printf("my id is %d\n",node_id);
  printf("geom is %d x %d x %d x %d x %d\n", sblock.shape.a,sblock.shape.b,sblock.shape.c,sblock.shape.d,sblock.shape.e);
  printf("my rel coord is %d x %d x %d x %d x %d\n", bg_map[node_id].a,bg_map[node_id].b,bg_map[node_id].c,bg_map[node_id].d,bg_map[node_id].e);

  printf("completing init_mapping_table\n" );
  #endif

  free(bg_map);
  return node_id;

};
//MUSPI_GIBarrier_t barrier;
  //MUSPI_GIBarrierInit( &barrier, 0 );
  //MUSPI_GIBarrierEnterAndWait( &barrier );
  //printf("GI Barrier done\n");

void printSpiDiagnostics(void)
{
   printf("node:%d startCalls %u waitCalls %u\n",
          staticId, nExecuteCalls, nCompleteCalls);
   exit(-1);
}


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

//      if ( myCoord[1] < Bc )
//      {
//        output=BGQ_CLASS_LINK_BP;
//        //am I in the edge?
//        if (bg_map[node_id].b>0) input=(BGQ_CLASS_LINK_BM | input); 
//      }
//      else if ( myCoord[1] > Bc)
//      {
//        output=BGQ_CLASS_LINK_BM;
//        if ( bg_map[node_id].b<sblock.shape.b-1) input=(BGQ_CLASS_LINK_BP | input); 
//      }
//      else if ( myCoord[1] == Bc )
//      {
//        if (bg_map[node_id].b>0) input=(BGQ_CLASS_LINK_BM | input); 
//        if (bg_map[node_id].b<sblock.shape.b-1) input=(BGQ_CLASS_LINK_BP | input); 
//     
//        if ( myCoord[2] < Cc )
//        {
//          output=BGQ_CLASS_LINK_CP;
//          //am I in the edge?
//          if (bg_map[node_id].c>0)
//            ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_CM | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//        }
//        else if ( myCoord[2] > Cc)
//        {
//          ((ClassRoute_t*)(spi_hdl->cr))->output=BGQ_CLASS_LINK_CM;
//          if ( bg_map[node_id].c<sblock.shape.c-1)
//            ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_CP | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//        }
//        else if ( myCoord[2] == Cc )
//        {
//          if (bg_map[node_id].c>0)
//            ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_CM | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//          if (bg_map[node_id].c<sblock.shape.c-1)
//            ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_CP | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//      
//          if ( myCoord[3] < Dc )
//          {
//            ((ClassRoute_t*)(spi_hdl->cr))->output=BGQ_CLASS_LINK_DP;
//            //am I in the edge?
//            if (bg_map[node_id].d>0)
//              ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_DM | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//          }
//          else if ( myCoord[3] > Dc)
//          {
//            ((ClassRoute_t*)(spi_hdl->cr))->output=BGQ_CLASS_LINK_DM;
//            if ( bg_map[node_id].d<sblock.shape.d-1)
//              ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_DP | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//          }
//          else if ( myCoord[3] == Dc )
//          {
//            if (bg_map[node_id].d>0)
//              ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_DM | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//            if (bg_map[node_id].d<sblock.shape.d-1)
//              ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_DP | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//       
//            if ( myCoord[4] < Ec )
//            {
//              ((ClassRoute_t*)(spi_hdl->cr))->output=BGQ_CLASS_LINK_EP;
//              //am I in the edge?
//              if (bg_map[node_id].e>0)
//                ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_EM | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//            }
//            else if ( myCoord[4] > Ec)
//            {
//              ((ClassRoute_t*)(spi_hdl->cr))->output=BGQ_CLASS_LINK_EM;
//              if ( bg_map[node_id].e<sblock.shape.e-1)
//                ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_EP | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//            }
//            else if ( myCoord[4] == Ec )
//            {
//              if (bg_map[node_id].e>0)
//                ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_EM | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//              if (bg_map[node_id].e<sblock.shape.e-1)
//                ((ClassRoute_t*)(spi_hdl->cr))->input=(BGQ_CLASS_LINK_EP | ((ClassRoute_t*)(spi_hdl->cr))->input); 
//            }
//          }
//        }
//      }
//    }

#endif
