#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <papi.h>

// -------------------------------------------------
// C or  C++ routines 
// -------------------------------------------------
// HPM_Init();    for serial jobs
// HPM_Start("label");
// HPM_Stop("label");
// HPM_Reset("label");
// HPM_Enable();
// HPM_Disable();
// HPM_Print();   for serial jobs
// -------------------------------------------------
// Fortran interface to control the sampled region
// -------------------------------------------------
// call hpm_init()   for serial jobs
// call hpm_start('label')
// call hpm_stop('label')
// call hpm_reset('label')
// call hpm_enable();
// call hpm_disable();
// call hpm_print()  for serial jobs
// -------------------------------------------------
#define NUM_COUNTERS 4
#define MAX_CODE_BLOCKS 100
#define LABEL_LEN 80
#define MAX_GROUPS 22

int index_from_label(char *);
#ifdef USE_MPI
int get_ranks_per_node(void);
#endif

#pragma weak hpm_init=HPM_Init
#pragma weak hpm_print=HPM_Print

static long long counter_in[MAX_CODE_BLOCKS][NUM_COUNTERS];
static long long counter_sum[MAX_CODE_BLOCKS][NUM_COUNTERS];
static int block_starts[MAX_CODE_BLOCKS];
static int block_stops[MAX_CODE_BLOCKS];
static double elapsed_in[MAX_CODE_BLOCKS];
static double elapsed_sum[MAX_CODE_BLOCKS];
static char code_block_label[MAX_CODE_BLOCKS][LABEL_LEN];

static int disabled = 0;
static int initialized = 0;
static int code_block = 0;
static int eventHandle = PAPI_NULL;

static int numcounters = NUM_COUNTERS;
static int group = 0;
static int group_list[MAX_GROUPS];
static int use_counter_list = 0;
static int num_groups = 0;
static int collect_all_counts = 1;
static int use_event_list = 0;
static char envname[NUM_COUNTERS][128];

static char event_name[][4][80] = { 
// group 0:
  {
    "perf::cycles",
    "perf::instructions",
    "perf::L1-dcache-loads",
    "perf::L1-dcache-load-misses"
  },  
// group 1:
  {
    "perf::ref-cycles",
    "perf::instructions",
    "perf::llc-prefetch-misses",
    "perf::llc-load-misses"
  },  
// group 2:
  {
    "perf::ref-cycles",
    "perf::instructions",
    "perf::L1-dcache-loads",
    "perf::L1-dcache-stores"
  },  
// group 3:
  {
    "PAPI_TOT_CYC",
    "PAPI_TOT_INS",
    "PAPI_REF_CYC",
    "PAPI_DP_OPS"
  },  
// group 4:
  {
    "PAPI_TOT_CYC",
    "PAPI_TOT_INS",
    "PAPI_REF_CYC",
    "PAPI_VEC_DP"
  },  
// group 5:
  {
    "PAPI_TOT_INS",
    "PAPI_LD_INS",
    "PAPI_SR_INS",
    "SIMD_FP_256:PACKED_DOUBLE"
  },
// group 6:  see likwid wiki for haswell avx instructions
  {
    "perf_raw::r01C6",
    "perf_raw::r02C6",
    "perf_raw::r04C6",
    "perf_raw::r003c"
  }
};


void HPM_Init(void)
{
   int i, j, myrank, nranks, rc, ic;
   char * ptr = NULL;
   char * list_ptr;
   char delimiters[] = {" ,"};
   char counter_name[128];

#ifdef USE_MPI
   PMPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   PMPI_Comm_size(MPI_COMM_WORLD, &nranks);
#else
   myrank = 0;
   nranks = 1;
#endif

   for (j=0; j<MAX_CODE_BLOCKS; j++) elapsed_sum[j] = 0.0;

   for (j=0; j<MAX_CODE_BLOCKS; j++) {
      block_starts[j] = 0;
      block_stops[j] = 0;
   }

   for (j=0; j<MAX_CODE_BLOCKS; j++) {
      for (ic=0; ic<numcounters; ic++) counter_sum[j][ic] = 0LL;
   }

   // set HPM_GROUP to make all ranks count using the same group
   ptr = getenv("HPM_GROUP");
   if (ptr == NULL) group = 0;
   else             group = atoi(ptr);

   use_event_list = 0;
   list_ptr = getenv("HPM_EVENT_LIST");
   if (list_ptr != NULL) {
      use_event_list = 1;
      ptr = strtok(list_ptr, delimiters);
      j = 0;
      while(ptr != NULL) {
         strcpy(envname[j], ptr);
         ptr = strtok(NULL, delimiters);
         j++;
         if (j == NUM_COUNTERS) break;
      }
      numcounters = j;
      group = 1000;
   }

   // set HPM_GROUP_LIST to make ranks round-robin over a list of groups
   // HPM_GROUP_LIST takes precedence if HPM_GROUP is also set
//   list_ptr = getenv("HPM_GROUP_LIST");
//   if (list_ptr != NULL) {
//      if (strncasecmp(list_ptr,"all",3) == 0) {
//         num_groups = MAX_GROUPS;
//         for (i=0; i<num_groups; i++) group_list[i] = i;
//      }
//      else {
//         num_groups = 0;
//         ptr = strtok(list_ptr, delimiters);
//         while(ptr != NULL) {
//             group_list[num_groups] = atoi(ptr);
//             ptr = strtok(NULL, delimiters);
//             num_groups++;
//         }
//
//         if ( (num_groups > nranks) &&  (myrank == 0) )
//             fprintf(stderr, "hpm warning: number of groups > nranks ... using the first nranks groups\n");
//      }
//
//      group = group_list[ myrank % num_groups ];
//
//      use_counter_list = 1;
//   }

   // optionally disable collection of all counts for the job summary
//   ptr = getenv("COLLECT_ALL_COUNTS");
//   if (ptr != NULL) {
//      if (strncasecmp(ptr,"no",2) == 0) collect_all_counts = 0;
//      else                              collect_all_counts = 1;
//   }

   // for now counting is done by just the master thread
   rc = PAPI_library_init(PAPI_VER_CURRENT);
   if (rc < 0) {
      fprintf(stderr, "PAPI initialization failed for group %d on rank %d ... counters not enabled\n", group, myrank);
      return;
   }

   rc = PAPI_create_eventset(&eventHandle);
   if (rc != PAPI_OK) {
      fprintf(stderr, "PAPI failed to create an event set for group %d on rank %d\n", group, myrank);
      return;
   }

   for (ic = 0; ic<numcounters; ic++) {
      memset(counter_name, '\0', sizeof(counter_name));
      if (use_event_list) sprintf(counter_name, "%s", envname[ic]);
      else                sprintf(counter_name, "%s", event_name[group][ic]);
      rc = PAPI_add_named_event(eventHandle, counter_name);
      if (rc != PAPI_OK) {
         fprintf(stderr, "PAPI failed to add event %s for group %d on rank %d ... counters not enabled\n",  counter_name, group, myrank);
         return;
      }
   }

   rc = PAPI_start(eventHandle);
   if (rc != PAPI_OK) {
      fprintf(stderr, "PAPI failed to start for group %d on rank %d ... counters not enabled\n", group, myrank);
      return;
   }

   if (myrank == 0) {
      if (use_counter_list) printf("HPM counting with a list of groups.\n");
      else                  printf("HPM counting with group = %d.\n", group);
   }

   initialized = 1;

}

void HPM_Disable(void)
{
   disabled = 1;
   return;
}

void HPM_Enable(void)
{
   disabled = 0;
   return;
}

#pragma weak hpm_disable_=hpm_disable
void hpm_disable(void)
{
   disabled = 1;
   return;
}

#pragma weak hpm_enable_=hpm_enable
void hpm_enable(void)
{
   disabled = 0;
   return;
}

void HPM_Reset(char * this_label)
{
   int ic, rc, j;

   j = index_from_label(this_label);

   block_starts[j] = 0;
   block_stops[j] = 0;
   elapsed_sum[j] = 0.0;

   for (ic=0; ic<NUM_COUNTERS; ic++) counter_sum[j][ic] = 0LL;
   
}

#pragma weak hpm_reset_=hpm_reset
void hpm_reset(char * f_label, int length)
{
   int ic, rc, j;
   char this_label[LABEL_LEN];
   strncpy(this_label, f_label, length);
   this_label[length] = '\0';
   j = index_from_label(this_label);

   block_starts[j] = 0;
   block_stops[j] = 0;
   elapsed_sum[j] = 0.0;

   for (ic=0; ic<NUM_COUNTERS; ic++) counter_sum[j][ic] = 0LL;
   
}


void HPM_Start(char * this_label)
{
   int ic, rc, j;
   struct timeval tv;

   if (disabled || !initialized) return;

   gettimeofday(&tv, NULL);

   j = index_from_label(this_label);

   rc = PAPI_read(eventHandle, counter_in[j]);
   if (rc != PAPI_OK) return;

   block_starts[j] += 1;

   elapsed_in[j] = ((double) tv.tv_sec) + 1.0e-6*((double) tv.tv_usec);

   return;

}

#pragma weak hpm_start_=hpm_start
void hpm_start(char * f_label, int length)
{
   int ic, rc, j;
   char this_label[LABEL_LEN];
   struct timeval tv;

   if (disabled || !initialized) return;

   strncpy(this_label, f_label, length);
   this_label[length] = '\0';

   gettimeofday(&tv, NULL);

   j = index_from_label(this_label);

   rc = PAPI_read(eventHandle, counter_in[j]);
   if (rc != PAPI_OK) return;

   block_starts[j] += 1;

   elapsed_in[j] = ((double) tv.tv_sec) + 1.0e-6*((double) tv.tv_usec);

   return;

}



void HPM_Stop(char * this_label)
{
   int ic, rc, j;
   struct timeval tv;
   long long values[NUM_COUNTERS];

   if (disabled || !initialized) return;

   if (code_block >= MAX_CODE_BLOCKS) return;

   gettimeofday(&tv, NULL);

   j = index_from_label(this_label);

   rc = PAPI_read(eventHandle, values);
   if (rc != PAPI_OK) return;

   for (ic=0; ic<numcounters; ic++) counter_sum[j][ic] += (values[ic] - counter_in[j][ic]);

   block_stops[j] += 1;

   elapsed_sum[j] += (((double) tv.tv_sec) + 1.0e-6*((double) tv.tv_usec)) - elapsed_in[j];

   return;
}


#pragma weak hpm_stop_=hpm_stop
void hpm_stop(char * f_label, int length)
{
   int ic, rc, j;
   struct timeval tv;
   long long values[NUM_COUNTERS];
   char this_label[LABEL_LEN];

   if (disabled || !initialized) return;

   if (code_block >= MAX_CODE_BLOCKS) return;

   strncpy(this_label, f_label, length);
   this_label[length] = '\0';

   gettimeofday(&tv, NULL);

   j = index_from_label(this_label);

   rc = PAPI_read(eventHandle, values);
   if (rc != PAPI_OK) return;

   for (ic=0; ic<numcounters; ic++) counter_sum[j][ic] += (values[ic] - counter_in[j][ic]);

   block_stops[j] += 1;

   elapsed_sum[j] += (((double) tv.tv_sec) + 1.0e-6*((double) tv.tv_usec)) - elapsed_in[j];

   return;
}


#ifdef USE_MPI
void HPM_Print(int jobid, int flag)
#else
void HPM_Print(void)
#endif
{
    int j, ic, myrank, nranks, rank, indx;
    int nblocks, max_nblocks, ranks_per_node, max_ranks_per_node;
    int * all_ranks;
    long long * all_counts;
    double * all_elapsed;
    FILE * fp;
    char filename[256], counter_label[256];
    struct countStruct {
                         double value;
                         int rank;
                       };
    struct countStruct my_val[NUM_COUNTERS], min_val[NUM_COUNTERS], max_val[NUM_COUNTERS];
    double dble_counts[NUM_COUNTERS], total_counts[NUM_COUNTERS];
    double avg_time, max_time, avg_value;
    double ld_bandwidth, st_bandwidth, mem_bandwidth;
    double max_ld_bw, max_st_bw, max_total_bw;
    int tag = 99;
    int report_rank, report_size, root_rank;
    
    if (!initialized) return;
    myrank = 0;
    nranks = 1;
    report_rank = 0;
    report_size = 1;

    if (code_block >= MAX_CODE_BLOCKS) nblocks = MAX_CODE_BLOCKS;
    else                               nblocks = code_block;

    ranks_per_node = 1;
    max_nblocks = nblocks;

    sprintf(filename, "hpm_group%d.txt", group);
    printf("filename= %s\n",filename);

    // skip the per-process files when using collect_all_counts
    //if (collect_all_counts) return;  

    //fp = fopen(filename, "w");
    //if (fp == NULL) fp = stderr;
    fp=stdout;

    fprintf(fp, "======================================================================\n");
    fprintf(fp, "Hardware counter data, group = %d.\n", group);
    fprintf(fp, "======================================================================\n");
    printf("nblocks= %d\n",nblocks);

    for (j=0; j<nblocks; j++) {
       if (block_starts[j] == block_stops[j]) {
          fprintf(fp, "----------------------------------------------------------------\n");
          fprintf(fp, "%s, call count = %d, elapsed time = %.3lf\n", code_block_label[j], block_starts[j], elapsed_sum[j]);
          for (ic=0; ic<numcounters; ic++) {
             memset(counter_label, '\0', sizeof(counter_label));
             if (use_event_list) sprintf(counter_label, "%s", envname[ic]);
             else                sprintf(counter_label, "%s", event_name[group][ic]);
             fprintf(fp, "%15lld : %s\n", counter_sum[j][ic], counter_label);
//           fprintf(fp, "%15lld : %s\n", counter_sum[j][ic], event_name[group][ic]);
          }
       }
       else {
          fprintf(fp, "mismatch in starts/stops for code block '%s'\n", code_block_label[j]);
          fprintf(fp, "  starts = %d\n", block_starts[j]);
          fprintf(fp, "  stops  = %d\n", block_stops[j]);

       }
       fprintf(fp, "\n");
    }

    if (fp != stderr) fclose(fp);
}

int index_from_label(char * this_label)
{
   int i, match;
   char * ptr;

   if (code_block < MAX_CODE_BLOCKS)
   {
       match = 0;
       for (i=code_block-1; i>=0; i--)
       {
           if (0 == strcmp(code_block_label[i], this_label))
           {
               match = 1;
               break;
           }
       }

       if (match == 0)
       {
           i = code_block;
           ptr = strcpy(code_block_label[i], this_label);
           if (ptr == NULL) code_block_label[i][0] = '\0';
           code_block ++;
       }
   }

   return i;

}


