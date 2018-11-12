#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cs_gettime.h"

static double gettimecost(void) {
  double t0,t1,twait = 0.1;
  long long int n = 0,i;

  if(0) {
    t0 = cs_gettime();
    t1 = t0;
    while(t1 - t0 < twait) {
      n++;
      t1 = cs_gettime();
    }
  } else {
    n = 0;
    t0 = cs_gettime();
    do {
      t1 = cs_gettime();
    } while(t1 == t0);
    while(cs_gettime() == t1) n++;
    
    if(n < 1000) n = 1000;
  }
  
  if(0) printf("Wait is %.3f s, n is %.3f M.\n",twait,1e-6*n);

  i = n;
  t0 = cs_gettime();
  while(i > 0) {
    (void) cs_gettime();
    i--;
  }
  t0 = cs_gettime() - t0;

  i = n;
  t1 = cs_gettime();
  while(i > 0) {
    (void) cs_gettime();
    (void) cs_gettime();
    i--;
  }
  t1 = cs_gettime() - t1;

  if(0) printf("The two times are t0=%.3e s, t1 = %.3e s\n",t0,t1);

  return (t1-t0)/n;
}

static int cmp_double(const void *ap,const void *bp) {
  double a = ((double *) ap)[0], b = ((double *) bp)[0];
  if(a < b) return -1;
  else if(a > b) return 1;
  else return 0;
}

double cs_clocksync(int np,int pid,int nmsg,
		 double gettime_cost[1],
		 double corr_std[1]) {
  const int nsamples = nmsg;
  double *tvec,*rvec,*corrvec;
  MPI_Status status;

  double tcost,tcorr,latavg,latstd,cavg,cstd;
  int shift;

  tvec = (double *) malloc(sizeof(double) * nsamples);
  rvec = (double *) malloc(sizeof(double) * nsamples);
  corrvec = (double *) malloc(sizeof(double) * nsamples);

  MPI_Barrier(MPI_COMM_WORLD);

  if(gettime_cost == NULL || gettime_cost[0] < 0.0)
    tcost = gettimecost();
  else
    tcost = gettime_cost[0];

  shift = 1;
  while((1<<shift) < np) shift = shift + 1;

  tcorr = 0.0;
  latavg = 0.0;
  latstd = -99.0;
  cavg = 0.0;
  cstd = -99.0;

  while(shift > 0) {
    if( ( pid & ((1<<(shift-1))-1) ) == 0) {
      int peer = pid ^ (1<<(shift-1));
      int i;

      if(0) printf("pid=%d  shift=%d  peer=%d\n",pid,peer,shift);

      tvec[0] = tcost;
      if(pid < peer) /* I am source, ready to initiate loop */
	MPI_Send(tvec,1,MPI_DOUBLE,peer,999,MPI_COMM_WORLD);
      for(i = 0; i<nsamples; i++) {
	MPI_Recv(rvec+i,1,MPI_DOUBLE,peer,999,MPI_COMM_WORLD,&status);
	tvec[i+1] = cs_gettime() + tcorr;
	MPI_Send(tvec+i+1,1,MPI_DOUBLE,peer,999,MPI_COMM_WORLD);
      }
      if(pid > peer) /* Receive last message from source */
	MPI_Recv(rvec+i,2,MPI_DOUBLE,peer,999,MPI_COMM_WORLD,&status);

      /* Now, compute time correction */
      if(pid > peer) {
	int n = nsamples - 2,i;
	double latsum = 0.0,latsum2 = 0.0;
	double csum,csum2;

	for(i = 0; i<n; i++) {
	  double lat = 0.5*(tvec[i+2] - tvec[i+1] - rvec[0] - tcost);
	  corrvec[i] = tvec[i+2] - rvec[i+1] - lat - 0.5*(tcost+rvec[0]);
	  latsum = latsum + lat;
	  latsum2 = latsum2 + lat*lat;
	}
	latavg = latsum/n;
	latstd = sqrt( (latsum2 - n*latavg*latavg)/(n-1) );

	qsort(corrvec,n,sizeof(double),cmp_double);
	csum = 0.0;
	csum2 = 0.0;
	for(i = n/4; i<n-n/4; i++) {
	  csum = csum + corrvec[i];
	  csum2 = csum2 + corrvec[i]*corrvec[i];
	}
	n = i-n/4;
	cavg = csum/n;	
	cstd = sqrt( (csum2 - n*cavg*cavg)/(n-1) );
	tcorr = -cavg;
      }

    }
    shift = shift - 1;
  }

  free(corrvec);
  free(rvec);
  free(tvec);

  MPI_Barrier(MPI_COMM_WORLD);

  if(0) { /* Print out a table of suynchronization data */
    int i,j;
    double data[5];
    double (*table)[5] = NULL;

    if(pid == 0)
      table = (double (*)[5]) malloc(sizeof(double[5]) * np);

    data[0] = tcorr;
    data[1] = tcost;
    data[2] = cstd;
    data[3] = latavg;
    data[4] = latstd;

    MPI_Gather(data,sizeof(data)/sizeof(double),MPI_DOUBLE,
	       table,sizeof(data)/sizeof(double),MPI_DOUBLE,
	       0,MPI_COMM_WORLD);
    if(pid == 0) {
      printf("%%%% TABLE OF LATACIES AND CORRECTION TERMS\n");
      printf("%%%% COLUMNS: tcorr tcost cstd latavg latstd\n");
      for(i = 0; i<np; i++) {
	printf("%6d",i);
	for(j = 0; j<sizeof(data)/sizeof(double); j++)
	  printf("  %20.10e",table[i][j]);
	printf("\n");
      }

      free(table);
    }
  }
  
  if(gettime_cost != NULL && gettime_cost[1] < 0.0)
    gettime_cost[0] = tcost;
  if(corr_std != NULL)
    corr_std[0] = cstd;
  return tcorr;
}

#if 0
double cs_pingpong(int peer,int nmesg,int sourceflag,
		   double tcorr,double tcost,double cstd_p[1]) {
  const int nsamples = nmesg;

  double *tvec,*rvec,*corrvec;
  MPI_Status status;

  double tcost,tcorr,latavg,latstd,cavg,cstd;
  int shift;

  tvec = (double *) malloc(sizeof(double) * nsamples);
  rvec = (double *) malloc(sizeof(double) * nsamples);
  corrvec = (double *) malloc(sizeof(double) * nsamples);

  latavg = 0.0;
  latstd = -99.0;
  cavg = 0.0;
  cstd = 0.0;

  tvec[0] = tcost;
  if(sourceflag != 0) /* I am source, ready to initiate loop */
    MPI_Send(tvec,1,MPI_DOUBLE,peer,999,MPI_COMM_WORLD);
  for(i = 0; i<nsamples; i++) {
    MPI_Recv(rvec+i,1,MPI_DOUBLE,peer,999,MPI_COMM_WORLD,&status);
    tvec[i+1] = cs_gettime() + tcorr;
    MPI_Send(tvec+i+1,1,MPI_DOUBLE,peer,999,MPI_COMM_WORLD);
  }
  if(pid > peer) /* Receive last message from source */
    MPI_Recv(rvec+i,2,MPI_DOUBLE,peer,999,MPI_COMM_WORLD,&status);
  
  /* Now, compute time correction */
  if(sourceflag == 0) {
    int n = nsamples - 2,i;
    double latsum = 0.0,latsum2 = 0.0;
    double csum,csum2;
    
    for(i = 0; i<n; i++) {
      double lat = 0.5*(tvec[i+2] - tvec[i+1] - rvec[0] - tcost);
      corrvec[i] = tvec[i+2] - rvec[i+1] - lat - 0.5*(tcost+rvec[0]);
      latsum = latsum + lat;
      latsum2 = latsum2 + lat*lat;
    }
    latavg = latsum/n;
    latstd = sqrt( (latsum2 - n*latavg*latavg)/(n-1) );
    
    qsort(corrvec,n,sizeof(double),cmp_double);
    csum = 0.0;
    csum2 = 0.0;
    for(i = n/4; i<n-n/4; i++) {
      csum = csum + corrvec[i];
      csum2 = csum2 + corrvec[i]*corrvec[i];
    }
    n = i-n/4;
    cavg = csum/n;	
    cstd = sqrt( (csum2 - n*cavg*cavg)/(n-1) );
  }

  free(corrvec);
  free(rvec);
  free(tvec);

  if(cstd_p != NULL)
    cstd_p[0] = cstd;
  return -cavg;
}

void cs_synctest(int np,int pid,int ntrial,int nmsg,
		 double tcorr,double tcost,
		 double corravg[1],double corrstd[1]) {

  int *list = (int *) malloc(sizeof(int) * np);

  srand(2136743); // So all processors generate the same sequence.
    
  for(itrial = 0; itrial<ntrial; itrial++) {
    for(i = 0; i<np; i++)
      list[i] = 0;

    for(i = 0; i<np-1; i++) { // Create random permutation.
      int t,j = rand()%(np-i);
      if(j < 0) j = -j;
      t = list[i];
      list[i] = list[i+j];
      list[i+j] = t;
    }

    dir = 1 - 2*(pid%2);
    peer = list[pid+dir];

  }

}
#endif
