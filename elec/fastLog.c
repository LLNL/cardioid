#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "fastLog.h" 

#define uint64_t long long unsigned
uint64_t getTime()
{
   struct timeval ptime;
   uint64_t  t = 0;
   gettimeofday(&ptime, (struct timezone *)NULL);
   t = ((uint64_t)1000000)*(uint64_t)ptime.tv_sec + (uint64_t)ptime.tv_usec;
   return t; 
}

int main()
{
    fastLogInit(); 
    double x ,y; 
    printf("byte_order=%d %d %d\n",BYTE_ORDER,LITTLE_ENDIAN,BIG_ENDIAN);
    for ( x=1.000;x<2.000001;x+=0.2) {y = fastLogNew(&x); printf("%e %e %e %e\n",x,y,log(x),y-log(x)); }
    for ( x=1.000;x<2.000001;x+=0.2) {y = fastLog(x); printf("%e %e %e %e\n",x,y,log(x),y-log(x)); }

    uint64_t t0,t1; 

    srand48(283458141); 
    double sum0=0.0; 
    t0 = getTime(); 
    int maxLoop = 100000000;
    for(int i=0;i<maxLoop;i++) 
    {
          double x = 100*(drand48()+0.00002); 
          sum0 += log(x) ;
    }
    t1 = getTime(); 
    printf("log       : sum=%e time=%e\n",sum0,(t1-t0)*1e-6); 

    srand48(283458141); 
    double error,sum=0.0; 
    t0 = getTime(); 
    for(int i=0;i<maxLoop;i++) 
    {
          double x = 100*(drand48()+0.00002); 
          sum += fastLog(x) ;
    }
    t1 = getTime(); 
    error=(sum-sum0)/sum0; 
    printf("fastLog   : sum=%e time=%e error=%e\n",sum,(t1-t0)*1e-6,error); 

    srand48(283458141); 
    sum =0.0; 
    t0 = getTime(); 
    for(int i=0;i<maxLoop;i++) 
    {
          double x = 100*(drand48()+0.00002); 
          double v = fastLogNew(&x) ;
          sum += v ;
    }
     t1 = getTime(); 
    error=(sum-sum0)/sum0; 
    printf("fastLogNew: sum=%e time=%e error=%e\n",sum,(t1-t0)*1e-6,error); 
}
