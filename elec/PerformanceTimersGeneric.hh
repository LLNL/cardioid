#ifndef MACHINE 
void genericInit() { }
void genericReadCounter(unsigned int counterIndex, int i, uint64_t *counter) { *counter=0; }
#define MAX_THREADS 1
unsigned int counterHandle[MAX_THREADS]; // store the event set handler for each thread
void (*machineSpecficInit)() = genericInit; 
void  (*readCounter)(unsigned int, int,  uint64_t*)=genericReadCounter;
int nCounters_=2; 
const char *counterNames_[] = {"#Calls" , "Time" };
uint64_t getTime()
{
   struct timeval ptime;
   uint64_t  t = 0;
   gettimeofday(&ptime, (struct timezone *)NULL);
   t = ((uint64_t)1000000)*(uint64_t)ptime.tv_sec + (uint64_t)ptime.tv_usec;
   return t; 
}
double getTick()
{
   double seconds_per_cycle = 1.0e-6;
   return seconds_per_cycle; 
}
int getNCores() { return 4;} 
#define MACHINE GENERIC
#endif
