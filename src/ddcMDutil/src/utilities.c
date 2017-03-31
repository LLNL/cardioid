/* $Id$ */
#include "utilities.h"
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/resource.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/times.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <wordexp.h>
#ifndef __APPLE__
#include <malloc.h>
#include "mpiUtils.h"

char* strptime(const char * restrict buf, const char * restrict format,
	       struct tm * restrict timeptr);

#endif
#include <assert.h>
#include "ddcMalloc.h"
#include "error.h"
#include "mpiUtils.h"

int  _wordexp(const char *restrict words , wordexp_t *restrict pwordexp, int flags )
{
	int rc = wordexp(words, pwordexp, flags);
#ifdef __APPLE__
	pwordexp->we_wordc=1; 
	pwordexp->we_wordv[0] = strdup(words); 
	rc =0; 
#endif
	return rc; 
}
void _wordfree(wordexp_t *pwordexp)
{
#ifndef __APPLE__
 wordfree(pwordexp);
#endif
}

/* Sorts integer data */
int intSortFunction(const void* av, const void* bv)
{
	const int* a = (const int*)av;
	const int* b = (const int*)bv;
	if (*a > *b) return 1;
	if (*a < *b) return -1;
	return 0;
}


char *location(char *string, char *file, int linenumber)
{
	static char buffer[256];
	sprintf(buffer, "%s in %s at linenumber %d", string, file, linenumber);
	return buffer;
}

/* MEMORY DEBUG --------------------------------------------------------------*/
/* PRINT OUT CHANGES IN MEMORY USAGE                                          */
/* usage: mem_debug("string indicating location in code");                    */
/* ---------------------------------------------------------------------------*/
/*   */
void mem_debug(char *mem_str)
{
  if (getRank(0) != 0) return; 
#ifdef __APPLE__
   static int flag =1; 
   if (flag) printf("No mallinfo on OS X.  Sorry no minfo for <%s>.\n",mem_str);
   flag =0; 
#else
   static struct mallinfo old_minfo, minfo;
   static int gstep = 0;
   {
      minfo = mallinfo();
      if (minfo.arena > old_minfo.arena)
      {
	 printf("MINFO step %d : %d %d kB %s\n", gstep, (minfo.arena - old_minfo.arena)/1024, minfo.arena/1024, mem_str);
	 fflush(stdout);
	 old_minfo = minfo;
      }
   }
#endif
}

/*
*/
/*
# arena is the total amount of space in the heap.
# ordblks is the number of chunks which are not in use.
# uordblks is the total amount of space allocated by malloc .
# fordblks is the total amount of space not in use.
# keepcost is the size of the top most memory block.
*/
void CHKMEM(char *string)
{
#ifndef __APPLE__
#include <sys/resource.h>
	static struct mallinfo minfo, minfosave;
	minfosave = minfo;
	minfo = mallinfo();
	if (minfo.arena > minfosave.arena)
	{
		printf("***********\n");
		printf("%s:%d %d %d %d\n", string, minfosave.arena, minfosave.ordblks, minfosave.uordblks, minfosave.fordblks);
		printf("%s:%d %d %d %d\n", string, minfo.arena, minfo.ordblks, minfo.uordblks, minfo.fordblks);
		printf("%s:%d %d %d %d\n", string, minfo.arena - minfosave.arena, minfo.ordblks - minfosave.ordblks, minfo.uordblks - minfosave.uordblks, 		       minfo.fordblks - minfosave.fordblks);
		printf("***********\n");
	}
#endif
}
clock_t clk()
{
	struct tms b;
	times(&b);
	return (b.tms_utime + b.tms_stime);
}

unsigned  long long getTick(void)
{
	struct timeval ptime;
	unsigned long long tick; 
//#ifdef BGL
//	unsigned long long rts_get_timebase(void);
//	tick = rts_get_timebase();
//#else
	gettimeofday(&ptime, (struct timezone *)NULL);
	unsigned long long sec = ptime.tv_sec  ;
	unsigned long long usec = ptime.tv_usec;
	tick = 1000000*sec + usec; 
//#endif
	return tick ;
}

double pelapsetime(void)
{
	struct timeval ptime;
	double t;
	unsigned long long rts_get_timebase(void);
#ifdef BGL
	double seconds_per_cycle = 1.4285714285714285714e-9;
	t = ((double)rts_get_timebase())*seconds_per_cycle;
#else
	gettimeofday(&ptime, (struct timezone *)NULL);
	t = ptime.tv_sec + 1e-6*ptime.tv_usec;
/*	t= MPI_Wtime();  */
/*	t = time(NULL); */
#endif
	return t;
}

double pcputime(void)
{
	double tusec, tssec, tumsec, tsmsec;
	double tsec, t0, t1;
	struct rusage ruval;

	
	getrusage(RUSAGE_SELF, &ruval);

	/*  Get user and system time in integer sec and microsec */
	/*  The user and system time field of the rusage structure */
	/*  are themselves structures. */
	tusec = (ruval.ru_utime).tv_sec;
	tumsec = (ruval.ru_utime).tv_usec;

	tssec = (ruval.ru_stime).tv_sec;
	tsmsec = (ruval.ru_stime).tv_usec;

	t0 = tusec + 1.e-6*tumsec;
	t1 = tssec + 1.e-6*tsmsec;

	/*  Return total time as the function value */
	tsec = t0 + t1;
	return (tsec);
}

/** Note that strptime only populates members of tm for which data are
 *  explicity provided in the string.  According to the strptime man
 *  page on OS X,
 *
 *    If the format string does not contain enough conversion
 *    specifications to completely specify the resulting struct tm, the
 *    unspecified members of timeptr are left untouched.  For example,
 *    if format is ``%H:%M:%S'', only tm_hour, tm_sec and tm_min will be
 *    modified.  If time relative to today is desired, initialize the
 *    timeptr structure with today's date before passing it to
 *    strptime().
 *
 *  Since we typically don't provide time zone or daylight saving
 *  information, we don't get the right value for tm.tm_isdst (and hence
 *  we get the wrong time) unless we properly initialize tm with today's
 *  date before calling strptime().
 */
time_t convert_timestring(char *time_string)
{
   time_t t0 = time(NULL);
   struct tm* local_tm = localtime(&t0);
   
   struct tm tm = *local_tm;
   char *end = strptime(time_string,"%m/%d/%Y-%X",&tm);
   time_t  t1; 
   if (end != NULL)
   {
      t1 = mktime(&tm);
   }
   else
   {
      if (getRank(0)==0) 
      {
         printf("time_string %s is the wrong format\nCorrect format is mm/dd/yyyy-hh:mm:ss\n",time_string); 
      }
      t1 = t0; 
   }
   return t1;
}


char *timestamp_string(void)
{
   char *c, *date;
   struct tm *local_tm ;
   time_t t;
   t = time(NULL);
   local_tm=localtime(&t);
   date = asctime(local_tm);
   c = strchr(date, '\n');
   if (c != NULL) *c = (char)'\0';
   return date; 
}

void timestamp_anyTask(const char* string)
{
   char *c, *date;
   struct tm *local_tm ;
   time_t t;
   t = time(NULL);
   local_tm=localtime(&t);
   char *dateString = asctime(local_tm); 
   date = strdup(dateString);

   c = index(date, '\n');
   if (c != NULL) *c = (char)'\0';
   printf("%s: %s\n", date, string);
   /*		mem_debug(string); */
   fflush(stdout);
   ddcFree(date);
}

void timestamp(const char *string)
{
   if (getRank(0) == 0) timestamp_anyTask(string);
}



/*
   __ctype_b()
   {
   }

*/

/** Make sure the builtin types are the sizes we think they are.
 * Inconsistencies will likely break something such as binary I/O. */
void checkLimits(void)
{
   assert(sizeof(char) == 1);
   assert(sizeof(int) == 4);
   assert(sizeof(unsigned) == 4);
   assert(sizeof(float) == 4);
   assert(sizeof(double) == 8);
   assert(sizeof(long long unsigned) == 8);
}

int removeDuplicates(unsigned* a, int na)
{
   if (na == 0) return 0;
   int newEnd = 0;
   int scan = 1;

   while (scan < na)
   {
      if (a[newEnd] != a[scan])
      {
         ++newEnd;
         if (newEnd != scan)
            a[newEnd] = a[scan];
      }
      ++scan;
   }
   return newEnd+1;
}
char *astrcat(char *dest, char *str)
{
   dest = realloc(dest,strlen(dest)+strlen(str)+1); 
   strcat(dest,str);
   return dest; 
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
