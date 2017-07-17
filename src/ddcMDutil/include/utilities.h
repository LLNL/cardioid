/* $Id$ */
#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <time.h>
#include <wordexp.h>

#ifdef __cplusplus
extern "C" {
#endif

void timestamp_anyTask(const char *string);
void timestamp(const char *string);
char *timestamp_string(void);
void mem_debug(char *mem_str);
void memoryinfo(FILE*file);
double pelapsetime(void);
double pcputime(void);
time_t convert_timestring(char *time_string);
void checkLimits(void);
int intSortFunction(const void* av, const void* bv);
int _wordexp(const char * words, wordexp_t * pwordexp, int flags);
void _wordfree(wordexp_t *pwordexp);
unsigned long long getTick(void);
int removeDuplicates(unsigned* a, int na);
char *astrcat(char *dest, char *str); 

#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
