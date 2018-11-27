/* $Id$ */
#include "match.h"
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include "ddcMalloc.h"
regex_t *compilepattern(char *pattern)
{
        regex_t *re;
	re = (regex_t *)ddcMalloc(sizeof(regex_t));
        if (regcomp(re, pattern, REG_EXTENDED) == 0)  return re; 
	return NULL ;             
}
char * findpattern1(char *string, regex_t  *re,char **end)
{
        int     status;
        regmatch_t pm[2];
	char *s;

	*end=string; 
        status = regexec(re, string, (size_t) 2, pm, 0);
        if (status == 0)   
	{
		s = string+pm[0].rm_so;
		*end = string+pm[0].rm_eo;
		return s;        
	}
        return NULL ;
}
char * findpattern(char *string, char *pattern,char **end)
{
        int     status;
        regex_t re;
        regmatch_t pm[2];
	char *s;

        if (regcomp(&re, pattern, REG_EXTENDED) != 0)  return NULL ;             /* report error */
	*end=string; 
        status = regexec(&re, string, (size_t) 2, pm, 0);
        if (status == 0)   
	{
		s = string+pm[0].rm_so;
		*end = string+pm[0].rm_eo;
        	regfree(&re);
		return s;        
	}
        regfree(&re);
        return NULL ;
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
