// $Id$

#ifndef MATCH_H
#define MATCH_H

#include <regex.h>

regex_t* compilepattern(char* pattern);
char* findpattern1(char* string, regex_t* re, char** end);
char* findpattern(char* string, char* pattern, char** end);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
