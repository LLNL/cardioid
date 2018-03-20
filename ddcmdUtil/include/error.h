/* $Id$ */
#ifndef ERROR_H
#define ERROR_H
#ifdef __cplusplus
extern "C" {
#endif

#define ERROR_IN(name,action)  NULL,name,__FILE__,__LINE__,action
#define LOCATION(name)  location(name,__FILE__,__LINE__)
char *location(char *string, char *file, int linenumber);
enum ACTION
{ CONTINUE, ABORT };
void error_action(char *start, ...);

#ifdef __cplusplus
}
#endif

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
