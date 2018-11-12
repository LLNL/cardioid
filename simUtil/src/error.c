/* $Id$ */
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include "error.h"
#include "mpiUtils.h"

static char message[1024];
static const int message_size = sizeof(message)/sizeof(char);

void error_action(char *start, ...)
{
	va_list ap;
	enum ACTION action;
	int line; 
	char *string;
	strncpy(message, start, message_size);
	va_start(ap, start);
	string = va_arg(ap, char *);
	while (string != NULL)
	{
	  strncat(message, " ",message_size - strlen(message) - 1);
	  strncat(message, string,message_size - strlen(message) - 1);
		string = va_arg(ap, char *);
	}
	fprintf(stderr,"\nMessage:%s\n\n", message);
	string = va_arg(ap, char *);
	fprintf(stderr, "Error in routine %s\n", string);
	string = va_arg(ap, char *);
	fprintf(stderr, "in file %s ", string);
	line = va_arg(ap, int);
	fprintf(stderr, "at line %d\n", line);
	action = va_arg(ap, enum ACTION);
	switch (action)
	{
	case CONTINUE:
		if (getRank(0)==0) fprintf(stderr, "Program will CONTINUE.\n");
		va_end(ap);
		return;
	case ABORT:
		if (getRank(0)==0) fprintf(stderr, "Program will ABORT.\n");
		va_end(ap);
		abortAll(-1);
		break;
	  default:
		if (getRank(0)==0) fprintf(stderr, "Program will ABORT because of unknown action.\n");
		va_end(ap);
		abortAll(-1); 
	}
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
