/* $Id$ */
#ifndef OBJECT_H
#define OBJECT_H
#include <stdio.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" 
{
#endif


#define TEST(A,B)  (B.mod != 0 && (A%B.mod )==B.offset)
#define TEST0(A,B)  (B != 0 && (A%B)==0)

enum OBJECTTYPES
{ DOUBLE, INT, SHORT, STRING, FILETYPE, FILELIST, LITERAL, I64, U64, WITH_UNITS };
enum MODE
{ ASCII, BINARY };
enum OBJECTACTION
{ ABORT_IF_NOT_FOUND, IGNORE_IF_NOT_FOUND, WARN_IF_NOT_FOUND };

typedef struct
{
	FILE *file;
	char *name;
} OBJECTFILE;

typedef struct object_st	/* object structure for parsing routine */
{
	char *name;
	char *objclass;
	char *value;
	char **valueptr;
	struct object_st *parent;
} OBJECT;


void object_set(const char *get, ... );
OBJECT *object_initialize(char *name, char *objclass, int size);
char *object_read(OBJECTFILE ofile);
int object_exists(const char *name, const char *objclass);
OBJECT *object_find(const char *name , const char *objclass);
OBJECT *object_find2(const char *, const char *, enum OBJECTACTION);
void object_free(OBJECT *);
/** Caller must free char* that are returned thru ptr. */
int object_get(OBJECT* obj,  const char *keyword, void *ptr, int type, int maxLength, const char* dvalue, ...);
/** Caller must free pointer that is returned through pptr.
 *  If pptr is a char** must also free individual strings (i.e. name[ii]) */
int object_getv(OBJECT* obj, const char *keyword, void *pptr, int type, int action_on_error);
void object_compilefilesubset(const char *filename, int first, int last);
void object_compilestring(char *string);
void object_replacekeyword(OBJECT *object,char *keyword,char* keywordvalue);
OBJECT* object_copy(OBJECT* );
int object_lineparse(char *, OBJECT *);
void object_compilefile(const char* filename);
void object_compile(void);
int object_testforkeyword(OBJECT*object, const char *name);
int object_get_with_units(OBJECT*object, char *name, void *ptr, int type,
			  int length, char *dvalue, char *unit_convert_from,
			  char *unit_convert_to);
char* object_keywordsAndValues(OBJECT* object);
OBJECT *object_longFind(const char *name, const char *objclass);
void object_list_purge();
OBJECTFILE object_fopen(const char *filename, char *mode);
void object_fclose(OBJECTFILE file);

void prune_spaces(char *string);
void zapChar(char zapped, char *string);
void trim(char *string);

void object_print_all(FILE* file);
void object_pretty_print(OBJECT* obj, FILE* file);

#ifdef WITH_MPI
void object_Bcast(int root, MPI_Comm comm);
#endif

#ifdef __cplusplus
}
#endif


#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
