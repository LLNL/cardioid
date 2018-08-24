/* $Id$ */
#include "object.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <assert.h>
#include <stdint.h> // uint64_t and int64_t
#include "error.h"
#include "match.h"
#include "mpiUtils.h"
#include "units.h"
#if 1 
#include "ddcMalloc.h"

typedef struct field_st		/* field structure for parsing routine */
{
	int n;
	int size;
	int element_size;
	void *v;
	char *string; 
} FIELD;

typedef struct FileList_st
{
	char *name;
	int start, end;
} FileList;

typedef struct pack_buf_st
{
	int n, nobject,mobject;
	char *buffer;
} PACKBUF;

static FIELD object_parse(OBJECT *, const char *, int, const char *, int);


#define Realloc ddcRealloc
#define Malloc ddcMalloc
#define Free ddcFree
#else 
#define Realloc realloc 
#define Malloc malloc 
#define Free free 
#endif

#ifndef ALPHA
#define EXTRA_ARGS ,...
#else
#define EXTRA_ARGS
#endif

#define MAXKEYWORDS 4096
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))


static char *filenames[] = { "object.data", NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
static char *filename_list=NULL;
static int nfiles = 1;
static int nobject = 0, mobject = 0;
static int niobject = 0, miobject = 0;
static OBJECT **object = NULL;
static OBJECT **object_list = NULL;
void object_set(const char *get, ... )
{
	va_list ap;
	char *what, *string, *list, *file, *what_ptr,*list_ptr;
   list_ptr = NULL;
   what_ptr = NULL;
	void **ptr;
	va_start(ap, get);
	what = strdup(get);
	string = strtok_r(what, " ",&what_ptr);
	while (string != NULL)
	{
		ptr = va_arg(ap, void **);
		if (strcmp(string, "files") == 0)
		{
			list = strdup((char *)ptr);
			nfiles = 0;
			file = strtok_r(list, " ",&list_ptr);
			while (file != NULL)
			{
				filenames[nfiles++] = strdup(file);
				if (filename_list == NULL) filename_list = strdup(file); 
				else 
				{
				   filename_list = (char*) Realloc(filename_list,strlen(filename_list)+strlen(file)+3);
					strcat(filename_list,", "); 
					strcat(filename_list,file); 
				}
				file = strtok_r(NULL, " ",&list_ptr);
			}
		}
		string = strtok_r(NULL, " ",&what_ptr);
	}
	va_end(ap);
}

//NOTE: this function deletes all initialize objects except routine manager
void object_list_purge()
{
	for(int i =0; i <niobject; i++)
	{ //ROUTINEMANAGER
	         if (strcmp(object_list[i]->objclass, "ROUTINEMANAGER"))
		{
			object_free(object_list[i]);
		}
		niobject=1;
	}
}

void prune_spaces(char *string)
{
	int i,l;
	if (string == NULL) return;
	l=strlen(string);
	if (l== 0) return;
	int j=0; 
	for (i=0;i<=l;i++) 
	{
		if (j < i ) string[j] = string[i] ;
		if (string[i] != ' ') j++; 
	}
}

/** Remove all instances of the zapped character from the string */
void zapChar(char zapped, char *string)
{
	int i,l;
	if (string == NULL) return;
	l=strlen(string);
	if (l== 0) return;
	int j=0; 
	for (i=0;i<=l;i++) 
	{
		if (j < i ) string[j] = string[i] ;
		if (string[i] != zapped) j++; 
	}
}

void trim(char *string)
{
	int i, j, k, l;
	if (string == NULL) return;
	if (strlen(string) == 0) return;
	i = -1;
	while (string[++i] == ' ');
	if(string[i] == '\0')
	  k = 0;
	else {
	  j = strlen(string);
	  while (string[--j] == ' ');
	  k = 0;
	  for (l = i; l <= j; l++)
		 string[k++] = string[l];
	}
	string[k] = 0x0;
}

/** The maxLength parameter specifies the maximum storage size of the
 *  ptr.  If fewer than maxLength values are supplied (either in the
 *  input or in the dvalue then the unspecified elements are left
 *  unitialized.
 *
 *  The dvalue is used only when the keyword is unspecified in input. */
int object_get_with_units(OBJECT*object, char *name, void *ptr, int type,
			  int maxLength, char *dvalue,
			  char *unit_convert_from, char *unit_convert_to)
{
	int i,l;
	FIELD f;
	f = object_parse(object, name, type, dvalue,ABORT_IF_NOT_FOUND);
	if (f.n == 0 || maxLength <= 0) return 0;
	l = MIN(f.n, maxLength);
	bcopy(f.v, ptr, l*f.element_size);

	double *d = (double*) ptr; 
	char *units=f.string;
	if (*units == (char)0) units = unit_convert_from; 
	double convert_factor = units_convert(1.0,units,unit_convert_to); 
	for (i=0;i<l;i++) d[i] *= convert_factor; 
	free(f.string); 
	return f.n;
}

/** When object_get is called for STRING or LITERAL types, object_parse
 *  allocates memory (via strdup) to hold a copy of the string.  It is
 *  this copy that is eventually returned to the caller.  To avoid
 *  memory leaks it is necesary to free all string ptrs that are passed
 *  to object_get.  For other types no memory is allocated.
 *
 *  Example:
 *  char* tmp;
 *  object_get(obj, "foo", &tmp, STRING, 1, "NONE");
 *  foo_init(tmp);
 *  ddcFree(tmp);
 *
 *  If an array of strings is populated, free must be called for each
 *  individual string.
 *
 *  Example:
 *  char* tmp[3];
 *  int n = object_get(obj, "foo", &tmp, STRING, 3, "NONE");
 *  for (int ii=0; ii<n; ++ii)
 *  {
 *      foo_init(tmp[ii]);
 *      ddcFree(tmp[ii]);
 *   }
 *
 *  The maxLength parameter specifies the maximum storage size of the
 *  ptr.  If fewer than maxLength values are supplied (either in the
 *  input or in the dvalue then the unspecified elements are left
 *  uninitialized.
 *
 *  The dvalue is used only when the keyword is unspecified in input. */
int object_get(OBJECT*object, const char *name, void *ptr, int type, int maxLength, const char *dvalue, ...)
{
	FIELD f = object_parse(object, name, type, dvalue,ABORT_IF_NOT_FOUND);
	if (f.n == 0 || maxLength <= 0) return 0;
	int l = MIN(f.n, maxLength);
	bcopy(f.v, ptr, l*f.element_size);
	if (type == WITH_UNITS) 
	{
		va_list ap; 
		va_start(ap,dvalue); 
		char *unit_convert_from = va_arg(ap,char *); 
		char *unit_convert_to   = va_arg(ap,char *); 
		va_end(ap); 

		double *d = (double*) ptr; 
		char *units=f.string;
		if (*units == (char)0) units = unit_convert_from; 
		int unitsBad = (units_check(units,unit_convert_to) ||  units_check(units,unit_convert_from)); 
		assert(unitsBad  == 0 ); 
		double convert_factor = units_convert(1.0,units,unit_convert_to); 
		for (int i=0;i<l;i++) d[i] *= convert_factor; 
		free(f.string); 
	}

	return f.n;
}

/** This function allocates storage for pptr.  The caller must free this
 *  pointer to avoid leaks.  This is true for all datatypes.
 *
 *  In the case of string arrays, you must also free each string in the
 *  array (since object_parse makes copies of the strings in the object
 *  database).
 *
 *  Example:
 *  char** tmp;
 *  int n = object_getv(obj, "foo", (void*)&tmp, STRING, IGNORE_IF_NOT_FOUND);
 *  for (int ii=0; ii<n; ++ii)
 *  {
 *      foo_init(tmp[ii]);
 *      ddcFree(tmp[ii]);
 *   }
 *   ddcFree(tmp);
 *
 *   This is just like object_get with the additional requirement that
 *   you must call free on the pointer passed in to object_getv.
 */
int object_getv(OBJECT*object, const char *name, void *pptr, int type, int  action_on_error)
{
	FIELD f;
	void *ptr;
	f = object_parse(object, name, type, NULL, action_on_error);
	if (f.n > 0)
	{
		ptr = Malloc(f.size);
		bcopy(f.v, ptr, f.size);
		*(void**)pptr = (void *)ptr;
	}
	else
	{
		*(void**)pptr = (void *)NULL;
	}
	return f.n;
}

int object_test(OBJECT*object, char *name, char **value_ptr, int *nvalue_ptr)
{
	static int lbuff = 0;
	static char *buffer = NULL, *line = NULL;
	static int nkeyword = 0;
	static char *keyword = NULL;
	static char *value;
	static int nvalue;
	int found=0;
	int nname;
	char *name_list[16], *name_save;
	char *string;
	int i, lstring;
	char *ptr;
	value = NULL;
	nvalue = 0;
	if (value_ptr != NULL) value = *value_ptr;
	if (nvalue != 0) nvalue = *nvalue_ptr;
	name_save = strdup(name);
	ptr = strtok(name_save, ";");
	nname = 0;
	while (ptr != NULL && nname < 16)
	{
		trim(ptr);
		name_list[nname] = ptr;
		ptr = strtok(NULL, ";");
		nname++;
	}
	string = object->value;
	lstring = strlen(string) + 1;
	if (lstring > lbuff)
	{
		lbuff = lstring;
		buffer = (char*) Realloc(buffer, lstring*sizeof(char));
		line = (char*) Realloc(line, lstring*sizeof(char));
	}
	strcpy(buffer, string);
	ptr = strtok(buffer, ";");
	while (ptr != NULL)
	{
		strcpy(line, ptr);
		trim(line);
		ptr = strchr(line, '=');
		*ptr = 0x0;
		while (nkeyword < (int)strlen(line) + 1)
		{
			nkeyword += 256;
			keyword = (char*) Realloc(keyword, nkeyword);
		}
		strncpy(keyword, line, nkeyword);
		trim(keyword);
		ptr++;
		while (nvalue < (int)strlen(ptr) + 1)
		{
			nvalue += 256;
			value = (char*) Realloc(value, nvalue);
		}
		strncpy(value, ptr, nvalue);
		trim(value);
		found = 0;
		for (i = 0; i < nname; i++)
		{
			if (strcmp(name_list[i], keyword) == 0)
			{
				found = 1;
				break;
			}
		}
		if (found) break;
		ptr = strtok(NULL, ";");
	}
	Free(name_save);
	if (value_ptr != NULL) *value_ptr = value;
	else Free(value);
	if (nvalue_ptr != NULL) *nvalue_ptr = nvalue;
	return found;
}

FILE *object_fileopen(char *filename)
{
	char *ptr, *name;
	FILE *file;
	int first, last;
	name = strdup(filename);
	ptr = index(name, '@');	/* filename */
	if (ptr != NULL)
	{
		*ptr = 0x0;	/* add a string termination to filename */
		file = fopen(name, "r");
		sscanf(ptr + 1, "%d-%d", &first, &last);
		fseek(file, first, SEEK_SET);
	}
	else
	{
		file = fopen(name, "r");
	}
	Free(name);
	return file;
}

void object_compilevalue(OBJECT*object)
{
	enum MODE
	{ FIRST, LAST };
	static int lbuff = 0;
	static char *buffer = NULL;
	char *keyword, *op, *value, *kvalue, *ptr;
	struct
	{
		char *keyword, *op, *value;
	} list[MAXKEYWORDS];
	int nlist;
	char *string, *opptr;
	int i, j, lstring;
	int mode = LAST;
	string = object->value;
	lstring = strlen(string) + 1;
	if (lstring > lbuff)	/* Create temp character  array */
	{
		lbuff = lstring;
		buffer = (char*) Realloc(buffer, lstring*sizeof(char));
	}
	strcpy(buffer, string);
	nlist = 0;
	ptr = strtok(buffer, ";");
	while (ptr != NULL)	/*  loop through   "keyword=value ; */
	{
		opptr = strchr(ptr, '=');
		keyword = ptr;
		value = opptr + 1;
		op = "=";
		*opptr = 0x0;
		if (*(opptr - 1) == '+')
		{
			op = "+=";
			*(opptr - 1) = 0x0;
		}
		trim(keyword);
		trim(value);
		list[nlist].keyword = keyword;
		list[nlist].op = op;
		list[nlist].value = value;
		nlist++;
		ptr = strtok(NULL, ";");
	}
	object->value[0] = 0x0;
	for (i = 0; i < nlist; i++)
	{
		keyword = list[i].keyword;
		if (keyword != NULL)
		{
			strcat(object->value, keyword);
			strcat(object->value, "=");
			kvalue = object->value + strlen(object->value);
			strcpy(kvalue, list[i].value);
			for (j = i + 1; j < nlist; j++)
			{
				if (list[j].keyword != NULL)
				{
					if (strcmp(list[j].keyword, keyword) == 0)
					{
						list[j].keyword = NULL;
						if (mode == LAST)
						{
							if (strcmp(list[j].op, "=") == 0) strcpy(kvalue, list[j].value);
						}
						if (strcmp(list[j].op, "+=") == 0)
						{
							strcat(kvalue, " ");
							strcat(kvalue, list[j].value);
						}
					}
				}
			}
			strcat(object->value, ";");
		}
	}
	if (object->valueptr) *object->valueptr = object->value;
}

/** Note that for STRING objects, this function calls strdup to make
 *  copies of the string(s).  Hence the FIELD returned by this function
 *  may contain newly allocated memory that someone will have to free
 *  (to avoid leaks).
 */
FIELD object_parse(OBJECT*object, const char *name, int type, const char *dvalue,int action_on_error)
{
	static int lbuff = 0, msize = 0;
	static char *buffer = NULL, *line = NULL, *tail, *sep;
	static union
	{
		double *d;
		int *i;
      uint64_t* u64;
      int64_t* i64;
		short *sh;
		char **s;
		FILE **file;
		FileList *filelist;
		void *v;
	} v;
	static int nvalue = 0;
	static int nkeyword = 0;
	static char *value = NULL;
	static char *keyword = NULL;
	FILE *file;
	int nname;
	char *name_list[16], *name_save;
	char *string;
	FIELD f;
	int i, lstring, nv, size, element_size=0, found, first, last;
	char *ptr, *vptr, *eptr;
	name_save = strdup(name);
	ptr = strtok(name_save, ";");
	nname = 0;
	found = 0;
	f.string = NULL; 
	while (ptr != NULL && nname < 16)	/* create list of names to parse */
	{
		trim(ptr);
		name_list[nname] = ptr;
		ptr = strtok(NULL, ";");
		nname++;
	}
	string = object->value;
	lstring = strlen(string) + 1;
	if (lstring > lbuff)	/* Create some temp character  array */
	{
		lbuff = lstring;
		if (lbuff < 512) lbuff = 512;
		buffer = (char*) Realloc(buffer, lbuff*sizeof(char));
		line = (char*) Realloc(line, lbuff*sizeof(char));
	}
	strcpy(buffer, string);
	ptr = strtok(buffer, ";");
	while (ptr != NULL)	/*  loop through   "keyword=value ; "  fields  and check if keyword matches one of the name fields */
	{
		strcpy(line, ptr);
		trim(line);
		ptr = strchr(line, '=');
		*ptr = 0x0;
		while (nkeyword < (int)strlen(line) + 1)
		{
			nkeyword += 256;
			keyword = (char*) Realloc(keyword, nkeyword);
		}
		strncpy(keyword, line, nkeyword);
		trim(keyword);
		ptr++;
		while (nvalue < (int)strlen(ptr) + 1)
		{
			nvalue += 256;
			value = (char*) Realloc(value, nvalue);
		}
		strncpy(value, ptr, nvalue);
		trim(value);
		found = 0;
		for (i = 0; i < nname; i++)
		{
			if (strcmp(name_list[i], keyword) == 0)
			{
				found = 1;
				break;
			}
		}
		if (found) break;
		ptr = strtok(NULL, ";");
	}
	if (found == 0)
	{
		if (action_on_error == ABORT_IF_NOT_FOUND && dvalue == NULL) error_action("Unable to locate ", name, " in object ", object->name, ERROR_IN("object_parse", ABORT));
		if (action_on_error == IGNORE_IF_NOT_FOUND)
		{
			f.size=0;
			f.element_size=0;
			f.n = 0;
			f.v = NULL;
			return f;
		}
		while (nvalue < (int)strlen(dvalue) + 1)
		{
			nvalue += 256;
			value = (char*) Realloc(value, nvalue);
		}
		strncpy(value, dvalue, nvalue);
	}

	trim(value);
	if (type == LITERAL)
	{
		vptr = value;
		v.s[0] = strdup(vptr);
		f.n = 1;
		f.v = v.v;
		f.size = f.element_size = sizeof(char *);
		return f;
	}
	tail = value;
	if (*tail == (char)'"') sep = "\"";
	else
		sep = " ";	/*Check if value start with "  */
	vptr = strtok_r(value, sep, &tail);
	nv = size = 0;
	if (vptr != NULL) if (!strncmp(vptr, "$B", 2))
		{
			void *ptr;
			sscanf(vptr + 2, "%d-%p", &size, &ptr);
			vptr = NULL;
			nv = 1;
		}
	while (vptr != NULL)
	{
		if (size + 256 > msize)
		{
			msize += 1024;
			v.v = Realloc(v.v, msize);
		}
		switch (type)
		{
		case DOUBLE:
			if (!strncmp(vptr, "0x", 2))
			{
				v.i[nv + 1] = strtol(vptr + 10, &eptr, 16);
				vptr[10] = '\0';
				v.i[nv] = strtol(vptr + 2, &eptr, 16);
			}
			else
			{
				v.d[nv] = strtod(vptr, &eptr);
			}
			element_size = sizeof(double);
			break;
		case WITH_UNITS:
			{
				double d = strtod(vptr, &eptr);
				v.d[nv] =d ; 
				element_size = sizeof(double);
			}
			break;
		case INT:
			v.i[nv] = strtol(vptr, &eptr, 0);
			element_size = sizeof(int);
			break;
		case U64:
			v.u64[nv] = strtoull(vptr, &eptr, 0);
			element_size = sizeof(uint64_t);
			break;
		case I64:
			v.i64[nv] = strtoll(vptr, &eptr, 0);
			element_size = sizeof(int64_t);
			break;
		case SHORT:
			v.sh[nv] = strtol(vptr, &eptr, 0);
			element_size = sizeof(short);
			break;
		case STRING:
			v.s[nv] = strdup(vptr);
			element_size = sizeof(char *);
			break;
		case FILETYPE:
			ptr = index(vptr, '@');	/* filename */
			if (ptr != NULL)
			{
				*ptr = 0x0;	/* add a string termination to filename */
				file = fopen(vptr, "r");
				sscanf(ptr + 1, "%d-%d", &first, &last);
				fseek(file, first, SEEK_SET);
			}
			else
			{
				file = fopen(vptr, "r");
			}
			element_size = sizeof(FILE *);
			v.file[nv] = file;
			break;
		case FILELIST:
			ptr = index(vptr, '@');	/* filename */
			if (ptr != NULL)
			{
				*ptr = 0x0;	/* add a string termination to filename */
				sscanf(ptr + 1, "%d-%d", &first, &last);
			}
			else
			{
				file = fopen(vptr, "r");
				first = 0;
				last = -1;
			}
			v.filelist[nv].name = strdup(vptr);
			v.filelist[nv].start = first;
			v.filelist[nv].end = last;
			element_size = sizeof(FileList);
			break;
		default:
			break;
		}
		nv++;
		size += element_size;
		if (tail == NULL || strlen(tail) == 0) break; 
		trim(tail);
		if (*tail == (char)'"') sep = "\"";
		else
			sep = " ";
		vptr = strtok_r(NULL, sep, &tail);
	}
	if (type == WITH_UNITS)
	{
		if (eptr == vptr) 
		{
			nv--; 
			size -= element_size; 
		}
		f.string = strdup(eptr); 
	}
	f.n = nv;
	f.v = v.v;
	f.size = size;
	f.element_size = element_size;
	return f;
}

void object_compilestring(char *string)
{
	OBJECT obj;
	char *line; 
	int rc, i, l;
	line = strdup(string); 
	rc = object_lineparse(line, &obj);
	if (rc < 2)
	{
		for (i = 0; i < nobject; i++)
		{
			if (!strcmp(object[i]->objclass, obj.objclass) && !strcmp(object[i]->name, obj.name)) break;
		}
		if (i == nobject)   //object is not in database
		{
			nobject++;
			if (mobject < nobject) 
			{
				mobject += 100;
				object = (OBJECT**) Realloc(object, mobject*sizeof(OBJECT*));
			}
         object[i] = (OBJECT *)ddcMalloc(sizeof(OBJECT)); 
			object[i]->name = strdup(obj.name);
			object[i]->objclass = strdup(obj.objclass);
			object[i]->value = strdup(obj.value);
			object[i]->valueptr = NULL;;
		}
		else
		{
			l = strlen(object[i]->value) + strlen(obj.value) + 1;
			object[i]->value = (char*) Realloc(object[i]->value, l);
			strcat(object[i]->value, obj.value);
		}
	}
	Free(obj.name);
	Free(obj.objclass);
	Free(obj.value);
	Free(line); 
}
void object_compilefilesubset(const char *filename, int first, int last)
{
	char *line;
	OBJECTFILE file;
	int i, n;
	file = object_fopen(filename, "r");
	n=0; 
	if (file.file != NULL)
	{
		while ((line = object_read(file)) != NULL )
		{
			if (n>=first)
			{
				object_compilestring(line);
			}
			if ( ++n > last ) break; 
		}
		object_fclose(file);
	}
	for (i = 0; i < nobject; i++) object_compilevalue(object[i]);
}
void object_compilefile(const char* filename) { object_compilefilesubset(filename,0,0x7fffffff); }

void object_replacekeyword(OBJECT *object,char *keyword,char* keywordvalue)
{
	char value[4096];
	int l; 
	sprintf(value,"%s=%s;",keyword,keywordvalue);
	l = strlen(object->value) + strlen(value) + 1;
	object->value = (char*) Realloc(object->value, l);
	strcat(object->value, value);
	object_compilevalue(object);
}
void object_compile(void)
{
	int k;
	if (mobject < nobject)
	{
		mobject += 100;
		object = (OBJECT**) Realloc(object, mobject*sizeof(OBJECT*));
	}
	for (k = 0; k < nfiles; k++)
	{
		object_compilefile(filenames[k]);
	}
}

void object_reset(OBJECT*target)
{
   OBJECT obj;
	char *line;
	OBJECTFILE file;
	int rc, l, k;
	if (target->name == NULL || target->objclass == NULL) return;
	*(target->value) = '\0';
	for (k = 0; k < nfiles; k++)
	{
		file = object_fopen(filenames[k], "r");
		if (file.file != NULL)
		{
			while ((line = object_read(file)) != NULL)
			{
				rc = object_lineparse(line, &obj);
				if (rc < 2)
				{
					if (!strcmp(target->objclass, obj.objclass) && !strcmp(target->name, obj.name))
					{
						l = strlen(target->value) + strlen(obj.value) + 1;
						target->value = (char*) Realloc(target->value, l);
						strcat(target->value, obj.value);
					}
				}
				Free(obj.name);
				Free(obj.objclass);
				Free(obj.value);
			}
			object_fclose(file);
		}
	}
}

int object_exists(const char *name, const char *objclass)
{
	int i;
	for (i = 0; i < nobject; i++)
	{
		if (!strcmp(object[i]->objclass, objclass) && !strcmp(object[i]->name, name)) return 1;
	}
	return 0; 
}

OBJECT *object_find(const char *name, const char *objclass)
{
	int i;
	for (i = 0; i < nobject; i++)
	{
		if (!strcmp(object[i]->objclass, objclass) && !strcmp(object[i]->name, name)) return object[i];
	}
	if (i == nobject) error_action("Unable to locate object <<", name,">> of class <<",objclass,">> in object files <<",filename_list ,">>", ERROR_IN("object_find", ABORT));
	return object[i];
}
OBJECT *object_longFind(const char *name, const char *objclass)
{
	int i;
	for (i = 0; i < niobject; i++)
	{
		if (!strcmp(object_list[i]->objclass, objclass) && !strcmp(object_list[i]->name, name)) return object_list[i];
	}
	if (i == niobject) error_action("Unable to locate object <<", name,">> of class <<",objclass,">> in object_list ", ERROR_IN("object_longFind", ABORT));
	return object_list[i];
}

OBJECT *object_find2(const char *name, const char *objclass, enum OBJECTACTION action)
{
	int i;
	for (i = 0; i < nobject; i++)
	{
		if (!strcmp(object[i]->objclass, objclass) && !strcmp(object[i]->name, name)) break;
	}
	if (i == nobject)
	{
		switch (action)
		{
		case ABORT_IF_NOT_FOUND:
			error_action("Unable to locate object ", name, "in data files", ERROR_IN("object_find2", ABORT));
			break;
		case WARN_IF_NOT_FOUND:
			error_action("Unable to locate object ", name, "in data files", ERROR_IN("object_find2", CONTINUE));
			return NULL;
			break;
		case IGNORE_IF_NOT_FOUND:
			return NULL;
			break;
		}
		return NULL;
	}
	return object[i];
}


OBJECT *object_initialize(char *name, char *objclass, int size)
{
	OBJECT *object_short, *object_long;
	object_long = (OBJECT*) Malloc(size);
	object_short = object_find(name, objclass);
	object_long->name = object_short->name;
	object_long->objclass = object_short->objclass;
	object_long->value = object_short->value;
	object_long->valueptr = &object_short->value;
	object_short->valueptr = &object_long->value;
	if (miobject <= niobject)
	{
		miobject += 100;
		object_list = (OBJECT**) Realloc(object_list, miobject*sizeof(OBJECT *));
	}
	object_list[niobject] = object_long;
	niobject++;
	return object_long;
}

OBJECT *object_find1(char *name, char *type)
{
	OBJECTFILE file;
	int rc, nfound, l, k;
	OBJECT *object = NULL;
	char *line, *value;
	if (object == NULL) object = (OBJECT*) Malloc(sizeof(OBJECT));
	value = (char*) Malloc(1);
	*value = 0x0;
	nfound = 0;
	for (k = 0; k < nfiles; k++)
	{
		file = object_fopen(filenames[k], "r");
		while ((line = object_read(file)) != NULL)
		{
			rc = object_lineparse(line, object);
			if (rc < 2 && (strcmp(object->objclass, type) == 0) && (strcmp(object->name, name) == 0))
			{
				nfound++;
				l = strlen(value) + 1;
				if (object->value != NULL)
				{
					l += strlen(object->value);
					value = (char*) Realloc(value, l);
					strcat(value, object->value);
				}
			}
			Free(object->name);
			Free(object->objclass);
			Free(object->value);
		}
		object_fclose(file);
	}
	if (nfound > 0)
	{
		object->name = strdup(name);
		object->objclass = strdup(type);
		object->value = strdup(value);
		return object;
	}
	error_action("Unable to locate object ", name, "in data files", ERROR_IN("object_find1", ABORT));
	return NULL;
}
char *object_read(OBJECTFILE ofile)
{
	FILE *file;
	static char *line = NULL;
	static int nline = 0;
	int n, c, len, nitems, offset;
	file = ofile.file;
	n = 0;
	c = getc(file);
	while (!feof(file))
	{
		while (nline < n + 5)
		{
			nline += 256;
			line = (char*) Realloc(line, nline);
		}
		switch (c)
		{
		case '"':
			line[n++] = (char)c;
			c = getc(file);
			int quoteStart = n;
			while (c != '"')
			{
				if (nline < n + 5 )
				{
					nline += 256;
					line = (char*) Realloc(line, nline);
				}
				line[n++] = (char)c;
				c = getc(file);
				if (c==EOF) 
				{
					int lEnd = MIN(quoteStart+1,n);
					line[lEnd]='\0';
					error_action("String termination  for string in line starting with:\n\n****" , line  ,"****\n\n not found in object file : ", ofile.name, ERROR_IN("object_read", ABORT));
				}
			}
			break;
		case '/':
			c = getc(file);
			int commentStart = n;
			int m=n; 
			switch (c)
			{
			  case '*':
				c = getc(file);
				char clast = c; 
				while (c != EOF)
				{
					if (c == '/'  && clast == '*') {c = ' '; break;}
					clast = c; 
					c = getc(file);
					m++; 
				}
				if (c==EOF) 
				{
					int lEnd = MIN(commentStart+2,m);
					line[lEnd]='\0';
					error_action("Comment termination of (/*) in line starting with:\n\n****",line,"****\n\n not found in object file : ", ofile.name, ERROR_IN("object_read", ABORT));
				}
				break;
			  case '/':
				clast = c = getc(file);
				while (c != EOF)
				{
					if (c == '\n'  && clast != '\\') {c = ' '; break;}
					clast = c; 
					c = getc(file);
					m++; 
				}
				if (c==EOF) 
				{
					int lEnd = MIN(commentStart+2,m);
					line[lEnd]='\0';
					error_action("Comment termination of (//) in line starting with:\n\n****",line,"****\n\n not found in object file : ", ofile.name, ERROR_IN("object_read", ABORT));
				}
				break;
			  default:
				ungetc(c, file);
				c = '/';
				break;
			}
			break;
		case '\n':
			c = ' ';
			break;
		case '\t':
			c = ' ';
			break;
		case ' ':
			break;
		case '{':
			break;
		case '}':
			break;
		case '\\':
			c = (char)getc(file);
			switch (c)
			{
			case 'n':
				c = '\n';
				break;
			case 't':
				c = '\t';
			}
			break;
		case '$':
			c = (char)getc(file);
			switch (c)
			{
			case 'S':
				nitems = fscanf(file, "%d", &offset);
				c = (char)getc(file);
				while (!isgraph(c))
				{
					if (feof(file)) return NULL;
					c = (char)getc(file);
				}
				long first = ftell(file) - 1;
				if (nitems > 0) fseek(file, offset, SEEK_CUR);
				while (c != ';' && c != '}')
				{
					if (feof(file)) return NULL;
					c = (char)getc(file);
				}
				long last = ftell(file) - 1;
				len = strlen(ofile.name);
				while (nline < n + len + 1 + 2 + 64)
				{
					nline += 256;
					line = (char*) Realloc(line, nline);
				}
				len = sprintf(line + n, "%s@%ld-%ld", ofile.name, first, last);
				n += len;
				break;
			case 'B':
				line[n++] = '$';
				break;
			}
			break;
		case ';':
			break;
		default:
			break;
		}
		line[n++] = (char)c;
		if (c == '}')
		{
			line[n++] = 0x0;
			trim(line);
			return line;
		}
		c = (char)getc(file);
	}
	return NULL;
}

void object_free(OBJECT*object)
{
   if (object)
   {
      Free(object->name);
      Free(object->objclass);
      Free(object->value);
      Free(object);
   }
}
void object_check(char *objectLine)
{
#define WORD  "[A-Za-z0-9_]+"
#define KEYWORD  "[A-Za-z0-9_-]+"
#define WS  " +" //white space
#define OWS " *" //optional white space
#define STRING "\"[^\"]+\"" //string
#define OP "(=|\\+=)"
	char *pattern = WORD WS WORD OWS "\\{ *(" KEYWORD OWS OP  OWS "([^;=}]+|" STRING ")+" OWS ";" OWS ")*}";
	char *end;
	char *start = findpattern(objectLine,pattern,&end);
	if (start != objectLine || strlen(end) >0  )
	{
		printf("Error in object:\n\n");
		printf("*** %s ***\n\n",objectLine);
		printf("Progam will terminate\n"); 
		abortAll(1); 
	}
	//printf("start =%s\n end=%s\n patten=%s\n",start,end,pattern); 
}

int object_lineparse(char *line, OBJECT*object)
{
	char *tok,*last;
   last = NULL;
	int rc, l;
	rc = 0;
	trim(line);
	object_check(line);
	tok = index(line,'{');
//	char *end = line+strlen(line); 
	*tok = '\0';
	tok++; 
	l = strlen(tok);
	if (tok[l-1] != '}' )  rc += 8; 
	tok[l-1] = '\0'; 
	object->value = strdup(tok);
	tok = strtok_r(line, " ",&last);
	object->name = strdup(tok);
	tok = strtok_r(NULL, " ",&last);
	object->objclass = strdup(tok);

	trim(object->name);
	trim(object->objclass);
	trim(object->value);
	if (object->value != NULL)
	{
		l = strlen(object->value);
		if (l>0 && object->value[l - 1] != ';')
		{
			object->value = (char*) Realloc(object->value, l + 2);
			strcat(object->value, ";");
		}
	}
	if (object->name == NULL) rc += 4;
	if (object->objclass == NULL) rc += 2;
	if (object->value == NULL) rc += 1;
	return rc;
}

char* object_keywordsAndValues(OBJECT* object)
{
	return object->value;
}

OBJECTFILE object_fopen(const char *filename, char *mode)
{
	OBJECTFILE file;
	file.file = fopen(filename, mode);
	if (file.file == NULL)
	{
		char *msg;
		int msg_length; 
		msg_length = strlen(filename)+256; 
		msg = (char*) Malloc(msg_length); 
		sprintf(msg, "%d: Error opening file=%s with mode %s from object_fopen", getRank(0), filename, mode);
		perror(msg);
		Free(msg); 
	}
	file.name = strdup(filename);
	return file;
}

void object_fclose(OBJECTFILE file)
{
	fclose(file.file);
	Free(file.name);
}

int object_register(char *name, char *type, int itype, void *address)
{
	return 0;
}

int modeindex(char *mode)
{
	if (!strcasecmp(mode, "FORMATTED")) return ASCII;
	if (!strcasecmp(mode, "ASCII")) return ASCII;
	if (!strcasecmp(mode, "BINARY")) return BINARY;
	if (!strcasecmp(mode, "UNFORMATTED")) return BINARY;
	return ASCII;
}
static void object_pack(PACKBUF *buf)
{
	int i,l,n;
	n=0;
	for (i=0;i<nobject;i++)
	{
		n+= l =  strlen(object[i]->name)+1 ; buf->buffer=(char*) Realloc(buf->buffer,n);memcpy(buf->buffer+n-l,object[i]->name,l); 
		n+= l =  strlen(object[i]->objclass)+1; buf->buffer=(char*) Realloc(buf->buffer,n);memcpy(buf->buffer+n-l,object[i]->objclass,l); 
		n+= l =  strlen(object[i]->value)+1; buf->buffer=(char*) Realloc(buf->buffer,n);memcpy(buf->buffer+n-l,object[i]->value,l);
	}
	buf->n =n; 
	buf->nobject =nobject; 
	buf->mobject =mobject; 
}
static void object_unpack(PACKBUF *buf)
{
	int i,l;
	char *ptr; 
	nobject=buf->nobject ; 
	mobject=buf->mobject ; 
	object = (OBJECT**) Realloc(object, mobject*sizeof(OBJECT*));
	ptr = buf->buffer; 
	for (i=0;i<nobject;i++)
	{
      object[i] = (OBJECT *)ddcMalloc(sizeof(OBJECT)); 
		l = strlen(ptr)+1 ; object[i]->name = strdup(ptr);ptr += l;
		l = strlen(ptr)+1 ; object[i]->objclass = strdup(ptr);ptr += l;
		l = strlen(ptr)+1 ; object[i]->value = strdup(ptr);ptr += l;
		object[i]->valueptr = NULL;;
	}
}

void object_print_all(FILE* file)
{
   for (int ii=0; ii<nobject; ++ii)
      object_pretty_print(object[ii], file);
}

void object_pretty_print(OBJECT* obj, FILE* file)
{
   fprintf(file, "%s %s\n", obj->name, obj->objclass);
   fprintf(file, "{\n");
   char* line = strdup(obj->value);
   char* lasts = NULL;
   char* pair = strtok_r(line, ";", &lasts);
   while (pair)
   {
      fprintf(file, "   %s;\n", pair);
      pair = strtok_r(NULL, ";", &lasts);
   }
   fprintf(file, "}\n\n");
   Free(line);
}



#ifdef WITH_MPI
void object_Bcast(int root, MPI_Comm comm)
{
	PACKBUF buf; 
	buf.buffer=NULL;
	int myRank;
	MPI_Comm_rank(comm, &myRank);
	if (myRank==root) object_pack(&buf);
	MPI_Bcast(&buf,3,MPI_INT,root,comm); 
	buf.buffer=Realloc(buf.buffer,buf.n);
	MPI_Bcast(buf.buffer,buf.n,MPI_CHAR,root,comm); 
	if (myRank!=root) object_unpack(&buf);
	Free(buf.buffer);
}
#endif 

OBJECT* object_copy(OBJECT* a)
{
   OBJECT* b = ddcMalloc(sizeof(OBJECT));
   b->name = strdup(a->name);
   b->objclass = strdup(a->objclass);
   b->value = strdup(a->value);
   b->valueptr = &b->value;
   b->parent = a->parent;
   return b;
}


int object_testforkeyword(OBJECT*object, const char *name)
{
	char *buffer;
	char *keyword;
	char *ptr;
	buffer = strdup(object->value); 
	ptr = strtok(buffer, ";");
	while (ptr != NULL)	/*  loop through   "keyword=value ; "  */
	{
		keyword = ptr; 
		ptr = strchr(ptr, '=');
		*ptr = 0x0;
		trim(keyword);
		if (strcmp(name, keyword) == 0) {Free(buffer); return 1; }
		ptr = strtok(NULL, ";");
	}
	Free(buffer); 
	return 0; 
}




/* Local Variables: */
/* tab-width: 3 */
/* End: */
