// $Id$

#include "ioUtils.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "ddcMalloc.h"
#include "utilities.h"

static int needsSwap=0;

int filetest(const char* filename, int type)
{
   struct stat buf;
   int rc = stat(filename, &buf);
   if (rc == -1 && errno == ENOENT)
   {
      printf("FILE: %s does not exist\n", filename);
      rc = 1;
   }
   if (rc == 0 && !(buf.st_mode & type))
   {
      printf("FILE: %s is not the specified type\n", filename);
      rc = 2;
   }
   return rc;
}


int DirTestCreate(const char *dirname)
{
	int mode = 0775;
	struct stat statbuf;
	char line[1024];
	int rc;
	rc = stat(dirname, &statbuf);
	if (rc == -1 && errno == ENOENT)
	{
		sprintf(line,"Creating Directory: %s\n", dirname);
		timestamp(line);
		rc = mkdir(dirname, mode);
		rc = stat(dirname, &statbuf);
	}
	if (rc != 0 || !(statbuf.st_mode & S_IFDIR))
	{
		printf("Can't Stat the Directory %s\n", dirname);
		printf("%d %x %x %x\n", rc, statbuf.st_mode, S_IFDIR, statbuf.st_mode & S_IFDIR);
	}
	return rc;
}

void
ioUtils_setSwap(unsigned endianKey)
{
   unsigned keylocal;
   memcpy(&keylocal,"1234",4);
   needsSwap = (endianKey != keylocal);
}

unsigned long long
mkInt(const unsigned char* data, const char* fieldType)
{
   unsigned long long result=0, i8;
   unsigned i4;
   int invalid = 1;

   if (fieldType[0] == 'u')
   {
      invalid = 0;
      long fieldLength = strtol(fieldType+1, NULL, 10);

      switch (fieldLength)
      {
	case 4:
	 copyBytes(&i4, data, 4);
	 if (needsSwap)  endianSwap(&i4, 4);
	 result = i4;
	 break;
	case 8:
	 copyBytes(&i8, data, 8);
	 if (needsSwap)  endianSwap(&i8, 8);
	 result = i8;
	 break;
	default:
	 invalid = 1;
      }
   }
   
   if (fieldType[0] == 'b')
   {
      invalid = 0;
      result = bFieldUnpack(data, fieldType);
   }

   if (invalid)
   {
      printf("ERROR: Invalid field type (%s) in mkInt\n", fieldType);
      exit(1);
   }

   return result;
}

double
mkDouble(const unsigned char* data, const char* fieldType)
{
   double result, f8;
   float f4;
   int invalid = 1;

   if (fieldType[0] == 'f')
   {
      invalid = 0;
      long fieldLength = strtol(fieldType+1, NULL, 10);

      switch (fieldLength)
      {
	case 4:
	 copyBytes(&f4, data, 4);
	 if (needsSwap)  endianSwap(&f4, 4);
	 result = f4;
	 break;
	case 8:
	 copyBytes(&f8, data, 8);
	 if (needsSwap)  endianSwap(&f8, 8);
	 result = f8;
	 break;
	default:
	 invalid = 1;
      }
   }
   
   if (invalid)
   {
      printf("ERROR: Invalid field type (%s) in mkDouble\n", fieldType);
      exit(1);
   }

   return result;
}

double
mkValue(const unsigned char* data, const char* fieldType)
{
   if (fieldType[0] == 'f') return mkDouble(data, fieldType);
   if (fieldType[0] == 'b' || fieldType[0] == 'u') return mkInt(data, fieldType);
   assert(1==0);
    return 0.0; 
}


void
copyBytes(void* to, const void* from, unsigned size)
{
   for (unsigned ii=0; ii<size; ++ii)
      *(((char*) to)+ii) = *(((char*) from)+ii);
}


unsigned long long
bFieldUnpack(const unsigned char* buf, const char* fieldType)
{
   assert(fieldType[0] == 'b');
   long fieldLength = strtol(fieldType+1, NULL, 10);
   unsigned long long result = 0;

   for (long ii=0; ii<fieldLength; ++ii)
   {
      result *= 256;
      result += *(buf+ii);
   }

   return result;
}

void
bFieldPack(unsigned char* buf, unsigned size, unsigned long long in)
{
   for (unsigned ii=0; ii<size; ++ii)
   {
      buf[(size-1) - ii] = in%256;
      in /= 256;
   }
}

unsigned bFieldSize(unsigned long long i)
{
   unsigned result=0;
   do
   {
      i /= 256;
      ++result;
   }
   while (i>0);
   return result;
      
}

   

unsigned* makeOffsets(char** fieldType, unsigned nFields)
{
   unsigned* offset = (unsigned *) ddcMalloc((nFields+1)*sizeof(unsigned));
   offset[0] = 0;
   for (unsigned ii=0; ii<nFields; ++ii)
   {
      unsigned fieldSize = strtoul( (fieldType[ii]+1), NULL, 10);
      if (fieldSize == 0)
      {
         printf("ERROR");
         exit(1);
      }
      offset[ii+1] = offset[ii] + (unsigned) fieldSize;
   }

   return offset;
}

enum IO_FIELDS* makeFieldMap(char** fieldNames, unsigned nFields)
{
   enum IO_FIELDS* map = (enum IO_FIELDS*) ddcMalloc(nFields*sizeof(enum IO_FIELDS));
   for (unsigned ii=0; ii<nFields; ++ii)
   {
      map[ii] = IOF_NONE;
      if      (strcasecmp(fieldNames[ii], "checksum") == 0) map[ii] = IOF_CHECKSUM;
      else if (strcasecmp(fieldNames[ii], "id")       == 0) map[ii] = IOF_ID;
      else if (strcasecmp(fieldNames[ii], "pinfo")    == 0) map[ii] = IOF_PINFO;
      else if (strcasecmp(fieldNames[ii], "rx")       == 0) map[ii] = IOF_RX;
      else if (strcasecmp(fieldNames[ii], "ry")       == 0) map[ii] = IOF_RY;
      else if (strcasecmp(fieldNames[ii], "rz")       == 0) map[ii] = IOF_RZ;
      else if (strcasecmp(fieldNames[ii], "vx")       == 0) map[ii] = IOF_VX;
      else if (strcasecmp(fieldNames[ii], "vy")       == 0) map[ii] = IOF_VY;
      else if (strcasecmp(fieldNames[ii], "vz")       == 0) map[ii] = IOF_VZ;
   }

   return map;


}

/**
 *  WARNING!! The offset array must have size nFields+1 where nFields is
 *  the number of fields in the field_types string.
 */
unsigned parseFieldSizes(unsigned* offset, const char* field_types)
{
   unsigned count = 0;
   char* string = strdup(field_types);
   char* word = strtok(string, " ");
   assert(word);
   offset[0] = 0;
   while (word != NULL)
   {
      assert(strlen(word) == 2);
      offset[count+1] = offset[count] + atoi(word+1);
      word = strtok(NULL, " ");
      ++count;
   }

   ddcFree(string);
   return count;
}

int makeOffsetsAndTypes(unsigned nInput, char** inputNames,
      char** inputTypes, const char* requestedNames,
      unsigned* offset, char** type)
{
   int nRequest = 0;
   int nWords = 0;
   unsigned* inputOffsets = makeOffsets(inputTypes, nInput);
   char* string = strdup(requestedNames);
   char* word = strtok(string, " ");
   if (!word)
      return -1;
   while (word != NULL)
   {
      ++nWords;
      for (unsigned ii=0; ii<nInput; ++ii)
      {
         if (strcmp(word, inputNames[ii]) == 0)
         {
            offset[nRequest] = inputOffsets[ii];
            strcpy(type[nRequest], inputTypes[ii]);
            ++nRequest;
            continue;
         }
      }
      if (nWords != nRequest)
         return -1;
      word = strtok(NULL, " ");
   }

   ddcFree(string);
   ddcFree(inputOffsets);
   return nRequest;
}



   void
endianSwap(void* data, int size)
{ 
   int i, j; 
   char* b, save; 
   b = (char*) (data) ; 
   j = size; 
   for (i=0; i<size/2; i++) 
   {
      --j; 
      save = b[i] ;  
      b[i] = b[j]; 
      b[j] = save; 
   }
} 

