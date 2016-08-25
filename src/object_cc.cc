#include "object_cc.hh"

#include "ddcMalloc.h"
using std::string;
using std::vector;


OBJECT* objectFind(const std::string& name, const std::string& objClass)
{
   return object_find(name.c_str(), objClass.c_str());
}

/** The object_get code doesn't handle empty strings for a default value
 *  very well.  The problem is that object_parse attempts to tokenize
 *  the default value and strtok_r returns a NULL pointer.  This results
 *  in a situation where the pointer passed to object_get (i.e., tmp)
 *  isn't set by object_get.  Fortunately, this can be detected by the
 *  fact that object_get will return 0 instead of 1.  We can catch this
 *  case and do the right thing with the value we return to the caller.
 */
void objectGet(OBJECT* obj, const string& name, string& value, const string& defVal)
{
   char* tmp;
   int nFound = object_get(obj, name.c_str(), &tmp, STRING, 1, defVal.c_str());
   if (nFound != 0)
   {
      value = tmp;
      ddcFree(tmp);
   }
   else
      value = "";
}

void objectGet(OBJECT* obj, const string& name, double& value, const string& defVal)
{
   double tmp;
   object_get(obj, name.c_str(), &tmp, DOUBLE, 1, defVal.c_str());
   value = tmp;
}

void objectGet(OBJECT* obj, const string& name, bool& value, const string& defVal)
{
   int tmp;
   object_get(obj, name.c_str(), &tmp, INT, 1, defVal.c_str());
   value = (tmp!=0);
}

void objectGet(OBJECT* obj, const string& name, int& value, const string& defVal)
{
   int tmp;
   object_get(obj, name.c_str(), &tmp, INT, 1, defVal.c_str());
   value = tmp;
}

void objectGet(OBJECT* obj, const string& name, unsigned& value, const string& defVal)
{
   int tmp;
   object_get(obj, name.c_str(), &tmp, INT, 1, defVal.c_str());
   value = tmp;
}

void objectGet(OBJECT* obj, const string& name, uint64_t& value, const string& defVal)
{
   uint64_t tmp;
   object_get(obj, name.c_str(), &tmp, U64, 1, defVal.c_str());
   value = tmp;
}
void objectGet(OBJECT* obj, const string& name, int64_t& value, const string& defVal)
{
   long long int  tmp;
   object_get(obj, name.c_str(), &tmp, I64, 1, defVal.c_str());
   value = tmp;
}

void objectGet(OBJECT* obj, const std::string& name, double& value, const std::string& defVal,
               const std::string& unitConvertTo)
{
   double tmp;
   object_get(obj, name.c_str(), &tmp, WITH_UNITS, 1, defVal.c_str(), NULL, unitConvertTo.c_str());
   value = tmp;
}


void objectGet(OBJECT* obj, const string& name, vector<string>& value)
{
   value.clear();
   char** tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, STRING, IGNORE_IF_NOT_FOUND);
   value.reserve(n);
   for (unsigned ii=0; ii<n; ++ii)
   {
      value.push_back(tmp[ii]);
      ddcFree(tmp[ii]);
   }
   ddcFree(tmp);
}

void objectGet(OBJECT* obj, const string& name, vector<int>& value)
{
   value.clear();
   int* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, INT, IGNORE_IF_NOT_FOUND);
   value.reserve(n);
   for (unsigned ii=0; ii<n; ++ii)
      value.push_back(tmp[ii]);
   ddcFree(tmp);
}
void objectGet(OBJECT* obj, const string& name, vector<unsigned>& value)
{
   value.clear();
   int* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, INT, IGNORE_IF_NOT_FOUND);
   value.reserve(n);
   for (unsigned ii=0; ii<n; ++ii)
      value.push_back(tmp[ii]);
   ddcFree(tmp);
}

void objectGet(OBJECT* obj, const string& name, vector<uint64_t>& value)
{
   value.clear();
   uint64_t* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, U64, IGNORE_IF_NOT_FOUND);
   value.reserve(n);
   for (unsigned ii=0; ii<n; ++ii)
      value.push_back(tmp[ii]);
   ddcFree(tmp);
}

void objectGet(OBJECT* obj, const string& name, vector<double>& value)
{
   value.clear();
   double* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, DOUBLE, IGNORE_IF_NOT_FOUND);
   value.reserve(n);
   for (unsigned ii=0; ii<n; ++ii)
      value.push_back(tmp[ii]);
   ddcFree(tmp);
}
