#include "object_cc.hh"

#include "ddcMalloc.h"

using std::string;
using std::vector;


OBJECT* objectFind(const std::string& name, const std::string& objClass)
{
   return object_find(name.c_str(), objClass.c_str());
}

void objectGet(OBJECT* obj, const string& name, string& value, const string& defVal)
{
   char* tmp;
   object_get(obj, name.c_str(), &tmp, STRING, 1, defVal.c_str());
   value = tmp;
   ddcFree(tmp);
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

void objectGet(OBJECT* obj, const string& name, vector<string>& value)
{
   value.clear();
   char** tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, STRING, IGNORE_IF_NOT_FOUND);
   for (unsigned ii=0; ii<n; ++ii)
   {
      value.push_back(tmp[ii]);
      ddcFree(tmp[ii]);
   }
   ddcFree(tmp);
}

void objectGet(OBJECT* obj, const string& name, vector<unsigned>& value)
{
   value.clear();
   int* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, INT, IGNORE_IF_NOT_FOUND);
   for (unsigned ii=0; ii<n; ++ii)
      value.push_back(tmp[ii]);
   ddcFree(tmp);
}

void objectGet(OBJECT* obj, const string& name, vector<double>& value)
{
   value.clear();
   double* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, DOUBLE, IGNORE_IF_NOT_FOUND);
   for (unsigned ii=0; ii<n; ++ii)
      value.push_back(tmp[ii]);
   ddcFree(tmp);
}
