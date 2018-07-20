#include "mfem.hpp"
#include "object.h"
#include "ddcMalloc.h"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cassert>
#include <memory>
#include <set>

#ifndef __ECGOBJUTIL_HPP__
#define __ECGOBJUTIL_HPP__

/** The object_get code doesn't handle empty strings for a default value
 *  very well.  The problem is that object_parse attempts to tokenize
 *  the default value and strtok_r returns a NULL pointer.  This results
 *  in a situation where the pointer passed to object_get (i.e., tmp)
 *  isn't set by object_get.  Fortunately, this can be detected by the
 *  fact that object_get will return 0 instead of 1.  We can catch this
 *  case and do the right thing with the value we return to the caller.
 */
void objectGet(OBJECT* obj, const std::string& name,
	       std::string& value, const std::string& defVal) {
  char* tmp;
  int nFound = object_get(obj, name.c_str(), &tmp, STRING, 1, defVal.c_str());
  if (nFound != 0) {
    value = tmp;
    ddcFree(tmp);
  }
  else
    value = "";
}

void objectGet(OBJECT* obj, const std::string& name, std::vector<int>& value) {
  value.clear();
  int* tmp;
  unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, INT, IGNORE_IF_NOT_FOUND);
  value.reserve(n);
  for (unsigned ii=0; ii<n; ++ii)
    value.push_back(tmp[ii]);
  ddcFree(tmp);
}

void objectGet(OBJECT* obj, const std::string& name, std::vector<double>& value)
{
  value.clear();
  double* tmp;
  unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, DOUBLE, IGNORE_IF_NOT_FOUND);
  value.reserve(n);
  for (unsigned ii=0; ii<n; ++ii)
    value.push_back(tmp[ii]);
  ddcFree(tmp);
}

std::set<int> readSet(std::istream& stream) {
  std::set<int> retval;
  while(1) {
    int value;
    stream >> value;
    if (!stream) {break;}
    retval.insert(value);
  }
  return retval;
}

std::set<int> ecg_readSet(OBJECT* obj, const std::string dataname) {
  std::string filename;
  objectGet(obj,dataname,filename,"");

  std::ifstream infile(filename);
  std::set<int> retval = readSet(infile);
#ifdef DEBUG
  std::cout << "Read " << retval.size() << " elements from "
	    << filename << std::endl;
#endif
  return retval;
}

void ecg_readGF(OBJECT* obj, const std::string dataname,
		mfem::Mesh* mesh, std::shared_ptr<mfem::GridFunction>& sp_gf) {
  std::string filename;
  objectGet(obj,dataname,filename,"");

  std::ifstream infile(filename);
  sp_gf = std::make_shared<mfem::GridFunction>(mesh, infile);
}

mfem::Mesh* ecg_readMeshptr(OBJECT* obj, const std::string dataname) {
  std::string filename;
  objectGet(obj, dataname, filename, "");

  // generate_edges=1, refine=1 (fix_orientation=true by default)
  return new mfem::Mesh(filename.c_str(), 1, 1);
}

void ecg_process_args(int argc, char* argv[]) {
  std::vector<std::string> objectFilenames;
  if (argc == 1)
    objectFilenames.push_back("ecg.data");

  for (int iargCursor=1; iargCursor<argc; iargCursor++)
    objectFilenames.push_back(argv[iargCursor]);

  for (int ii=0; ii<objectFilenames.size(); ii++)
    object_compilefile(objectFilenames[ii].c_str());
}

#endif
