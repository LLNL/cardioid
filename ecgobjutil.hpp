#include "mfem.hpp"
#include "object.h"
#include "ddcMalloc.h"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cassert>
#include <memory>
#include <set>
#include <ctime>

#pragma once

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

template <typename T> struct PioTypeFromCType;

template <> struct PioTypeFromCType<double> {
   enum { result = DOUBLE };
};

template <> struct PioTypeFromCType<int> {
   enum { result = INT };
};

template <> struct PioTypeFromCType<std::string> {
   enum { result = STRING };
};

template <typename T>
void objectGetv(OBJECT* obj, const std::string& name, std::vector<T>& value) {
   value.clear();
   T* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, PioTypeFromCType<T>::result, IGNORE_IF_NOT_FOUND);
   value.reserve(n);
   for (unsigned ii=0; ii<n; ++ii)
   {
      value.push_back(tmp[ii]);
   }
   ddcFree(tmp);
}

template <>
void objectGetv<std::string>(OBJECT* obj, const std::string& name, std::vector<std::string>& value)
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

std::unordered_map<int,int> readMap(std::istream& stream) {
   std::unordered_map<int,int> retval;
   int i = 0;
   while(1) {
      int value;
      stream >> value;
      if (!stream) {break;}
      retval[i++] = value;
   }
   return retval;
}

std::unordered_map<int,int> ecg_readMap(OBJECT* obj, const std::string dataname) {
   std::string filename;
   objectGet(obj,dataname,filename,"");

   std::ifstream infile(filename);
   std::unordered_map<int,int> retval = readMap(infile);
#ifdef DEBUG
   std::cout << "Read " << retval.size() << " elements from "
	     << filename << std::endl;
#endif
   return retval;
}

void ecg_readGF(OBJECT* obj, const std::string dataname,
		mfem::ParMesh* pmesh, std::shared_ptr<mfem::ParGridFunction>& sp_gf) {
   std::string filename;
   objectGet(obj,dataname,filename,"");

   std::ifstream infile(filename);
   sp_gf = std::make_shared<mfem::ParGridFunction>(pmesh, infile);
}

mfem::Mesh* ecg_readMeshptr(OBJECT* obj, const std::string dataname) {
   std::string filename;
   objectGet(obj, dataname, filename, "");

   // generate_edges=1, refine=1 (fix_orientation=true by default)
   return new mfem::Mesh(filename.c_str(), 1, 1);
}

// If rank==0, then object_bcast()

void ecg_process_args(int argc, char* argv[]) {
   std::vector<std::string> objectFilenames;
   if (argc == 1)
      objectFilenames.push_back("ecg.data");

   for (int iargCursor=1; iargCursor<argc; iargCursor++)
      objectFilenames.push_back(argv[iargCursor]);

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank == 0) {
      for (int ii=0; ii<objectFilenames.size(); ii++)
	 object_compilefile(objectFilenames[ii].c_str());
   }
   object_Bcast(0,MPI_COMM_WORLD);
}

// Read sensors
std::unordered_map<int,int> ecg_readInverseMap(OBJECT* obj, const std::string dataname) {
   std::string filename;
   objectGet(obj,dataname,filename,"");
   std::ifstream sensorStream(filename);

   std::unordered_map<int,int> gfFromGid;
   //read in the sensor file
   {
      int token = 0;
      do {
	 //read in from the sensor file
	 int sensorGid;
	 sensorStream >> sensorGid;
	 if (!sensorStream) break;
	 gfFromGid[sensorGid] = token++;
      } while(1);
   }

   return gfFromGid;
}

// This is really awful.  Too specific to be general, but presented as if it could be.
std::map<int,std::string> ecg_readAssignments(OBJECT* obj, const std::string dataname) {
   std::string filename;
   objectGet(obj,dataname,filename,"");
   std::ifstream electrodeStream(filename);

   std::map<int,std::string> electrodeFromGid;
   do {
      int gid;
      std::string label;
      electrodeStream >> label; // First column becomes the electrode name
      electrodeStream >> gid;   // Second column gets discarded
      electrodeStream >> gid;   // Third column becomes the GID
      if (!electrodeStream) break;
      electrodeFromGid[gid] = label;
   } while(1);

   return electrodeFromGid;
}

// Basic inline progress reporting

time_t delta_t() {
   static time_t last_time;
   time_t curr_time = time(nullptr);
   time_t retval = curr_time - last_time;
   last_time = curr_time;
   return retval;
}

void StartTimer(std::string msg) {
#ifdef DEBUG
   std::cout << std::endl << msg << "... ";
   time_t discard = delta_t();  // Update timestamp but don't report it
#endif
}

void EndTimer() {
#ifdef DEBUG
   std::cout << "Done in " << delta_t() << "s." << std::endl;
#endif
}

// electrodes.txt is a bit different from the others.
// <label> = <gid>
// We could overload the template but let's not for now since we don't even know 

// Instead of trying to use BucketOfBits, just look for where it uses calls which occur above,
// and modify the above to use similar patterns.

// Flatten part(s) of readPioFile

// readAscii

// pfgets distributes all data among ranks, but only some are
// getting it from the file(s)

// While you have no way of knowing which data you will get from a pfgets,
// You can know nobody else is getting it.

// Reassemble as required

// Verify whether/how mfem sends data when assigned into a ParGridFunction
// (look for any PGF assemble mechanisms for clues)

// string::find_first_...
