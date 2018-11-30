#pragma once

#include "mfem.hpp"
#include "object.h"
#include <iostream>
#include <unordered_map>
#include <memory>
#include <set>
#include "object_cc.hh"
#include "pioFixedRecordHelper.h"

std::set<int> readSet(std::istream& stream);
std::set<int> ecg_readSet(OBJECT* obj, const std::string dataname);
std::unordered_map<int,int> readMap(std::istream& stream);
std::unordered_map<int,int> ecg_readMap(OBJECT* obj, const std::string dataname);
void ecg_readGF(OBJECT* obj, const std::string dataname,
		mfem::Mesh* mesh, std::shared_ptr<mfem::GridFunction>& sp_gf);
mfem::Mesh* ecg_readMeshptr(OBJECT* obj, const std::string dataname);
void ecg_process_args(int argc, char* argv[]);
std::unordered_map<int,int> ecg_readInverseMap(OBJECT* obj, const std::string dataname);
// This is really awful.  Too specific to be general, but presented as if it could be.
std::map<int,std::string> ecg_readAssignments(OBJECT* obj, const std::string dataname);
// Basic inline progress reporting

double ecg_readParGF(OBJECT* obj, const std::string VmFilename, const int global_size, const std::unordered_map<int,int> &gfFromGid, std::vector<double>& final_values);
