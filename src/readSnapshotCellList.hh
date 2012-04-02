#ifndef READ_SNAPSHOT_CELL_LIST_HH
#define READ_SNAPSHOT_CELL_LIST_HH

#include <string>
#include <set>
#include "Long64.hh"
#include "object_cc.hh"

class Simulate;

bool readSnapshotCellList(std::string filename, Simulate& sim, OBJECT* obj);

#endif
