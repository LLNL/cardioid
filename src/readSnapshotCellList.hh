#ifndef READ_SNAPSHOT_CELL_LIST_HH
#define READ_SNAPSHOT_CELL_LIST_HH

#include <string>
#include <set>
#include "Long64.hh"

class Simulate;

bool readSnapshotCellList(std::string filename, Simulate& sim);

#endif
