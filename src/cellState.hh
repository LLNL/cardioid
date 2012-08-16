#ifndef CELL_STATE_HH
#define CELL_STATE_HH

#include "Long64.hh"
#include <string>

class cellState
{
  public:
    Long64 gid_;
    double value_;
    int cellType_;
    int dest_;
    Long64 sortind_;
    static bool destLessThan(const cellState &a, const cellState &b)
        { return a.dest_ < b.dest_; }
    static bool indLessThan(const cellState &a, const cellState &b)
        { return a.sortind_ < b.sortind_; }
};

#endif
