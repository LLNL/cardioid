#ifndef ANATOMY_CELL_HH
#define ANATOMY_CELL_HH

#include "Long64.hh"
#include "SymmetricTensor.hh"

class AnatomyCell
{
  public:
   Long64 gid_;
   int cellType_;
   SymmetricTensor sigma_;
   int dest_;
   int sortind_;
   static bool destLessThan(const AnatomyCell &a, const AnatomyCell &b)
        { return a.dest_ < b.dest_; }
   static bool indLessThan(const AnatomyCell &a, const AnatomyCell &b)
        { return a.sortind_ < b.sortind_; }
};


class AnatomyCellDestSort
{
 public:
   bool operator()(const AnatomyCell& a, const AnatomyCell& b)
   {
      return a.dest_ < b.dest_;
   }
};

class AnatomyCellGidSort
{
 public:
   bool operator()(const AnatomyCell& a, const AnatomyCell& b)
   {
      return a.gid_ < b.gid_;
   }
};

#endif
