#ifndef ANATOMY_CELL_HH
#define ANATOMY_CELL_HH

#include "Long64.hh"

class AnatomyCell
{
  public:
   Long64 gid_;
   int cellType_;
   int theta_;
   int phi_;
   int dest_;
};

#endif
