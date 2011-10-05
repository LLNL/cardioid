#ifndef ANATOMY_CELL_HH
#define ANATOMY_CELL_HH

#include "Long64.hh"

class AnatomyCell
{
  public:
   Long64 _gid;
   int _cellType;
   int _theta;
   int _phi;
   int _dest;
};

#endif
