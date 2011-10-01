#ifndef MECHMODEL_H
#define MECHMODEL_H
//#include "ddc.h"
#include <vector>
using namespace std;

class MechModel {

  protected:

  int _nloc;
  
  public:

  //  DDC* ddc;
  vector<double> dataLoc;  // local mechanical data
  int np(void) { return _nloc; };

  MechModel() {}

  // virtual destructor needed to ensure proper deallocation
  virtual ~MechModel() {}

  virtual void calcReact(void) = 0;
};
#endif
