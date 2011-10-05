#ifndef MECHMODEL_H
#define MECHMODEL_H

#include <vector>

class MechModel {

  protected:

  int _nloc;
  
  public:

  std::vector<double> dataLoc;  // local mechanical data
  int np(void) { return _nloc; };

  MechModel() {}

  // virtual destructor needed to ensure proper deallocation
  virtual ~MechModel() {}

  virtual void calcReact(void) = 0;
};
#endif
