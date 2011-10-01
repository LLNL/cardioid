#ifndef EPHYSMODEL_H
#define EPHYSMODEL_H
//#include "ddc.h"
#include <vector>
#include <string>
using namespace std;

class EPhysModel {

  protected:

  int _nloc;
  
  public:

  //  DDC* ddc;
  vector<double> potLoc;  // local potential grid 
  int np(void) { return _nloc; };

  EPhysModel() {}

  // virtual destructor needed to ensure proper deallocation
  virtual ~EPhysModel() {}
  virtual void calcReact(void) = 0;
  virtual void read(string filebase) = 0;
  virtual void write(string filebase) = 0;
};
#endif
