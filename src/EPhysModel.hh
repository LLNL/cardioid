#ifndef EPHYSMODEL_H
#define EPHYSMODEL_H

#include <vector>
#include <string>


class EPhysModel {

  protected:

  int _nloc;
  
  public:

  std::vector<double> potLoc;  // local potential grid 
  int np(void) { return _nloc; };

  EPhysModel() {}

  // virtual destructor needed to ensure proper deallocation
  virtual ~EPhysModel() {}
  virtual void calcReact(void) = 0;
  virtual void read(const std::string& filebase) = 0;
  virtual void write(const std::string& filebase) = 0;
};
#endif
