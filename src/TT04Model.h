#ifndef TT04MODEL_H
#define TT04MODEL_H

#include <vector>
#include "EPhysModel.h"
using namespace std;

class TT04Model : public EPhysModel
{

  public:
  
  TT04Model();
  void calcReact();
  void read(string filebase);
  void write(string filebase);
  
};
#endif

