#ifndef TT04MODEL_H
#define TT04MODEL_H

#include <string>
#include "EPhysModel.hh"

class TT04Model : public EPhysModel
{

 public:
  
   TT04Model();
   void calcReact();
   void read(const std::string& filebase);
   void write(const std::string& filebase);
  
};
#endif

