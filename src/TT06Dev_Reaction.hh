#ifndef TT06DEV_REACTION_HH
#define TT06DEV_REACTION_HH
#include "Reaction.hh"
class Anatomy;
class TT06DevState;

class TT06Dev_Reaction : public Reaction
{
 public:
   
   TT06Dev_Reaction(const Anatomy& anatomy, double tolerance,int mod);
   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT06Dev_Reaction(const TT06Dev_Reaction&);
   TT06Dev_Reaction& operator=(const TT06Dev_Reaction&);
   ~TT06Dev_Reaction();

   void calc(double dt, const std::vector<double>& Vm, const std::vector<double>& iStim, std::vector<double>& dVm);

 private:

   unsigned nCells_;
   double dtForFit_; 
   double tolerance_; 
   int mod_; 
   
   std::vector<int>              ttType_; // maps cellType to ttType
   std::vector<TT06DevState> s_;
};

#endif
