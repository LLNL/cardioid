#ifndef TT06DEV_REACTION_HH
#define TT06DEV_REACTION_HH
#include "Reaction.hh"
#include "TT06Func.hh"

class ThreadTeam;
class Anatomy;

class TT06Dev_Reaction : public Reaction
{
 public:
   
   TT06Dev_Reaction(const Anatomy& anatomy, double tolerance, int mod, const ThreadTeam& group);
   std::string methodName() const {return "TT06_Dev";}
   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT06Dev_Reaction(const TT06Dev_Reaction&);
   TT06Dev_Reaction& operator=(const TT06Dev_Reaction&);
   ~TT06Dev_Reaction();

   void updateNonGate(double dt, const std::vector<double>&Vm, std::vector<double>&dVR);
   void updateGate   (double dt, const std::vector<double>&Vm) ;
   void calc(double dt, const std::vector<double>& Vm, const std::vector<double>& iStim, std::vector<double>& dVm);
   void initializeMembraneVoltage(std::vector<double>& Vm);
   void writeStateDev(int loop); 

 private:


   int nonGateWorkPartition(int& offset);
   int idPEStart_; 
   int nPEs_; 
   unsigned nCells_;
   double dtForFit_; 
   double tolerance_; 
   int mod_; 
   const ThreadTeam& group_;
   
   std::vector<int>              ttType_; // maps cellType to ttType
   std::vector<TT06Func::TT06DevState> s_;
   double **gates_; 
};

#endif
