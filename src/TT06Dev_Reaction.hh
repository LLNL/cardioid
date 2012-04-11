#ifndef TT06DEV_REACTION_HH
#define TT06DEV_REACTION_HH
#include "Reaction.hh"
#include "TT06Func.hh"
using namespace std;
using namespace TT06Func ;

class ThreadTeam;
class Anatomy;

class TT06Dev_Reaction : public Reaction
{
 public:
   
   TT06Dev_Reaction(Anatomy& anatomy,  map<string,TT06Func::CellTypeParmsFull>cellTypeParms, vector<string>cellTypeNames, double tolerance, int mod, const ThreadTeam& group);
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


   unsigned nCellTypes_;
   vector<TT06Func::CellTypeParms> cellTypeParms_; 
   vector<double> initialVm_; 
   unsigned nCells_;
   vector<unsigned> nCellsOfType_; 

   int nonGateWorkPartition(int& offset);
   int idPEStart_; 
   int nPEs_; 
   int nC_[2]; 
   double dtForFit_; 
   double tolerance_; 
   int mod_; 
   PADE **fit_;
   double *mhu_[nGateVar];
   double *tauR_[nGateVar];
   const ThreadTeam& group_;
   
   std::vector<int>              ttType_; // maps cellType to ttType
   int nCellBuffer_; 
   double *stateBuffer_; 
   std::vector<double*>state_; 
   std::vector<int>cellTypeVector_; 
};

#endif
