#ifndef TT06DEV_REACTION_HH
#define TT06DEV_REACTION_HH
#include "Reaction.hh"
#include "TT06Func.hh"
#include "TT06NonGates.h"
using namespace std;
using namespace TT06Func ;

class ThreadTeam;
class Anatomy;

class TT06Dev_Reaction : public Reaction
{
 public:
   
   TT06Dev_Reaction(Anatomy& anatomy,  map<string,TT06Func::CellTypeParmsFull>cellTypeParms, vector<string>cellTypeNames, double tolerance, int mod, int  fastReaction, const ThreadTeam& group);
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
   vector<CellTypeParms> cellTypeParms_; 
   vector<double> initialVm_; 
   unsigned nCells_;
   vector<unsigned> nCellsOfType_; 

    void (*update_gate_)   (double dt,                                      int nCells, int *cellType, double *Vm, int offset, double **state, WORK& work);
    void (*update_nonGate_)(double dt, struct CellTypeParms *cellTypeParms, int nCells, int *cellType, double *Vm, int offset, double **state, double *dVdt);
   int nonGateWorkPartition(int& offset);
   double dtForFit_; 
   double tolerance_; 
   int mod_; 
   int fastReaction_; 
   PADE **fit_;
   double *mhu_[nGateVar];
   double *tauR_[nGateVar];
   const ThreadTeam& group_;
   vector<WORK> gateWork_; 
   LogParms logParms_[64]; 
   
   std::vector<int>              ttType_; // maps cellType to ttType
   int nCellBuffer_; 
   double *stateBuffer_; 
   std::vector<double*>state_; 
   std::vector<int>cellTypeVector_; 
};

#endif
