#ifndef TT06DEV_REACTION_HH
#define TT06DEV_REACTION_HH
#include <map>
#include <vector>
#include <string>
#include "Reaction.hh"
#include "TT06Func.hh"
#include "TT06NonGates.h"
#include "CheckpointVarInfo.hh"
#include "pade.hh"

typedef double (*OVF)(double x, void *parms) ;

OVF fitFuncMap(std::string name) ;

class ThreadTeam;
class Anatomy;

struct TT06Dev_ReactionParms
{
   std::map<std::string, TT06Func::CellTypeParmsFull> cellTypeParms;
   std::vector<std::string> cellTypeNames;
   std::string fitFile; 
   double tolerance;
   std::vector<int> gateThreadMap; 
   int jhTauSmooth;
   PADE *fit; 
   int fastReaction;
};


class TT06Dev_Reaction : public Reaction
{
 public:
   
   TT06Dev_Reaction(double dt, Anatomy& anatomy, TT06Dev_ReactionParms& parms,
      const ThreadTeam& group);
   std::string methodName() const {return "TT06_Dev";}
   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT06Dev_Reaction(const TT06Dev_Reaction&);
   TT06Dev_Reaction& operator=(const TT06Dev_Reaction&);
   ~TT06Dev_Reaction();

   void updateNonGate(double dt, const VectorDouble32&Vm, VectorDouble32&dVR);
   void updateGate   (double dt, const VectorDouble32&Vm) ;
   void calc(double dt, const VectorDouble32& Vm, const std::vector<double>& iStim, VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);
   void writeStateDev(int loop);

   /** Support for checkpoint/restart */
   void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                  std::vector<std::string>& fieldUnits) const;
   int getVarHandle(const std::string& varName) const;
   void setValue(int iCell, int varHandle, double value);
   double getValue(int iCell, int varHandle) const;
   void getValue(int iCell,
                 const std::vector<int>& handle,
                 std::vector<double>& value) const;
   const std::string getUnit(const std::string& varName) const;

   

 private:

   HandleMap& getHandleMap() const;

   unsigned nCellTypes_;
   std::vector<CellTypeParms> cellTypeParms_; 
   std::vector<double> initialVm_; 
   unsigned nCells_;
   std::vector<unsigned> nCellsOfType_; 

   void (*update_gate_)   (double dt,                                      int nCells, int *cellType, double *Vm, int offset, double **state, PADE* xfit, TT06Func::WORK& work);
   void (*update_nonGate_)(void *fit, double dt, struct CellTypeParms *cellTypeParms, int nCells, int *cellType, double *Vm, int offset, double **state, double *dVdt);
   int nonGateWorkPartition(int& offset);
   
   double tolerance_; 
   int mod_; 
   int fastReaction_; 
   PADE *fit_;
   double *mhu_[nGateVar];
   double *tauR_[nGateVar];
   const ThreadTeam& group_;
   std::vector<TT06Func::WORK> gateWork_; 
   typedef struct { int offsetCell,nCell; } ngwork;
   std::vector<ngwork> nonGateWork_;
   LogParms logParms_[64]; 
   
   std::vector<int>              ttType_; // maps cellType to ttType
   int nCellBuffer_; 
   double *stateBuffer_; 
   std::vector<double*>state_; 
   std::vector<int>cellTypeVector_; 
};

#endif
