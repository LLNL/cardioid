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
   TT06Func::CellTypeParmsFull cellTypeParm;
   std::vector<std::string> currentNames;
   std::string fitFile; 
   double tolerance;
   std::vector<int> gateThreadMap; 
   int jhTauSmooth;
   PADE *fit; 
   int fastReaction;
};

typedef void  (*UPDATEGATE)(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a) ; 

class TT06Dev_Reaction : public Reaction
{
 public:
   
   TT06Dev_Reaction(const double dt, const int numPoints, TT06Dev_ReactionParms& parms, const ThreadTeam& group);

   std::string methodName() const {return "TT06_Dev";}
   // copy constructor and assignment operator intentionally left unimplemented.
   TT06Dev_Reaction(const TT06Dev_Reaction&);
   TT06Dev_Reaction& operator=(const TT06Dev_Reaction&);
   ~TT06Dev_Reaction();

   void calc(double dt,
             ro_mgarray_ptr<double> Vm,
             ro_mgarray_ptr<double> iStim,
             wo_mgarray_ptr<double> dVm);
   void updateNonGate(double dt, ro_mgarray_ptr<double> Vm, wo_mgarray_ptr<double> dVR);
   void updateGate   (double dt, ro_mgarray_ptr<double> Vm);
   void initializeMembraneVoltage(wo_mgarray_ptr<double> Vm);

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

   CellTypeParms cellTypeParm_;
   double XXXinitialVm_;
   int nCells_;

   std::vector<double> XXXstateInitial_;

   void (*update_gate_)   (double dt,                                      int nCells, int s_switch, double *Vm, int offset, double **state, PADE* xfit, TT06Func::WORK& work);
   void (*update_nonGate_)(void *fit, double dt, struct CellTypeParms *cellTypeParms, int nCells, double *Vm, int offset, double **state, double *dVdt);
   int nonGateWorkPartition_(int& offset);
   void mkCellTypeParms_(TT06Dev_ReactionParms& parms);
   void mkState_(TT06Dev_ReactionParms& parms);
   void mkFitParms_(TT06Dev_ReactionParms& parms);
   void mkWorkBundles_(TT06Dev_ReactionParms& parms); 
   
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
   
   int nCellBuffer_; 
   double *stateBuffer_; 
   std::vector<double*>state_; 
   double* gateX_[13];
   double* mhuX_[13];
   double* tauRX_[13];
   UPDATEGATE gateEqX_[13];
   int gateThreadMap_[12];

};

#endif
