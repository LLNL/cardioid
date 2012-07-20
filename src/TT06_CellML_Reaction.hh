#ifndef TT06_CELLML_REACTION_HH
#define TT06_CELLML_REACTION_HH

#include "Reaction.hh"
#include <map>
#include "CheckpointVarInfo.hh"

class Anatomy;
class TT06_CellML;
class TT06_CellMLState;

class TT06_CellML_Reaction : public Reaction
{
 public:
   enum IntegratorType {forwardEuler, rushLarsen};
   
   // There is a map of strings to VarHandle in getHandleMap that must
   // be maintained in parallel with this enum.  The value undefinedName
   // must be -1 and nVars must be last in this enum.
   enum VarHandle{undefinedName = -1,
                  // These are the state variables:
                  Vm, K_i, Na_i, Ca_i, Xr1_gate, Xr2_gate, Xs_gate, m_gate,
                  h_gate, j_gate, Ca_ss, d_gate, f_gate, f2_gate, fCass_gate,
                  s_gate, r_gate, Ca_SR, R_prime,
                  // end marker
                  nVars};
   
   TT06_CellML_Reaction(const Anatomy& anatomy, IntegratorType integrator);
   ~TT06_CellML_Reaction();
   std::string methodName() const {return "TT06_CellML";}

   void calc(double dt,
             const VectorDouble32& Vm,
             const std::vector<double>& iStim,
             VectorDouble32& dVm);
   void initializeMembraneVoltage(VectorDouble32& Vm);

   /** Functions needed for checkpoint/restart */
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
   void forwardEulerIntegrator(double dt, const VectorDouble32& Vm,
      const std::vector<double>& iStim, VectorDouble32& dVm);
   void rushLarsenIntegrator(double dt, const VectorDouble32& Vm,
      const std::vector<double>& iStim, VectorDouble32& dVm);

   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT06_CellML_Reaction(const TT06_CellML_Reaction&);
   TT06_CellML_Reaction& operator=(const TT06_CellML_Reaction&);

   int nCells_;
   IntegratorType integrator_;
   
   std::vector<int>              ttType_; // maps cellType to ttType
   std::vector<TT06_CellML*>     cellModel_;
   std::vector<TT06_CellMLState> s_;
};

#endif
