#ifndef TT06_RRG_HH
#define TT06_RRG_HH

#include <string>
#include <vector>
#include <map>
#include "CheckpointVarInfo.hh"

class TT06_RRG
{
 public:
   // There is a map of strings to VarHandle in getHandleMap that must
   // be maintained in parallel with this enum.  The value undefinedName
   // must be -1 and nVars must be last in this enum.
   enum VarHandle{undefinedName = -1,
                  // These are the cell-wise parameters:
                  s_switch, g_Ks, g_to, P_NaK, g_NaL,
                  // These are the state variables:
                  Vm, K_i, Na_i, Ca_i, Xr1_gate, Xr2_gate, Xs_gate, m_gate,
                  h_gate, j_gate, Ca_ss, d_gate, f_gate, f2_gate, fCass_gate,
                  s_gate, r_gate, Ca_SR, R_prime, jL_gate,
                  // end marker
                  nVars};
      

   TT06_RRG(int cellType);
   double calc(double dt, double Vm, double iStim);
   double defaultVoltage();
   static void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                 std::vector<std::string>& fieldUnits);

   static int getVarHandle(const std::string& varName);
   void setValue(int varHandle, double value);

   double getValue(int handle) const;
   void getValue(const std::vector<int>& handle,
                 std::vector<double>& value) const;
   static const std::string getUnit(const std::string& varName);


   
 private:

   static HandleMap& getHandleMap();
   void initConsts(int cellType);
   void initStates(int cellType);

   double computeRates(double dt, double iStim);

   static double constants_[53];
   int s_switch_;
   double defaultVoltage_;
   double g_Ks_;  // formerly CONSTANTS[15]
   double g_to_;  // formerly CONSTANTS[20]
   double P_NaK_; // formerly CONSTANTS[21]
   double g_NaL_;
   double states_[20];
};

#endif
