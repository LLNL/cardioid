#ifndef TT06_RRG_HH
#define TT06_RRG_HH

#include <string>

class TT06_RRG
{
 public:
   // There is a map of strings to VarHandle in getVarHandle that must
   // be maintained in parallel with this enum.  The value undefinedName
   // must be first and nVars must be last in this enum.
   enum VarHandle{undefinedName, switchTauS, g_Ks, g_to, P_NaK, g_NaL, nVars};
      


   TT06_RRG(int cellType);
   double calc(double dt, double Vm, double iStim);
   double defaultVoltage();
   static VarHandle getVarHandle(const std::string& varName);
   void setVariable(VarHandle varHandle, double value);

 private:

   void initConsts(int cellType);
   void initStates(int cellType);

   double computeRates(double dt, double iStim);

   static double constants_[53];
   int switchTauS_;
   double defaultVoltage_;
   double g_Ks_;  // formerly CONSTANTS[15]
   double g_to_;  // formerly CONSTANTS[20]
   double P_NaK_; // formerly CONSTANTS[21]
   double g_NaL_;
   double states_[20];
};

#endif
