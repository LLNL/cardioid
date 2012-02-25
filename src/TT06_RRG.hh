#ifndef TT06_RRG_HH
#define TT06_RRG_HH

#include <string>
#include <vector>
#include <map>

class TT06_RRG
{
 public:
   // There is a map of strings to VarHandle in getVarHandle that must
   // be maintained in parallel with this enum.  The value undefinedName
   // must be first and nVars must be last in this enum.
   enum VarHandle{undefinedName,
                  // These are the cell-wise parameters:
                  s_switch, g_Ks, g_to, P_NaK, g_NaL,
                  // These are the state variables:
                  Vm, K_i, Na_i, Ca_i, Xr1, Xr2, Xs, m, h, j,
                  Ca_ss, d, f, f2, fCass, s, r, Ca_SR, R_prime, NaL_i,
                  // end marker
                  nVars};
      
   class VarInfo
   {
    public:
      VarInfo()
      : handle_(undefinedName), checkpoint_(false), unit_("1")
      {};
      VarInfo(VarHandle handle, bool checkpoint, std::string unit)
      : handle_(handle), checkpoint_(checkpoint), unit_(unit)
      {};
      
      VarHandle   handle_;
      bool        checkpoint_;
      std::string unit_; // output unit
   };

   typedef std::map<std::string, VarInfo> HandleMap; 

   TT06_RRG(int cellType);
   double calc(double dt, double Vm, double iStim);
   double defaultVoltage();
   static VarHandle getVarHandle(const std::string& varName);
   static std::vector<int> getVarHandle(const std::vector<std::string>& varName);
   void setVariable(VarHandle varHandle, double value);

   double getValue(VarHandle handle) const;
   void getValue(const std::vector<int>& handle,
                 std::vector<double>& value) const;

   static const std::string& getUnit(const std::string& varName);

   static void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                 std::vector<std::string>& fieldUnits);


   
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
