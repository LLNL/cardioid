#ifndef TT06_RRG_HH
#define TT06_RRG_HH


class TT06_RRG
{
 public:
   TT06_RRG(int cellType);
   double calc(double dt, double Vm, double iStim);
   double defaultVoltage();
   
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
