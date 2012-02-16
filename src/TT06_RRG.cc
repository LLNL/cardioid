#include "TT06_RRG.hh"
#include <cmath>
#include <cassert>
#include <map>
//ddt #include <iostream>

using namespace std;

double TT06_RRG::constants_[53];


TT06_RRG::TT06_RRG(int cellType)
{
   initConsts(cellType);
   initStates(cellType);
   defaultVoltage_ = states_[0];
}

double TT06_RRG::calc(double dt, double Vm, double iStim)
{
   states_[0] = Vm;
   double dVmdt = computeRates(dt, iStim);
   return dVmdt;
}

double TT06_RRG::defaultVoltage()
{
   return defaultVoltage_;
}

/** This function maps the string representation of a variable name to
 * the handle representation.  Returns the value undefinedName for
 * unrecognized varName. */
TT06_RRG::VarHandle TT06_RRG::getVarHandle(const string& varName)
{
   typedef map<string, VarHandle> HandleMap; 
   static HandleMap handleMap;
   if (handleMap.size() == 0)
   {
      handleMap["s_switch"] = s_switch;
      handleMap["g_Ks"]     = g_Ks;
      handleMap["g_to"]     = g_to;
      handleMap["P_NaK"]    = P_NaK;
      handleMap["g_NaL"]    = g_NaL;
      handleMap["Vm"]       = Vm;
      handleMap["K_i"]      = K_i;
      handleMap["Na_i"]     = Na_i;
      handleMap["Ca_i"]     = Ca_i;
      handleMap["Xr1"]      = Xr1;
      handleMap["Xr2"]      = Xr2;
      handleMap["Xs"]       = Xs;
      handleMap["m"]        = m;
      handleMap["h"]        = h;
      handleMap["j"]        = j;
      handleMap["Ca_ss"]    = Ca_ss;
      handleMap["d"]        = d;
      handleMap["f"]        = f;
      handleMap["f2"]       = f2;
      handleMap["fCass"]    = fCass;
      handleMap["s"]        = s;
      handleMap["r"]        = r;
      handleMap["Ca_SR"]    = Ca_SR;
      handleMap["R_prime"]  = R_prime;
      handleMap["NaL_i"]    = NaL_i;
      assert(handleMap.size() == nVars-1);
   }
   return handleMap[varName];
}


void TT06_RRG::setVariable(VarHandle varHandle, double value)
{
//ddt   cout << "Setting var "<< varHandle << "=" << value<<endl;
   
   switch (varHandle)
   {
     case undefinedName:
      assert(false);
      break;
     case s_switch:   s_switch_ = int(value);  break;
     case g_Ks:       g_Ks_ = value;           break;
     case g_to:       g_to_ = value;           break;
     case P_NaK:      P_NaK_ = value;          break;
     case g_NaL:      g_NaL_ = value;          break;
     case Vm:         states_[0] = value;      break;
     case K_i:        states_[1] = value;      break;
     case Na_i:       states_[2] = value;      break;
     case Ca_i:       states_[3] = value;      break;
     case Xr1:        states_[4] = value;      break;
     case Xr2:        states_[5] = value;      break;
     case Xs:         states_[6] = value;      break;
     case m:          states_[7] = value;      break;
     case h:          states_[8] = value;      break;
     case j:          states_[9] = value;      break;
     case Ca_ss:      states_[10] = value;     break;
     case d:          states_[11] = value;     break;
     case f:          states_[12] = value;     break;
     case f2:         states_[13] = value;     break;
     case fCass:      states_[14] = value;     break;
     case s:          states_[15] = value;     break;
     case r:          states_[16] = value;     break;
     case Ca_SR:      states_[17] = value;     break;
     case R_prime:    states_[18] = value;     break;
     case NaL_i:      states_[19] = value;     break;
     case nVars:
      assert(false);
      break;
   }
}



/** Everything below this point started life as the TT06 mid model
 *  generated code from the CellML web site (retrieved 8-Nov-2011).
 *  This code has been altered as follows:
 *
 *  1.  initConsts was split into two functions, initConsts and
 *  initStates.  Both are member functions, depend on cell type and
 *  initialize the appropriate class members directly.  initConsts 
 *  initializes the constants that vary from cell to cell as well as the
 *  static constants_ array.
 *
 *  2.  computeRates was converted to a member function.  The argument
 *  list was simplified.  The rates and algebraics arrays are created in
 *  the scope of the function.  Constants 15, 20, and 21, are not taken
 *  from the constants_ arrary, but rather from cellwise values g_Ks_,
 *  g_to_, and P_NaK_, respectively.  The calculation of ALGEBRAIC[12]
 *  (stimulus current) is replaced by the value of iStim that is passed
 *  in.  This obsoletes constants 5-8
 *
 *  3.  In computeRates we have modified the calculation of RATES[0] so
 *  that it no longer includes ALGEBRAIC[12].  This is for consistency
 *  with our conventions for dVm/dt (i.e., it is for the reaction part
 *  only).
 *
 *  4.  Integration of the internal state variables is included in
 *  computeRates.
 *
 *  5.  Additional state variable added for a slow Na gate.
 */ 


void TT06_RRG::initConsts(int cellType)
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      constants_[0] = 8314.472;
      constants_[1] = 310;
      constants_[2] = 96485.3415;
      constants_[3] = 0.185;
      constants_[4] = 0.016404;
//unused      constants_[5] = 10;
//unused      constants_[6] = 1000;
//unused      constants_[7] = 1;
//unused      constants_[8] = 52;
      constants_[9] = 0.03;
      constants_[10] = 5.4;
      constants_[11] = 140;
      constants_[12] = 2;
      constants_[13] = 5.405;
      constants_[14] = 0.153;
// unused      constants_[15] = 0.098;
      constants_[16] = 14.838;
      constants_[17] = 0.00029;
      constants_[18] = 0.0000398;
      constants_[19] = 0.000592;
//unused      constants_[20] = 0.294;
//unused      constants_[21] = 2.724;
      constants_[22] = 1;
      constants_[23] = 40;
      constants_[24] = 1000;
      constants_[25] = 0.1;
      constants_[26] = 2.5;
      constants_[27] = 0.35;
      constants_[28] = 1.38;
      constants_[29] = 87.5;
      constants_[30] = 0.1238;
      constants_[31] = 0.0005;
      constants_[32] = 0.0146;
      constants_[33] = 0.15;
      constants_[34] = 0.045;
      constants_[35] = 0.06;
      constants_[36] = 0.005;
      constants_[37] = 1.5;
      constants_[38] = 2.5;
      constants_[39] = 1;
      constants_[40] = 0.102;
      constants_[41] = 0.0038;
      constants_[42] = 0.00025;
      constants_[43] = 0.00036;
      constants_[44] = 0.006375;
      constants_[45] = 0.2;
      constants_[46] = 0.001;
      constants_[47] = 10;
      constants_[48] = 0.3;
      constants_[49] = 0.4;
      constants_[50] = 0.00025;
      constants_[51] = 0.001094;
      constants_[52] = 0.00005468;
   }

   switch (cellType)
   {
     case 0: // Endo
      g_Ks_ = 0.392;
      g_to_ = 0.073;
      P_NaK_ = 3.0;
      g_NaL_ = 0.15;
      s_switch_ = 0;
      break;
     case 1: // Mid
      g_Ks_ = 0.098;
      g_to_ = 0.294;
      P_NaK_ = 3.1;
      g_NaL_ = 0.3;
      s_switch_ = 1;
      break;
     case 2: // Epi
      g_Ks_ = 0.392;
      g_to_ = 0.294;
      P_NaK_ = 3.0;
      g_NaL_ = 0.15;
      s_switch_ = 1;
      break;
     default:
      assert(false);
   }
}

void TT06_RRG::initStates(int cellType)
{
   switch (cellType)
   {
     case 0: // endo
      states_[0] = -86.709;
      states_[1] = 138.4;
      states_[2] = 10.355;
      states_[3] = 0.00013;
      states_[4] = 0.00448;
      states_[5] = 0.476;
      states_[6] = 0.0087;
      states_[7] = 0.00155;
      states_[8] = 0.7573;
      states_[9] = 0.7225;
      states_[10] = 0.00036;
      states_[11] = 3.164e-5;
      states_[12] = 0.8009;
      states_[13] = 0.9778;
      states_[14] = 0.9953;
      states_[15] = 0.3212;
      states_[16] = 2.235e-8;
      states_[17] = 3.715;
      states_[18] = 0.9068;
      states_[19] = 0.066;
      break;
     case 1: // mid
      states_[0] = -85.423;
      states_[1] = 138.52;
      states_[2] = 10.132;
      states_[3] = 0.000153;
      states_[4] = 0.0165;
      states_[5] = 0.473;
      states_[6] = 0.0174;
      states_[7] = 0.00165;
      states_[8] = 0.749;
      states_[9] = 0.6788;
      states_[10] = 0.00042;
      states_[11] = 3.288e-5;
      states_[12] = 0.7026;
      states_[13] = 0.9526;
      states_[14] = 0.9942;
      states_[15] = 0.999998;
      states_[16] = 2.347e-8;
      states_[17] = 4.272;
      states_[18] = 0.8978;
      states_[19] = 0.066;
      break;
     case 2: // epi
      states_[0] = -85.23;
      states_[1] = 136.89;
      states_[2] = 8.604;
      states_[3] = 0.000126;
      states_[4] = 0.00621;
      states_[5] = 0.4712;
      states_[6] = 0.0095;
      states_[7] = 0.00172;
      states_[8] = 0.7444;
      states_[9] = 0.7045;
      states_[10] = 0.00036;
      states_[11] = 3.373e-5;
      states_[12] = 0.7888;
      states_[13] = 0.9755;
      states_[14] = 0.9953;
      states_[15] = 0.999998;
      states_[16] = 2.42e-8;
      states_[17] = 3.64;
      states_[18] = 0.9073;
      states_[19] = 0.066;
      break;
     default:
      assert(false);
   }
   
}

double TT06_RRG::computeRates(double dt, double iStim)
{
   double algebraic[70];
   double rates[20];
   
   algebraic[12] = iStim;
   algebraic[7] = 1.00000/(1.00000+(exp(((states_[0]+20.0000)/7.00000))));
   algebraic[20] =  1102.50*(exp((- (pow((states_[0]+27.0000), 2.00000))/225.000)))+200.000/(1.00000+(exp(((13.0000 - states_[0])/10.0000))))+180.000/(1.00000+(exp(((states_[0]+30.0000)/10.0000))))+20.0000;
   rates[12] = (algebraic[7] - states_[12])/algebraic[20];
   algebraic[8] = 0.670000/(1.00000+(exp(((states_[0]+35.0000)/7.00000))))+0.330000;
   algebraic[21] =  562.000*(exp((- (pow((states_[0]+27.0000), 2.00000))/240.000)))+31.0000/(1.00000+(exp(((25.0000 - states_[0])/10.0000))))+80.0000/(1.00000+(exp(((states_[0]+30.0000)/10.0000))));
   rates[13] = (algebraic[8] - states_[13])/algebraic[21];
   algebraic[9] = 0.600000/(1.00000+(pow((states_[10]/0.0500000), 2.00000)))+0.400000;
   algebraic[22] = 80.0000/(1.00000+(pow((states_[10]/0.0500000), 2.00000)))+2.00000;
   rates[14] = (algebraic[9] - states_[14])/algebraic[22];
   switch (s_switch_)
   {
     case 0:
      algebraic[10] = 1.00000/(1.00000+(exp(((states_[0]+28.0000)/5.00000))));
      algebraic[23] =  1000.00*(exp((- (pow((states_[0]+67.0000), 2.00000))/1000.00)))+8.00000;
      break;
     case 1:
      algebraic[10] = 1.00000/(1.00000+(exp(((states_[0]+20.0000)/5.00000))));
      algebraic[23] =  85.0000*(exp((- (pow((states_[0]+45.0000), 2.00000))/320.000)))+5.00000/(1.00000+(exp(((states_[0] - 20.0000)/5.00000))))+3.00000;
      break;
     default:
      assert(false);
   }

   rates[15] = (algebraic[10] - states_[15])/algebraic[23];
   algebraic[11] = 1.00000/(1.00000+(exp(((20.0000 - states_[0])/6.00000))));
   algebraic[24] =  9.50000*(exp((- (pow((states_[0]+40.0000), 2.00000))/1800.00)))+0.800000;
   rates[16] = (algebraic[11] - states_[16])/algebraic[24];
   algebraic[0] = 1.00000/(1.00000+(exp(((- 26.0000 - states_[0])/7.00000))));
   algebraic[13] = 450.000/(1.00000+(exp(((- 45.0000 - states_[0])/10.0000))));
   algebraic[26] = 6.00000/(1.00000+(exp(((states_[0]+30.0000)/11.5000))));
   algebraic[34] =  1.00000*algebraic[13]*algebraic[26];
   rates[4] = (algebraic[0] - states_[4])/algebraic[34];
   algebraic[1] = 1.00000/(1.00000+(exp(((states_[0]+88.0000)/24.0000))));
   algebraic[14] = 3.00000/(1.00000+(exp(((- 60.0000 - states_[0])/20.0000))));
   algebraic[27] = 1.12000/(1.00000+(exp(((states_[0] - 60.0000)/20.0000))));
   algebraic[35] =  1.00000*algebraic[14]*algebraic[27];
   rates[5] = (algebraic[1] - states_[5])/algebraic[35];
   algebraic[2] = 1.00000/(1.00000+(exp(((- 5.00000 - states_[0])/14.0000))));
   algebraic[15] = 1400.00/ pow((1.00000+(exp(((5.00000 - states_[0])/6.00000)))), 1.0 / 2);
   algebraic[28] = 1.00000/(1.00000+(exp(((states_[0] - 35.0000)/15.0000))));
   algebraic[36] =  1.00000*algebraic[15]*algebraic[28]+80.0000;
   rates[6] = (algebraic[2] - states_[6])/algebraic[36];
   algebraic[3] = 1.00000/(pow((1.00000+(exp(((- 56.8600 - states_[0])/9.03000)))), 2.00000));
   algebraic[16] = 1.00000/(1.00000+(exp(((- 60.0000 - states_[0])/5.00000))));
   algebraic[29] = 0.100000/(1.00000+(exp(((states_[0]+35.0000)/5.00000))))+0.100000/(1.00000+(exp(((states_[0] - 50.0000)/200.000))));
   algebraic[37] =  1.00000*algebraic[16]*algebraic[29];
   rates[7] = (algebraic[3] - states_[7])/algebraic[37];
   algebraic[4] = 1.00000/(pow((1.00000+(exp(((states_[0]+71.5500)/7.43000)))), 2.00000));
   algebraic[17] = (states_[0]<- 40.0000 ?  0.0570000*(exp((- (states_[0]+80.0000)/6.80000))) : 0.00000);
   algebraic[30] = (states_[0]<- 40.0000 ?  2.70000*(exp(( 0.0790000*states_[0])))+ 310000.*(exp(( 0.348500*states_[0]))) : 0.770000/( 0.130000*(1.00000+(exp(((states_[0]+10.6600)/- 11.1000))))));
   algebraic[38] = 1.00000/(algebraic[17]+algebraic[30]);
   rates[8] = (algebraic[4] - states_[8])/algebraic[38];
   algebraic[5] = 1.00000/(pow((1.00000+(exp(((states_[0]+71.5500)/7.43000)))), 2.00000));
   algebraic[18] = (states_[0]<- 40.0000 ? (( ( - 25428.0*(exp(( 0.244400*states_[0]))) -  6.94800e-06*(exp(( - 0.0439100*states_[0]))))*(states_[0]+37.7800))/1.00000)/(1.00000+(exp(( 0.311000*(states_[0]+79.2300))))) : 0.00000);
   algebraic[31] = (states_[0]<- 40.0000 ? ( 0.0242400*(exp(( - 0.0105200*states_[0]))))/(1.00000+(exp(( - 0.137800*(states_[0]+40.1400))))) : ( 0.600000*(exp(( 0.0570000*states_[0]))))/(1.00000+(exp(( - 0.100000*(states_[0]+32.0000))))));
   algebraic[39] = 1.00000/(algebraic[18]+algebraic[31]);
   rates[9] = (algebraic[5] - states_[9])/algebraic[39];
   algebraic[6] = 1.00000/(1.00000+(exp(((- 8.00000 - states_[0])/7.50000))));
   algebraic[19] = 1.40000/(1.00000+(exp(((- 35.0000 - states_[0])/13.0000))))+0.250000;
   algebraic[32] = 1.40000/(1.00000+(exp(((states_[0]+5.00000)/5.00000))));
   algebraic[40] = 1.00000/(1.00000+(exp(((50.0000 - states_[0])/20.0000))));
   algebraic[42] =  1.00000*algebraic[19]*algebraic[32]+algebraic[40];
   rates[11] = (algebraic[6] - states_[11])/algebraic[42];
   algebraic[55] = (( (( P_NaK_*constants_[10])/(constants_[10]+constants_[22]))*states_[2])/(states_[2]+constants_[23]))/(1.00000+ 0.124500*(exp((( - 0.100000*states_[0]*constants_[2])/( constants_[0]*constants_[1]))))+ 0.0353000*(exp((( - states_[0]*constants_[2])/( constants_[0]*constants_[1])))));
   algebraic[25] =  (( constants_[0]*constants_[1])/constants_[2])*(log((constants_[11]/states_[2])));
   algebraic[50] =  constants_[16]*(pow(states_[7], 3.00000))*states_[8]*states_[9]*(states_[0] - algebraic[25]);
   algebraic[51] =  constants_[17]*(states_[0] - algebraic[25]);
   algebraic[56] = ( constants_[24]*( (exp((( constants_[27]*states_[0]*constants_[2])/( constants_[0]*constants_[1]))))*(pow(states_[2], 3.00000))*constants_[12] -  (exp((( (constants_[27] - 1.00000)*states_[0]*constants_[2])/( constants_[0]*constants_[1]))))*(pow(constants_[11], 3.00000))*states_[3]*constants_[26]))/( ((pow(constants_[29], 3.00000))+(pow(constants_[11], 3.00000)))*(constants_[28]+constants_[12])*(1.00000+ constants_[25]*(exp((( (constants_[27] - 1.00000)*states_[0]*constants_[2])/( constants_[0]*constants_[1]))))));

   double jLinf = 1.0/exp((states_[0]+91.0)/6.1);
   jLinf *= jLinf;
   rates[19] = (jLinf - states_[19])/670.0;

   double iNaL = g_NaL_*states_[7]*states_[7]*states_[7]*states_[19]*(states_[0]-algebraic[25]);
   
   rates[2] =  (( - 1.00000*(algebraic[50]+iNaL+algebraic[51]+ 3.00000*algebraic[55]+ 3.00000*algebraic[56]))/( 1.00000*constants_[4]*constants_[2]))*constants_[3];
   algebraic[33] =  (( constants_[0]*constants_[1])/constants_[2])*(log((constants_[10]/states_[1])));
   algebraic[44] = 0.100000/(1.00000+(exp(( 0.0600000*((states_[0] - algebraic[33]) - 200.000)))));
   algebraic[45] = ( 3.00000*(exp(( 0.000200000*((states_[0] - algebraic[33])+100.000))))+(exp(( 0.100000*((states_[0] - algebraic[33]) - 10.0000)))))/(1.00000+(exp(( - 0.500000*(states_[0] - algebraic[33])))));
   algebraic[46] = algebraic[44]/(algebraic[44]+algebraic[45]);
   algebraic[47] =  constants_[13]*algebraic[46]* pow((constants_[10]/5.40000), 1.0 / 2)*(states_[0] - algebraic[33]);
   algebraic[54] =  g_to_*states_[16]*states_[15]*(states_[0] - algebraic[33]);
   algebraic[48] =  constants_[14]* pow((constants_[10]/5.40000), 1.0 / 2)*states_[4]*states_[5]*(states_[0] - algebraic[33]);
   algebraic[41] =  (( constants_[0]*constants_[1])/constants_[2])*(log(((constants_[10]+ constants_[9]*constants_[11])/(states_[1]+ constants_[9]*states_[2]))));
   algebraic[49] =  g_Ks_*(pow(states_[6], 2.00000))*(states_[0] - algebraic[41]);
   algebraic[52] = ( (( constants_[18]*states_[11]*states_[12]*states_[13]*states_[14]*4.00000*(states_[0] - 15.0000)*(pow(constants_[2], 2.00000)))/( constants_[0]*constants_[1]))*( 0.250000*states_[10]*(exp((( 2.00000*(states_[0] - 15.0000)*constants_[2])/( constants_[0]*constants_[1])))) - constants_[12]))/((exp((( 2.00000*(states_[0] - 15.0000)*constants_[2])/( constants_[0]*constants_[1])))) - 1.00000);
   algebraic[43] =  (( 0.500000*constants_[0]*constants_[1])/constants_[2])*(log((constants_[12]/states_[3])));
   algebraic[53] =  constants_[19]*(states_[0] - algebraic[43]);
   algebraic[58] = ( constants_[32]*(states_[0] - algebraic[33]))/(1.00000+(exp(((25.0000 - states_[0])/5.98000))));
   algebraic[57] = ( constants_[30]*states_[3])/(states_[3]+constants_[31]);
   rates[0] = - (algebraic[47]+algebraic[54]+algebraic[48]+algebraic[49]+algebraic[52]+algebraic[55]+algebraic[50]+algebraic[51]+algebraic[56]+algebraic[53]+algebraic[58]+algebraic[57]+iNaL);
   rates[1] =  (( - 1.00000*((algebraic[47]+algebraic[54]+algebraic[48]+algebraic[49]+algebraic[58]+algebraic[12]) -  2.00000*algebraic[55]))/( 1.00000*constants_[4]*constants_[2]))*constants_[3];
   algebraic[59] = constants_[44]/(1.00000+(pow(constants_[42], 2.00000))/(pow(states_[3], 2.00000)));
   algebraic[60] =  constants_[43]*(states_[17] - states_[3]);
   algebraic[61] =  constants_[41]*(states_[10] - states_[3]);
   algebraic[63] = 1.00000/(1.00000+( constants_[45]*constants_[46])/(pow((states_[3]+constants_[46]), 2.00000)));
   rates[3] =  algebraic[63]*((( (algebraic[60] - algebraic[59])*constants_[51])/constants_[4]+algebraic[61]) - ( 1.00000*((algebraic[53]+algebraic[57]) -  2.00000*algebraic[56])*constants_[3])/( 2.00000*1.00000*constants_[4]*constants_[2]));
   algebraic[62] = constants_[38] - (constants_[38] - constants_[39])/(1.00000+(pow((constants_[37]/states_[17]), 2.00000)));
   algebraic[65] =  constants_[34]*algebraic[62];
   rates[18] =  - algebraic[65]*states_[10]*states_[18]+ constants_[36]*(1.00000 - states_[18]);
   algebraic[64] = constants_[33]/algebraic[62];
   algebraic[66] = ( algebraic[64]*(pow(states_[10], 2.00000))*states_[18])/(constants_[35]+ algebraic[64]*(pow(states_[10], 2.00000)));
   algebraic[67] =  constants_[40]*algebraic[66]*(states_[17] - states_[10]);
   algebraic[68] = 1.00000/(1.00000+( constants_[47]*constants_[48])/(pow((states_[17]+constants_[48]), 2.00000)));
   rates[17] =  algebraic[68]*(algebraic[59] - (algebraic[67]+algebraic[60]));
   algebraic[69] = 1.00000/(1.00000+( constants_[49]*constants_[50])/(pow((states_[10]+constants_[50]), 2.00000)));
   rates[10] =  algebraic[69]*((( - 1.00000*algebraic[52]*constants_[3])/( 2.00000*1.00000*constants_[52]*constants_[2])+( algebraic[67]*constants_[51])/constants_[52]) - ( algebraic[61]*constants_[4])/constants_[52]);

   // forward euler for all states except rushLarsen for fast sodium m gate.
   for (unsigned jj=1; jj<7; ++jj)
      states_[jj] += rates[jj] * dt;
   states_[7] = algebraic[3] - (algebraic[3]-states_[7])*exp(-dt/algebraic[37]);
   for (unsigned jj=8; jj<20; ++jj)
      states_[jj] += rates[jj] * dt;

   return rates[0];
}


// This is the original variable key from cellML.
/*
   There are a total of 70 entries in the algebraic variable array.
   There are a total of 19 entries in each of the rate and state variable arrays.
   There are a total of 53 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (joule_per_mole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_millimole).
 * CONSTANTS[3] is Cm in component membrane (microF).
 * CONSTANTS[4] is V_c in component membrane (micrometre3).
 * ALGEBRAIC[47] is i_K1 in component inward_rectifier_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[54] is i_to in component transient_outward_current (picoA_per_picoF).
 * ALGEBRAIC[48] is i_Kr in component rapid_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[49] is i_Ks in component slow_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[52] is i_CaL in component L_type_Ca_current (picoA_per_picoF).
 * ALGEBRAIC[55] is i_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[50] is i_Na in component fast_sodium_current (picoA_per_picoF).
 * ALGEBRAIC[51] is i_b_Na in component sodium_background_current (picoA_per_picoF).
 * ALGEBRAIC[56] is i_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * ALGEBRAIC[53] is i_b_Ca in component calcium_background_current (picoA_per_picoF).
 * ALGEBRAIC[58] is i_p_K in component potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[57] is i_p_Ca in component calcium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[12] is i_Stim in component membrane (picoA_per_picoF).
 * CONSTANTS[5] is stim_start in component membrane (millisecond).
 * CONSTANTS[6] is stim_period in component membrane (millisecond).
 * CONSTANTS[7] is stim_duration in component membrane (millisecond).
 * CONSTANTS[8] is stim_amplitude in component membrane (picoA_per_picoF).
 * ALGEBRAIC[25] is E_Na in component reversal_potentials (millivolt).
 * ALGEBRAIC[33] is E_K in component reversal_potentials (millivolt).
 * ALGEBRAIC[41] is E_Ks in component reversal_potentials (millivolt).
 * ALGEBRAIC[43] is E_Ca in component reversal_potentials (millivolt).
 * CONSTANTS[9] is P_kna in component reversal_potentials (dimensionless).
 * CONSTANTS[10] is K_o in component potassium_dynamics (millimolar).
 * CONSTANTS[11] is Na_o in component sodium_dynamics (millimolar).
 * STATES[1] is K_i in component potassium_dynamics (millimolar).
 * STATES[2] is Na_i in component sodium_dynamics (millimolar).
 * CONSTANTS[12] is Ca_o in component calcium_dynamics (millimolar).
 * STATES[3] is Ca_i in component calcium_dynamics (millimolar).
 * CONSTANTS[13] is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
 * ALGEBRAIC[46] is xK1_inf in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[44] is alpha_K1 in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[45] is beta_K1 in component inward_rectifier_potassium_current (dimensionless).
 * CONSTANTS[14] is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[4] is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * STATES[5] is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[0] is xr1_inf in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[13] is alpha_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[26] is beta_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[34] is tau_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (millisecond).
 * ALGEBRAIC[1] is xr2_inf in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[14] is alpha_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[27] is beta_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[35] is tau_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (millisecond).
 * CONSTANTS[15] is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[6] is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[2] is xs_inf in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[15] is alpha_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[28] is beta_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[36] is tau_xs in component slow_time_dependent_potassium_current_Xs_gate (millisecond).
 * CONSTANTS[16] is g_Na in component fast_sodium_current (nanoS_per_picoF).
 * STATES[7] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[8] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[9] is j in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[3] is m_inf in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[16] is alpha_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[29] is beta_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[37] is tau_m in component fast_sodium_current_m_gate (millisecond).
 * ALGEBRAIC[4] is h_inf in component fast_sodium_current_h_gate (dimensionless).
 * ALGEBRAIC[17] is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[30] is beta_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[38] is tau_h in component fast_sodium_current_h_gate (millisecond).
 * ALGEBRAIC[5] is j_inf in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[18] is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[31] is beta_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[39] is tau_j in component fast_sodium_current_j_gate (millisecond).
 * CONSTANTS[17] is g_bna in component sodium_background_current (nanoS_per_picoF).
 * CONSTANTS[18] is g_CaL in component L_type_Ca_current (litre_per_farad_second).
 * STATES[10] is Ca_ss in component calcium_dynamics (millimolar).
 * STATES[11] is d in component L_type_Ca_current_d_gate (dimensionless).
 * STATES[12] is f in component L_type_Ca_current_f_gate (dimensionless).
 * STATES[13] is f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * STATES[14] is fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * ALGEBRAIC[6] is d_inf in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[19] is alpha_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[32] is beta_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[40] is gamma_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[42] is tau_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[7] is f_inf in component L_type_Ca_current_f_gate (dimensionless).
 * ALGEBRAIC[20] is tau_f in component L_type_Ca_current_f_gate (millisecond).
 * ALGEBRAIC[8] is f2_inf in component L_type_Ca_current_f2_gate (dimensionless).
 * ALGEBRAIC[21] is tau_f2 in component L_type_Ca_current_f2_gate (millisecond).
 * ALGEBRAIC[9] is fCass_inf in component L_type_Ca_current_fCass_gate (dimensionless).
 * ALGEBRAIC[22] is tau_fCass in component L_type_Ca_current_fCass_gate (millisecond).
 * CONSTANTS[19] is g_bca in component calcium_background_current (nanoS_per_picoF).
 * CONSTANTS[20] is g_to in component transient_outward_current (nanoS_per_picoF).
 * STATES[15] is s in component transient_outward_current_s_gate (dimensionless).
 * STATES[16] is r in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[10] is s_inf in component transient_outward_current_s_gate (dimensionless).
 * ALGEBRAIC[23] is tau_s in component transient_outward_current_s_gate (millisecond).
 * ALGEBRAIC[11] is r_inf in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[24] is tau_r in component transient_outward_current_r_gate (millisecond).
 * CONSTANTS[21] is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * CONSTANTS[22] is K_mk in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[23] is K_mNa in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[24] is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * CONSTANTS[25] is K_sat in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[26] is alpha in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[27] is gamma in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[28] is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[29] is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[30] is g_pCa in component calcium_pump_current (picoA_per_picoF).
 * CONSTANTS[31] is K_pCa in component calcium_pump_current (millimolar).
 * CONSTANTS[32] is g_pK in component potassium_pump_current (nanoS_per_picoF).
 * STATES[17] is Ca_SR in component calcium_dynamics (millimolar).
 * ALGEBRAIC[67] is i_rel in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[59] is i_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[60] is i_leak in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[61] is i_xfer in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[66] is O in component calcium_dynamics (dimensionless).
 * STATES[18] is R_prime in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[64] is k1 in component calcium_dynamics (per_millimolar2_per_millisecond).
 * ALGEBRAIC[65] is k2 in component calcium_dynamics (per_millimolar_per_millisecond).
 * CONSTANTS[33] is k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond).
 * CONSTANTS[34] is k2_prime in component calcium_dynamics (per_millimolar_per_millisecond).
 * CONSTANTS[35] is k3 in component calcium_dynamics (per_millisecond).
 * CONSTANTS[36] is k4 in component calcium_dynamics (per_millisecond).
 * CONSTANTS[37] is EC in component calcium_dynamics (millimolar).
 * CONSTANTS[38] is max_sr in component calcium_dynamics (dimensionless).
 * CONSTANTS[39] is min_sr in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[62] is kcasr in component calcium_dynamics (dimensionless).
 * CONSTANTS[40] is V_rel in component calcium_dynamics (per_millisecond).
 * CONSTANTS[41] is V_xfer in component calcium_dynamics (per_millisecond).
 * CONSTANTS[42] is K_up in component calcium_dynamics (millimolar).
 * CONSTANTS[43] is V_leak in component calcium_dynamics (per_millisecond).
 * CONSTANTS[44] is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[63] is Ca_i_bufc in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[68] is Ca_sr_bufsr in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[69] is Ca_ss_bufss in component calcium_dynamics (dimensionless).
 * CONSTANTS[45] is Buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[46] is K_buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[47] is Buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[48] is K_buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[49] is Buf_ss in component calcium_dynamics (millimolar).
 * CONSTANTS[50] is K_buf_ss in component calcium_dynamics (millimolar).
 * CONSTANTS[51] is V_sr in component calcium_dynamics (micrometre3).
 * CONSTANTS[52] is V_ss in component calcium_dynamics (micrometre3).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[4] is d/dt Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * RATES[5] is d/dt Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * RATES[6] is d/dt Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * RATES[7] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[8] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[9] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
 * RATES[11] is d/dt d in component L_type_Ca_current_d_gate (dimensionless).
 * RATES[12] is d/dt f in component L_type_Ca_current_f_gate (dimensionless).
 * RATES[13] is d/dt f2 in component L_type_Ca_current_f2_gate (dimensionless).
 * RATES[14] is d/dt fCass in component L_type_Ca_current_fCass_gate (dimensionless).
 * RATES[15] is d/dt s in component transient_outward_current_s_gate (dimensionless).
 * RATES[16] is d/dt r in component transient_outward_current_r_gate (dimensionless).
 * RATES[18] is d/dt R_prime in component calcium_dynamics (dimensionless).
 * RATES[3] is d/dt Ca_i in component calcium_dynamics (millimolar).
 * RATES[17] is d/dt Ca_SR in component calcium_dynamics (millimolar).
 * RATES[10] is d/dt Ca_ss in component calcium_dynamics (millimolar).
 * RATES[2] is d/dt Na_i in component sodium_dynamics (millimolar).
 * RATES[1] is d/dt K_i in component potassium_dynamics (millimolar).
 */

