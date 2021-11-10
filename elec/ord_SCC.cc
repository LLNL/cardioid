/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "ord_SCC.hh"
#include "object_cc.hh"
#include "mpiUtils.h"
#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>

using namespace std;

#define setDefault(name, value) objectGet(obj, #name, name, TO_STRING(value))

template<typename TTT>
static inline string TO_STRING(const TTT x)
{
   stringstream ss;
   ss << x;
   return ss.str();
}

static const char* interpName[] = {
    NULL
};


   REACTION_FACTORY(ord_SCC)(OBJECT* obj, const double _dt, const int numPoints, const ThreadTeam&)
   {
      ord_SCC::ThisReaction* reaction = new ord_SCC::ThisReaction(numPoints, _dt);

      //override the defaults
      //EDIT_PARAMETERS
      double JrelStiffConst;
      double celltype;
      double g_Na;
      setDefault(celltype, 0);
      setDefault(JrelStiffConst, 0.0050000000000000001);
      setDefault(g_Na, 14.837999999999999);
      reaction->JrelStiffConst = JrelStiffConst;
      reaction->celltype = celltype;
      reaction->g_Na = g_Na;
      bool reusingInterpolants = false;
      string fitName;
      objectGet(obj, "fit", fitName, "");
      int funcCount = sizeof(reaction->_interpolant)/sizeof(reaction->_interpolant[0])-1; //BGQ_HACKFIX, compiler bug with zero length arrays
      if (fitName != "")
      {
         OBJECT* fitObj = objectFind(fitName, "FIT");
         if (1
         )
         {
            vector<string> functions;
            objectGet(fitObj, "functions", functions);
            OBJECT* funcObj;
            if (functions.size() == funcCount)
            {
               for (int _ii=0; _ii<functions.size(); _ii++)
               {
                  OBJECT* funcObj = objectFind(functions[_ii], "FUNCTION");
                  objectGet(funcObj, "numer", reaction->_interpolant[_ii].numNumer_, "-1");
                  objectGet(funcObj, "denom", reaction->_interpolant[_ii].numDenom_, "-1");
                  objectGet(funcObj, "coeff", reaction->_interpolant[_ii].coeff_);
               }
               reusingInterpolants = true;
            }
         }
      }

      if (!reusingInterpolants)
      {
      reaction->createInterpolants(_dt);

      //save the interpolants
      if (funcCount > 0 && getRank(0) == 0)
      {
         ofstream outfile((string(obj->name) +".fit.data").c_str());
         outfile.precision(16);
         fitName = string(obj->name) + "_fit";
         outfile << obj->name << " REACTION { fit=" << fitName << "; }\n";
         outfile << fitName << " FIT {\n";
         outfile << "   functions = ";
         for (int _ii=0; _ii<funcCount; _ii++) {
            outfile << obj->name << "_interpFunc" << _ii << "_" << interpName[_ii] << " ";
         }
         outfile << ";\n";
         outfile << "}\n";

         for (int _ii=0; _ii<funcCount; _ii++)
         {
            outfile << obj->name << "_interpFunc" << _ii << "_" << interpName[_ii] << " FUNCTION { "
                    << "numer=" << reaction->_interpolant[_ii].numNumer_ << "; "
                    << "denom=" << reaction->_interpolant[_ii].numDenom_ << "; "
                    << "coeff=";
            for (int _jj=0; _jj<reaction->_interpolant[_ii].coeff_.size(); _jj++)
            {
               outfile << reaction->_interpolant[_ii].coeff_[_jj] << " ";
            }
            outfile << "; }\n";
         }
         outfile.close();
      }
      }
#ifdef USE_CUDA
      reaction->constructKernel();
#endif
      return reaction;
   }
#undef setDefault

namespace ord_SCC 
{

void ThisReaction::createInterpolants(const double _dt) {

}

#ifdef USE_CUDA

void generateInterpString(stringstream& ss, const Interpolation& interp, const char* interpVar)
{
   ss <<
   "{\n"
   "   const double _numerCoeff[]={";
    for (int _ii=interp.numNumer_-1; _ii>=0; _ii--)
   {
      if (_ii != interp.numNumer_-1) { ss << ", "; }
      ss << interp.coeff_[_ii];
   }
   ss<< "};\n"
   "   const double _denomCoeff[]={";
   for (int _ii=interp.numDenom_+interp.numNumer_-2; _ii>=interp.numNumer_; _ii--)
   {
      ss << interp.coeff_[_ii] << ", ";
   }
   ss<< "1};\n"
   "   double _inVal = " << interpVar << ";\n"
   "   double _numerator=_numerCoeff[0];\n"
   "   for (int _jj=1; _jj<sizeof(_numerCoeff)/sizeof(_numerCoeff[0]); _jj++)\n"
   "   {\n"
   "      _numerator = _numerCoeff[_jj] + _inVal*_numerator;\n"
   "   }\n"
   "   if (sizeof(_denomCoeff)/sizeof(_denomCoeff[0]) == 1)\n"
   "   {\n"
   "      _ratPoly = _numerator;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      double _denominator=_denomCoeff[0];\n"
   "      for (int _jj=1; _jj<sizeof(_denomCoeff)/sizeof(_denomCoeff[0]); _jj++)\n"
   "      {\n"
   "         _denominator = _denomCoeff[_jj] + _inVal*_denominator;\n"
   "      }\n"
   "      _ratPoly = _numerator/_denominator;\n"
   "   }\n"
   "}"
      ;
}

void ThisReaction::constructKernel()
{

   stringstream ss;
   ss.precision(16);
   ss <<
   "enum StateOffset {\n"

   "   CaMKt_off,\n"
   "   Jrelnp_off,\n"
   "   Jrelp_off,\n"
   "   a_off,\n"
   "   ap_off,\n"
   "   cai_off,\n"
   "   cajsr_off,\n"
   "   cansr_off,\n"
   "   cass_off,\n"
   "   d_off,\n"
   "   fcaf_off,\n"
   "   fcafp_off,\n"
   "   fcas_off,\n"
   "   ff_off,\n"
   "   ffp_off,\n"
   "   fs_off,\n"
   "   h_off,\n"
   "   hL_off,\n"
   "   hLp_off,\n"
   "   iF_off,\n"
   "   iFp_off,\n"
   "   iS_off,\n"
   "   iSp_off,\n"
   "   j_off,\n"
   "   jca_off,\n"
   "   ki_off,\n"
   "   kss_off,\n"
   "   m_off,\n"
   "   mL_off,\n"
   "   nai_off,\n"
   "   nass_off,\n"
   "   nca_off,\n"
   "   xk1_off,\n"
   "   xrf_off,\n"
   "   xrs_off,\n"
   "   xs1_off,\n"
   "   xs2_off,\n"
   "   NUMSTATES\n"
   "};\n"
   "extern \"C\"\n"
   "__global__ void ord_SCC_kernel(const int* _indexArray, const double* _Vm, const double* _iStim, double* _dVm, double* _state) {\n"
   "const double _dt = " << __cachedDt << ";\n"
   "const int _nCells = " << nCells_ << ";\n"

   "const double JrelStiffConst = " << JrelStiffConst << ";\n"
   "const double celltype = " << celltype << ";\n"
   "const double g_Na = " << g_Na << ";\n"
   "const int _ii = threadIdx.x + blockIdx.x*blockDim.x;\n"
   "if (_ii >= _nCells) { return; }\n"
   "const double V = _Vm[_indexArray[_ii]];\n"
   "double _ratPoly;\n"

   "const double CaMKt = _state[_ii+CaMKt_off*_nCells];\n"
   "const double Jrelnp = _state[_ii+Jrelnp_off*_nCells];\n"
   "const double Jrelp = _state[_ii+Jrelp_off*_nCells];\n"
   "const double a = _state[_ii+a_off*_nCells];\n"
   "const double ap = _state[_ii+ap_off*_nCells];\n"
   "const double cai = _state[_ii+cai_off*_nCells];\n"
   "const double cajsr = _state[_ii+cajsr_off*_nCells];\n"
   "const double cansr = _state[_ii+cansr_off*_nCells];\n"
   "const double cass = _state[_ii+cass_off*_nCells];\n"
   "const double d = _state[_ii+d_off*_nCells];\n"
   "const double fcaf = _state[_ii+fcaf_off*_nCells];\n"
   "const double fcafp = _state[_ii+fcafp_off*_nCells];\n"
   "const double fcas = _state[_ii+fcas_off*_nCells];\n"
   "const double ff = _state[_ii+ff_off*_nCells];\n"
   "const double ffp = _state[_ii+ffp_off*_nCells];\n"
   "const double fs = _state[_ii+fs_off*_nCells];\n"
   "const double h = _state[_ii+h_off*_nCells];\n"
   "const double hL = _state[_ii+hL_off*_nCells];\n"
   "const double hLp = _state[_ii+hLp_off*_nCells];\n"
   "const double iF = _state[_ii+iF_off*_nCells];\n"
   "const double iFp = _state[_ii+iFp_off*_nCells];\n"
   "const double iS = _state[_ii+iS_off*_nCells];\n"
   "const double iSp = _state[_ii+iSp_off*_nCells];\n"
   "const double j = _state[_ii+j_off*_nCells];\n"
   "const double jca = _state[_ii+jca_off*_nCells];\n"
   "const double ki = _state[_ii+ki_off*_nCells];\n"
   "const double kss = _state[_ii+kss_off*_nCells];\n"
   "const double m = _state[_ii+m_off*_nCells];\n"
   "const double mL = _state[_ii+mL_off*_nCells];\n"
   "const double nai = _state[_ii+nai_off*_nCells];\n"
   "const double nass = _state[_ii+nass_off*_nCells];\n"
   "const double nca = _state[_ii+nca_off*_nCells];\n"
   "const double xk1 = _state[_ii+xk1_off*_nCells];\n"
   "const double xrf = _state[_ii+xrf_off*_nCells];\n"
   "const double xrs = _state[_ii+xrs_off*_nCells];\n"
   "const double xs1 = _state[_ii+xs1_off*_nCells];\n"
   "const double xs2 = _state[_ii+xs2_off*_nCells];\n"
   "//get the gate updates (diagonalized exponential integrator)\n"
   "double v = V;\n"
   "double ko = 5.4000000000000004;\n"
   "double _expensive_functions_003 = exp(-0.18996960486322187*v);\n"
   "double mLss = 1.0/(0.00029157958563553099*_expensive_functions_003 + 1.0);\n"
   "double _expensive_functions_004 = exp(0.13354700854700854*v);\n"
   "double hLss = 1.0/(120578.15595522427*_expensive_functions_004 + 1.0);\n"
   "double thL = 200.0;\n"
   "double _expensive_functions_005 = exp(0.13354700854700854*v);\n"
   "double hLpss = 1.0/(275969.29038698709*_expensive_functions_005 + 1.0);\n"
   "double thLp = 3.0*thL;\n"
   "double _expensive_functions_006 = exp(-0.067476383265856948*v);\n"
   "double ass = 1.0/(2.6316508161673635*_expensive_functions_006 + 1.0);\n"
   "double _expensive_functions_007 = exp(-0.034035137876343539*v);\n"
   "double _expensive_functions_008 = exp(0.034035137876343539*v);\n"
   "double ta = 1.0515000000000001/(3.5/(30.069572727397507*_expensive_functions_008 + 1.0) + 1.0/(2.2621017070578837*_expensive_functions_007 + 1.2089000000000001));\n"
   "double _expensive_functions_009 = exp(0.17510068289266328*v);\n"
   "double iss = 1.0/(2194.970764538301*_expensive_functions_009 + 1.0);\n"
   "double __melodee_temp_007 = celltype == 1;\n"
   "double delta_epi;\n"
   "if (__melodee_temp_007)\n"
   "{\n"
   "   double _expensive_functions_010 = exp(0.20000000000000001*v);\n"
   "   delta_epi = 1.0 - 0.94999999999999996/(1202604.2841647768*_expensive_functions_010 + 1.0);\n"
   "}\n"
   "else\n"
   "{\n"
   "   delta_epi = 1.0;\n"
   "}\n"
   "double _expensive_functions_010 = exp(-0.01*v);\n"
   "double _expensive_functions_011 = exp(0.060277275467148887*v);\n"
   "double tiF = 4.5620000000000003 + (1.0/(0.14468698421272827*_expensive_functions_010 + 1.6300896349780942*_expensive_functions_011));\n"
   "double _expensive_functions_012 = exp(-0.016934801016088061*v);\n"
   "double _expensive_functions_013 = exp(0.12377769525931426*v);\n"
   "double tiS = 23.620000000000001 + (1.0/(0.00027617763953377436*_expensive_functions_012 + 0.024208962804604526*_expensive_functions_013));\n"
   "double ORd_tiF = delta_epi*tiF;\n"
   "double ORd_tiS = delta_epi*tiS;\n"
   "double _expensive_functions_015 = exp(-0.067476383265856948*v);\n"
   "double apss = 1.0/(5.1674284622306663*_expensive_functions_015 + 1.0);\n"
   "double _expensive_functions_016 = exp(0.062932662051604776*v);\n"
   "double _expensive_functions_017 = exp(-4.6425255338904359*v);\n"
   "double dti_develop = 1.3540000000000001 + 0.0001/(2.6591269045230603e-5*_expensive_functions_016 + 4.5541779737128264e+24*_expensive_functions_017);\n"
   "double _expensive_functions_018 = exp(0.050000000000000003*v);\n"
   "double dti_recover = 1.0 - 0.5/(33.115451958692312*_expensive_functions_018 + 1.0);\n"
   "double tiFp = ORd_tiF*dti_develop*dti_recover;\n"
   "double tiSp = ORd_tiS*dti_develop*dti_recover;\n"
   "double _expensive_functions_019 = exp(-0.23640661938534277*v);\n"
   "double dss = 1.0/(0.39398514226669484*_expensive_functions_019 + 1.0);\n"
   "double _expensive_functions_020 = exp(0.089999999999999997*v);\n"
   "double _expensive_functions_021 = exp(-0.050000000000000003*v);\n"
   "double td = 0.59999999999999998 + 1.0/(3.5254214873653824*_expensive_functions_020 + 0.74081822068171788*_expensive_functions_021);\n"
   "double _expensive_functions_022 = exp(0.27056277056277056*v);\n"
   "double fss = 1.0/(199.86038496778565*_expensive_functions_022 + 1.0);\n"
   "double _expensive_functions_023 = exp(0.10000000000000001*v);\n"
   "double _expensive_functions_024 = exp(-0.10000000000000001*v);\n"
   "double tff = 7.0 + 1.0/(0.033250752445187923*_expensive_functions_023 + 0.00060900877456475714*_expensive_functions_024);\n"
   "double _expensive_functions_025 = exp(-0.25*v);\n"
   "double _expensive_functions_026 = exp(0.16666666666666666*v);\n"
   "double tfs = 1000.0 + 1.0/(1.0027667890106652e-5*_expensive_functions_025 + 8.0534156181248854e-5*_expensive_functions_026);\n"
   "double fcass = fss;\n"
   "double _expensive_functions_027 = exp(-0.14285714285714285*v);\n"
   "double _expensive_functions_028 = exp(0.14285714285714285*v);\n"
   "double tfcaf = 7.0 + 1.0/(0.070831798097406196*_expensive_functions_027 + 0.022588724880310371*_expensive_functions_028);\n"
   "double _expensive_functions_029 = exp(0.14285714285714285*v);\n"
   "double _expensive_functions_030 = exp(-0.33333333333333331*v);\n"
   "double tfcas = 100.0 + 1.0/(0.00012*_expensive_functions_029 + 0.00012*_expensive_functions_030);\n"
   "double tjca = 75.0;\n"
   "double tffp = 2.5*tff;\n"
   "double tfcafp = 2.5*tfcaf;\n"
   "double Kmn = 0.002;\n"
   "double k2n = 1000.0;\n"
   "double km2n = 1.0*jca;\n"
   "double anca = 1.0/(k2n/km2n + ((Kmn/cass + 1.0)*(Kmn/cass + 1.0)*(Kmn/cass + 1.0)*(Kmn/cass + 1.0)));\n"
   "double _expensive_functions_038 = exp(-0.14729709824716453*v);\n"
   "double xrss = 1.0/(0.29287308872377504*_expensive_functions_038 + 1.0);\n"
   "double _expensive_functions_039 = exp(0.25846471956577927*v);\n"
   "double _expensive_functions_040 = exp(-0.049067713444553483*v);\n"
   "double txrf = 12.98 + 1.0/(0.0001020239312894894*_expensive_functions_039 + 0.00042992960891929087*_expensive_functions_040);\n"
   "double _expensive_functions_041 = exp(0.13596193065941536*v);\n"
   "double _expensive_functions_042 = exp(-0.038550501156515031*v);\n"
   "double txrs = 1.865 + 1.0/(0.00059224200368093944*_expensive_functions_041 + 3.5499661118024631e-5*_expensive_functions_042);\n"
   "double _expensive_functions_047 = exp(-0.11195700850873264*v);\n"
   "double xs1ss = 1.0/(0.27288596035656526*_expensive_functions_047 + 1.0);\n"
   "double _expensive_functions_048 = exp(0.056179775280898875*v);\n"
   "double _expensive_functions_049 = exp(-0.0043478260869565218*v);\n"
   "double txs1 = 817.29999999999995 + 1.0/(0.0035040677630748581*_expensive_functions_048 + 0.0005184809083581659*_expensive_functions_049);\n"
   "double xs2ss = xs1ss;\n"
   "double _expensive_functions_050 = exp(-0.032258064516129031*v);\n"
   "double _expensive_functions_051 = exp(0.050000000000000003*v);\n"
   "double txs2 = 1.0/(0.0022561357010639103*_expensive_functions_050 + 0.00082084998623898806*_expensive_functions_051);\n"
   "double _expensive_functions_053 = exp((-2.5537999999999998*ko - v - 144.59)/(1.5691999999999999*ko + 3.8115000000000001));\n"
   "double xk1ss = 1.0/(_expensive_functions_053 + 1.0);\n"
   "double _expensive_functions_054 = exp(-0.049115913555992145*v);\n"
   "double _expensive_functions_055 = exp(0.014423770373575654*v);\n"
   "double txk1 = 122.2/(0.0019352007631390235*_expensive_functions_054 + 30.433647575249029*_expensive_functions_055);\n"
   "double __melodee_temp_000 = V < -40;\n"
   "double alpha_h;\n"
   "double beta_h;\n"
   "if (__melodee_temp_000)\n"
   "{\n"
   "   double _expensive_functions_065 = exp(-0.14705882352941177*V);\n"
   "   alpha_h = 4.4312679295805147e-7*_expensive_functions_065;\n"
   "   double _expensive_functions_066 = exp(0.34849999999999998*V);\n"
   "   double _expensive_functions_067 = exp(0.079000000000000001*V);\n"
   "   beta_h = 310000*_expensive_functions_066 + 2.7000000000000002*_expensive_functions_067;\n"
   "}\n"
   "else\n"
   "{\n"
   "   alpha_h = 0;\n"
   "   double _expensive_functions_065 = exp(-0.0900900900900901*V);\n"
   "   beta_h = 0.77000000000000002/(0.049758141083938695*_expensive_functions_065 + 0.13);\n"
   "}\n"
   "double _expensive_functions_065 = exp(0.13458950201884254*V);\n"
   "double h_inf = 4.3210917837689708e-9/((_expensive_functions_065 + 6.5735011856460265e-5)*(_expensive_functions_065 + 6.5735011856460265e-5));\n"
   "double tau_h = (1.0/(alpha_h + beta_h));\n"
   "double __melodee_temp_001 = V < -40;\n"
   "double alpha_j;\n"
   "double beta_j;\n"
   "if (__melodee_temp_001)\n"
   "{\n"
   "   double _expensive_functions_066 = exp(0.311*V);\n"
   "   double _expensive_functions_067 = exp(0.24440000000000001*V);\n"
   "   double _expensive_functions_068 = exp(-0.043909999999999998*V);\n"
   "   alpha_j = (-25428*_expensive_functions_067 - 6.9480000000000002e-6*_expensive_functions_068)*(V + 37.780000000000001)/(50262745825.953987*_expensive_functions_066 + 1);\n"
   "   double _expensive_functions_069 = exp(-0.13780000000000001*V);\n"
   "   double _expensive_functions_070 = exp(-0.01052*V);\n"
   "   beta_j = 0.024240000000000001*_expensive_functions_070/(0.003960868339904256*_expensive_functions_069 + 1);\n"
   "}\n"
   "else\n"
   "{\n"
   "   alpha_j = 0;\n"
   "   double _expensive_functions_066 = exp(-0.10000000000000001*V);\n"
   "   double _expensive_functions_067 = exp(0.057000000000000002*V);\n"
   "   beta_j = 0.59999999999999998*_expensive_functions_067/(0.040762203978366204*_expensive_functions_066 + 1);\n"
   "}\n"
   "double _expensive_functions_066 = exp(0.13458950201884254*V);\n"
   "double j_inf = 4.3210917837689708e-9/((_expensive_functions_066 + 6.5735011856460265e-5)*(_expensive_functions_066 + 6.5735011856460265e-5));\n"
   "double tau_j = (1.0/(alpha_j + beta_j));\n"
   "double _expensive_functions_067 = exp(-1.0/5.0*V - 12);\n"
   "double alpha_m = (1.0/(_expensive_functions_067 + 1));\n"
   "double _expensive_functions_068 = exp((1.0/5.0)*V + 7);\n"
   "double _expensive_functions_069 = exp((1.0/200.0)*V - 1.0/4.0);\n"
   "double beta_m = 0.10000000000000001/(_expensive_functions_069 + 1) + 0.10000000000000001/(_expensive_functions_068 + 1);\n"
   "double _expensive_functions_070 = exp(-0.11074197120708749*V);\n"
   "double m_inf = (1.0/(0.0018422115811651339*_expensive_functions_070 + 1)/(0.0018422115811651339*_expensive_functions_070 + 1));\n"
   "double tm = alpha_m*beta_m;\n"
   "double tmL = tm;\n"
   "double _expensive_functions_071 = exp(-_dt/ta);\n"
   "double _a_RLA = _expensive_functions_071 - 1;\n"
   "double _a_RLB = -ass;\n"
   "double _expensive_functions_072 = exp(-_dt/ta);\n"
   "double _ap_RLA = _expensive_functions_072 - 1;\n"
   "double _ap_RLB = -apss;\n"
   "double _expensive_functions_073 = exp(-_dt/td);\n"
   "double _d_RLA = _expensive_functions_073 - 1;\n"
   "double _d_RLB = -dss;\n"
   "double _expensive_functions_074 = exp(-_dt/tfcaf);\n"
   "double _fcaf_RLA = _expensive_functions_074 - 1;\n"
   "double _fcaf_RLB = -fcass;\n"
   "double _expensive_functions_075 = exp(-_dt/tfcafp);\n"
   "double _fcafp_RLA = _expensive_functions_075 - 1;\n"
   "double _fcafp_RLB = -fcass;\n"
   "double _expensive_functions_076 = exp(-_dt/tfcas);\n"
   "double _fcas_RLA = _expensive_functions_076 - 1;\n"
   "double _fcas_RLB = -fcass;\n"
   "double _expensive_functions_077 = exp(-_dt/tff);\n"
   "double _ff_RLA = _expensive_functions_077 - 1;\n"
   "double _ff_RLB = -fss;\n"
   "double _expensive_functions_078 = exp(-_dt/tffp);\n"
   "double _ffp_RLA = _expensive_functions_078 - 1;\n"
   "double _ffp_RLB = -fss;\n"
   "double _expensive_functions_079 = exp(-_dt/tfs);\n"
   "double _fs_RLA = _expensive_functions_079 - 1;\n"
   "double _fs_RLB = -fss;\n"
   "double _expensive_functions_080 = exp(-_dt/tau_h);\n"
   "double _h_RLA = _expensive_functions_080 - 1;\n"
   "double _h_RLB = -h_inf;\n"
   "double _expensive_functions_081 = exp(-_dt/thL);\n"
   "double _hL_RLA = _expensive_functions_081 - 1;\n"
   "double _hL_RLB = -hLss;\n"
   "double _expensive_functions_082 = exp(-_dt/thLp);\n"
   "double _hLp_RLA = _expensive_functions_082 - 1;\n"
   "double _hLp_RLB = -hLpss;\n"
   "double _expensive_functions_083 = exp(-_dt/ORd_tiF);\n"
   "double _iF_RLA = _expensive_functions_083 - 1;\n"
   "double _iF_RLB = -iss;\n"
   "double _expensive_functions_084 = exp(-_dt/tiFp);\n"
   "double _iFp_RLA = _expensive_functions_084 - 1;\n"
   "double _iFp_RLB = -iss;\n"
   "double _expensive_functions_085 = exp(-_dt/ORd_tiS);\n"
   "double _iS_RLA = _expensive_functions_085 - 1;\n"
   "double _iS_RLB = -iss;\n"
   "double _expensive_functions_086 = exp(-_dt/tiSp);\n"
   "double _iSp_RLA = _expensive_functions_086 - 1;\n"
   "double _iSp_RLB = -iss;\n"
   "double _expensive_functions_087 = exp(-_dt/tau_j);\n"
   "double _j_RLA = _expensive_functions_087 - 1;\n"
   "double _j_RLB = -j_inf;\n"
   "double _expensive_functions_088 = exp(-_dt/tjca);\n"
   "double _jca_RLA = _expensive_functions_088 - 1;\n"
   "double _jca_RLB = -fcass;\n"
   "double _expensive_functions_089 = exp(-_dt/tm);\n"
   "double _m_RLA = _expensive_functions_089 - 1;\n"
   "double _m_RLB = -m_inf;\n"
   "double _expensive_functions_090 = exp(-_dt/tmL);\n"
   "double _mL_RLA = _expensive_functions_090 - 1;\n"
   "double _mL_RLB = -mLss;\n"
   "double _expensive_functions_091 = exp(-_dt*km2n);\n"
   "double _nca_RLA = _expensive_functions_091 - 1;\n"
   "double _nca_RLB = -anca*k2n/km2n;\n"
   "double _expensive_functions_092 = exp(-_dt/txk1);\n"
   "double _xk1_RLA = _expensive_functions_092 - 1;\n"
   "double _xk1_RLB = -xk1ss;\n"
   "double _expensive_functions_093 = exp(-_dt/txrf);\n"
   "double _xrf_RLA = _expensive_functions_093 - 1;\n"
   "double _xrf_RLB = -xrss;\n"
   "double _expensive_functions_094 = exp(-_dt/txrs);\n"
   "double _xrs_RLA = _expensive_functions_094 - 1;\n"
   "double _xrs_RLB = -xrss;\n"
   "double _expensive_functions_095 = exp(-_dt/txs1);\n"
   "double _xs1_RLA = _expensive_functions_095 - 1;\n"
   "double _xs1_RLB = -xs1ss;\n"
   "double _expensive_functions_096 = exp(-_dt/txs2);\n"
   "double _xs2_RLA = _expensive_functions_096 - 1;\n"
   "double _xs2_RLB = -xs2ss;\n"
   "//get the other differential updates\n"
   "double nao = 140.0;\n"
   "double cao = 1.8;\n"
   "double R = 8314.0;\n"
   "double T = 310.0;\n"
   "double F = 96485.0;\n"
   "double L = 0.01;\n"
   "double rad = 0.0011000000000000001;\n"
   "double vcell = 3140.0*(rad*rad)*L;\n"
   "double Ageo = 6.2800000000000002*L*rad + 6.2800000000000002*(rad*rad);\n"
   "double Acap = 2*Ageo;\n"
   "double vmyo = 0.68000000000000005*vcell;\n"
   "double vnsr = 0.055199999999999999*vcell;\n"
   "double vjsr = 0.0047999999999999996*vcell;\n"
   "double vss = 0.02*vcell;\n"
   "double _expensive_functions = log(nao/nai);\n"
   "double ENa = R*T*_expensive_functions/F;\n"
   "double _expensive_functions_001 = log(ko/ki);\n"
   "double EK = R*T*_expensive_functions_001/F;\n"
   "double PKNa = 0.018329999999999999;\n"
   "double _expensive_functions_002 = log((PKNa*nao + ko)/(PKNa*nai + ki));\n"
   "double EKs = R*T*_expensive_functions_002/F;\n"
   "double vffrt = (F*F)*v/(R*T);\n"
   "double vfrt = F*v/(R*T);\n"
   "double KmCaMK = 0.14999999999999999;\n"
   "double aCaMK = 0.050000000000000003;\n"
   "double bCaMK = 0.00068000000000000005;\n"
   "double CaMKo = 0.050000000000000003;\n"
   "double KmCaM = 0.0015;\n"
   "double CaMKb = CaMKo*(1.0 - CaMKt)/(KmCaM/cass + 1.0);\n"
   "double CaMKa = CaMKb + CaMKt;\n"
   "double CaMKt_diff = CaMKb*aCaMK*(CaMKb + CaMKt) - CaMKt*bCaMK;\n"
   "double GNaL = 0.0074999999999999997;\n"
   "double __melodee_temp_006 = celltype == 1;\n"
   "double ORd_GNaL;\n"
   "if (__melodee_temp_006)\n"
   "{\n"
   "   ORd_GNaL = 0.59999999999999998*GNaL;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_GNaL = GNaL;\n"
   "}\n"
   "double fINaLp = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double INaL = ORd_GNaL*mL*(-ENa + v)*(fINaLp*hLp + hL*(1.0 - fINaLp));\n"
   "double _expensive_functions_014 = exp(0.0066137566137566143*v);\n"
   "double AiF = 1.0/(0.24348537187522867*_expensive_functions_014 + 1.0);\n"
   "double AiS = 1.0 - AiF;\n"
   "double i = AiF*iF + AiS*iS;\n"
   "double ip = AiF*iFp + AiS*iSp;\n"
   "double Gto = 0.02;\n"
   "double __melodee_temp_009 = celltype == 1;\n"
   "double ORd_Gto;\n"
   "if (__melodee_temp_009)\n"
   "{\n"
   "   ORd_Gto = 4.0*Gto;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_008 = celltype == 2;\n"
   "   if (__melodee_temp_008)\n"
   "   {\n"
   "      ORd_Gto = 4.0*Gto;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      ORd_Gto = Gto;\n"
   "   }\n"
   "}\n"
   "double fItop = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double Ito = ORd_Gto*(-EK + v)*(a*i*(1.0 - fItop) + ap*fItop*ip);\n"
   "double Aff = 0.59999999999999998;\n"
   "double Afs = 1.0 - Aff;\n"
   "double f = Aff*ff + Afs*fs;\n"
   "double _expensive_functions_031 = exp(0.10000000000000001*v);\n"
   "double Afcaf = 0.29999999999999999 + 0.59999999999999998/(0.36787944117144233*_expensive_functions_031 + 1.0);\n"
   "double Afcas = 1.0 - Afcaf;\n"
   "double fca = Afcaf*fcaf + Afcas*fcas;\n"
   "double fp = Aff*ffp + Afs*fs;\n"
   "double fcap = Afcaf*fcafp + Afcas*fcas;\n"
   "double _expensive_functions_032 = exp(2.0*vfrt);\n"
   "double _expensive_functions_033 = exp(2.0*vfrt);\n"
   "double PhiCaL = 4.0*vffrt*(_expensive_functions_033*cass - 0.34100000000000003*cao)/(_expensive_functions_032 - 1.0);\n"
   "double _expensive_functions_034 = exp(1.0*vfrt);\n"
   "double _expensive_functions_035 = exp(1.0*vfrt);\n"
   "double PhiCaNa = 1.0*vffrt*(0.75*_expensive_functions_035*nass - 0.75*nao)/(_expensive_functions_034 - 1.0);\n"
   "double _expensive_functions_036 = exp(1.0*vfrt);\n"
   "double _expensive_functions_037 = exp(1.0*vfrt);\n"
   "double PhiCaK = 1.0*vffrt*(0.75*_expensive_functions_037*kss - 0.75*ko)/(_expensive_functions_036 - 1.0);\n"
   "double zca = 2.0;\n"
   "double PCa = 0.0001;\n"
   "double __melodee_temp_011 = celltype == 1;\n"
   "double ORd_PCa;\n"
   "if (__melodee_temp_011)\n"
   "{\n"
   "   ORd_PCa = 1.2*PCa;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_010 = celltype == 2;\n"
   "   if (__melodee_temp_010)\n"
   "   {\n"
   "      ORd_PCa = 2.5*PCa;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      ORd_PCa = PCa;\n"
   "   }\n"
   "}\n"
   "double PCap = 1.1000000000000001*ORd_PCa;\n"
   "double PCaNa = 0.00125*ORd_PCa;\n"
   "double PCaK = 0.00035740000000000001*ORd_PCa;\n"
   "double PCaNap = 0.00125*PCap;\n"
   "double PCaKp = 0.00035740000000000001*PCap;\n"
   "double fICaLp = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double ICaL = ORd_PCa*PhiCaL*d*(1.0 - fICaLp)*(f*(1.0 - nca) + fca*jca*nca) + PCap*PhiCaL*d*fICaLp*(fcap*jca*nca + fp*(1.0 - nca));\n"
   "double ICaNa = PCaNa*PhiCaNa*d*(1.0 - fICaLp)*(f*(1.0 - nca) + fca*jca*nca) + PCaNap*PhiCaNa*d*fICaLp*(fcap*jca*nca + fp*(1.0 - nca));\n"
   "double ICaK = PCaK*PhiCaK*d*(1.0 - fICaLp)*(f*(1.0 - nca) + fca*jca*nca) + PCaKp*PhiCaK*d*fICaLp*(fcap*jca*nca + fp*(1.0 - nca));\n"
   "double _expensive_functions_043 = exp(0.026171159382360639*v);\n"
   "double Axrf = 1.0/(4.197299094734718*_expensive_functions_043 + 1.0);\n"
   "double Axrs = 1.0 - Axrf;\n"
   "double xr = Axrf*xrf + Axrs*xrs;\n"
   "double _expensive_functions_044 = exp(0.013333333333333334*v);\n"
   "double _expensive_functions_045 = exp(0.033333333333333333*v);\n"
   "double rkr = 1.0/((2.0820090840784555*_expensive_functions_044 + 1.0)*(0.71653131057378927*_expensive_functions_045 + 1.0));\n"
   "double GKr = 0.045999999999999999;\n"
   "double __melodee_temp_013 = celltype == 1;\n"
   "double ORd_GKr;\n"
   "if (__melodee_temp_013)\n"
   "{\n"
   "   ORd_GKr = 1.3*GKr;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_012 = celltype == 2;\n"
   "   if (__melodee_temp_012)\n"
   "   {\n"
   "      ORd_GKr = 0.80000000000000004*GKr;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      ORd_GKr = GKr;\n"
   "   }\n"
   "}\n"
   "double _expensive_functions_046 = sqrt(ko);\n"
   "double IKr = 0.43033148291193518*ORd_GKr*_expensive_functions_046*rkr*xr*(-EK + v);\n"
   "double _expensive_functions_052 = pow((1.0/cai), 1.3999999999999999);\n"
   "double KsCa = 1.0 + 0.59999999999999998/(6.4818210260626455e-7*_expensive_functions_052 + 1.0);\n"
   "double GKs = 0.0033999999999999998;\n"
   "double __melodee_temp_014 = celltype == 1;\n"
   "double ORd_GKs;\n"
   "if (__melodee_temp_014)\n"
   "{\n"
   "   ORd_GKs = 1.3999999999999999*GKs;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_GKs = GKs;\n"
   "}\n"
   "double IKs = KsCa*ORd_GKs*xs1*xs2*(-EKs + v);\n"
   "double _expensive_functions_056 = exp(-0.27388602127883704*ko + 0.10534077741493732*v);\n"
   "double rk1 = 1.0/(69220.632210676704*_expensive_functions_056 + 1.0);\n"
   "double GK1 = 0.1908;\n"
   "double __melodee_temp_016 = celltype == 1;\n"
   "double ORd_GK1;\n"
   "if (__melodee_temp_016)\n"
   "{\n"
   "   ORd_GK1 = 1.2*GK1;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_015 = celltype == 2;\n"
   "   if (__melodee_temp_015)\n"
   "   {\n"
   "      ORd_GK1 = 1.3*GK1;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      ORd_GK1 = GK1;\n"
   "   }\n"
   "}\n"
   "double _expensive_functions_057 = sqrt(ko);\n"
   "double IK1 = ORd_GK1*_expensive_functions_057*rk1*xk1*(-EK + v);\n"
   "double kna1 = 15.0;\n"
   "double kna2 = 5.0;\n"
   "double kna3 = 88.120000000000005;\n"
   "double kasymm = 12.5;\n"
   "double wna = 60000.0;\n"
   "double wca = 60000.0;\n"
   "double wnaca = 5000.0;\n"
   "double kcaon = 1500000.0;\n"
   "double kcaoff = 5000.0;\n"
   "double qna = 0.52239999999999998;\n"
   "double qca = 0.16700000000000001;\n"
   "double hca = exp(F*qca*v/(R*T));\n"
   "double hna = exp(F*qna*v/(R*T));\n"
   "double h1 = 1 + nai*(hna + 1)/kna3;\n"
   "double h2 = hna*nai/(h1*kna3);\n"
   "double h3 = 1.0/h1;\n"
   "double h4 = 1.0 + nai*(1 + nai/kna2)/kna1;\n"
   "double h5 = (nai*nai)/(h4*kna1*kna2);\n"
   "double h6 = 1.0/h4;\n"
   "double h7 = 1.0 + nao*(1.0 + 1.0/hna)/kna3;\n"
   "double h8 = nao/(h7*hna*kna3);\n"
   "double h9 = 1.0/h7;\n"
   "double h10 = kasymm + 1.0 + nao*(1.0 + nao/kna2)/kna1;\n"
   "double h11 = (nao*nao)/(h10*kna1*kna2);\n"
   "double h12 = 1.0/h10;\n"
   "double k1 = cao*h12*kcaon;\n"
   "double k2 = kcaoff;\n"
   "double k3p = h9*wca;\n"
   "double k3pp = h8*wnaca;\n"
   "double k3 = k3p + k3pp;\n"
   "double k4p = h3*wca/hca;\n"
   "double k4pp = h2*wnaca;\n"
   "double k4 = k4p + k4pp;\n"
   "double k5 = kcaoff;\n"
   "double k6 = cai*h6*kcaon;\n"
   "double k7 = h2*h5*wna;\n"
   "double k8 = h11*h8*wna;\n"
   "double x1 = k2*k4*(k6 + k7) + k5*k7*(k2 + k3);\n"
   "double x2 = k1*k7*(k4 + k5) + k4*k6*(k1 + k8);\n"
   "double x3 = k1*k3*(k6 + k7) + k6*k8*(k2 + k3);\n"
   "double x4 = k2*k8*(k4 + k5) + k3*k5*(k1 + k8);\n"
   "double E1 = x1/(x1 + x2 + x3 + x4);\n"
   "double E2 = x2/(x1 + x2 + x3 + x4);\n"
   "double E3 = x3/(x1 + x2 + x3 + x4);\n"
   "double E4 = x4/(x1 + x2 + x3 + x4);\n"
   "double KmCaAct = 0.00014999999999999999;\n"
   "double allo = 1.0/((KmCaAct*KmCaAct)/(cai*cai) + 1.0);\n"
   "double zna = 1.0;\n"
   "double JncxNa = -3.0*E1*k8 - E2*k3pp + E3*k4pp + 3.0*E4*k7;\n"
   "double JncxCa = -E1*k1 + E2*k2;\n"
   "double Gncx = 0.00080000000000000004;\n"
   "double __melodee_temp_018 = celltype == 1;\n"
   "double ORd_Gncx;\n"
   "if (__melodee_temp_018)\n"
   "{\n"
   "   ORd_Gncx = 1.1000000000000001*Gncx;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_017 = celltype == 2;\n"
   "   if (__melodee_temp_017)\n"
   "   {\n"
   "      ORd_Gncx = 1.3999999999999999*Gncx;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      ORd_Gncx = Gncx;\n"
   "   }\n"
   "}\n"
   "double INaCa_i = 0.80000000000000004*ORd_Gncx*allo*(JncxCa*zca + JncxNa*zna);\n"
   "double ORd_h1 = 1 + nass*(hna + 1)/kna3;\n"
   "double ORd_h2 = hna*nass/(ORd_h1*kna3);\n"
   "double ORd_h3 = 1.0/ORd_h1;\n"
   "double ORd_h4 = 1.0 + nass*(1 + nass/kna2)/kna1;\n"
   "double ORd_h5 = (nass*nass)/(ORd_h4*kna1*kna2);\n"
   "double ORd_h6 = 1.0/ORd_h4;\n"
   "double ORd_h7 = 1.0 + nao*(1.0 + 1.0/hna)/kna3;\n"
   "double ORd_h8 = nao/(ORd_h7*hna*kna3);\n"
   "double ORd_h9 = 1.0/ORd_h7;\n"
   "double ORd_h10 = kasymm + 1.0 + nao*(1 + nao/kna2)/kna1;\n"
   "double ORd_h11 = (nao*nao)/(ORd_h10*kna1*kna2);\n"
   "double ORd_h12 = 1.0/ORd_h10;\n"
   "double ORd_k1 = ORd_h12*cao*kcaon;\n"
   "double ORd_k2 = kcaoff;\n"
   "double ORd_k3p = ORd_h9*wca;\n"
   "double ORd_k3pp = ORd_h8*wnaca;\n"
   "double ORd_k3 = ORd_k3p + ORd_k3pp;\n"
   "double ORd_k4p = ORd_h3*wca/hca;\n"
   "double ORd_k4pp = ORd_h2*wnaca;\n"
   "double ORd_k4 = ORd_k4p + ORd_k4pp;\n"
   "double ORd_k5 = kcaoff;\n"
   "double ORd_k6 = ORd_h6*cass*kcaon;\n"
   "double ORd_k7 = ORd_h2*ORd_h5*wna;\n"
   "double ORd_k8 = ORd_h11*ORd_h8*wna;\n"
   "double ORd_x1 = ORd_k2*ORd_k4*(ORd_k6 + ORd_k7) + ORd_k5*ORd_k7*(ORd_k2 + ORd_k3);\n"
   "double ORd_x2 = ORd_k1*ORd_k7*(ORd_k4 + ORd_k5) + ORd_k4*ORd_k6*(ORd_k1 + ORd_k8);\n"
   "double ORd_x3 = ORd_k1*ORd_k3*(ORd_k6 + ORd_k7) + ORd_k6*ORd_k8*(ORd_k2 + ORd_k3);\n"
   "double ORd_x4 = ORd_k2*ORd_k8*(ORd_k4 + ORd_k5) + ORd_k3*ORd_k5*(ORd_k1 + ORd_k8);\n"
   "double ORd_E1 = ORd_x1/(ORd_x1 + ORd_x2 + ORd_x3 + ORd_x4);\n"
   "double ORd_E2 = ORd_x2/(ORd_x1 + ORd_x2 + ORd_x3 + ORd_x4);\n"
   "double ORd_E3 = ORd_x3/(ORd_x1 + ORd_x2 + ORd_x3 + ORd_x4);\n"
   "double ORd_E4 = ORd_x4/(ORd_x1 + ORd_x2 + ORd_x3 + ORd_x4);\n"
   "double ORd_KmCaAct = 0.00014999999999999999;\n"
   "double ORd_allo = 1.0/((ORd_KmCaAct*ORd_KmCaAct)/(cass*cass) + 1.0);\n"
   "double ORd_JncxNa = -3.0*ORd_E1*ORd_k8 - ORd_E2*ORd_k3pp + ORd_E3*ORd_k4pp + 3.0*ORd_E4*ORd_k7;\n"
   "double ORd_JncxCa = -ORd_E1*ORd_k1 + ORd_E2*ORd_k2;\n"
   "double INaCa_ss = 0.20000000000000001*ORd_Gncx*ORd_allo*(ORd_JncxCa*zca + ORd_JncxNa*zna);\n"
   "double k1p = 949.5;\n"
   "double k1m = 182.40000000000001;\n"
   "double k2p = 687.20000000000005;\n"
   "double k2m = 39.399999999999999;\n"
   "double ORd_k3p_001 = 1899.0;\n"
   "double k3m = 79300.0;\n"
   "double ORd_k4p_001 = 639.0;\n"
   "double k4m = 40.0;\n"
   "double Knai0 = 9.0730000000000004;\n"
   "double Knao0 = 27.780000000000001;\n"
   "double delta = -0.155;\n"
   "double _expensive_functions_058 = exp(0.33333333333333331*F*delta*v/(R*T));\n"
   "double Knai = Knai0*_expensive_functions_058;\n"
   "double _expensive_functions_059 = exp(0.33333333333333331*F*v*(1.0 - delta)/(R*T));\n"
   "double Knao = Knao0*_expensive_functions_059;\n"
   "double Kki = 0.5;\n"
   "double Kko = 0.35820000000000002;\n"
   "double MgADP = 0.050000000000000003;\n"
   "double MgATP = 9.8000000000000007;\n"
   "double Kmgatp = 1.698e-7;\n"
   "double H = 9.9999999999999995e-8;\n"
   "double eP = 4.2000000000000002;\n"
   "double Khp = 1.698e-7;\n"
   "double Knap = 224.0;\n"
   "double Kxkur = 292.0;\n"
   "double P = eP/(H/Khp + 1.0 + ki/Kxkur + nai/Knap);\n"
   "double a1 = k1p*((nai/Knai)*(nai/Knai)*(nai/Knai))/(((1.0 + ki/Kki)*(1.0 + ki/Kki)) + ((1.0 + nai/Knai)*(1.0 + nai/Knai)*(1.0 + nai/Knai)) - 1.0);\n"
   "double b1 = MgADP*k1m;\n"
   "double a2 = k2p;\n"
   "double b2 = k2m*((nao/Knao)*(nao/Knao)*(nao/Knao))/(((1.0 + ko/Kko)*(1.0 + ko/Kko)) + ((1.0 + nao/Knao)*(1.0 + nao/Knao)*(1.0 + nao/Knao)) - 1.0);\n"
   "double a3 = ORd_k3p_001*((ko/Kko)*(ko/Kko))/(((1.0 + ko/Kko)*(1.0 + ko/Kko)) + ((1.0 + nao/Knao)*(1.0 + nao/Knao)*(1.0 + nao/Knao)) - 1.0);\n"
   "double b3 = H*P*k3m/(1.0 + MgATP/Kmgatp);\n"
   "double a4 = MgATP*ORd_k4p_001/(Kmgatp*(1.0 + MgATP/Kmgatp));\n"
   "double b4 = k4m*((ki/Kki)*(ki/Kki))/(((1.0 + ki/Kki)*(1.0 + ki/Kki)) + ((1.0 + nai/Knai)*(1.0 + nai/Knai)*(1.0 + nai/Knai)) - 1.0);\n"
   "double ORd_x1_001 = a1*a2*a4 + a1*a2*b3 + a2*b3*b4 + b2*b3*b4;\n"
   "double ORd_x2_001 = a1*a2*a3 + a2*a3*b4 + a3*b1*b4 + b1*b2*b4;\n"
   "double ORd_x3_001 = a2*a3*a4 + a3*a4*b1 + a4*b1*b2 + b1*b2*b3;\n"
   "double ORd_x4_001 = a1*a3*a4 + a1*a4*b2 + a1*b2*b3 + b2*b3*b4;\n"
   "double ORd_E1_001 = ORd_x1_001/(ORd_x1_001 + ORd_x2_001 + ORd_x3_001 + ORd_x4_001);\n"
   "double ORd_E2_001 = ORd_x2_001/(ORd_x1_001 + ORd_x2_001 + ORd_x3_001 + ORd_x4_001);\n"
   "double ORd_E3_001 = ORd_x3_001/(ORd_x1_001 + ORd_x2_001 + ORd_x3_001 + ORd_x4_001);\n"
   "double ORd_E4_001 = ORd_x4_001/(ORd_x1_001 + ORd_x2_001 + ORd_x3_001 + ORd_x4_001);\n"
   "double zk = 1.0;\n"
   "double JnakNa = 3.0*ORd_E1_001*a3 - 3.0*ORd_E2_001*b3;\n"
   "double JnakK = -2.0*ORd_E3_001*a1 + 2.0*ORd_E4_001*b1;\n"
   "double Pnak = 30;\n"
   "double __melodee_temp_020 = celltype == 1;\n"
   "double ORd_Pnak;\n"
   "if (__melodee_temp_020)\n"
   "{\n"
   "   ORd_Pnak = 0.90000000000000002*Pnak;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_019 = celltype == 2;\n"
   "   if (__melodee_temp_019)\n"
   "   {\n"
   "      ORd_Pnak = 0.69999999999999996*Pnak;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      ORd_Pnak = Pnak;\n"
   "   }\n"
   "}\n"
   "double INaK = ORd_Pnak*(JnakK*zk + JnakNa*zna);\n"
   "double _expensive_functions_060 = exp(-0.054525627044711013*v);\n"
   "double xkb = 1.0/(2.2023634509492389*_expensive_functions_060 + 1.0);\n"
   "double GKb = 0.0030000000000000001;\n"
   "double __melodee_temp_021 = celltype == 1;\n"
   "double ORd_GKb;\n"
   "if (__melodee_temp_021)\n"
   "{\n"
   "   ORd_GKb = 0.59999999999999998*GKb;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_GKb = GKb;\n"
   "}\n"
   "double IKb = ORd_GKb*xkb*(-EK + v);\n"
   "double PNab = 3.75e-10;\n"
   "double _expensive_functions_061 = exp(vfrt);\n"
   "double _expensive_functions_062 = exp(vfrt);\n"
   "double INab = PNab*vffrt*(_expensive_functions_062*nai - nao)/(_expensive_functions_061 - 1.0);\n"
   "double PCab = 2.4999999999999999e-8;\n"
   "double _expensive_functions_063 = exp(2.0*vfrt);\n"
   "double _expensive_functions_064 = exp(2.0*vfrt);\n"
   "double ICab = 4.0*PCab*vffrt*(_expensive_functions_064*cai - 0.34100000000000003*cao)/(_expensive_functions_063 - 1.0);\n"
   "double GpCa = 0.00050000000000000001;\n"
   "double IpCa = GpCa*cai/(cai + 0.00050000000000000001);\n"
   "double JdiffNa = -0.5*nai + 0.5*nass;\n"
   "double JdiffK = -0.5*ki + 0.5*kss;\n"
   "double Jdiff = -5.0*cai + 5.0*cass;\n"
   "double bt = 4.75;\n"
   "double a_rel = 0.5*bt;\n"
   "double Jrel_inf = -ICaL*a_rel/(25.62890625*(((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))) + 1.0);\n"
   "double __melodee_temp_022 = celltype == 2;\n"
   "double ORd_Jrel_inf;\n"
   "if (__melodee_temp_022)\n"
   "{\n"
   "   ORd_Jrel_inf = 1.7*Jrel_inf;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_Jrel_inf = Jrel_inf;\n"
   "}\n"
   "double tau_rel = bt/(1.0 + 0.0123/cajsr);\n"
   "double __melodee_temp_023 = tau_rel < 0.001;\n"
   "double ORd_tau_rel;\n"
   "if (__melodee_temp_023)\n"
   "{\n"
   "   ORd_tau_rel = 0.001;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_tau_rel = tau_rel;\n"
   "}\n"
   "double Jrelnp_diff = (-Jrelnp + ORd_Jrel_inf)/ORd_tau_rel;\n"
   "double btp = 1.25*bt;\n"
   "double a_relp = 0.5*btp;\n"
   "double Jrelp_inf = -ICaL*a_relp/(25.62890625*(((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))) + 1.0);\n"
   "double __melodee_temp_024 = celltype == 2;\n"
   "double ORd_Jrelp_inf;\n"
   "if (__melodee_temp_024)\n"
   "{\n"
   "   ORd_Jrelp_inf = 1.7*Jrelp_inf;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_Jrelp_inf = Jrelp_inf;\n"
   "}\n"
   "double tau_relp = btp/(1.0 + 0.0123/cajsr);\n"
   "double __melodee_temp_025 = tau_relp < 0.001;\n"
   "double ORd_tau_relp;\n"
   "if (__melodee_temp_025)\n"
   "{\n"
   "   ORd_tau_relp = 0.001;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_tau_relp = tau_relp;\n"
   "}\n"
   "double Jrelp_diff = (-Jrelp + ORd_Jrelp_inf)/ORd_tau_relp;\n"
   "double fJrelp = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double Jrel = Jrelnp*(1.0 - fJrelp) + Jrelp*fJrelp;\n"
   "double __melodee_temp_026 = Jrel*JrelStiffConst > cajsr;\n"
   "double ORd_Jrel;\n"
   "if (__melodee_temp_026)\n"
   "{\n"
   "   ORd_Jrel = cajsr/JrelStiffConst;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_Jrel = Jrel;\n"
   "}\n"
   "double Jupnp = 0.0043750000000000004*cai/(cai + 0.00092000000000000003);\n"
   "double Jupp = 0.01203125*cai/(cai + 0.00075000000000000002);\n"
   "double __melodee_temp_027 = celltype == 1;\n"
   "double ORd_Jupnp;\n"
   "double ORd_Jupp;\n"
   "if (__melodee_temp_027)\n"
   "{\n"
   "   ORd_Jupnp = 1.3*Jupnp;\n"
   "   ORd_Jupp = 1.3*Jupp;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_Jupnp = Jupnp;\n"
   "   ORd_Jupp = Jupp;\n"
   "}\n"
   "double fJupp = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double Jleak = 0.00026249999999999998*cansr;\n"
   "double Jup = -Jleak + ORd_Jupnp*(1.0 - fJupp) + ORd_Jupp*fJupp;\n"
   "double Jtr = -0.01*cajsr + 0.01*cansr;\n"
   "double cmdnmax = 0.050000000000000003;\n"
   "double __melodee_temp_028 = celltype == 1;\n"
   "double ORd_cmdnmax;\n"
   "if (__melodee_temp_028)\n"
   "{\n"
   "   ORd_cmdnmax = 1.3*cmdnmax;\n"
   "}\n"
   "else\n"
   "{\n"
   "   ORd_cmdnmax = cmdnmax;\n"
   "}\n"
   "double kmcmdn = 0.0023800000000000002;\n"
   "double trpnmax = 0.070000000000000007;\n"
   "double kmtrpn = 0.00050000000000000001;\n"
   "double BSRmax = 0.047;\n"
   "double KmBSR = 0.00087000000000000001;\n"
   "double BSLmax = 1.1240000000000001;\n"
   "double KmBSL = 0.0086999999999999994;\n"
   "double csqnmax = 10.0;\n"
   "double kmcsqn = 0.80000000000000004;\n"
   "double nass_diff = Acap*(-ICaNa - 3.0*INaCa_ss)/(F*vss) - JdiffNa;\n"
   "double ki_diff = Acap*(-IK1 - IKb - IKr - IKs + 2.0*INaK - Ito)/(F*vmyo) + JdiffK*vss/vmyo;\n"
   "double kss_diff = -Acap*ICaK/(F*vss) - JdiffK;\n"
   "double Bcai = 1.0/(ORd_cmdnmax*kmcmdn*(1.0/(cai + kmcmdn)/(cai + kmcmdn)) + kmtrpn*trpnmax*(1.0/(cai + kmtrpn)/(cai + kmtrpn)) + 1.0);\n"
   "double cai_diff = Bcai*(0.5*Acap*(-ICab + 2.0*INaCa_i - IpCa)/(F*vmyo) + Jdiff*vss/vmyo - Jup*vnsr/vmyo);\n"
   "double Bcass = 1.0/(BSLmax*KmBSL*(1.0/(KmBSL + cass)/(KmBSL + cass)) + BSRmax*KmBSR*(1.0/(KmBSR + cass)/(KmBSR + cass)) + 1.0);\n"
   "double cass_diff = Bcass*(0.5*Acap*(-ICaL + 2.0*INaCa_ss)/(F*vss) - Jdiff + ORd_Jrel*vjsr/vss);\n"
   "double cansr_diff = -Jtr*vjsr/vnsr + Jup;\n"
   "double Bcajsr = 1.0/(csqnmax*kmcsqn*(1.0/(cajsr + kmcsqn)/(cajsr + kmcsqn)) + 1.0);\n"
   "double cajsr_diff = Bcajsr*(Jtr - ORd_Jrel);\n"
   "double i_Na = (m*m*m)*g_Na*h*j*(-ENa + V);\n"
   "double i_Naitot = i_Na;\n"
   "double INa = i_Naitot;\n"
   "double nai_diff = Acap*(-INa - 3.0*INaCa_i - 3.0*INaK - INaL - INab)/(F*vmyo) + JdiffNa*vss/vmyo;\n"
   "//get Iion\n"
   "double Iion = ICaK + ICaL + ICaNa + ICab + IK1 + IKb + IKr + IKs + INa + INaCa_i + INaCa_ss + INaK + INaL + INab + IpCa + Ito;\n"
   "double Iion_001 = Iion;\n"
   "//Do the markov update (1 step rosenbrock with gauss siedel)\n"
   "int _count=0;\n"
   "do\n"
   "{\n"
   "   _count++;\n"
   "} while (_count<50);\n"
   "//EDIT_STATE\n"
   "_state[_ii+CaMKt_off*_nCells] += _dt*CaMKt_diff;\n"
   "_state[_ii+Jrelnp_off*_nCells] += _dt*Jrelnp_diff;\n"
   "_state[_ii+Jrelp_off*_nCells] += _dt*Jrelp_diff;\n"
   "_state[_ii+cai_off*_nCells] += _dt*cai_diff;\n"
   "_state[_ii+cajsr_off*_nCells] += _dt*cajsr_diff;\n"
   "_state[_ii+cansr_off*_nCells] += _dt*cansr_diff;\n"
   "_state[_ii+cass_off*_nCells] += _dt*cass_diff;\n"
   "_state[_ii+ki_off*_nCells] += _dt*ki_diff;\n"
   "_state[_ii+kss_off*_nCells] += _dt*kss_diff;\n"
   "_state[_ii+nai_off*_nCells] += _dt*nai_diff;\n"
   "_state[_ii+nass_off*_nCells] += _dt*nass_diff;\n"
   "_state[_ii+a_off*_nCells] += _a_RLA*(a+_a_RLB);\n"
   "_state[_ii+ap_off*_nCells] += _ap_RLA*(ap+_ap_RLB);\n"
   "_state[_ii+d_off*_nCells] += _d_RLA*(d+_d_RLB);\n"
   "_state[_ii+fcaf_off*_nCells] += _fcaf_RLA*(fcaf+_fcaf_RLB);\n"
   "_state[_ii+fcafp_off*_nCells] += _fcafp_RLA*(fcafp+_fcafp_RLB);\n"
   "_state[_ii+fcas_off*_nCells] += _fcas_RLA*(fcas+_fcas_RLB);\n"
   "_state[_ii+ff_off*_nCells] += _ff_RLA*(ff+_ff_RLB);\n"
   "_state[_ii+ffp_off*_nCells] += _ffp_RLA*(ffp+_ffp_RLB);\n"
   "_state[_ii+fs_off*_nCells] += _fs_RLA*(fs+_fs_RLB);\n"
   "_state[_ii+h_off*_nCells] += _h_RLA*(h+_h_RLB);\n"
   "_state[_ii+hL_off*_nCells] += _hL_RLA*(hL+_hL_RLB);\n"
   "_state[_ii+hLp_off*_nCells] += _hLp_RLA*(hLp+_hLp_RLB);\n"
   "_state[_ii+iF_off*_nCells] += _iF_RLA*(iF+_iF_RLB);\n"
   "_state[_ii+iFp_off*_nCells] += _iFp_RLA*(iFp+_iFp_RLB);\n"
   "_state[_ii+iS_off*_nCells] += _iS_RLA*(iS+_iS_RLB);\n"
   "_state[_ii+iSp_off*_nCells] += _iSp_RLA*(iSp+_iSp_RLB);\n"
   "_state[_ii+j_off*_nCells] += _j_RLA*(j+_j_RLB);\n"
   "_state[_ii+jca_off*_nCells] += _jca_RLA*(jca+_jca_RLB);\n"
   "_state[_ii+m_off*_nCells] += _m_RLA*(m+_m_RLB);\n"
   "_state[_ii+mL_off*_nCells] += _mL_RLA*(mL+_mL_RLB);\n"
   "_state[_ii+nca_off*_nCells] += _nca_RLA*(nca+_nca_RLB);\n"
   "_state[_ii+xk1_off*_nCells] += _xk1_RLA*(xk1+_xk1_RLB);\n"
   "_state[_ii+xrf_off*_nCells] += _xrf_RLA*(xrf+_xrf_RLB);\n"
   "_state[_ii+xrs_off*_nCells] += _xrs_RLA*(xrs+_xrs_RLB);\n"
   "_state[_ii+xs1_off*_nCells] += _xs1_RLA*(xs1+_xs1_RLB);\n"
   "_state[_ii+xs2_off*_nCells] += _xs2_RLA*(xs2+_xs2_RLB);\n"
   "_dVm[_indexArray[_ii]] = -Iion_001;\n"
   "}\n";

   _program_code = ss.str();
   //cout << ss.str();
   nvrtcCreateProgram(&_program,
                      _program_code.c_str(),
                      "ord_SCC_program",
                      0,
                      NULL,
                      NULL);
   nvrtcCompileProgram(_program,
                       0,
                       NULL);
   std::size_t size;
   nvrtcGetPTXSize(_program, &size);
   _ptx.resize(size);
   nvrtcGetPTX(_program, &_ptx[0]);

   cuModuleLoadDataEx(&_module, &_ptx[0], 0, 0, 0);
   cuModuleGetFunction(&_kernel, _module, "ord_SCC_kernel");
}

void ThisReaction::calc(double dt,
                ro_mgarray_ptr<int> indexArray_m,
                ro_mgarray_ptr<double> Vm_m,
                ro_mgarray_ptr<double> iStim_m,
                wo_mgarray_ptr<double> dVm_m)
{
   if (nCells_ == 0) { return; }

   {
      int errorCode=-1;
      if (blockSize_ == -1) { blockSize_ = 1024; }
      while(1)
      {
         const int* indexArrayRaw = indexArray_m.useOn(GPU).raw();
         const double* VmRaw = Vm_m.useOn(GPU).raw();
         const double* iStimRaw = iStim_m.useOn(GPU).raw();
         double* dVmRaw = dVm_m.useOn(GPU).raw();
         double* stateRaw= stateTransport_.readwrite(GPU).raw();
         void* args[] = { &indexArrayRaw,
                          &VmRaw,
                          &iStimRaw,
                          &dVmRaw,
                          &stateRaw};
         int errorCode = cuLaunchKernel(_kernel,
                                        (nCells_+blockSize_-1)/blockSize_, 1, 1,
                                        blockSize_,1,1,
                                        0, NULL,
                                        args, 0);
         if (errorCode == CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES && blockSize_ > 0)
         {
            blockSize_ /= 2;
            continue;
         }
         else if (errorCode != CUDA_SUCCESS)
         {
            printf("Cuda return %d;\n", errorCode);
            assert(0 && "Launch failed");
         }
         else
         {
            break;
         }
         //cuCtxSynchronize();
      }
   }
}

enum StateOffset {
   _CaMKt_off,
   _Jrelnp_off,
   _Jrelp_off,
   _a_off,
   _ap_off,
   _cai_off,
   _cajsr_off,
   _cansr_off,
   _cass_off,
   _d_off,
   _fcaf_off,
   _fcafp_off,
   _fcas_off,
   _ff_off,
   _ffp_off,
   _fs_off,
   _h_off,
   _hL_off,
   _hLp_off,
   _iF_off,
   _iFp_off,
   _iS_off,
   _iSp_off,
   _j_off,
   _jca_off,
   _ki_off,
   _kss_off,
   _m_off,
   _mL_off,
   _nai_off,
   _nass_off,
   _nca_off,
   _xk1_off,
   _xrf_off,
   _xrs_off,
   _xs1_off,
   _xs2_off,
   NUMSTATES
};

ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   stateTransport_.resize(nCells_*NUMSTATES);
   __cachedDt = __dt;
   blockSize_ = -1;
   _program = NULL;
   _module = NULL;
}

ThisReaction::~ThisReaction() {
   if (_program)
   {
      nvrtcDestroyProgram(&_program);
      _program = NULL;
   }
   if (_module)
   {
      cuModuleUnload(_module);
      _module = NULL;
   }
}

#else //USE_CUDA

#define width SIMDOPS_FLOAT64V_WIDTH
#define real simdops::float64v
#define load simdops::load

ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   state_.resize((nCells_+width-1)/width);
   __cachedDt = __dt;
}

ThisReaction::~ThisReaction() {}

void ThisReaction::calc(double _dt,
                ro_mgarray_ptr<int> ___indexArray,
                ro_mgarray_ptr<double> ___Vm,
                ro_mgarray_ptr<double> ,
                wo_mgarray_ptr<double> ___dVm)
{
   ro_array_ptr<int>    __indexArray = ___indexArray.useOn(CPU);
   ro_array_ptr<double> __Vm = ___Vm.useOn(CPU);
   wo_array_ptr<double> __dVm = ___dVm.useOn(CPU);

   //define the constants
   double nao = 140.0;
   double cao = 1.8;
   double ko = 5.4000000000000004;
   double R = 8314.0;
   double T = 310.0;
   double F = 96485.0;
   double L = 0.01;
   double rad = 0.0011000000000000001;
   double vcell = 3140.0*(rad*rad)*L;
   double Ageo = 6.2800000000000002*L*rad + 6.2800000000000002*(rad*rad);
   double Acap = 2*Ageo;
   double vmyo = 0.68000000000000005*vcell;
   double vnsr = 0.055199999999999999*vcell;
   double vjsr = 0.0047999999999999996*vcell;
   double vss = 0.02*vcell;
   double PKNa = 0.018329999999999999;
   double KmCaMK = 0.14999999999999999;
   double aCaMK = 0.050000000000000003;
   double bCaMK = 0.00068000000000000005;
   double CaMKo = 0.050000000000000003;
   double KmCaM = 0.0015;
   double thL = 200.0;
   double thLp = 3.0*thL;
   double GNaL = 0.0074999999999999997;
   double __melodee_temp_006 = celltype == 1;
   double ORd_GNaL;
   if (__melodee_temp_006)
   {
      ORd_GNaL = 0.59999999999999998*GNaL;
   }
   else
   {
      ORd_GNaL = GNaL;
   }
   double __melodee_temp_007 = celltype == 1;
   double delta_epi;
   if (__melodee_temp_007)
   {
   }
   else
   {
      delta_epi = 1.0;
   }
   double Gto = 0.02;
   double __melodee_temp_009 = celltype == 1;
   double ORd_Gto;
   if (__melodee_temp_009)
   {
      ORd_Gto = 4.0*Gto;
   }
   else
   {
      double __melodee_temp_008 = celltype == 2;
      if (__melodee_temp_008)
      {
         ORd_Gto = 4.0*Gto;
      }
      else
      {
         ORd_Gto = Gto;
      }
   }
   double Aff = 0.59999999999999998;
   double Afs = 1.0 - Aff;
   double tjca = 75.0;
   double Kmn = 0.002;
   double k2n = 1000.0;
   double zca = 2.0;
   double PCa = 0.0001;
   double __melodee_temp_011 = celltype == 1;
   double ORd_PCa;
   if (__melodee_temp_011)
   {
      ORd_PCa = 1.2*PCa;
   }
   else
   {
      double __melodee_temp_010 = celltype == 2;
      if (__melodee_temp_010)
      {
         ORd_PCa = 2.5*PCa;
      }
      else
      {
         ORd_PCa = PCa;
      }
   }
   double PCap = 1.1000000000000001*ORd_PCa;
   double PCaNa = 0.00125*ORd_PCa;
   double PCaK = 0.00035740000000000001*ORd_PCa;
   double PCaNap = 0.00125*PCap;
   double PCaKp = 0.00035740000000000001*PCap;
   double GKr = 0.045999999999999999;
   double __melodee_temp_013 = celltype == 1;
   double ORd_GKr;
   if (__melodee_temp_013)
   {
      ORd_GKr = 1.3*GKr;
   }
   else
   {
      double __melodee_temp_012 = celltype == 2;
      if (__melodee_temp_012)
      {
         ORd_GKr = 0.80000000000000004*GKr;
      }
      else
      {
         ORd_GKr = GKr;
      }
   }
   double _expensive_functions_046 = sqrt(ko);
   double GKs = 0.0033999999999999998;
   double __melodee_temp_014 = celltype == 1;
   double ORd_GKs;
   if (__melodee_temp_014)
   {
      ORd_GKs = 1.3999999999999999*GKs;
   }
   else
   {
      ORd_GKs = GKs;
   }
   double GK1 = 0.1908;
   double __melodee_temp_016 = celltype == 1;
   double ORd_GK1;
   if (__melodee_temp_016)
   {
      ORd_GK1 = 1.2*GK1;
   }
   else
   {
      double __melodee_temp_015 = celltype == 2;
      if (__melodee_temp_015)
      {
         ORd_GK1 = 1.3*GK1;
      }
      else
      {
         ORd_GK1 = GK1;
      }
   }
   double _expensive_functions_057 = sqrt(ko);
   double kna1 = 15.0;
   double kna2 = 5.0;
   double kna3 = 88.120000000000005;
   double kasymm = 12.5;
   double wna = 60000.0;
   double wca = 60000.0;
   double wnaca = 5000.0;
   double kcaon = 1500000.0;
   double kcaoff = 5000.0;
   double qna = 0.52239999999999998;
   double qca = 0.16700000000000001;
   double h10 = kasymm + 1.0 + nao*(1.0 + nao/kna2)/kna1;
   double h11 = (nao*nao)/(h10*kna1*kna2);
   double h12 = 1.0/h10;
   double k1 = cao*h12*kcaon;
   double k2 = kcaoff;
   double k5 = kcaoff;
   double KmCaAct = 0.00014999999999999999;
   double zna = 1.0;
   double Gncx = 0.00080000000000000004;
   double __melodee_temp_018 = celltype == 1;
   double ORd_Gncx;
   if (__melodee_temp_018)
   {
      ORd_Gncx = 1.1000000000000001*Gncx;
   }
   else
   {
      double __melodee_temp_017 = celltype == 2;
      if (__melodee_temp_017)
      {
         ORd_Gncx = 1.3999999999999999*Gncx;
      }
      else
      {
         ORd_Gncx = Gncx;
      }
   }
   double ORd_h10 = kasymm + 1.0 + nao*(1 + nao/kna2)/kna1;
   double ORd_h11 = (nao*nao)/(ORd_h10*kna1*kna2);
   double ORd_h12 = 1.0/ORd_h10;
   double ORd_k1 = ORd_h12*cao*kcaon;
   double ORd_k2 = kcaoff;
   double ORd_k5 = kcaoff;
   double ORd_KmCaAct = 0.00014999999999999999;
   double k1p = 949.5;
   double k1m = 182.40000000000001;
   double k2p = 687.20000000000005;
   double k2m = 39.399999999999999;
   double ORd_k3p_001 = 1899.0;
   double k3m = 79300.0;
   double ORd_k4p_001 = 639.0;
   double k4m = 40.0;
   double Knai0 = 9.0730000000000004;
   double Knao0 = 27.780000000000001;
   double delta = -0.155;
   double Kki = 0.5;
   double Kko = 0.35820000000000002;
   double MgADP = 0.050000000000000003;
   double MgATP = 9.8000000000000007;
   double Kmgatp = 1.698e-7;
   double H = 9.9999999999999995e-8;
   double eP = 4.2000000000000002;
   double Khp = 1.698e-7;
   double Knap = 224.0;
   double Kxkur = 292.0;
   double b1 = MgADP*k1m;
   double a2 = k2p;
   double a4 = MgATP*ORd_k4p_001/(Kmgatp*(1.0 + MgATP/Kmgatp));
   double zk = 1.0;
   double Pnak = 30;
   double __melodee_temp_020 = celltype == 1;
   double ORd_Pnak;
   if (__melodee_temp_020)
   {
      ORd_Pnak = 0.90000000000000002*Pnak;
   }
   else
   {
      double __melodee_temp_019 = celltype == 2;
      if (__melodee_temp_019)
      {
         ORd_Pnak = 0.69999999999999996*Pnak;
      }
      else
      {
         ORd_Pnak = Pnak;
      }
   }
   double GKb = 0.0030000000000000001;
   double __melodee_temp_021 = celltype == 1;
   double ORd_GKb;
   if (__melodee_temp_021)
   {
      ORd_GKb = 0.59999999999999998*GKb;
   }
   else
   {
      ORd_GKb = GKb;
   }
   double PNab = 3.75e-10;
   double PCab = 2.4999999999999999e-8;
   double GpCa = 0.00050000000000000001;
   double bt = 4.75;
   double a_rel = 0.5*bt;
   double __melodee_temp_022 = celltype == 2;
   double btp = 1.25*bt;
   double a_relp = 0.5*btp;
   double __melodee_temp_024 = celltype == 2;
   double __melodee_temp_027 = celltype == 1;
   double cmdnmax = 0.050000000000000003;
   double __melodee_temp_028 = celltype == 1;
   double ORd_cmdnmax;
   if (__melodee_temp_028)
   {
      ORd_cmdnmax = 1.3*cmdnmax;
   }
   else
   {
      ORd_cmdnmax = cmdnmax;
   }
   double kmcmdn = 0.0023800000000000002;
   double trpnmax = 0.070000000000000007;
   double kmtrpn = 0.00050000000000000001;
   double BSRmax = 0.047;
   double KmBSR = 0.00087000000000000001;
   double BSLmax = 1.1240000000000001;
   double KmBSL = 0.0086999999999999994;
   double csqnmax = 10.0;
   double kmcsqn = 0.80000000000000004;
   double _expensive_functions_081 = exp(-_dt/thL);
   double _hL_RLA = _expensive_functions_081 - 1;
   double _expensive_functions_082 = exp(-_dt/thLp);
   double _hLp_RLA = _expensive_functions_082 - 1;
   double _expensive_functions_088 = exp(-_dt/tjca);
   double _jca_RLA = _expensive_functions_088 - 1;
   for (unsigned __jj=0; __jj<(nCells_+width-1)/width; __jj++)
   {
      const int __ii = __jj*width;
      //set Vm
      double __Vm_local[width];
      {
         int cursor = 0;
         for (int __kk=0; __kk<width; __kk++)
         {
            __Vm_local[__kk] = __Vm[__indexArray[__ii+cursor]];
            if (__ii+__kk < nCells_) { cursor++; }
         }
      }
      const real V = load(&__Vm_local[0]);
      //set all state variables
      real CaMKt=load(state_[__jj].CaMKt);
      real Jrelnp=load(state_[__jj].Jrelnp);
      real Jrelp=load(state_[__jj].Jrelp);
      real a=load(state_[__jj].a);
      real ap=load(state_[__jj].ap);
      real cai=load(state_[__jj].cai);
      real cajsr=load(state_[__jj].cajsr);
      real cansr=load(state_[__jj].cansr);
      real cass=load(state_[__jj].cass);
      real d=load(state_[__jj].d);
      real fcaf=load(state_[__jj].fcaf);
      real fcafp=load(state_[__jj].fcafp);
      real fcas=load(state_[__jj].fcas);
      real ff=load(state_[__jj].ff);
      real ffp=load(state_[__jj].ffp);
      real fs=load(state_[__jj].fs);
      real h=load(state_[__jj].h);
      real hL=load(state_[__jj].hL);
      real hLp=load(state_[__jj].hLp);
      real iF=load(state_[__jj].iF);
      real iFp=load(state_[__jj].iFp);
      real iS=load(state_[__jj].iS);
      real iSp=load(state_[__jj].iSp);
      real j=load(state_[__jj].j);
      real jca=load(state_[__jj].jca);
      real ki=load(state_[__jj].ki);
      real kss=load(state_[__jj].kss);
      real m=load(state_[__jj].m);
      real mL=load(state_[__jj].mL);
      real nai=load(state_[__jj].nai);
      real nass=load(state_[__jj].nass);
      real nca=load(state_[__jj].nca);
      real xk1=load(state_[__jj].xk1);
      real xrf=load(state_[__jj].xrf);
      real xrs=load(state_[__jj].xrs);
      real xs1=load(state_[__jj].xs1);
      real xs2=load(state_[__jj].xs2);
      //get the gate updates (diagonalized exponential integrator)
      real v = V;
      real _expensive_functions_003 = exp(-0.18996960486322187*v);
      real mLss = 1.0/(0.00029157958563553099*_expensive_functions_003 + 1.0);
      real _expensive_functions_004 = exp(0.13354700854700854*v);
      real hLss = 1.0/(120578.15595522427*_expensive_functions_004 + 1.0);
      real _expensive_functions_005 = exp(0.13354700854700854*v);
      real hLpss = 1.0/(275969.29038698709*_expensive_functions_005 + 1.0);
      real _expensive_functions_006 = exp(-0.067476383265856948*v);
      real ass = 1.0/(2.6316508161673635*_expensive_functions_006 + 1.0);
      real _expensive_functions_007 = exp(-0.034035137876343539*v);
      real _expensive_functions_008 = exp(0.034035137876343539*v);
      real ta = 1.0515000000000001/(3.5/(30.069572727397507*_expensive_functions_008 + 1.0) + 1.0/(2.2621017070578837*_expensive_functions_007 + 1.2089000000000001));
      real _expensive_functions_009 = exp(0.17510068289266328*v);
      real iss = 1.0/(2194.970764538301*_expensive_functions_009 + 1.0);
      if (__melodee_temp_007)
      {
         real _expensive_functions_010 = exp(0.20000000000000001*v);
         delta_epi = 1.0 - 0.94999999999999996/(1202604.2841647768*_expensive_functions_010 + 1.0);
      }
      else
      {
      }
      real _expensive_functions_010 = exp(-0.01*v);
      real _expensive_functions_011 = exp(0.060277275467148887*v);
      real tiF = 4.5620000000000003 + (1.0/(0.14468698421272827*_expensive_functions_010 + 1.6300896349780942*_expensive_functions_011));
      real _expensive_functions_012 = exp(-0.016934801016088061*v);
      real _expensive_functions_013 = exp(0.12377769525931426*v);
      real tiS = 23.620000000000001 + (1.0/(0.00027617763953377436*_expensive_functions_012 + 0.024208962804604526*_expensive_functions_013));
      real ORd_tiF = delta_epi*tiF;
      real ORd_tiS = delta_epi*tiS;
      real _expensive_functions_015 = exp(-0.067476383265856948*v);
      real apss = 1.0/(5.1674284622306663*_expensive_functions_015 + 1.0);
      real _expensive_functions_016 = exp(0.062932662051604776*v);
      real _expensive_functions_017 = exp(-4.6425255338904359*v);
      real dti_develop = 1.3540000000000001 + 0.0001/(2.6591269045230603e-5*_expensive_functions_016 + 4.5541779737128264e+24*_expensive_functions_017);
      real _expensive_functions_018 = exp(0.050000000000000003*v);
      real dti_recover = 1.0 - 0.5/(33.115451958692312*_expensive_functions_018 + 1.0);
      real tiFp = ORd_tiF*dti_develop*dti_recover;
      real tiSp = ORd_tiS*dti_develop*dti_recover;
      real _expensive_functions_019 = exp(-0.23640661938534277*v);
      real dss = 1.0/(0.39398514226669484*_expensive_functions_019 + 1.0);
      real _expensive_functions_020 = exp(0.089999999999999997*v);
      real _expensive_functions_021 = exp(-0.050000000000000003*v);
      real td = 0.59999999999999998 + 1.0/(3.5254214873653824*_expensive_functions_020 + 0.74081822068171788*_expensive_functions_021);
      real _expensive_functions_022 = exp(0.27056277056277056*v);
      real fss = 1.0/(199.86038496778565*_expensive_functions_022 + 1.0);
      real _expensive_functions_023 = exp(0.10000000000000001*v);
      real _expensive_functions_024 = exp(-0.10000000000000001*v);
      real tff = 7.0 + 1.0/(0.033250752445187923*_expensive_functions_023 + 0.00060900877456475714*_expensive_functions_024);
      real _expensive_functions_025 = exp(-0.25*v);
      real _expensive_functions_026 = exp(0.16666666666666666*v);
      real tfs = 1000.0 + 1.0/(1.0027667890106652e-5*_expensive_functions_025 + 8.0534156181248854e-5*_expensive_functions_026);
      real fcass = fss;
      real _expensive_functions_027 = exp(-0.14285714285714285*v);
      real _expensive_functions_028 = exp(0.14285714285714285*v);
      real tfcaf = 7.0 + 1.0/(0.070831798097406196*_expensive_functions_027 + 0.022588724880310371*_expensive_functions_028);
      real _expensive_functions_029 = exp(0.14285714285714285*v);
      real _expensive_functions_030 = exp(-0.33333333333333331*v);
      real tfcas = 100.0 + 1.0/(0.00012*_expensive_functions_029 + 0.00012*_expensive_functions_030);
      real tffp = 2.5*tff;
      real tfcafp = 2.5*tfcaf;
      real km2n = 1.0*jca;
      real anca = 1.0/(k2n/km2n + ((Kmn/cass + 1.0)*(Kmn/cass + 1.0)*(Kmn/cass + 1.0)*(Kmn/cass + 1.0)));
      real _expensive_functions_038 = exp(-0.14729709824716453*v);
      real xrss = 1.0/(0.29287308872377504*_expensive_functions_038 + 1.0);
      real _expensive_functions_039 = exp(0.25846471956577927*v);
      real _expensive_functions_040 = exp(-0.049067713444553483*v);
      real txrf = 12.98 + 1.0/(0.0001020239312894894*_expensive_functions_039 + 0.00042992960891929087*_expensive_functions_040);
      real _expensive_functions_041 = exp(0.13596193065941536*v);
      real _expensive_functions_042 = exp(-0.038550501156515031*v);
      real txrs = 1.865 + 1.0/(0.00059224200368093944*_expensive_functions_041 + 3.5499661118024631e-5*_expensive_functions_042);
      real _expensive_functions_047 = exp(-0.11195700850873264*v);
      real xs1ss = 1.0/(0.27288596035656526*_expensive_functions_047 + 1.0);
      real _expensive_functions_048 = exp(0.056179775280898875*v);
      real _expensive_functions_049 = exp(-0.0043478260869565218*v);
      real txs1 = 817.29999999999995 + 1.0/(0.0035040677630748581*_expensive_functions_048 + 0.0005184809083581659*_expensive_functions_049);
      real xs2ss = xs1ss;
      real _expensive_functions_050 = exp(-0.032258064516129031*v);
      real _expensive_functions_051 = exp(0.050000000000000003*v);
      real txs2 = 1.0/(0.0022561357010639103*_expensive_functions_050 + 0.00082084998623898806*_expensive_functions_051);
      real _expensive_functions_053 = exp((-2.5537999999999998*ko - v - 144.59)/(1.5691999999999999*ko + 3.8115000000000001));
      real xk1ss = 1.0/(_expensive_functions_053 + 1.0);
      real _expensive_functions_054 = exp(-0.049115913555992145*v);
      real _expensive_functions_055 = exp(0.014423770373575654*v);
      real txk1 = 122.2/(0.0019352007631390235*_expensive_functions_054 + 30.433647575249029*_expensive_functions_055);
      real __melodee_temp_000 = V < -40;
      real alpha_h;
      real beta_h;
      if (__melodee_temp_000)
      {
         real _expensive_functions_065 = exp(-0.14705882352941177*V);
         alpha_h = 4.4312679295805147e-7*_expensive_functions_065;
         real _expensive_functions_066 = exp(0.34849999999999998*V);
         real _expensive_functions_067 = exp(0.079000000000000001*V);
         beta_h = 310000*_expensive_functions_066 + 2.7000000000000002*_expensive_functions_067;
      }
      else
      {
         alpha_h = 0;
         real _expensive_functions_065 = exp(-0.0900900900900901*V);
         beta_h = 0.77000000000000002/(0.049758141083938695*_expensive_functions_065 + 0.13);
      }
      real _expensive_functions_065 = exp(0.13458950201884254*V);
      real h_inf = 4.3210917837689708e-9/((_expensive_functions_065 + 6.5735011856460265e-5)*(_expensive_functions_065 + 6.5735011856460265e-5));
      real tau_h = (1.0/(alpha_h + beta_h));
      real __melodee_temp_001 = V < -40;
      real alpha_j;
      real beta_j;
      if (__melodee_temp_001)
      {
         real _expensive_functions_066 = exp(0.311*V);
         real _expensive_functions_067 = exp(0.24440000000000001*V);
         real _expensive_functions_068 = exp(-0.043909999999999998*V);
         alpha_j = (-25428*_expensive_functions_067 - 6.9480000000000002e-6*_expensive_functions_068)*(V + 37.780000000000001)/(50262745825.953987*_expensive_functions_066 + 1);
         real _expensive_functions_069 = exp(-0.13780000000000001*V);
         real _expensive_functions_070 = exp(-0.01052*V);
         beta_j = 0.024240000000000001*_expensive_functions_070/(0.003960868339904256*_expensive_functions_069 + 1);
      }
      else
      {
         alpha_j = 0;
         real _expensive_functions_066 = exp(-0.10000000000000001*V);
         real _expensive_functions_067 = exp(0.057000000000000002*V);
         beta_j = 0.59999999999999998*_expensive_functions_067/(0.040762203978366204*_expensive_functions_066 + 1);
      }
      real _expensive_functions_066 = exp(0.13458950201884254*V);
      real j_inf = 4.3210917837689708e-9/((_expensive_functions_066 + 6.5735011856460265e-5)*(_expensive_functions_066 + 6.5735011856460265e-5));
      real tau_j = (1.0/(alpha_j + beta_j));
      real _expensive_functions_067 = exp(-1.0/5.0*V - 12);
      real alpha_m = (1.0/(_expensive_functions_067 + 1));
      real _expensive_functions_068 = exp((1.0/5.0)*V + 7);
      real _expensive_functions_069 = exp((1.0/200.0)*V - 1.0/4.0);
      real beta_m = 0.10000000000000001/(_expensive_functions_069 + 1) + 0.10000000000000001/(_expensive_functions_068 + 1);
      real _expensive_functions_070 = exp(-0.11074197120708749*V);
      real m_inf = (1.0/(0.0018422115811651339*_expensive_functions_070 + 1)/(0.0018422115811651339*_expensive_functions_070 + 1));
      real tm = alpha_m*beta_m;
      real tmL = tm;
      real _expensive_functions_071 = exp(-_dt/ta);
      real _a_RLA = _expensive_functions_071 - 1;
      real _a_RLB = -ass;
      real _expensive_functions_072 = exp(-_dt/ta);
      real _ap_RLA = _expensive_functions_072 - 1;
      real _ap_RLB = -apss;
      real _expensive_functions_073 = exp(-_dt/td);
      real _d_RLA = _expensive_functions_073 - 1;
      real _d_RLB = -dss;
      real _expensive_functions_074 = exp(-_dt/tfcaf);
      real _fcaf_RLA = _expensive_functions_074 - 1;
      real _fcaf_RLB = -fcass;
      real _expensive_functions_075 = exp(-_dt/tfcafp);
      real _fcafp_RLA = _expensive_functions_075 - 1;
      real _fcafp_RLB = -fcass;
      real _expensive_functions_076 = exp(-_dt/tfcas);
      real _fcas_RLA = _expensive_functions_076 - 1;
      real _fcas_RLB = -fcass;
      real _expensive_functions_077 = exp(-_dt/tff);
      real _ff_RLA = _expensive_functions_077 - 1;
      real _ff_RLB = -fss;
      real _expensive_functions_078 = exp(-_dt/tffp);
      real _ffp_RLA = _expensive_functions_078 - 1;
      real _ffp_RLB = -fss;
      real _expensive_functions_079 = exp(-_dt/tfs);
      real _fs_RLA = _expensive_functions_079 - 1;
      real _fs_RLB = -fss;
      real _expensive_functions_080 = exp(-_dt/tau_h);
      real _h_RLA = _expensive_functions_080 - 1;
      real _h_RLB = -h_inf;
      real _hL_RLB = -hLss;
      real _hLp_RLB = -hLpss;
      real _expensive_functions_083 = exp(-_dt/ORd_tiF);
      real _iF_RLA = _expensive_functions_083 - 1;
      real _iF_RLB = -iss;
      real _expensive_functions_084 = exp(-_dt/tiFp);
      real _iFp_RLA = _expensive_functions_084 - 1;
      real _iFp_RLB = -iss;
      real _expensive_functions_085 = exp(-_dt/ORd_tiS);
      real _iS_RLA = _expensive_functions_085 - 1;
      real _iS_RLB = -iss;
      real _expensive_functions_086 = exp(-_dt/tiSp);
      real _iSp_RLA = _expensive_functions_086 - 1;
      real _iSp_RLB = -iss;
      real _expensive_functions_087 = exp(-_dt/tau_j);
      real _j_RLA = _expensive_functions_087 - 1;
      real _j_RLB = -j_inf;
      real _jca_RLB = -fcass;
      real _expensive_functions_089 = exp(-_dt/tm);
      real _m_RLA = _expensive_functions_089 - 1;
      real _m_RLB = -m_inf;
      real _expensive_functions_090 = exp(-_dt/tmL);
      real _mL_RLA = _expensive_functions_090 - 1;
      real _mL_RLB = -mLss;
      real _expensive_functions_091 = exp(-_dt*km2n);
      real _nca_RLA = _expensive_functions_091 - 1;
      real _nca_RLB = -anca*k2n/km2n;
      real _expensive_functions_092 = exp(-_dt/txk1);
      real _xk1_RLA = _expensive_functions_092 - 1;
      real _xk1_RLB = -xk1ss;
      real _expensive_functions_093 = exp(-_dt/txrf);
      real _xrf_RLA = _expensive_functions_093 - 1;
      real _xrf_RLB = -xrss;
      real _expensive_functions_094 = exp(-_dt/txrs);
      real _xrs_RLA = _expensive_functions_094 - 1;
      real _xrs_RLB = -xrss;
      real _expensive_functions_095 = exp(-_dt/txs1);
      real _xs1_RLA = _expensive_functions_095 - 1;
      real _xs1_RLB = -xs1ss;
      real _expensive_functions_096 = exp(-_dt/txs2);
      real _xs2_RLA = _expensive_functions_096 - 1;
      real _xs2_RLB = -xs2ss;
      //get the other differential updates
      real _expensive_functions = log(nao/nai);
      real ENa = R*T*_expensive_functions/F;
      real _expensive_functions_001 = log(ko/ki);
      real EK = R*T*_expensive_functions_001/F;
      real _expensive_functions_002 = log((PKNa*nao + ko)/(PKNa*nai + ki));
      real EKs = R*T*_expensive_functions_002/F;
      real vffrt = (F*F)*v/(R*T);
      real vfrt = F*v/(R*T);
      real CaMKb = CaMKo*(1.0 - CaMKt)/(KmCaM/cass + 1.0);
      real CaMKa = CaMKb + CaMKt;
      real CaMKt_diff = CaMKb*aCaMK*(CaMKb + CaMKt) - CaMKt*bCaMK;
      real fINaLp = 1.0/(1.0 + KmCaMK/CaMKa);
      real INaL = ORd_GNaL*mL*(-ENa + v)*(fINaLp*hLp + hL*(1.0 - fINaLp));
      real _expensive_functions_014 = exp(0.0066137566137566143*v);
      real AiF = 1.0/(0.24348537187522867*_expensive_functions_014 + 1.0);
      real AiS = 1.0 - AiF;
      real i = AiF*iF + AiS*iS;
      real ip = AiF*iFp + AiS*iSp;
      real fItop = 1.0/(1.0 + KmCaMK/CaMKa);
      real Ito = ORd_Gto*(-EK + v)*(a*i*(1.0 - fItop) + ap*fItop*ip);
      real f = Aff*ff + Afs*fs;
      real _expensive_functions_031 = exp(0.10000000000000001*v);
      real Afcaf = 0.29999999999999999 + 0.59999999999999998/(0.36787944117144233*_expensive_functions_031 + 1.0);
      real Afcas = 1.0 - Afcaf;
      real fca = Afcaf*fcaf + Afcas*fcas;
      real fp = Aff*ffp + Afs*fs;
      real fcap = Afcaf*fcafp + Afcas*fcas;
      real _expensive_functions_032 = exp(2.0*vfrt);
      real _expensive_functions_033 = exp(2.0*vfrt);
      real PhiCaL = 4.0*vffrt*(_expensive_functions_033*cass - 0.34100000000000003*cao)/(_expensive_functions_032 - 1.0);
      real _expensive_functions_034 = exp(1.0*vfrt);
      real _expensive_functions_035 = exp(1.0*vfrt);
      real PhiCaNa = 1.0*vffrt*(0.75*_expensive_functions_035*nass - 0.75*nao)/(_expensive_functions_034 - 1.0);
      real _expensive_functions_036 = exp(1.0*vfrt);
      real _expensive_functions_037 = exp(1.0*vfrt);
      real PhiCaK = 1.0*vffrt*(0.75*_expensive_functions_037*kss - 0.75*ko)/(_expensive_functions_036 - 1.0);
      real fICaLp = 1.0/(1.0 + KmCaMK/CaMKa);
      real ICaL = ORd_PCa*PhiCaL*d*(1.0 - fICaLp)*(f*(1.0 - nca) + fca*jca*nca) + PCap*PhiCaL*d*fICaLp*(fcap*jca*nca + fp*(1.0 - nca));
      real ICaNa = PCaNa*PhiCaNa*d*(1.0 - fICaLp)*(f*(1.0 - nca) + fca*jca*nca) + PCaNap*PhiCaNa*d*fICaLp*(fcap*jca*nca + fp*(1.0 - nca));
      real ICaK = PCaK*PhiCaK*d*(1.0 - fICaLp)*(f*(1.0 - nca) + fca*jca*nca) + PCaKp*PhiCaK*d*fICaLp*(fcap*jca*nca + fp*(1.0 - nca));
      real _expensive_functions_043 = exp(0.026171159382360639*v);
      real Axrf = 1.0/(4.197299094734718*_expensive_functions_043 + 1.0);
      real Axrs = 1.0 - Axrf;
      real xr = Axrf*xrf + Axrs*xrs;
      real _expensive_functions_044 = exp(0.013333333333333334*v);
      real _expensive_functions_045 = exp(0.033333333333333333*v);
      real rkr = 1.0/((2.0820090840784555*_expensive_functions_044 + 1.0)*(0.71653131057378927*_expensive_functions_045 + 1.0));
      real IKr = 0.43033148291193518*ORd_GKr*_expensive_functions_046*rkr*xr*(-EK + v);
      real _expensive_functions_052 = pow((1.0/cai), 1.3999999999999999);
      real KsCa = 1.0 + 0.59999999999999998/(6.4818210260626455e-7*_expensive_functions_052 + 1.0);
      real IKs = KsCa*ORd_GKs*xs1*xs2*(-EKs + v);
      real _expensive_functions_056 = exp(-0.27388602127883704*ko + 0.10534077741493732*v);
      real rk1 = 1.0/(69220.632210676704*_expensive_functions_056 + 1.0);
      real IK1 = ORd_GK1*_expensive_functions_057*rk1*xk1*(-EK + v);
      real hca = exp(F*qca*v/(R*T));
      real hna = exp(F*qna*v/(R*T));
      real h1 = 1 + nai*(hna + 1)/kna3;
      real h2 = hna*nai/(h1*kna3);
      real h3 = 1.0/h1;
      real h4 = 1.0 + nai*(1 + nai/kna2)/kna1;
      real h5 = (nai*nai)/(h4*kna1*kna2);
      real h6 = 1.0/h4;
      real h7 = 1.0 + nao*(1.0 + 1.0/hna)/kna3;
      real h8 = nao/(h7*hna*kna3);
      real h9 = 1.0/h7;
      real k3p = h9*wca;
      real k3pp = h8*wnaca;
      real k3 = k3p + k3pp;
      real k4p = h3*wca/hca;
      real k4pp = h2*wnaca;
      real k4 = k4p + k4pp;
      real k6 = cai*h6*kcaon;
      real k7 = h2*h5*wna;
      real k8 = h11*h8*wna;
      real x1 = k2*k4*(k6 + k7) + k5*k7*(k2 + k3);
      real x2 = k1*k7*(k4 + k5) + k4*k6*(k1 + k8);
      real x3 = k1*k3*(k6 + k7) + k6*k8*(k2 + k3);
      real x4 = k2*k8*(k4 + k5) + k3*k5*(k1 + k8);
      real E1 = x1/(x1 + x2 + x3 + x4);
      real E2 = x2/(x1 + x2 + x3 + x4);
      real E3 = x3/(x1 + x2 + x3 + x4);
      real E4 = x4/(x1 + x2 + x3 + x4);
      real allo = 1.0/((KmCaAct*KmCaAct)/(cai*cai) + 1.0);
      real JncxNa = -3.0*E1*k8 - E2*k3pp + E3*k4pp + 3.0*E4*k7;
      real JncxCa = -E1*k1 + E2*k2;
      real INaCa_i = 0.80000000000000004*ORd_Gncx*allo*(JncxCa*zca + JncxNa*zna);
      real ORd_h1 = 1 + nass*(hna + 1)/kna3;
      real ORd_h2 = hna*nass/(ORd_h1*kna3);
      real ORd_h3 = 1.0/ORd_h1;
      real ORd_h4 = 1.0 + nass*(1 + nass/kna2)/kna1;
      real ORd_h5 = (nass*nass)/(ORd_h4*kna1*kna2);
      real ORd_h6 = 1.0/ORd_h4;
      real ORd_h7 = 1.0 + nao*(1.0 + 1.0/hna)/kna3;
      real ORd_h8 = nao/(ORd_h7*hna*kna3);
      real ORd_h9 = 1.0/ORd_h7;
      real ORd_k3p = ORd_h9*wca;
      real ORd_k3pp = ORd_h8*wnaca;
      real ORd_k3 = ORd_k3p + ORd_k3pp;
      real ORd_k4p = ORd_h3*wca/hca;
      real ORd_k4pp = ORd_h2*wnaca;
      real ORd_k4 = ORd_k4p + ORd_k4pp;
      real ORd_k6 = ORd_h6*cass*kcaon;
      real ORd_k7 = ORd_h2*ORd_h5*wna;
      real ORd_k8 = ORd_h11*ORd_h8*wna;
      real ORd_x1 = ORd_k2*ORd_k4*(ORd_k6 + ORd_k7) + ORd_k5*ORd_k7*(ORd_k2 + ORd_k3);
      real ORd_x2 = ORd_k1*ORd_k7*(ORd_k4 + ORd_k5) + ORd_k4*ORd_k6*(ORd_k1 + ORd_k8);
      real ORd_x3 = ORd_k1*ORd_k3*(ORd_k6 + ORd_k7) + ORd_k6*ORd_k8*(ORd_k2 + ORd_k3);
      real ORd_x4 = ORd_k2*ORd_k8*(ORd_k4 + ORd_k5) + ORd_k3*ORd_k5*(ORd_k1 + ORd_k8);
      real ORd_E1 = ORd_x1/(ORd_x1 + ORd_x2 + ORd_x3 + ORd_x4);
      real ORd_E2 = ORd_x2/(ORd_x1 + ORd_x2 + ORd_x3 + ORd_x4);
      real ORd_E3 = ORd_x3/(ORd_x1 + ORd_x2 + ORd_x3 + ORd_x4);
      real ORd_E4 = ORd_x4/(ORd_x1 + ORd_x2 + ORd_x3 + ORd_x4);
      real ORd_allo = 1.0/((ORd_KmCaAct*ORd_KmCaAct)/(cass*cass) + 1.0);
      real ORd_JncxNa = -3.0*ORd_E1*ORd_k8 - ORd_E2*ORd_k3pp + ORd_E3*ORd_k4pp + 3.0*ORd_E4*ORd_k7;
      real ORd_JncxCa = -ORd_E1*ORd_k1 + ORd_E2*ORd_k2;
      real INaCa_ss = 0.20000000000000001*ORd_Gncx*ORd_allo*(ORd_JncxCa*zca + ORd_JncxNa*zna);
      real _expensive_functions_058 = exp(0.33333333333333331*F*delta*v/(R*T));
      real Knai = Knai0*_expensive_functions_058;
      real _expensive_functions_059 = exp(0.33333333333333331*F*v*(1.0 - delta)/(R*T));
      real Knao = Knao0*_expensive_functions_059;
      real P = eP/(H/Khp + 1.0 + ki/Kxkur + nai/Knap);
      real a1 = k1p*((nai/Knai)*(nai/Knai)*(nai/Knai))/(((1.0 + ki/Kki)*(1.0 + ki/Kki)) + ((1.0 + nai/Knai)*(1.0 + nai/Knai)*(1.0 + nai/Knai)) - 1.0);
      real b2 = k2m*((nao/Knao)*(nao/Knao)*(nao/Knao))/(((1.0 + ko/Kko)*(1.0 + ko/Kko)) + ((1.0 + nao/Knao)*(1.0 + nao/Knao)*(1.0 + nao/Knao)) - 1.0);
      real a3 = ORd_k3p_001*((ko/Kko)*(ko/Kko))/(((1.0 + ko/Kko)*(1.0 + ko/Kko)) + ((1.0 + nao/Knao)*(1.0 + nao/Knao)*(1.0 + nao/Knao)) - 1.0);
      real b3 = H*P*k3m/(1.0 + MgATP/Kmgatp);
      real b4 = k4m*((ki/Kki)*(ki/Kki))/(((1.0 + ki/Kki)*(1.0 + ki/Kki)) + ((1.0 + nai/Knai)*(1.0 + nai/Knai)*(1.0 + nai/Knai)) - 1.0);
      real ORd_x1_001 = a1*a2*a4 + a1*a2*b3 + a2*b3*b4 + b2*b3*b4;
      real ORd_x2_001 = a1*a2*a3 + a2*a3*b4 + a3*b1*b4 + b1*b2*b4;
      real ORd_x3_001 = a2*a3*a4 + a3*a4*b1 + a4*b1*b2 + b1*b2*b3;
      real ORd_x4_001 = a1*a3*a4 + a1*a4*b2 + a1*b2*b3 + b2*b3*b4;
      real ORd_E1_001 = ORd_x1_001/(ORd_x1_001 + ORd_x2_001 + ORd_x3_001 + ORd_x4_001);
      real ORd_E2_001 = ORd_x2_001/(ORd_x1_001 + ORd_x2_001 + ORd_x3_001 + ORd_x4_001);
      real ORd_E3_001 = ORd_x3_001/(ORd_x1_001 + ORd_x2_001 + ORd_x3_001 + ORd_x4_001);
      real ORd_E4_001 = ORd_x4_001/(ORd_x1_001 + ORd_x2_001 + ORd_x3_001 + ORd_x4_001);
      real JnakNa = 3.0*ORd_E1_001*a3 - 3.0*ORd_E2_001*b3;
      real JnakK = -2.0*ORd_E3_001*a1 + 2.0*ORd_E4_001*b1;
      real INaK = ORd_Pnak*(JnakK*zk + JnakNa*zna);
      real _expensive_functions_060 = exp(-0.054525627044711013*v);
      real xkb = 1.0/(2.2023634509492389*_expensive_functions_060 + 1.0);
      real IKb = ORd_GKb*xkb*(-EK + v);
      real _expensive_functions_061 = exp(vfrt);
      real _expensive_functions_062 = exp(vfrt);
      real INab = PNab*vffrt*(_expensive_functions_062*nai - nao)/(_expensive_functions_061 - 1.0);
      real _expensive_functions_063 = exp(2.0*vfrt);
      real _expensive_functions_064 = exp(2.0*vfrt);
      real ICab = 4.0*PCab*vffrt*(_expensive_functions_064*cai - 0.34100000000000003*cao)/(_expensive_functions_063 - 1.0);
      real IpCa = GpCa*cai/(cai + 0.00050000000000000001);
      real JdiffNa = -0.5*nai + 0.5*nass;
      real JdiffK = -0.5*ki + 0.5*kss;
      real Jdiff = -5.0*cai + 5.0*cass;
      real Jrel_inf = -ICaL*a_rel/(25.62890625*(((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))) + 1.0);
      real ORd_Jrel_inf;
      if (__melodee_temp_022)
      {
         ORd_Jrel_inf = 1.7*Jrel_inf;
      }
      else
      {
         ORd_Jrel_inf = Jrel_inf;
      }
      real tau_rel = bt/(1.0 + 0.0123/cajsr);
      real __melodee_temp_023 = tau_rel < 0.001;
      real ORd_tau_rel;
      if (__melodee_temp_023)
      {
         ORd_tau_rel = 0.001;
      }
      else
      {
         ORd_tau_rel = tau_rel;
      }
      real Jrelnp_diff = (-Jrelnp + ORd_Jrel_inf)/ORd_tau_rel;
      real Jrelp_inf = -ICaL*a_relp/(25.62890625*(((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))*((1.0/cajsr))) + 1.0);
      real ORd_Jrelp_inf;
      if (__melodee_temp_024)
      {
         ORd_Jrelp_inf = 1.7*Jrelp_inf;
      }
      else
      {
         ORd_Jrelp_inf = Jrelp_inf;
      }
      real tau_relp = btp/(1.0 + 0.0123/cajsr);
      real __melodee_temp_025 = tau_relp < 0.001;
      real ORd_tau_relp;
      if (__melodee_temp_025)
      {
         ORd_tau_relp = 0.001;
      }
      else
      {
         ORd_tau_relp = tau_relp;
      }
      real Jrelp_diff = (-Jrelp + ORd_Jrelp_inf)/ORd_tau_relp;
      real fJrelp = 1.0/(1.0 + KmCaMK/CaMKa);
      real Jrel = Jrelnp*(1.0 - fJrelp) + Jrelp*fJrelp;
      real __melodee_temp_026 = Jrel*JrelStiffConst > cajsr;
      real ORd_Jrel;
      if (__melodee_temp_026)
      {
         ORd_Jrel = cajsr/JrelStiffConst;
      }
      else
      {
         ORd_Jrel = Jrel;
      }
      real Jupnp = 0.0043750000000000004*cai/(cai + 0.00092000000000000003);
      real Jupp = 0.01203125*cai/(cai + 0.00075000000000000002);
      real ORd_Jupnp;
      real ORd_Jupp;
      if (__melodee_temp_027)
      {
         ORd_Jupnp = 1.3*Jupnp;
         ORd_Jupp = 1.3*Jupp;
      }
      else
      {
         ORd_Jupnp = Jupnp;
         ORd_Jupp = Jupp;
      }
      real fJupp = 1.0/(1.0 + KmCaMK/CaMKa);
      real Jleak = 0.00026249999999999998*cansr;
      real Jup = -Jleak + ORd_Jupnp*(1.0 - fJupp) + ORd_Jupp*fJupp;
      real Jtr = -0.01*cajsr + 0.01*cansr;
      real nass_diff = Acap*(-ICaNa - 3.0*INaCa_ss)/(F*vss) - JdiffNa;
      real ki_diff = Acap*(-IK1 - IKb - IKr - IKs + 2.0*INaK - Ito)/(F*vmyo) + JdiffK*vss/vmyo;
      real kss_diff = -Acap*ICaK/(F*vss) - JdiffK;
      real Bcai = 1.0/(ORd_cmdnmax*kmcmdn*(1.0/(cai + kmcmdn)/(cai + kmcmdn)) + kmtrpn*trpnmax*(1.0/(cai + kmtrpn)/(cai + kmtrpn)) + 1.0);
      real cai_diff = Bcai*(0.5*Acap*(-ICab + 2.0*INaCa_i - IpCa)/(F*vmyo) + Jdiff*vss/vmyo - Jup*vnsr/vmyo);
      real Bcass = 1.0/(BSLmax*KmBSL*(1.0/(KmBSL + cass)/(KmBSL + cass)) + BSRmax*KmBSR*(1.0/(KmBSR + cass)/(KmBSR + cass)) + 1.0);
      real cass_diff = Bcass*(0.5*Acap*(-ICaL + 2.0*INaCa_ss)/(F*vss) - Jdiff + ORd_Jrel*vjsr/vss);
      real cansr_diff = -Jtr*vjsr/vnsr + Jup;
      real Bcajsr = 1.0/(csqnmax*kmcsqn*(1.0/(cajsr + kmcsqn)/(cajsr + kmcsqn)) + 1.0);
      real cajsr_diff = Bcajsr*(Jtr - ORd_Jrel);
      real i_Na = (m*m*m)*g_Na*h*j*(-ENa + V);
      real i_Naitot = i_Na;
      real INa = i_Naitot;
      real nai_diff = Acap*(-INa - 3.0*INaCa_i - 3.0*INaK - INaL - INab)/(F*vmyo) + JdiffNa*vss/vmyo;
      //get Iion
      real Iion = ICaK + ICaL + ICaNa + ICab + IK1 + IKb + IKr + IKs + INa + INaCa_i + INaCa_ss + INaK + INaL + INab + IpCa + Ito;
      real Iion_001 = Iion;
      //Do the markov update (1 step rosenbrock with gauss siedel)
      //EDIT_STATE
      CaMKt += _dt*CaMKt_diff;
      Jrelnp += _dt*Jrelnp_diff;
      Jrelp += _dt*Jrelp_diff;
      cai += _dt*cai_diff;
      cajsr += _dt*cajsr_diff;
      cansr += _dt*cansr_diff;
      cass += _dt*cass_diff;
      ki += _dt*ki_diff;
      kss += _dt*kss_diff;
      nai += _dt*nai_diff;
      nass += _dt*nass_diff;
      a += _a_RLA*(a+_a_RLB);
      ap += _ap_RLA*(ap+_ap_RLB);
      d += _d_RLA*(d+_d_RLB);
      fcaf += _fcaf_RLA*(fcaf+_fcaf_RLB);
      fcafp += _fcafp_RLA*(fcafp+_fcafp_RLB);
      fcas += _fcas_RLA*(fcas+_fcas_RLB);
      ff += _ff_RLA*(ff+_ff_RLB);
      ffp += _ffp_RLA*(ffp+_ffp_RLB);
      fs += _fs_RLA*(fs+_fs_RLB);
      h += _h_RLA*(h+_h_RLB);
      hL += _hL_RLA*(hL+_hL_RLB);
      hLp += _hLp_RLA*(hLp+_hLp_RLB);
      iF += _iF_RLA*(iF+_iF_RLB);
      iFp += _iFp_RLA*(iFp+_iFp_RLB);
      iS += _iS_RLA*(iS+_iS_RLB);
      iSp += _iSp_RLA*(iSp+_iSp_RLB);
      j += _j_RLA*(j+_j_RLB);
      jca += _jca_RLA*(jca+_jca_RLB);
      m += _m_RLA*(m+_m_RLB);
      mL += _mL_RLA*(mL+_mL_RLB);
      nca += _nca_RLA*(nca+_nca_RLB);
      xk1 += _xk1_RLA*(xk1+_xk1_RLB);
      xrf += _xrf_RLA*(xrf+_xrf_RLB);
      xrs += _xrs_RLA*(xrs+_xrs_RLB);
      xs1 += _xs1_RLA*(xs1+_xs1_RLB);
      xs2 += _xs2_RLA*(xs2+_xs2_RLB);
      store(state_[__jj].CaMKt, CaMKt);
      store(state_[__jj].Jrelnp, Jrelnp);
      store(state_[__jj].Jrelp, Jrelp);
      store(state_[__jj].a, a);
      store(state_[__jj].ap, ap);
      store(state_[__jj].cai, cai);
      store(state_[__jj].cajsr, cajsr);
      store(state_[__jj].cansr, cansr);
      store(state_[__jj].cass, cass);
      store(state_[__jj].d, d);
      store(state_[__jj].fcaf, fcaf);
      store(state_[__jj].fcafp, fcafp);
      store(state_[__jj].fcas, fcas);
      store(state_[__jj].ff, ff);
      store(state_[__jj].ffp, ffp);
      store(state_[__jj].fs, fs);
      store(state_[__jj].h, h);
      store(state_[__jj].hL, hL);
      store(state_[__jj].hLp, hLp);
      store(state_[__jj].iF, iF);
      store(state_[__jj].iFp, iFp);
      store(state_[__jj].iS, iS);
      store(state_[__jj].iSp, iSp);
      store(state_[__jj].j, j);
      store(state_[__jj].jca, jca);
      store(state_[__jj].ki, ki);
      store(state_[__jj].kss, kss);
      store(state_[__jj].m, m);
      store(state_[__jj].mL, mL);
      store(state_[__jj].nai, nai);
      store(state_[__jj].nass, nass);
      store(state_[__jj].nca, nca);
      store(state_[__jj].xk1, xk1);
      store(state_[__jj].xrf, xrf);
      store(state_[__jj].xrs, xrs);
      store(state_[__jj].xs1, xs1);
      store(state_[__jj].xs2, xs2);
      double __dVm_local[width];
      simdops::store(&__dVm_local[0],-Iion_001);
      {
         int cursor = 0;
         for (int __kk=0; __kk<width && __ii+__kk<nCells_; __kk++)
         {
            __dVm[__indexArray[__ii+__kk]] = __dVm_local[__kk];
         }
      }
   }
}
#endif //USE_CUDA
   
string ThisReaction::methodName() const
{
   return "ord_SCC";
}

void ThisReaction::initializeMembraneVoltage(ro_mgarray_ptr<int> __indexArray_m, wo_mgarray_ptr<double> __Vm_m)
{
   assert(__Vm_m.size() >= nCells_);

   wo_array_ptr<double> __Vm = __Vm_m.useOn(CPU);
   ro_array_ptr<int> __indexArray = __indexArray_m.useOn(CPU);
#ifdef USE_CUDA
#define READ_STATE(state,index) (stateData[_##state##_off*nCells_+index])
   wo_array_ptr<double> stateData = stateTransport_.useOn(CPU);
#else //USE_CUDA
#define READ_STATE(state,index) (state_[index/width].state[index % width])
   state_.resize((nCells_+width-1)/width);
#endif //USE_CUDA


   double V_init = -85;
   double V = V_init;
   double nai_init = 7;
   double nai = nai_init;
   double nass_init = 7;
   double nass = nass_init;
   double ki_init = 145;
   double ki = ki_init;
   double kss_init = 145;
   double kss = kss_init;
   double cai_init = 0.0001;
   double cai = cai_init;
   double cass_init = 0.0001;
   double cass = cass_init;
   double cansr_init = 1.2;
   double cansr = cansr_init;
   double cajsr_init = 1.2;
   double cajsr = cajsr_init;
   double Jrelnp_init = 0;
   double Jrelnp = Jrelnp_init;
   double Jrelp_init = 0;
   double Jrelp = Jrelp_init;
   double CaMKt_init = 0;
   double CaMKt = CaMKt_init;
   double mL_init = 0;
   double mL = mL_init;
   double hL_init = 1;
   double hL = hL_init;
   double hLp_init = 1;
   double hLp = hLp_init;
   double a_init = 0;
   double a = a_init;
   double iF_init = 1;
   double iF = iF_init;
   double iS_init = 1;
   double iS = iS_init;
   double ap_init = 0;
   double ap = ap_init;
   double iFp_init = 1;
   double iFp = iFp_init;
   double iSp_init = 1;
   double iSp = iSp_init;
   double d_init = 0;
   double d = d_init;
   double ff_init = 1;
   double ff = ff_init;
   double fs_init = 1;
   double fs = fs_init;
   double fcaf_init = 1;
   double fcaf = fcaf_init;
   double fcas_init = 1;
   double fcas = fcas_init;
   double jca_init = 1;
   double jca = jca_init;
   double nca_init = 0;
   double nca = nca_init;
   double ffp_init = 1;
   double ffp = ffp_init;
   double fcafp_init = 1;
   double fcafp = fcafp_init;
   double xrf_init = 0;
   double xrf = xrf_init;
   double xrs_init = 0;
   double xrs = xrs_init;
   double xs1_init = 0;
   double xs1 = xs1_init;
   double xs2_init = 0;
   double xs2 = xs2_init;
   double xk1_init = 1;
   double xk1 = xk1_init;
   double h_init = 0.55000000000000004;
   double h = h_init;
   double j_init = 0.66000000000000003;
   double j = j_init;
   double m_init = 0.0015499999999999999;
   double m = m_init;
   for (int iCell=0; iCell<nCells_; iCell++)
   {
      READ_STATE(CaMKt,iCell) = CaMKt;
      READ_STATE(Jrelnp,iCell) = Jrelnp;
      READ_STATE(Jrelp,iCell) = Jrelp;
      READ_STATE(a,iCell) = a;
      READ_STATE(ap,iCell) = ap;
      READ_STATE(cai,iCell) = cai;
      READ_STATE(cajsr,iCell) = cajsr;
      READ_STATE(cansr,iCell) = cansr;
      READ_STATE(cass,iCell) = cass;
      READ_STATE(d,iCell) = d;
      READ_STATE(fcaf,iCell) = fcaf;
      READ_STATE(fcafp,iCell) = fcafp;
      READ_STATE(fcas,iCell) = fcas;
      READ_STATE(ff,iCell) = ff;
      READ_STATE(ffp,iCell) = ffp;
      READ_STATE(fs,iCell) = fs;
      READ_STATE(h,iCell) = h;
      READ_STATE(hL,iCell) = hL;
      READ_STATE(hLp,iCell) = hLp;
      READ_STATE(iF,iCell) = iF;
      READ_STATE(iFp,iCell) = iFp;
      READ_STATE(iS,iCell) = iS;
      READ_STATE(iSp,iCell) = iSp;
      READ_STATE(j,iCell) = j;
      READ_STATE(jca,iCell) = jca;
      READ_STATE(ki,iCell) = ki;
      READ_STATE(kss,iCell) = kss;
      READ_STATE(m,iCell) = m;
      READ_STATE(mL,iCell) = mL;
      READ_STATE(nai,iCell) = nai;
      READ_STATE(nass,iCell) = nass;
      READ_STATE(nca,iCell) = nca;
      READ_STATE(xk1,iCell) = xk1;
      READ_STATE(xrf,iCell) = xrf;
      READ_STATE(xrs,iCell) = xrs;
      READ_STATE(xs1,iCell) = xs1;
      READ_STATE(xs2,iCell) = xs2;
      __Vm[__indexArray[iCell]] = V_init;
   }
}

enum varHandles
{
   CaMKt_handle,
   Jrelnp_handle,
   Jrelp_handle,
   a_handle,
   ap_handle,
   cai_handle,
   cajsr_handle,
   cansr_handle,
   cass_handle,
   d_handle,
   fcaf_handle,
   fcafp_handle,
   fcas_handle,
   ff_handle,
   ffp_handle,
   fs_handle,
   h_handle,
   hL_handle,
   hLp_handle,
   iF_handle,
   iFp_handle,
   iS_handle,
   iSp_handle,
   j_handle,
   jca_handle,
   ki_handle,
   kss_handle,
   m_handle,
   mL_handle,
   nai_handle,
   nass_handle,
   nca_handle,
   xk1_handle,
   xrf_handle,
   xrs_handle,
   xs1_handle,
   xs2_handle,
   i_Na_handle,
   NUMHANDLES
};

const string ThisReaction::getUnit(const std::string& varName) const
{
   if(0) {}
   else if (varName == "CaMKt") { return "mM"; }
   else if (varName == "Jrelnp") { return "M/s"; }
   else if (varName == "Jrelp") { return "M/s"; }
   else if (varName == "a") { return "1"; }
   else if (varName == "ap") { return "1"; }
   else if (varName == "cai") { return "mM"; }
   else if (varName == "cajsr") { return "mM"; }
   else if (varName == "cansr") { return "mM"; }
   else if (varName == "cass") { return "mM"; }
   else if (varName == "d") { return "1"; }
   else if (varName == "fcaf") { return "1"; }
   else if (varName == "fcafp") { return "1"; }
   else if (varName == "fcas") { return "1"; }
   else if (varName == "ff") { return "1"; }
   else if (varName == "ffp") { return "1"; }
   else if (varName == "fs") { return "1"; }
   else if (varName == "h") { return "1"; }
   else if (varName == "hL") { return "1"; }
   else if (varName == "hLp") { return "1"; }
   else if (varName == "iF") { return "1"; }
   else if (varName == "iFp") { return "1"; }
   else if (varName == "iS") { return "1"; }
   else if (varName == "iSp") { return "1"; }
   else if (varName == "i_Na") { return "INVALID"; }
   else if (varName == "j") { return "1"; }
   else if (varName == "jca") { return "1"; }
   else if (varName == "ki") { return "mM"; }
   else if (varName == "kss") { return "mM"; }
   else if (varName == "m") { return "1"; }
   else if (varName == "mL") { return "1"; }
   else if (varName == "nai") { return "mM"; }
   else if (varName == "nass") { return "mM"; }
   else if (varName == "nca") { return "1"; }
   else if (varName == "xk1") { return "1"; }
   else if (varName == "xrf") { return "1"; }
   else if (varName == "xrs") { return "1"; }
   else if (varName == "xs1") { return "1"; }
   else if (varName == "xs2") { return "1"; }
   return "INVALID";
}

int ThisReaction::getVarHandle(const std::string& varName) const
{
   if (0) {}
   else if (varName == "CaMKt") { return CaMKt_handle; }
   else if (varName == "Jrelnp") { return Jrelnp_handle; }
   else if (varName == "Jrelp") { return Jrelp_handle; }
   else if (varName == "a") { return a_handle; }
   else if (varName == "ap") { return ap_handle; }
   else if (varName == "cai") { return cai_handle; }
   else if (varName == "cajsr") { return cajsr_handle; }
   else if (varName == "cansr") { return cansr_handle; }
   else if (varName == "cass") { return cass_handle; }
   else if (varName == "d") { return d_handle; }
   else if (varName == "fcaf") { return fcaf_handle; }
   else if (varName == "fcafp") { return fcafp_handle; }
   else if (varName == "fcas") { return fcas_handle; }
   else if (varName == "ff") { return ff_handle; }
   else if (varName == "ffp") { return ffp_handle; }
   else if (varName == "fs") { return fs_handle; }
   else if (varName == "h") { return h_handle; }
   else if (varName == "hL") { return hL_handle; }
   else if (varName == "hLp") { return hLp_handle; }
   else if (varName == "iF") { return iF_handle; }
   else if (varName == "iFp") { return iFp_handle; }
   else if (varName == "iS") { return iS_handle; }
   else if (varName == "iSp") { return iSp_handle; }
   else if (varName == "i_Na") { return i_Na_handle; }
   else if (varName == "j") { return j_handle; }
   else if (varName == "jca") { return jca_handle; }
   else if (varName == "ki") { return ki_handle; }
   else if (varName == "kss") { return kss_handle; }
   else if (varName == "m") { return m_handle; }
   else if (varName == "mL") { return mL_handle; }
   else if (varName == "nai") { return nai_handle; }
   else if (varName == "nass") { return nass_handle; }
   else if (varName == "nca") { return nca_handle; }
   else if (varName == "xk1") { return xk1_handle; }
   else if (varName == "xrf") { return xrf_handle; }
   else if (varName == "xrs") { return xrs_handle; }
   else if (varName == "xs1") { return xs1_handle; }
   else if (varName == "xs2") { return xs2_handle; }
   return -1;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readwrite(CPU);
#endif //USE_CUDA



   if (0) {}
   else if (varHandle == CaMKt_handle) { READ_STATE(CaMKt,iCell) = value; }
   else if (varHandle == Jrelnp_handle) { READ_STATE(Jrelnp,iCell) = value; }
   else if (varHandle == Jrelp_handle) { READ_STATE(Jrelp,iCell) = value; }
   else if (varHandle == a_handle) { READ_STATE(a,iCell) = value; }
   else if (varHandle == ap_handle) { READ_STATE(ap,iCell) = value; }
   else if (varHandle == cai_handle) { READ_STATE(cai,iCell) = value; }
   else if (varHandle == cajsr_handle) { READ_STATE(cajsr,iCell) = value; }
   else if (varHandle == cansr_handle) { READ_STATE(cansr,iCell) = value; }
   else if (varHandle == cass_handle) { READ_STATE(cass,iCell) = value; }
   else if (varHandle == d_handle) { READ_STATE(d,iCell) = value; }
   else if (varHandle == fcaf_handle) { READ_STATE(fcaf,iCell) = value; }
   else if (varHandle == fcafp_handle) { READ_STATE(fcafp,iCell) = value; }
   else if (varHandle == fcas_handle) { READ_STATE(fcas,iCell) = value; }
   else if (varHandle == ff_handle) { READ_STATE(ff,iCell) = value; }
   else if (varHandle == ffp_handle) { READ_STATE(ffp,iCell) = value; }
   else if (varHandle == fs_handle) { READ_STATE(fs,iCell) = value; }
   else if (varHandle == h_handle) { READ_STATE(h,iCell) = value; }
   else if (varHandle == hL_handle) { READ_STATE(hL,iCell) = value; }
   else if (varHandle == hLp_handle) { READ_STATE(hLp,iCell) = value; }
   else if (varHandle == iF_handle) { READ_STATE(iF,iCell) = value; }
   else if (varHandle == iFp_handle) { READ_STATE(iFp,iCell) = value; }
   else if (varHandle == iS_handle) { READ_STATE(iS,iCell) = value; }
   else if (varHandle == iSp_handle) { READ_STATE(iSp,iCell) = value; }
   else if (varHandle == j_handle) { READ_STATE(j,iCell) = value; }
   else if (varHandle == jca_handle) { READ_STATE(jca,iCell) = value; }
   else if (varHandle == ki_handle) { READ_STATE(ki,iCell) = value; }
   else if (varHandle == kss_handle) { READ_STATE(kss,iCell) = value; }
   else if (varHandle == m_handle) { READ_STATE(m,iCell) = value; }
   else if (varHandle == mL_handle) { READ_STATE(mL,iCell) = value; }
   else if (varHandle == nai_handle) { READ_STATE(nai,iCell) = value; }
   else if (varHandle == nass_handle) { READ_STATE(nass,iCell) = value; }
   else if (varHandle == nca_handle) { READ_STATE(nca,iCell) = value; }
   else if (varHandle == xk1_handle) { READ_STATE(xk1,iCell) = value; }
   else if (varHandle == xrf_handle) { READ_STATE(xrf,iCell) = value; }
   else if (varHandle == xrs_handle) { READ_STATE(xrs,iCell) = value; }
   else if (varHandle == xs1_handle) { READ_STATE(xs1,iCell) = value; }
   else if (varHandle == xs2_handle) { READ_STATE(xs2,iCell) = value; }
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   if (0) {}
   else if (varHandle == CaMKt_handle) { return READ_STATE(CaMKt,iCell); }
   else if (varHandle == Jrelnp_handle) { return READ_STATE(Jrelnp,iCell); }
   else if (varHandle == Jrelp_handle) { return READ_STATE(Jrelp,iCell); }
   else if (varHandle == a_handle) { return READ_STATE(a,iCell); }
   else if (varHandle == ap_handle) { return READ_STATE(ap,iCell); }
   else if (varHandle == cai_handle) { return READ_STATE(cai,iCell); }
   else if (varHandle == cajsr_handle) { return READ_STATE(cajsr,iCell); }
   else if (varHandle == cansr_handle) { return READ_STATE(cansr,iCell); }
   else if (varHandle == cass_handle) { return READ_STATE(cass,iCell); }
   else if (varHandle == d_handle) { return READ_STATE(d,iCell); }
   else if (varHandle == fcaf_handle) { return READ_STATE(fcaf,iCell); }
   else if (varHandle == fcafp_handle) { return READ_STATE(fcafp,iCell); }
   else if (varHandle == fcas_handle) { return READ_STATE(fcas,iCell); }
   else if (varHandle == ff_handle) { return READ_STATE(ff,iCell); }
   else if (varHandle == ffp_handle) { return READ_STATE(ffp,iCell); }
   else if (varHandle == fs_handle) { return READ_STATE(fs,iCell); }
   else if (varHandle == h_handle) { return READ_STATE(h,iCell); }
   else if (varHandle == hL_handle) { return READ_STATE(hL,iCell); }
   else if (varHandle == hLp_handle) { return READ_STATE(hLp,iCell); }
   else if (varHandle == iF_handle) { return READ_STATE(iF,iCell); }
   else if (varHandle == iFp_handle) { return READ_STATE(iFp,iCell); }
   else if (varHandle == iS_handle) { return READ_STATE(iS,iCell); }
   else if (varHandle == iSp_handle) { return READ_STATE(iSp,iCell); }
   else if (varHandle == j_handle) { return READ_STATE(j,iCell); }
   else if (varHandle == jca_handle) { return READ_STATE(jca,iCell); }
   else if (varHandle == ki_handle) { return READ_STATE(ki,iCell); }
   else if (varHandle == kss_handle) { return READ_STATE(kss,iCell); }
   else if (varHandle == m_handle) { return READ_STATE(m,iCell); }
   else if (varHandle == mL_handle) { return READ_STATE(mL,iCell); }
   else if (varHandle == nai_handle) { return READ_STATE(nai,iCell); }
   else if (varHandle == nass_handle) { return READ_STATE(nass,iCell); }
   else if (varHandle == nca_handle) { return READ_STATE(nca,iCell); }
   else if (varHandle == xk1_handle) { return READ_STATE(xk1,iCell); }
   else if (varHandle == xrf_handle) { return READ_STATE(xrf,iCell); }
   else if (varHandle == xrs_handle) { return READ_STATE(xrs,iCell); }
   else if (varHandle == xs1_handle) { return READ_STATE(xs1,iCell); }
   else if (varHandle == xs2_handle) { return READ_STATE(xs2,iCell); }
   return NAN;
}

double ThisReaction::getValue(int iCell, int varHandle, double V) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   const double CaMKt=READ_STATE(CaMKt,iCell);
   const double Jrelnp=READ_STATE(Jrelnp,iCell);
   const double Jrelp=READ_STATE(Jrelp,iCell);
   const double a=READ_STATE(a,iCell);
   const double ap=READ_STATE(ap,iCell);
   const double cai=READ_STATE(cai,iCell);
   const double cajsr=READ_STATE(cajsr,iCell);
   const double cansr=READ_STATE(cansr,iCell);
   const double cass=READ_STATE(cass,iCell);
   const double d=READ_STATE(d,iCell);
   const double fcaf=READ_STATE(fcaf,iCell);
   const double fcafp=READ_STATE(fcafp,iCell);
   const double fcas=READ_STATE(fcas,iCell);
   const double ff=READ_STATE(ff,iCell);
   const double ffp=READ_STATE(ffp,iCell);
   const double fs=READ_STATE(fs,iCell);
   const double h=READ_STATE(h,iCell);
   const double hL=READ_STATE(hL,iCell);
   const double hLp=READ_STATE(hLp,iCell);
   const double iF=READ_STATE(iF,iCell);
   const double iFp=READ_STATE(iFp,iCell);
   const double iS=READ_STATE(iS,iCell);
   const double iSp=READ_STATE(iSp,iCell);
   const double j=READ_STATE(j,iCell);
   const double jca=READ_STATE(jca,iCell);
   const double ki=READ_STATE(ki,iCell);
   const double kss=READ_STATE(kss,iCell);
   const double m=READ_STATE(m,iCell);
   const double mL=READ_STATE(mL,iCell);
   const double nai=READ_STATE(nai,iCell);
   const double nass=READ_STATE(nass,iCell);
   const double nca=READ_STATE(nca,iCell);
   const double xk1=READ_STATE(xk1,iCell);
   const double xrf=READ_STATE(xrf,iCell);
   const double xrs=READ_STATE(xrs,iCell);
   const double xs1=READ_STATE(xs1,iCell);
   const double xs2=READ_STATE(xs2,iCell);
   if (0) {}
   else if (varHandle == CaMKt_handle)
   {
      return CaMKt;
   }
   else if (varHandle == Jrelnp_handle)
   {
      return Jrelnp;
   }
   else if (varHandle == Jrelp_handle)
   {
      return Jrelp;
   }
   else if (varHandle == a_handle)
   {
      return a;
   }
   else if (varHandle == ap_handle)
   {
      return ap;
   }
   else if (varHandle == cai_handle)
   {
      return cai;
   }
   else if (varHandle == cajsr_handle)
   {
      return cajsr;
   }
   else if (varHandle == cansr_handle)
   {
      return cansr;
   }
   else if (varHandle == cass_handle)
   {
      return cass;
   }
   else if (varHandle == d_handle)
   {
      return d;
   }
   else if (varHandle == fcaf_handle)
   {
      return fcaf;
   }
   else if (varHandle == fcafp_handle)
   {
      return fcafp;
   }
   else if (varHandle == fcas_handle)
   {
      return fcas;
   }
   else if (varHandle == ff_handle)
   {
      return ff;
   }
   else if (varHandle == ffp_handle)
   {
      return ffp;
   }
   else if (varHandle == fs_handle)
   {
      return fs;
   }
   else if (varHandle == h_handle)
   {
      return h;
   }
   else if (varHandle == hL_handle)
   {
      return hL;
   }
   else if (varHandle == hLp_handle)
   {
      return hLp;
   }
   else if (varHandle == iF_handle)
   {
      return iF;
   }
   else if (varHandle == iFp_handle)
   {
      return iFp;
   }
   else if (varHandle == iS_handle)
   {
      return iS;
   }
   else if (varHandle == iSp_handle)
   {
      return iSp;
   }
   else if (varHandle == i_Na_handle)
   {
      double nao = 140.0;
      double R = 8314.0;
      double T = 310.0;
      double F = 96485.0;
      double _expensive_functions = log(nao/nai);
      double ENa = R*T*_expensive_functions/F;
      double i_Na = (m*m*m)*g_Na*h*j*(-ENa + V);
      return i_Na;
   }
   else if (varHandle == j_handle)
   {
      return j;
   }
   else if (varHandle == jca_handle)
   {
      return jca;
   }
   else if (varHandle == ki_handle)
   {
      return ki;
   }
   else if (varHandle == kss_handle)
   {
      return kss;
   }
   else if (varHandle == m_handle)
   {
      return m;
   }
   else if (varHandle == mL_handle)
   {
      return mL;
   }
   else if (varHandle == nai_handle)
   {
      return nai;
   }
   else if (varHandle == nass_handle)
   {
      return nass;
   }
   else if (varHandle == nca_handle)
   {
      return nca;
   }
   else if (varHandle == xk1_handle)
   {
      return xk1;
   }
   else if (varHandle == xrf_handle)
   {
      return xrf;
   }
   else if (varHandle == xrs_handle)
   {
      return xrs;
   }
   else if (varHandle == xs1_handle)
   {
      return xs1;
   }
   else if (varHandle == xs2_handle)
   {
      return xs2;
   }
   return NAN;
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.clear();
   fieldUnits.clear();
   fieldNames.push_back("CaMKt");
   fieldUnits.push_back(getUnit("CaMKt"));
   fieldNames.push_back("Jrelnp");
   fieldUnits.push_back(getUnit("Jrelnp"));
   fieldNames.push_back("Jrelp");
   fieldUnits.push_back(getUnit("Jrelp"));
   fieldNames.push_back("a");
   fieldUnits.push_back(getUnit("a"));
   fieldNames.push_back("ap");
   fieldUnits.push_back(getUnit("ap"));
   fieldNames.push_back("cai");
   fieldUnits.push_back(getUnit("cai"));
   fieldNames.push_back("cajsr");
   fieldUnits.push_back(getUnit("cajsr"));
   fieldNames.push_back("cansr");
   fieldUnits.push_back(getUnit("cansr"));
   fieldNames.push_back("cass");
   fieldUnits.push_back(getUnit("cass"));
   fieldNames.push_back("d");
   fieldUnits.push_back(getUnit("d"));
   fieldNames.push_back("fcaf");
   fieldUnits.push_back(getUnit("fcaf"));
   fieldNames.push_back("fcafp");
   fieldUnits.push_back(getUnit("fcafp"));
   fieldNames.push_back("fcas");
   fieldUnits.push_back(getUnit("fcas"));
   fieldNames.push_back("ff");
   fieldUnits.push_back(getUnit("ff"));
   fieldNames.push_back("ffp");
   fieldUnits.push_back(getUnit("ffp"));
   fieldNames.push_back("fs");
   fieldUnits.push_back(getUnit("fs"));
   fieldNames.push_back("h");
   fieldUnits.push_back(getUnit("h"));
   fieldNames.push_back("hL");
   fieldUnits.push_back(getUnit("hL"));
   fieldNames.push_back("hLp");
   fieldUnits.push_back(getUnit("hLp"));
   fieldNames.push_back("iF");
   fieldUnits.push_back(getUnit("iF"));
   fieldNames.push_back("iFp");
   fieldUnits.push_back(getUnit("iFp"));
   fieldNames.push_back("iS");
   fieldUnits.push_back(getUnit("iS"));
   fieldNames.push_back("iSp");
   fieldUnits.push_back(getUnit("iSp"));
   fieldNames.push_back("j");
   fieldUnits.push_back(getUnit("j"));
   fieldNames.push_back("jca");
   fieldUnits.push_back(getUnit("jca"));
   fieldNames.push_back("ki");
   fieldUnits.push_back(getUnit("ki"));
   fieldNames.push_back("kss");
   fieldUnits.push_back(getUnit("kss"));
   fieldNames.push_back("m");
   fieldUnits.push_back(getUnit("m"));
   fieldNames.push_back("mL");
   fieldUnits.push_back(getUnit("mL"));
   fieldNames.push_back("nai");
   fieldUnits.push_back(getUnit("nai"));
   fieldNames.push_back("nass");
   fieldUnits.push_back(getUnit("nass"));
   fieldNames.push_back("nca");
   fieldUnits.push_back(getUnit("nca"));
   fieldNames.push_back("xk1");
   fieldUnits.push_back(getUnit("xk1"));
   fieldNames.push_back("xrf");
   fieldUnits.push_back(getUnit("xrf"));
   fieldNames.push_back("xrs");
   fieldUnits.push_back(getUnit("xrs"));
   fieldNames.push_back("xs1");
   fieldUnits.push_back(getUnit("xs1"));
   fieldNames.push_back("xs2");
   fieldUnits.push_back(getUnit("xs2"));
}

}
