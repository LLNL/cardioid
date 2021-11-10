/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "ohara_cipa_SCC.hh"
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


   REACTION_FACTORY(ohara_cipa_SCC)(OBJECT* obj, const double _dt, const int numPoints, const ThreadTeam&)
   {
      ohara_cipa_SCC::ThisReaction* reaction = new ohara_cipa_SCC::ThisReaction(numPoints, _dt);

      //override the defaults
      //EDIT_PARAMETERS
      double JrelStiffConst;
      double celltype;
      setDefault(celltype, 0);
      setDefault(JrelStiffConst, 0.0050000000000000001);
      reaction->JrelStiffConst = JrelStiffConst;
      reaction->celltype = celltype;
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

namespace ohara_cipa_SCC 
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

   "   C1_off,\n"
   "   C2_off,\n"
   "   CaMKt_off,\n"
   "   Cbound_off,\n"
   "   D_off,\n"
   "   IC1_off,\n"
   "   IC2_off,\n"
   "   IO_off,\n"
   "   IObound_off,\n"
   "   Jrelnp_off,\n"
   "   Jrelp_off,\n"
   "   O_off,\n"
   "   Obound_off,\n"
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
   "   hL_off,\n"
   "   hLp_off,\n"
   "   hf_off,\n"
   "   hs_off,\n"
   "   hsp_off,\n"
   "   iF_off,\n"
   "   iFp_off,\n"
   "   iS_off,\n"
   "   iSp_off,\n"
   "   j_off,\n"
   "   jca_off,\n"
   "   jp_off,\n"
   "   ki_off,\n"
   "   kss_off,\n"
   "   m_off,\n"
   "   mL_off,\n"
   "   nai_off,\n"
   "   nass_off,\n"
   "   nca_off,\n"
   "   xk1_off,\n"
   "   xs1_off,\n"
   "   xs2_off,\n"
   "   NUMSTATES\n"
   "};\n"
   "extern \"C\"\n"
   "__global__ void ohara_cipa_SCC_kernel(const int* _indexArray, const double* _Vm, const double* _iStim, double* _dVm, double* _state) {\n"
   "const double _dt = " << __cachedDt << ";\n"
   "const int _nCells = " << nCells_ << ";\n"

   "const double JrelStiffConst = " << JrelStiffConst << ";\n"
   "const double celltype = " << celltype << ";\n"
   "const int _ii = threadIdx.x + blockIdx.x*blockDim.x;\n"
   "if (_ii >= _nCells) { return; }\n"
   "const double V = _Vm[_indexArray[_ii]];\n"
   "double _ratPoly;\n"

   "const double C1 = _state[_ii+C1_off*_nCells];\n"
   "const double C2 = _state[_ii+C2_off*_nCells];\n"
   "const double CaMKt = _state[_ii+CaMKt_off*_nCells];\n"
   "const double Cbound = _state[_ii+Cbound_off*_nCells];\n"
   "const double D = _state[_ii+D_off*_nCells];\n"
   "const double IC1 = _state[_ii+IC1_off*_nCells];\n"
   "const double IC2 = _state[_ii+IC2_off*_nCells];\n"
   "const double IO = _state[_ii+IO_off*_nCells];\n"
   "const double IObound = _state[_ii+IObound_off*_nCells];\n"
   "const double Jrelnp = _state[_ii+Jrelnp_off*_nCells];\n"
   "const double Jrelp = _state[_ii+Jrelp_off*_nCells];\n"
   "const double O = _state[_ii+O_off*_nCells];\n"
   "const double Obound = _state[_ii+Obound_off*_nCells];\n"
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
   "const double hL = _state[_ii+hL_off*_nCells];\n"
   "const double hLp = _state[_ii+hLp_off*_nCells];\n"
   "const double hf = _state[_ii+hf_off*_nCells];\n"
   "const double hs = _state[_ii+hs_off*_nCells];\n"
   "const double hsp = _state[_ii+hsp_off*_nCells];\n"
   "const double iF = _state[_ii+iF_off*_nCells];\n"
   "const double iFp = _state[_ii+iFp_off*_nCells];\n"
   "const double iS = _state[_ii+iS_off*_nCells];\n"
   "const double iSp = _state[_ii+iSp_off*_nCells];\n"
   "const double j = _state[_ii+j_off*_nCells];\n"
   "const double jca = _state[_ii+jca_off*_nCells];\n"
   "const double jp = _state[_ii+jp_off*_nCells];\n"
   "const double ki = _state[_ii+ki_off*_nCells];\n"
   "const double kss = _state[_ii+kss_off*_nCells];\n"
   "const double m = _state[_ii+m_off*_nCells];\n"
   "const double mL = _state[_ii+mL_off*_nCells];\n"
   "const double nai = _state[_ii+nai_off*_nCells];\n"
   "const double nass = _state[_ii+nass_off*_nCells];\n"
   "const double nca = _state[_ii+nca_off*_nCells];\n"
   "const double xk1 = _state[_ii+xk1_off*_nCells];\n"
   "const double xs1 = _state[_ii+xs1_off*_nCells];\n"
   "const double xs2 = _state[_ii+xs2_off*_nCells];\n"
   "//get the gate updates (diagonalized exponential integrator)\n"
   "double v = V;\n"
   "double cao = 1.8;\n"
   "double ko = 5.4000000000000004;\n"
   "double F = 96485;\n"
   "double R = 8314;\n"
   "double T = 310;\n"
   "double frt = F/(R*T);\n"
   "double ffrt = F*frt;\n"
   "double vfrt = frt*v;\n"
   "double KmCaMK = 0.14999999999999999;\n"
   "double CaMKo = 0.050000000000000003;\n"
   "double KmCaM = 0.0015;\n"
   "double _expensive_functions = exp((1.0/10.0)*v - 1);\n"
   "double Afcaf = 0.29999999999999999 + 0.59999999999999998/(_expensive_functions + 1);\n"
   "double Aff = 0.59999999999999998;\n"
   "double B_1 = 2*frt;\n"
   "double PCa_b = 0.00010069999999999999;\n"
   "double v0 = 0;\n"
   "double Afcas = 1 - Afcaf;\n"
   "double Afs = 1 - Aff;\n"
   "double __melodee_temp_030 = celltype == 2;\n"
   "double __melodee_temp_031;\n"
   "if (__melodee_temp_030)\n"
   "{\n"
   "   __melodee_temp_031 = 2.5*PCa_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_031 = PCa_b;\n"
   "}\n"
   "double __melodee_temp_032 = celltype == 1;\n"
   "double __melodee_temp_033;\n"
   "if (__melodee_temp_032)\n"
   "{\n"
   "   __melodee_temp_033 = 1.2*PCa_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_033 = __melodee_temp_031;\n"
   "}\n"
   "double PCa = __melodee_temp_033;\n"
   "double U_1 = B_1*(v - v0);\n"
   "double PCap = 1.1000000000000001*PCa;\n"
   "double __melodee_temp_036 = -9.9999999999999995e-8 >= U_1 && U_1 >= 9.9999999999999995e-8;\n"
   "double f = Aff*ff + Afs*fs;\n"
   "double fp = Aff*ffp + Afs*fs;\n"
   "double _expensive_functions_014 = exp(-0.049115913555992145*v);\n"
   "double _expensive_functions_015 = exp(0.014423770373575654*v);\n"
   "double txk1 = 122.2/(0.0019352007631390235*_expensive_functions_014 + 30.433647575249029*_expensive_functions_015);\n"
   "double _expensive_functions_016 = exp((-2.5537999999999998*ko - v - 144.59)/(1.5691999999999999*ko + 3.8115000000000001));\n"
   "double xk1ss = (1.0/(_expensive_functions_016 + 1));\n"
   "double txs1_max = 817.29999999999995;\n"
   "double _expensive_functions_034 = exp(-1.0/31.0*v);\n"
   "double _expensive_functions_035 = exp((1.0/20.0)*v - 5.0/2.0);\n"
   "double txs2 = (1.0/(0.0022561357010639103*_expensive_functions_034 + 0.01*_expensive_functions_035));\n"
   "double _expensive_functions_036 = exp(-0.11195700850873264*v);\n"
   "double xs1ss = (1.0/(0.27288596035656526*_expensive_functions_036 + 1));\n"
   "double _expensive_functions_037 = exp(0.056179775280898875*v);\n"
   "double _expensive_functions_038 = exp(-1.0/230.0*v - 21.0/23.0);\n"
   "double txs1 = txs1_max + (1.0/(0.0035040677630748581*_expensive_functions_037 + 0.001292*_expensive_functions_038));\n"
   "double xs2ss = xs1ss;\n"
   "double hssV1 = 82.900000000000006;\n"
   "double hssV2 = 6.0860000000000003;\n"
   "double mssV1 = 39.57;\n"
   "double mssV2 = 9.8710000000000004;\n"
   "double mtD1 = 6.7649999999999997;\n"
   "double mtD2 = 8.5519999999999996;\n"
   "double mtV1 = 11.640000000000001;\n"
   "double mtV2 = 34.770000000000003;\n"
   "double mtV3 = 77.420000000000002;\n"
   "double mtV4 = 5.9550000000000001;\n"
   "double shift_INa_inact = 0;\n"
   "double _expensive_functions_039 = exp((hssV1 - shift_INa_inact + v)/hssV2);\n"
   "double hss = (1.0/(_expensive_functions_039 + 1));\n"
   "double _expensive_functions_040 = exp(-0.16431153466973381*shift_INa_inact + 0.16431153466973381*v);\n"
   "double hssp = (1.0/(2281075.8166971938*_expensive_functions_040 + 1));\n"
   "double _expensive_functions_041 = exp((-mssV1 - v)/mssV2);\n"
   "double mss = (1.0/(_expensive_functions_041 + 1));\n"
   "double _expensive_functions_042 = exp(0.15910898965791567*shift_INa_inact - 0.15910898965791567*v);\n"
   "double _expensive_functions_043 = exp(-0.049333991119881598*shift_INa_inact + 0.049333991119881598*v);\n"
   "double thf = (1.0/(1.183856958289087e-5*_expensive_functions_042 + 6.3055491858172754*_expensive_functions_043));\n"
   "double _expensive_functions_044 = exp(0.035650623885918005*shift_INa_inact - 0.035650623885918005*v);\n"
   "double _expensive_functions_045 = exp(-0.017649135192375574*shift_INa_inact + 0.017649135192375574*v);\n"
   "double ths = (1.0/(0.0051646702353817919*_expensive_functions_044 + 0.36987619372096325*_expensive_functions_045));\n"
   "double _expensive_functions_046 = exp(-0.02600780234070221*shift_INa_inact + 0.02600780234070221*v);\n"
   "double _expensive_functions_047 = exp(0.12075836251660427*shift_INa_inact - 0.12075836251660427*v);\n"
   "double tj = 2.0379999999999998 + (1.0/(0.31319363947387729*_expensive_functions_046 + 1.1315282095590072e-7*_expensive_functions_047));\n"
   "double _expensive_functions_048 = exp((mtV1 + v)/mtV2);\n"
   "double _expensive_functions_049 = exp((-mtV3 - v)/mtV4);\n"
   "double tm = (1.0/(_expensive_functions_048*mtD1 + _expensive_functions_049*mtD2));\n"
   "double jss = hss;\n"
   "double thsp = 3*ths;\n"
   "double tjp = 1.46*tj;\n"
   "double _expensive_functions_052 = exp(0.13354700854700854*v);\n"
   "double hLss = (1.0/(120578.15595522427*_expensive_functions_052 + 1));\n"
   "double _expensive_functions_053 = exp(0.13354700854700854*v);\n"
   "double hLssp = (1.0/(275969.29038698709*_expensive_functions_053 + 1));\n"
   "double _expensive_functions_054 = exp(-0.18996960486322187*v);\n"
   "double mLss = (1.0/(0.00029157958563553099*_expensive_functions_054 + 1));\n"
   "double thL = 200;\n"
   "double tmL = tm;\n"
   "double thLp = 3*thL;\n"
   "double _expensive_functions_056 = exp(-0.067476383265856948*v);\n"
   "double ass = (1.0/(2.6316508161673635*_expensive_functions_056 + 1));\n"
   "double _expensive_functions_057 = exp(-0.067476383265856948*v);\n"
   "double assp = (1.0/(5.1674284622306663*_expensive_functions_057 + 1));\n"
   "double __melodee_temp_020 = celltype == 1;\n"
   "double __melodee_temp_021;\n"
   "if (__melodee_temp_020)\n"
   "{\n"
   "   double _expensive_functions_058 = exp((1.0/5.0)*v + 14);\n"
   "   __melodee_temp_021 = 1 - 0.94999999999999996/(_expensive_functions_058 + 1);\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_021 = 1;\n"
   "}\n"
   "double delta_epi = __melodee_temp_021;\n"
   "double _expensive_functions_058 = exp(0.062932662051604776*v);\n"
   "double _expensive_functions_059 = exp(-4.6425255338904359*v);\n"
   "double dti_develop = 1.3540000000000001 + 0.0001/(2.6591269045230603e-5*_expensive_functions_058 + 4.5541779737128264e+24*_expensive_functions_059);\n"
   "double _expensive_functions_060 = exp((1.0/20.0)*v + 7.0/2.0);\n"
   "double dti_recover = 1 - 0.5/(_expensive_functions_060 + 1);\n"
   "double _expensive_functions_061 = exp(0.17510068289266328*v);\n"
   "double iss = (1.0/(2194.970764538301*_expensive_functions_061 + 1));\n"
   "double _expensive_functions_062 = exp(-0.034035137876343539*v);\n"
   "double _expensive_functions_063 = exp(0.034035137876343539*v);\n"
   "double ta = 1.0515000000000001/(3.5/(30.069572727397507*_expensive_functions_063 + 1) + (1.0/(2.2621017070578837*_expensive_functions_062 + 1.2089000000000001)));\n"
   "double _expensive_functions_064 = exp(-1.0/100.0*v - 1);\n"
   "double _expensive_functions_065 = exp(0.060277275467148887*v);\n"
   "double tiF_b = 4.5620000000000003 + (1.0/(0.39329999999999998*_expensive_functions_064 + 1.6300896349780942*_expensive_functions_065));\n"
   "double _expensive_functions_066 = exp(-0.016934801016088061*v);\n"
   "double _expensive_functions_067 = exp(0.12377769525931426*v);\n"
   "double tiS_b = 23.620000000000001 + (1.0/(0.00027617763953377436*_expensive_functions_066 + 0.024208962804604526*_expensive_functions_067));\n"
   "double tiF = delta_epi*tiF_b;\n"
   "double tiS = delta_epi*tiS_b;\n"
   "double tiFp = dti_develop*dti_recover*tiF;\n"
   "double tiSp = dti_develop*dti_recover*tiS;\n"
   "double bt = 4.75;\n"
   "double a_rel = 0.5*bt;\n"
   "double btp = 1.25*bt;\n"
   "double tau_rel_temp = bt/(1 + 0.0123/cajsr);\n"
   "double a_relp = 0.5*btp;\n"
   "double __melodee_temp_047 = tau_rel_temp < 0.001;\n"
   "double __melodee_temp_048;\n"
   "if (__melodee_temp_047)\n"
   "{\n"
   "   __melodee_temp_048 = 0.001;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_048 = tau_rel_temp;\n"
   "}\n"
   "double tau_rel = __melodee_temp_048;\n"
   "double tau_relp_temp = btp/(1 + 0.0123/cajsr);\n"
   "double __melodee_temp_049 = celltype == 2;\n"
   "double __melodee_temp_051 = tau_relp_temp < 0.001;\n"
   "double __melodee_temp_052;\n"
   "if (__melodee_temp_051)\n"
   "{\n"
   "   __melodee_temp_052 = 0.001;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_052 = tau_relp_temp;\n"
   "}\n"
   "double tau_relp = __melodee_temp_052;\n"
   "double __melodee_temp_053 = celltype == 2;\n"
   "double CaMKb = CaMKo*(1 - CaMKt)/(KmCaM/cass + 1);\n"
   "double CaMKa = CaMKb + CaMKt;\n"
   "double fICaLp = (1.0/(1 + KmCaMK/CaMKa));\n"
   "double _expensive_functions_071 = exp(2*vfrt);\n"
   "double A_1 = 4*ffrt*(_expensive_functions_071*cass - 0.34100000000000003*cao)/B_1;\n"
   "double __melodee_temp_037;\n"
   "if (__melodee_temp_036)\n"
   "{\n"
   "   __melodee_temp_037 = A_1*(1 - 0.5*U_1);\n"
   "}\n"
   "else\n"
   "{\n"
   "   double _expensive_functions_074 = exp(U_1);\n"
   "   __melodee_temp_037 = A_1*U_1/(_expensive_functions_074 - 1);\n"
   "}\n"
   "double PhiCaL = __melodee_temp_037;\n"
   "double fca = Afcaf*fcaf + Afcas*fcas;\n"
   "double fcap = Afcaf*fcafp + Afcas*fcas;\n"
   "double ICaL = PCa*PhiCaL*d*(1 - fICaLp)*(f*(1 - nca) + fca*jca*nca) + PCap*PhiCaL*d*fICaLp*(fcap*jca*nca + fp*(1 - nca));\n"
   "double Jrel_inf_temp = -ICaL*a_rel/(1 + 25.62890625/(cajsr*cajsr*cajsr*cajsr*cajsr*cajsr*cajsr*cajsr));\n"
   "double __melodee_temp_050;\n"
   "if (__melodee_temp_049)\n"
   "{\n"
   "   __melodee_temp_050 = 1.7*Jrel_inf_temp;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_050 = Jrel_inf_temp;\n"
   "}\n"
   "double Jrel_inf = __melodee_temp_050;\n"
   "double Jrel_temp = -ICaL*a_relp/(1 + 25.62890625/(cajsr*cajsr*cajsr*cajsr*cajsr*cajsr*cajsr*cajsr));\n"
   "double __melodee_temp_054;\n"
   "if (__melodee_temp_053)\n"
   "{\n"
   "   __melodee_temp_054 = 1.7*Jrel_temp;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_054 = Jrel_temp;\n"
   "}\n"
   "double Jrel_infp = __melodee_temp_054;\n"
   "double _expensive_functions_178 = exp(-_dt/tau_rel);\n"
   "double _Jrelnp_RLA = _expensive_functions_178 - 1;\n"
   "double _Jrelnp_RLB = -Jrel_inf;\n"
   "double _expensive_functions_179 = exp(-_dt/tau_relp);\n"
   "double _Jrelp_RLA = _expensive_functions_179 - 1;\n"
   "double _Jrelp_RLB = -Jrel_infp;\n"
   "double _expensive_functions_180 = exp(-_dt/ta);\n"
   "double _a_RLA = _expensive_functions_180 - 1;\n"
   "double _a_RLB = -ass;\n"
   "double _expensive_functions_181 = exp(-_dt/ta);\n"
   "double _ap_RLA = _expensive_functions_181 - 1;\n"
   "double _ap_RLB = -assp;\n"
   "double _expensive_functions_182 = exp(-_dt/thL);\n"
   "double _hL_RLA = _expensive_functions_182 - 1;\n"
   "double _hL_RLB = -hLss;\n"
   "double _expensive_functions_183 = exp(-_dt/thLp);\n"
   "double _hLp_RLA = _expensive_functions_183 - 1;\n"
   "double _hLp_RLB = -hLssp;\n"
   "double _expensive_functions_184 = exp(-_dt/thf);\n"
   "double _hf_RLA = _expensive_functions_184 - 1;\n"
   "double _hf_RLB = -hss;\n"
   "double _expensive_functions_185 = exp(-_dt/ths);\n"
   "double _hs_RLA = _expensive_functions_185 - 1;\n"
   "double _hs_RLB = -hss;\n"
   "double _expensive_functions_186 = exp(-_dt/thsp);\n"
   "double _hsp_RLA = _expensive_functions_186 - 1;\n"
   "double _hsp_RLB = -hssp;\n"
   "double _expensive_functions_187 = exp(-_dt/tiF);\n"
   "double _iF_RLA = _expensive_functions_187 - 1;\n"
   "double _iF_RLB = -iss;\n"
   "double _expensive_functions_188 = exp(-_dt/tiFp);\n"
   "double _iFp_RLA = _expensive_functions_188 - 1;\n"
   "double _iFp_RLB = -iss;\n"
   "double _expensive_functions_189 = exp(-_dt/tiS);\n"
   "double _iS_RLA = _expensive_functions_189 - 1;\n"
   "double _iS_RLB = -iss;\n"
   "double _expensive_functions_190 = exp(-_dt/tiSp);\n"
   "double _iSp_RLA = _expensive_functions_190 - 1;\n"
   "double _iSp_RLB = -iss;\n"
   "double _expensive_functions_191 = exp(-_dt/tj);\n"
   "double _j_RLA = _expensive_functions_191 - 1;\n"
   "double _j_RLB = -jss;\n"
   "double _expensive_functions_192 = exp(-_dt/tjp);\n"
   "double _jp_RLA = _expensive_functions_192 - 1;\n"
   "double _jp_RLB = -jss;\n"
   "double _expensive_functions_193 = exp(-_dt/tm);\n"
   "double _m_RLA = _expensive_functions_193 - 1;\n"
   "double _m_RLB = -mss;\n"
   "double _expensive_functions_194 = exp(-_dt/tmL);\n"
   "double _mL_RLA = _expensive_functions_194 - 1;\n"
   "double _mL_RLB = -mLss;\n"
   "double _expensive_functions_195 = exp(-_dt/txk1);\n"
   "double _xk1_RLA = _expensive_functions_195 - 1;\n"
   "double _xk1_RLB = -xk1ss;\n"
   "double _expensive_functions_196 = exp(-_dt/txs1);\n"
   "double _xs1_RLA = _expensive_functions_196 - 1;\n"
   "double _xs1_RLB = -xs1ss;\n"
   "double _expensive_functions_197 = exp(-_dt/txs2);\n"
   "double _xs2_RLA = _expensive_functions_197 - 1;\n"
   "double _xs2_RLB = -xs2ss;\n"
   "//get the other differential updates\n"
   "double nao = 140;\n"
   "double L = 0.01;\n"
   "double rad = 0.0011000000000000001;\n"
   "double Ageo = 6.2800000000000002*L*rad + 6.2800000000000002*(rad*rad);\n"
   "double vcell = 3140.0*(rad*rad)*L;\n"
   "double Acap = 2*Ageo;\n"
   "double vjsr = 0.0047999999999999996*vcell;\n"
   "double vmyo = 0.68000000000000005*vcell;\n"
   "double vnsr = 0.055199999999999999*vcell;\n"
   "double vss = 0.02*vcell;\n"
   "double zca = 2;\n"
   "double zk = 1;\n"
   "double zna = 1;\n"
   "double PKNa = 0.018329999999999999;\n"
   "double cm = 1;\n"
   "double aCaMK = 0.050000000000000003;\n"
   "double bCaMK = 0.00068000000000000005;\n"
   "double B_2 = frt;\n"
   "double B_3 = frt;\n"
   "double Kmn = 0.002;\n"
   "double _expensive_functions_001 = exp(-0.23640661938534277*v);\n"
   "double dss = (1.0/(0.39398514226669484*_expensive_functions_001 + 1));\n"
   "double _expensive_functions_002 = exp(0.27056277056277056*v);\n"
   "double fss = (1.0/(199.86038496778565*_expensive_functions_002 + 1));\n"
   "double k2n = 1000;\n"
   "double _expensive_functions_003 = exp(0.089999999999999997*v);\n"
   "double _expensive_functions_004 = exp(-0.050000000000000003*v);\n"
   "double td = 0.59999999999999998 + (1.0/(3.5254214873653824*_expensive_functions_003 + 0.74081822068171788*_expensive_functions_004));\n"
   "double _expensive_functions_005 = exp((1.0/7.0)*v - 4.0/7.0);\n"
   "double _expensive_functions_006 = exp(4.0/7.0 - 1.0/7.0*v);\n"
   "double tfcaf = 7 + (1.0/(0.040000000000000001*_expensive_functions_005 + 0.040000000000000001*_expensive_functions_006));\n"
   "double _expensive_functions_007 = exp(-1.0/3.0*v);\n"
   "double _expensive_functions_008 = exp((1.0/7.0)*v);\n"
   "double tfcas = 100 + (1.0/(0.00012*_expensive_functions_007 + 0.00012*_expensive_functions_008));\n"
   "double _expensive_functions_009 = exp(-1.0/10.0*v - 2);\n"
   "double _expensive_functions_010 = exp((1.0/10.0)*v + 2);\n"
   "double tff = 7 + (1.0/(0.0044999999999999997*_expensive_functions_009 + 0.0044999999999999997*_expensive_functions_010));\n"
   "double _expensive_functions_011 = exp(-1.0/4.0*v - 5.0/4.0);\n"
   "double _expensive_functions_012 = exp((1.0/6.0)*v + 5.0/6.0);\n"
   "double tfs = 1000 + (1.0/(3.4999999999999997e-5*_expensive_functions_011 + 3.4999999999999997e-5*_expensive_functions_012));\n"
   "double tjca = 75;\n"
   "double U_2 = B_2*(v - v0);\n"
   "double U_3 = B_3*(v - v0);\n"
   "double d_diff = (-d + dss)/td;\n"
   "double fcass = fss;\n"
   "double ff_diff = (-ff + fss)/tff;\n"
   "double fs_diff = (-fs + fss)/tfs;\n"
   "double tfcafp = 2.5*tfcaf;\n"
   "double tffp = 2.5*tff;\n"
   "double PCaK = 0.00035740000000000001*PCa;\n"
   "double PCaNa = 0.00125*PCa;\n"
   "double __melodee_temp_034 = -9.9999999999999995e-8 >= U_3 && U_3 >= 9.9999999999999995e-8;\n"
   "double __melodee_temp_038 = -9.9999999999999995e-8 >= U_2 && U_2 >= 9.9999999999999995e-8;\n"
   "double fcaf_diff = (-fcaf + fcass)/tfcaf;\n"
   "double fcafp_diff = (-fcafp + fcass)/tfcafp;\n"
   "double fcas_diff = (-fcas + fcass)/tfcas;\n"
   "double ffp_diff = (-ffp + fss)/tffp;\n"
   "double jca_diff = (fcass - jca)/tjca;\n"
   "double PCaKp = 0.00035740000000000001*PCap;\n"
   "double PCaNap = 0.00125*PCap;\n"
   "double B = 2*frt;\n"
   "double PCab = 2.4999999999999999e-8;\n"
   "double ICab_v0 = 0;\n"
   "double U = B*(-ICab_v0 + v);\n"
   "double __melodee_temp_040 = -9.9999999999999995e-8 >= U && U >= 9.9999999999999995e-8;\n"
   "double GK1_b = 0.3239783999999998;\n"
   "double _expensive_functions_013 = exp(-0.27388602127883704*ko + 0.10534077741493732*v);\n"
   "double rk1 = (1.0/(69220.632210676704*_expensive_functions_013 + 1));\n"
   "double __melodee_temp_000 = celltype == 2;\n"
   "double __melodee_temp_001;\n"
   "if (__melodee_temp_000)\n"
   "{\n"
   "   __melodee_temp_001 = 1.3*GK1_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_001 = GK1_b;\n"
   "}\n"
   "double __melodee_temp_002 = celltype == 1;\n"
   "double __melodee_temp_003;\n"
   "if (__melodee_temp_002)\n"
   "{\n"
   "   __melodee_temp_003 = 1.2*GK1_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_003 = __melodee_temp_001;\n"
   "}\n"
   "double GK1 = __melodee_temp_003;\n"
   "double GKb_b = 0.0030000000000000001;\n"
   "double _expensive_functions_017 = exp(-0.054525627044711013*v);\n"
   "double xkb = (1.0/(2.2023634509492389*_expensive_functions_017 + 1));\n"
   "double __melodee_temp_004 = celltype == 1;\n"
   "double __melodee_temp_005;\n"
   "if (__melodee_temp_004)\n"
   "{\n"
   "   __melodee_temp_005 = 0.59999999999999998*GKb_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_005 = GKb_b;\n"
   "}\n"
   "double GKb = __melodee_temp_005;\n"
   "double A1 = 0.0264;\n"
   "double A11 = 0.00078680000000000004;\n"
   "double A2 = 4.9860000000000002e-6;\n"
   "double A21 = 5.4550000000000003e-6;\n"
   "double A3 = 0.001214;\n"
   "double A31 = 0.005509;\n"
   "double A4 = 1.8539999999999999e-5;\n"
   "double A41 = 0.0014159999999999999;\n"
   "double A51 = 0.44919999999999999;\n"
   "double A52 = 0.31809999999999999;\n"
   "double A53 = 0.14899999999999999;\n"
   "double A61 = 0.012409999999999999;\n"
   "double A62 = 0.3226;\n"
   "double A63 = 0.0089779999999999999;\n"
   "double B1 = 4.6310000000000002e-5;\n"
   "double B11 = 1.5349999999999998e-8;\n"
   "double B2 = -0.0042259999999999997;\n"
   "double B21 = -0.16880000000000001;\n"
   "double B3 = 0.0085159999999999993;\n"
   "double B31 = 7.7710000000000006e-9;\n"
   "double B4 = -0.04641;\n"
   "double B41 = -0.02877;\n"
   "double B51 = 0.0085950000000000002;\n"
   "double B52 = 3.613e-8;\n"
   "double B53 = 0.0046680000000000003;\n"
   "double B61 = 0.17249999999999999;\n"
   "double B62 = -0.00065749999999999999;\n"
   "double B63 = -0.02215;\n"
   "double D_diff = 0;\n"
   "double GKr_b = 0.04658545454545456;\n"
   "double Kmax = 0;\n"
   "double Kt = 0;\n"
   "double Ku = 0;\n"
   "double Temp = 37;\n"
   "double Vhalf = 1;\n"
   "double halfmax = 1;\n"
   "double n = 1;\n"
   "double q1 = 4.843;\n"
   "double q11 = 4.9420000000000002;\n"
   "double q2 = 4.2300000000000004;\n"
   "double q21 = 4.1559999999999997;\n"
   "double q3 = 4.9619999999999997;\n"
   "double q31 = 4.2199999999999998;\n"
   "double q4 = 3.7690000000000001;\n"
   "double q41 = 1.4590000000000001;\n"
   "double q51 = 5;\n"
   "double q52 = 4.6630000000000003;\n"
   "double q53 = 2.4119999999999999;\n"
   "double q61 = 5.5679999999999996;\n"
   "double q62 = 5;\n"
   "double q63 = 5.6820000000000004;\n"
   "double __melodee_temp_006 = celltype == 2;\n"
   "double __melodee_temp_007;\n"
   "if (__melodee_temp_006)\n"
   "{\n"
   "   __melodee_temp_007 = 0.80000000000000004*GKr_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_007 = GKr_b;\n"
   "}\n"
   "double __melodee_temp_008 = celltype == 1;\n"
   "double __melodee_temp_009;\n"
   "if (__melodee_temp_008)\n"
   "{\n"
   "   __melodee_temp_009 = 1.3*GKr_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_009 = __melodee_temp_007;\n"
   "}\n"
   "double GKr = __melodee_temp_009;\n"
   "double _expensive_functions_018 = exp(0.14729709824716453*Vhalf - 0.14729709824716453*v);\n"
   "double _expensive_functions_019 = log(D);\n"
   "double _expensive_functions_020 = exp(_expensive_functions_019*n);\n"
   "double _expensive_functions_021 = log(D);\n"
   "double _expensive_functions_022 = exp(_expensive_functions_021*n);\n"
   "double _expensive_functions_023 = exp(B53*v);\n"
   "double _expensive_functions_024 = exp(-B63*v);\n"
   "double _expensive_functions_025 = log(q63);\n"
   "double _expensive_functions_026 = exp(-1.0/10.0*_expensive_functions_025*(Temp - 20));\n"
   "double _expensive_functions_027 = log(q53);\n"
   "double _expensive_functions_028 = exp((1.0/10.0)*_expensive_functions_027*(Temp - 20));\n"
   "double IObound_diff = -A53*IObound*Ku*_expensive_functions_023*_expensive_functions_024*_expensive_functions_026*_expensive_functions_028/A63 + Cbound*Kt/(_expensive_functions_018 + 1) + IO*Kmax*Ku*_expensive_functions_022/(_expensive_functions_020 + halfmax) - IObound*Kt;\n"
   "double _expensive_functions_029 = exp(0.14729709824716453*Vhalf - 0.14729709824716453*v);\n"
   "double _expensive_functions_030 = log(D);\n"
   "double _expensive_functions_031 = exp(_expensive_functions_030*n);\n"
   "double _expensive_functions_032 = log(D);\n"
   "double _expensive_functions_033 = exp(_expensive_functions_032*n);\n"
   "double Obound_diff = Cbound*Kt/(_expensive_functions_029 + 1) + Kmax*Ku*O*_expensive_functions_033/(_expensive_functions_031 + halfmax) - Kt*Obound - Ku*Obound;\n"
   "double GKs_b = 0.006358000000000001;\n"
   "double __melodee_temp_010 = celltype == 1;\n"
   "double __melodee_temp_011;\n"
   "if (__melodee_temp_010)\n"
   "{\n"
   "   __melodee_temp_011 = 1.3999999999999999*GKs_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_011 = GKs_b;\n"
   "}\n"
   "double GKs = __melodee_temp_011;\n"
   "double Ahf = 0.98999999999999999;\n"
   "double GNa = 75;\n"
   "double Ahs = 1 - Ahf;\n"
   "double Gncx_b = 0.00080000000000000004;\n"
   "double KmCaAct = 0.00014999999999999999;\n"
   "double kasymm = 12.5;\n"
   "double kcaoff = 5000.0;\n"
   "double kcaon = 1500000.0;\n"
   "double kna1 = 15;\n"
   "double kna2 = 5;\n"
   "double kna3 = 88.120000000000005;\n"
   "double qca = 0.16700000000000001;\n"
   "double qna = 0.52239999999999998;\n"
   "double wca = 60000.0;\n"
   "double wna = 60000.0;\n"
   "double wnaca = 5000.0;\n"
   "double __melodee_temp_012 = celltype == 2;\n"
   "double __melodee_temp_013;\n"
   "if (__melodee_temp_012)\n"
   "{\n"
   "   __melodee_temp_013 = 1.3999999999999999*Gncx_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_013 = Gncx_b;\n"
   "}\n"
   "double __melodee_temp_014 = celltype == 1;\n"
   "double __melodee_temp_015;\n"
   "if (__melodee_temp_014)\n"
   "{\n"
   "   __melodee_temp_015 = 1.1000000000000001*Gncx_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_015 = __melodee_temp_013;\n"
   "}\n"
   "double Gncx = __melodee_temp_015;\n"
   "double h10_i = kasymm + 1 + nao*(1 + nao/kna2)/kna1;\n"
   "double h10_ss = kasymm + 1 + nao*(1 + nao/kna2)/kna1;\n"
   "double hca = exp(F*qca*v/(R*T));\n"
   "double hna = exp(F*qna*v/(R*T));\n"
   "double k2_i = kcaoff;\n"
   "double k2_ss = kcaoff;\n"
   "double k5_i = kcaoff;\n"
   "double k5_ss = kcaoff;\n"
   "double h11_i = (nao*nao)/(h10_i*kna1*kna2);\n"
   "double h11_ss = (nao*nao)/(h10_ss*kna1*kna2);\n"
   "double h12_i = (1.0/h10_i);\n"
   "double h12_ss = (1.0/h10_ss);\n"
   "double h7_i = 1 + nao*(1 + (1.0/hna))/kna3;\n"
   "double h7_ss = 1 + nao*(1 + (1.0/hna))/kna3;\n"
   "double h8_i = nao/(h7_i*hna*kna3);\n"
   "double h8_ss = nao/(h7_ss*hna*kna3);\n"
   "double h9_i = (1.0/h7_i);\n"
   "double h9_ss = (1.0/h7_ss);\n"
   "double k1_i = cao*h12_i*kcaon;\n"
   "double k1_ss = cao*h12_ss*kcaon;\n"
   "double k3p_i = h9_i*wca;\n"
   "double k3p_ss = h9_ss*wca;\n"
   "double k3pp_i = h8_i*wnaca;\n"
   "double k3pp_ss = h8_ss*wnaca;\n"
   "double k8_i = h11_i*h8_i*wna;\n"
   "double k8_ss = h11_ss*h8_ss*wna;\n"
   "double k3_i = k3p_i + k3pp_i;\n"
   "double k3_ss = k3p_ss + k3pp_ss;\n"
   "double H = 9.9999999999999995e-8;\n"
   "double Khp = 1.698e-7;\n"
   "double Kki = 0.5;\n"
   "double Kko = 0.35820000000000002;\n"
   "double Kmgatp = 1.698e-7;\n"
   "double Knai0 = 9.0730000000000004;\n"
   "double Knao0 = 27.780000000000001;\n"
   "double Knap = 224;\n"
   "double Kxkur = 292;\n"
   "double MgADP = 0.050000000000000003;\n"
   "double MgATP = 9.8000000000000007;\n"
   "double Pnak_b = 30;\n"
   "double delta = -0.155;\n"
   "double eP = 4.2000000000000002;\n"
   "double k1m = 182.40000000000001;\n"
   "double k1p = 949.5;\n"
   "double k2m = 39.399999999999999;\n"
   "double k2p = 687.20000000000005;\n"
   "double k3m = 79300;\n"
   "double k3p = 1899;\n"
   "double k4m = 40;\n"
   "double k4p = 639;\n"
   "double _expensive_functions_050 = exp((1.0/3.0)*F*delta*v/(R*T));\n"
   "double Knai = Knai0*_expensive_functions_050;\n"
   "double _expensive_functions_051 = exp((1.0/3.0)*F*v*(1 - delta)/(R*T));\n"
   "double Knao = Knao0*_expensive_functions_051;\n"
   "double __melodee_temp_016 = celltype == 2;\n"
   "double __melodee_temp_017;\n"
   "if (__melodee_temp_016)\n"
   "{\n"
   "   __melodee_temp_017 = 0.69999999999999996*Pnak_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_017 = Pnak_b;\n"
   "}\n"
   "double __melodee_temp_018 = celltype == 1;\n"
   "double __melodee_temp_019;\n"
   "if (__melodee_temp_018)\n"
   "{\n"
   "   __melodee_temp_019 = 0.90000000000000002*Pnak_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_019 = __melodee_temp_017;\n"
   "}\n"
   "double Pnak = __melodee_temp_019;\n"
   "double a2 = k2p;\n"
   "double a4 = MgATP*k4p/(Kmgatp*(1 + MgATP/Kmgatp));\n"
   "double b1 = MgADP*k1m;\n"
   "double a3 = (ko*ko)*k3p/((Kko*Kko)*(((1 + ko/Kko)*(1 + ko/Kko)) + ((1 + nao/Knao)*(1 + nao/Knao)*(1 + nao/Knao)) - 1));\n"
   "double b2 = (nao*nao*nao)*k2m/((Knao*Knao*Knao)*(((1 + ko/Kko)*(1 + ko/Kko)) + ((1 + nao/Knao)*(1 + nao/Knao)*(1 + nao/Knao)) - 1));\n"
   "double GNaL_b = 0.019957499999999975;\n"
   "double __melodee_temp_028 = celltype == 1;\n"
   "double __melodee_temp_029;\n"
   "if (__melodee_temp_028)\n"
   "{\n"
   "   __melodee_temp_029 = 0.59999999999999998*GNaL_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_029 = GNaL_b;\n"
   "}\n"
   "double GNaL = __melodee_temp_029;\n"
   "double INab_B = frt;\n"
   "double PNab = 3.75e-10;\n"
   "double INab_v0 = 0;\n"
   "double INab_U = INab_B*(-INab_v0 + v);\n"
   "double __melodee_temp_042 = -9.9999999999999995e-8 >= INab_U && INab_U >= 9.9999999999999995e-8;\n"
   "double GpCa = 0.00050000000000000001;\n"
   "double KmCap = 0.00050000000000000001;\n"
   "double _expensive_functions_055 = exp(0.0066137566137566143*v);\n"
   "double AiF = (1.0/(0.24348537187522867*_expensive_functions_055 + 1));\n"
   "double Gto_b = 0.02;\n"
   "double AiS = 1 - AiF;\n"
   "double __melodee_temp_022 = celltype == 2;\n"
   "double __melodee_temp_023;\n"
   "if (__melodee_temp_022)\n"
   "{\n"
   "   __melodee_temp_023 = 4*Gto_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_023 = Gto_b;\n"
   "}\n"
   "double __melodee_temp_024 = celltype == 1;\n"
   "double __melodee_temp_025;\n"
   "if (__melodee_temp_024)\n"
   "{\n"
   "   __melodee_temp_025 = 4*Gto_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_025 = __melodee_temp_023;\n"
   "}\n"
   "double Gto = __melodee_temp_025;\n"
   "double Jup_b = 1.0;\n"
   "double __melodee_temp_026 = celltype == 1;\n"
   "double __melodee_temp_027;\n"
   "if (__melodee_temp_026)\n"
   "{\n"
   "   __melodee_temp_027 = 1.3;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_027 = 1;\n"
   "}\n"
   "double upScale = __melodee_temp_027;\n"
   "double BSLmax = 1.1240000000000001;\n"
   "double BSRmax = 0.047;\n"
   "double KmBSL = 0.0086999999999999994;\n"
   "double KmBSR = 0.00087000000000000001;\n"
   "double cmdnmax_b = 0.050000000000000003;\n"
   "double csqnmax = 10;\n"
   "double kmcmdn = 0.0023800000000000002;\n"
   "double kmcsqn = 0.80000000000000004;\n"
   "double kmtrpn = 0.00050000000000000001;\n"
   "double trpnmax = 0.070000000000000007;\n"
   "double __melodee_temp_044 = celltype == 1;\n"
   "double __melodee_temp_045;\n"
   "if (__melodee_temp_044)\n"
   "{\n"
   "   __melodee_temp_045 = 1.3*cmdnmax_b;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_045 = cmdnmax_b;\n"
   "}\n"
   "double cmdnmax = __melodee_temp_045;\n"
   "double Jrel_scaling_factor = 1.0;\n"
   "double Jtr = -1.0/100.0*cajsr + (1.0/100.0)*cansr;\n"
   "double _expensive_functions_068 = log(ko/ki);\n"
   "double EK = R*T*_expensive_functions_068/F;\n"
   "double _expensive_functions_069 = log(nao/nai);\n"
   "double ENa = R*T*_expensive_functions_069/F;\n"
   "double _expensive_functions_070 = log((PKNa*nao + ko)/(PKNa*nai + ki));\n"
   "double EKs = R*T*_expensive_functions_070/F;\n"
   "double Jdiff = -5.0*cai + 5.0*cass;\n"
   "double JdiffK = -1.0/2.0*ki + (1.0/2.0)*kss;\n"
   "double JdiffNa = -1.0/2.0*nai + (1.0/2.0)*nass;\n"
   "double CaMKt_diff = CaMKb*aCaMK*(CaMKb + CaMKt) - CaMKt*bCaMK;\n"
   "double km2n = jca;\n"
   "double _expensive_functions_072 = exp(vfrt);\n"
   "double A_2 = 0.75*ffrt*(_expensive_functions_072*nass - nao)/B_2;\n"
   "double _expensive_functions_073 = exp(vfrt);\n"
   "double A_3 = 0.75*ffrt*(_expensive_functions_073*kss - ko)/B_3;\n"
   "double anca = (1.0/(k2n/km2n + ((Kmn/cass + 1)*(Kmn/cass + 1)*(Kmn/cass + 1)*(Kmn/cass + 1))));\n"
   "double __melodee_temp_035;\n"
   "if (__melodee_temp_034)\n"
   "{\n"
   "   __melodee_temp_035 = A_3*(1 - 0.5*U_3);\n"
   "}\n"
   "else\n"
   "{\n"
   "   double _expensive_functions_074 = exp(U_3);\n"
   "   __melodee_temp_035 = A_3*U_3/(_expensive_functions_074 - 1);\n"
   "}\n"
   "double PhiCaK = __melodee_temp_035;\n"
   "double __melodee_temp_039;\n"
   "if (__melodee_temp_038)\n"
   "{\n"
   "   __melodee_temp_039 = A_2*(1 - 0.5*U_2);\n"
   "}\n"
   "else\n"
   "{\n"
   "   double _expensive_functions_074 = exp(U_2);\n"
   "   __melodee_temp_039 = A_2*U_2/(_expensive_functions_074 - 1);\n"
   "}\n"
   "double PhiCaNa = __melodee_temp_039;\n"
   "double nca_diff = anca*k2n - km2n*nca;\n"
   "double ICaK = PCaK*PhiCaK*d*(1 - fICaLp)*(f*(1 - nca) + fca*jca*nca) + PCaKp*PhiCaK*d*fICaLp*(fcap*jca*nca + fp*(1 - nca));\n"
   "double ICaNa = PCaNa*PhiCaNa*d*(1 - fICaLp)*(f*(1 - nca) + fca*jca*nca) + PCaNap*PhiCaNa*d*fICaLp*(fcap*jca*nca + fp*(1 - nca));\n"
   "double _expensive_functions_074 = exp(2*vfrt);\n"
   "double A = 4*PCab*ffrt*(_expensive_functions_074*cai - 0.34100000000000003*cao)/B;\n"
   "double __melodee_temp_041;\n"
   "if (__melodee_temp_040)\n"
   "{\n"
   "   __melodee_temp_041 = A*(1 - 0.5*U);\n"
   "}\n"
   "else\n"
   "{\n"
   "   double _expensive_functions_075 = exp(U);\n"
   "   __melodee_temp_041 = A*U/(_expensive_functions_075 - 1);\n"
   "}\n"
   "double ICab = __melodee_temp_041;\n"
   "double _expensive_functions_075 = sqrt(ko);\n"
   "double IK1 = GK1*_expensive_functions_075*rk1*xk1*(-EK + v);\n"
   "double IKb = GKb*xkb*(-EK + v);\n"
   "double _expensive_functions_076 = exp(B2*v);\n"
   "double _expensive_functions_077 = log(q2);\n"
   "double _expensive_functions_078 = exp((1.0/10.0)*_expensive_functions_077*(Temp - 20));\n"
   "double _expensive_functions_079 = exp(B61*v);\n"
   "double _expensive_functions_080 = log(q61);\n"
   "double _expensive_functions_081 = exp((1.0/10.0)*_expensive_functions_080*(Temp - 20));\n"
   "double _expensive_functions_082 = exp(B1*v);\n"
   "double _expensive_functions_083 = log(q1);\n"
   "double _expensive_functions_084 = exp((1.0/10.0)*_expensive_functions_083*(Temp - 20));\n"
   "double _expensive_functions_085 = exp(B51*v);\n"
   "double _expensive_functions_086 = log(q51);\n"
   "double _expensive_functions_087 = exp((1.0/10.0)*_expensive_functions_086*(Temp - 20));\n"
   "double C1_diff = -A1*C1*_expensive_functions_082*_expensive_functions_084 + A2*C2*_expensive_functions_076*_expensive_functions_078 - A51*C1*_expensive_functions_085*_expensive_functions_087 + A61*IC1*_expensive_functions_079*_expensive_functions_081;\n"
   "double _expensive_functions_088 = exp(B1*v);\n"
   "double _expensive_functions_089 = log(q1);\n"
   "double _expensive_functions_090 = exp((1.0/10.0)*_expensive_functions_089*(Temp - 20));\n"
   "double _expensive_functions_091 = exp(B41*v);\n"
   "double _expensive_functions_092 = log(q41);\n"
   "double _expensive_functions_093 = exp((1.0/10.0)*_expensive_functions_092*(Temp - 20));\n"
   "double _expensive_functions_094 = exp(B62*v);\n"
   "double _expensive_functions_095 = log(q62);\n"
   "double _expensive_functions_096 = exp((1.0/10.0)*_expensive_functions_095*(Temp - 20));\n"
   "double _expensive_functions_097 = exp(B2*v);\n"
   "double _expensive_functions_098 = log(q2);\n"
   "double _expensive_functions_099 = exp((1.0/10.0)*_expensive_functions_098*(Temp - 20));\n"
   "double _expensive_functions_100 = exp(B31*v);\n"
   "double _expensive_functions_101 = log(q31);\n"
   "double _expensive_functions_102 = exp((1.0/10.0)*_expensive_functions_101*(Temp - 20));\n"
   "double _expensive_functions_103 = exp(B52*v);\n"
   "double _expensive_functions_104 = log(q52);\n"
   "double _expensive_functions_105 = exp((1.0/10.0)*_expensive_functions_104*(Temp - 20));\n"
   "double C2_diff = A1*C1*_expensive_functions_088*_expensive_functions_090 - A2*C2*_expensive_functions_097*_expensive_functions_099 - A31*C2*_expensive_functions_100*_expensive_functions_102 + A41*O*_expensive_functions_091*_expensive_functions_093 - A52*C2*_expensive_functions_103*_expensive_functions_105 + A62*IC2*_expensive_functions_094*_expensive_functions_096;\n"
   "double _expensive_functions_106 = exp(0.14729709824716453*Vhalf - 0.14729709824716453*v);\n"
   "double Cbound_diff = -2*Cbound*Kt/(_expensive_functions_106 + 1) + IObound*Kt + Kt*Obound;\n"
   "double _expensive_functions_107 = exp(B21*v);\n"
   "double _expensive_functions_108 = log(q21);\n"
   "double _expensive_functions_109 = exp((1.0/10.0)*_expensive_functions_108*(Temp - 20));\n"
   "double _expensive_functions_110 = exp(B51*v);\n"
   "double _expensive_functions_111 = log(q51);\n"
   "double _expensive_functions_112 = exp((1.0/10.0)*_expensive_functions_111*(Temp - 20));\n"
   "double _expensive_functions_113 = exp(B11*v);\n"
   "double _expensive_functions_114 = log(q11);\n"
   "double _expensive_functions_115 = exp((1.0/10.0)*_expensive_functions_114*(Temp - 20));\n"
   "double _expensive_functions_116 = exp(B61*v);\n"
   "double _expensive_functions_117 = log(q61);\n"
   "double _expensive_functions_118 = exp((1.0/10.0)*_expensive_functions_117*(Temp - 20));\n"
   "double IC1_diff = -A11*IC1*_expensive_functions_113*_expensive_functions_115 + A21*IC2*_expensive_functions_107*_expensive_functions_109 + A51*C1*_expensive_functions_110*_expensive_functions_112 - A61*IC1*_expensive_functions_116*_expensive_functions_118;\n"
   "double _expensive_functions_119 = exp(B11*v);\n"
   "double _expensive_functions_120 = log(q11);\n"
   "double _expensive_functions_121 = exp((1.0/10.0)*_expensive_functions_120*(Temp - 20));\n"
   "double _expensive_functions_122 = exp(B4*v);\n"
   "double _expensive_functions_123 = log(q4);\n"
   "double _expensive_functions_124 = exp((1.0/10.0)*_expensive_functions_123*(Temp - 20));\n"
   "double _expensive_functions_125 = exp(B52*v);\n"
   "double _expensive_functions_126 = log(q52);\n"
   "double _expensive_functions_127 = exp((1.0/10.0)*_expensive_functions_126*(Temp - 20));\n"
   "double _expensive_functions_128 = exp(B21*v);\n"
   "double _expensive_functions_129 = log(q21);\n"
   "double _expensive_functions_130 = exp((1.0/10.0)*_expensive_functions_129*(Temp - 20));\n"
   "double _expensive_functions_131 = exp(B3*v);\n"
   "double _expensive_functions_132 = log(q3);\n"
   "double _expensive_functions_133 = exp((1.0/10.0)*_expensive_functions_132*(Temp - 20));\n"
   "double _expensive_functions_134 = exp(B62*v);\n"
   "double _expensive_functions_135 = log(q62);\n"
   "double _expensive_functions_136 = exp((1.0/10.0)*_expensive_functions_135*(Temp - 20));\n"
   "double IC2_diff = A11*IC1*_expensive_functions_119*_expensive_functions_121 - A21*IC2*_expensive_functions_128*_expensive_functions_130 - A3*IC2*_expensive_functions_131*_expensive_functions_133 + A4*IO*_expensive_functions_122*_expensive_functions_124 + A52*C2*_expensive_functions_125*_expensive_functions_127 - A62*IC2*_expensive_functions_134*_expensive_functions_136;\n"
   "double _expensive_functions_137 = exp(B3*v);\n"
   "double _expensive_functions_138 = log(q3);\n"
   "double _expensive_functions_139 = exp((1.0/10.0)*_expensive_functions_138*(Temp - 20));\n"
   "double _expensive_functions_140 = exp(B53*v);\n"
   "double _expensive_functions_141 = log(q53);\n"
   "double _expensive_functions_142 = exp((1.0/10.0)*_expensive_functions_141*(Temp - 20));\n"
   "double _expensive_functions_143 = exp(B4*v);\n"
   "double _expensive_functions_144 = log(q4);\n"
   "double _expensive_functions_145 = exp((1.0/10.0)*_expensive_functions_144*(Temp - 20));\n"
   "double _expensive_functions_146 = exp(B63*v);\n"
   "double _expensive_functions_147 = log(q63);\n"
   "double _expensive_functions_148 = exp((1.0/10.0)*_expensive_functions_147*(Temp - 20));\n"
   "double _expensive_functions_149 = log(D);\n"
   "double _expensive_functions_150 = exp(_expensive_functions_149*n);\n"
   "double _expensive_functions_151 = log(D);\n"
   "double _expensive_functions_152 = exp(_expensive_functions_151*n);\n"
   "double _expensive_functions_153 = exp(B53*v);\n"
   "double _expensive_functions_154 = exp(-B63*v);\n"
   "double _expensive_functions_155 = log(q63);\n"
   "double _expensive_functions_156 = exp(-1.0/10.0*_expensive_functions_155*(Temp - 20));\n"
   "double _expensive_functions_157 = log(q53);\n"
   "double _expensive_functions_158 = exp((1.0/10.0)*_expensive_functions_157*(Temp - 20));\n"
   "double IO_diff = A3*IC2*_expensive_functions_137*_expensive_functions_139 - A4*IO*_expensive_functions_143*_expensive_functions_145 + A53*O*_expensive_functions_140*_expensive_functions_142 + A53*IObound*Ku*_expensive_functions_153*_expensive_functions_154*_expensive_functions_156*_expensive_functions_158/A63 - A63*IO*_expensive_functions_146*_expensive_functions_148 - IO*Kmax*Ku*_expensive_functions_152/(_expensive_functions_150 + halfmax);\n"
   "double _expensive_functions_159 = exp(B31*v);\n"
   "double _expensive_functions_160 = log(q31);\n"
   "double _expensive_functions_161 = exp((1.0/10.0)*_expensive_functions_160*(Temp - 20));\n"
   "double _expensive_functions_162 = exp(B63*v);\n"
   "double _expensive_functions_163 = log(q63);\n"
   "double _expensive_functions_164 = exp((1.0/10.0)*_expensive_functions_163*(Temp - 20));\n"
   "double _expensive_functions_165 = exp(B41*v);\n"
   "double _expensive_functions_166 = log(q41);\n"
   "double _expensive_functions_167 = exp((1.0/10.0)*_expensive_functions_166*(Temp - 20));\n"
   "double _expensive_functions_168 = exp(B53*v);\n"
   "double _expensive_functions_169 = log(q53);\n"
   "double _expensive_functions_170 = exp((1.0/10.0)*_expensive_functions_169*(Temp - 20));\n"
   "double _expensive_functions_171 = log(D);\n"
   "double _expensive_functions_172 = exp(_expensive_functions_171*n);\n"
   "double _expensive_functions_173 = log(D);\n"
   "double _expensive_functions_174 = exp(_expensive_functions_173*n);\n"
   "double O_diff = A31*C2*_expensive_functions_159*_expensive_functions_161 - A41*O*_expensive_functions_165*_expensive_functions_167 - A53*O*_expensive_functions_168*_expensive_functions_170 + A63*IO*_expensive_functions_162*_expensive_functions_164 - Kmax*Ku*O*_expensive_functions_174/(_expensive_functions_172 + halfmax) + Ku*Obound;\n"
   "double _expensive_functions_175 = sqrt(ko);\n"
   "double IKr = 0.43033148291193518*GKr*O*_expensive_functions_175*(-EK + v);\n"
   "double _expensive_functions_176 = pow((1.0/cai), 1.3999999999999999);\n"
   "double KsCa = 1 + 0.59999999999999998/(6.4818210260626455e-7*_expensive_functions_176 + 1);\n"
   "double IKs = GKs*KsCa*xs1*xs2*(-EKs + v);\n"
   "double fINap = (1.0/(1 + KmCaMK/CaMKa));\n"
   "double h = Ahf*hf + Ahs*hs;\n"
   "double hp = Ahf*hf + Ahs*hsp;\n"
   "double INa = (m*m*m)*GNa*(-ENa + v)*(fINap*hp*jp + h*j*(1 - fINap));\n"
   "double allo_i = (1.0/((KmCaAct*KmCaAct)/(cai*cai) + 1));\n"
   "double allo_ss = (1.0/((KmCaAct*KmCaAct)/(cass*cass) + 1));\n"
   "double h4_i = 1 + nai*(1 + nai/kna2)/kna1;\n"
   "double h4_ss = 1 + nass*(1 + nass/kna2)/kna1;\n"
   "double h1_i = 1 + nai*(hna + 1)/kna3;\n"
   "double h1_ss = 1 + nass*(hna + 1)/kna3;\n"
   "double h5_i = (nai*nai)/(h4_i*kna1*kna2);\n"
   "double h5_ss = (nass*nass)/(h4_ss*kna1*kna2);\n"
   "double h6_i = (1.0/h4_i);\n"
   "double h6_ss = (1.0/h4_ss);\n"
   "double h2_i = hna*nai/(h1_i*kna3);\n"
   "double h2_ss = hna*nass/(h1_ss*kna3);\n"
   "double h3_i = (1.0/h1_i);\n"
   "double h3_ss = (1.0/h1_ss);\n"
   "double k6_i = cai*h6_i*kcaon;\n"
   "double k6_ss = cass*h6_ss*kcaon;\n"
   "double k4p_i = h3_i*wca/hca;\n"
   "double k4p_ss = h3_ss*wca/hca;\n"
   "double k4pp_i = h2_i*wnaca;\n"
   "double k4pp_ss = h2_ss*wnaca;\n"
   "double k7_i = h2_i*h5_i*wna;\n"
   "double k7_ss = h2_ss*h5_ss*wna;\n"
   "double k4_i = k4p_i + k4pp_i;\n"
   "double k4_ss = k4p_ss + k4pp_ss;\n"
   "double x1_i = k2_i*k4_i*(k6_i + k7_i) + k5_i*k7_i*(k2_i + k3_i);\n"
   "double x1_ss = k2_ss*k4_ss*(k6_ss + k7_ss) + k5_ss*k7_ss*(k2_ss + k3_ss);\n"
   "double x2_i = k1_i*k7_i*(k4_i + k5_i) + k4_i*k6_i*(k1_i + k8_i);\n"
   "double x2_ss = k1_ss*k7_ss*(k4_ss + k5_ss) + k4_ss*k6_ss*(k1_ss + k8_ss);\n"
   "double x3_i = k1_i*k3_i*(k6_i + k7_i) + k6_i*k8_i*(k2_i + k3_i);\n"
   "double x3_ss = k1_ss*k3_ss*(k6_ss + k7_ss) + k6_ss*k8_ss*(k2_ss + k3_ss);\n"
   "double x4_i = k2_i*k8_i*(k4_i + k5_i) + k3_i*k5_i*(k1_i + k8_i);\n"
   "double x4_ss = k2_ss*k8_ss*(k4_ss + k5_ss) + k3_ss*k5_ss*(k1_ss + k8_ss);\n"
   "double E1_i = x1_i/(x1_i + x2_i + x3_i + x4_i);\n"
   "double E1_ss = x1_ss/(x1_ss + x2_ss + x3_ss + x4_ss);\n"
   "double E2_i = x2_i/(x1_i + x2_i + x3_i + x4_i);\n"
   "double E2_ss = x2_ss/(x1_ss + x2_ss + x3_ss + x4_ss);\n"
   "double E3_i = x3_i/(x1_i + x2_i + x3_i + x4_i);\n"
   "double E3_ss = x3_ss/(x1_ss + x2_ss + x3_ss + x4_ss);\n"
   "double E4_i = x4_i/(x1_i + x2_i + x3_i + x4_i);\n"
   "double E4_ss = x4_ss/(x1_ss + x2_ss + x3_ss + x4_ss);\n"
   "double JncxCa_i = -E1_i*k1_i + E2_i*k2_i;\n"
   "double JncxCa_ss = -E1_ss*k1_ss + E2_ss*k2_ss;\n"
   "double JncxNa_i = -3*E1_i*k8_i - E2_i*k3pp_i + E3_i*k4pp_i + 3*E4_i*k7_i;\n"
   "double JncxNa_ss = -3*E1_ss*k8_ss - E2_ss*k3pp_ss + E3_ss*k4pp_ss + 3*E4_ss*k7_ss;\n"
   "double INaCa_i = 0.80000000000000004*Gncx*allo_i*(JncxCa_i*zca + JncxNa_i*zna);\n"
   "double INaCa_ss = 0.20000000000000001*Gncx*allo_ss*(JncxCa_ss*zca + JncxNa_ss*zna);\n"
   "double P = eP/(H/Khp + 1 + ki/Kxkur + nai/Knap);\n"
   "double a1 = (nai*nai*nai)*k1p/((Knai*Knai*Knai)*(((1 + ki/Kki)*(1 + ki/Kki)) + ((1 + nai/Knai)*(1 + nai/Knai)*(1 + nai/Knai)) - 1));\n"
   "double b3 = H*P*k3m/(1 + MgATP/Kmgatp);\n"
   "double b4 = (ki*ki)*k4m/((Kki*Kki)*(((1 + ki/Kki)*(1 + ki/Kki)) + ((1 + nai/Knai)*(1 + nai/Knai)*(1 + nai/Knai)) - 1));\n"
   "double x1 = a1*a2*a4 + a1*a2*b3 + a2*b3*b4 + b2*b3*b4;\n"
   "double x2 = a1*a2*a3 + a2*a3*b4 + a3*b1*b4 + b1*b2*b4;\n"
   "double x3 = a2*a3*a4 + a3*a4*b1 + a4*b1*b2 + b1*b2*b3;\n"
   "double x4 = a1*a3*a4 + a1*a4*b2 + a1*b2*b3 + b2*b3*b4;\n"
   "double E1 = x1/(x1 + x2 + x3 + x4);\n"
   "double E2 = x2/(x1 + x2 + x3 + x4);\n"
   "double E3 = x3/(x1 + x2 + x3 + x4);\n"
   "double E4 = x4/(x1 + x2 + x3 + x4);\n"
   "double JnakK = -2*E3*a1 + 2*E4*b1;\n"
   "double JnakNa = 3*E1*a3 - 3*E2*b3;\n"
   "double INaK = Pnak*(JnakK*zk + JnakNa*zna);\n"
   "double fINaLp = (1.0/(1 + KmCaMK/CaMKa));\n"
   "double INaL = GNaL*mL*(-ENa + v)*(fINaLp*hLp + hL*(1 - fINaLp));\n"
   "double _expensive_functions_177 = exp(vfrt);\n"
   "double INab_A = PNab*ffrt*(_expensive_functions_177*nai - nao)/INab_B;\n"
   "double __melodee_temp_043;\n"
   "if (__melodee_temp_042)\n"
   "{\n"
   "   __melodee_temp_043 = INab_A*(1 - 0.5*INab_U);\n"
   "}\n"
   "else\n"
   "{\n"
   "   _expensive_functions_178 = exp(INab_U);\n"
   "   __melodee_temp_043 = INab_A*INab_U/(_expensive_functions_178 - 1);\n"
   "}\n"
   "double INab = __melodee_temp_043;\n"
   "double IpCa = GpCa*cai/(KmCap + cai);\n"
   "double fItop = (1.0/(1 + KmCaMK/CaMKa));\n"
   "double i = AiF*iF + AiS*iS;\n"
   "double ip = AiF*iFp + AiS*iSp;\n"
   "double Ito = Gto*(-EK + v)*(a*i*(1 - fItop) + ap*fItop*ip);\n"
   "double Jleak = 0.00026249999999999998*cansr;\n"
   "double fJupp = (1.0/(1 + KmCaMK/CaMKa));\n"
   "double Jupnp = 0.0043750000000000004*cai*upScale/(cai + 0.00092000000000000003);\n"
   "double Jupp = 0.01203125*cai*upScale/(cai + 0.00075000000000000002);\n"
   "double Jup = Jup_b*(-Jleak + Jupnp*(1 - fJupp) + Jupp*fJupp);\n"
   "double cansr_diff = -Jtr*vjsr/vnsr + Jup;\n"
   "double Bcajsr = (1.0/(csqnmax*kmcsqn/((cajsr + kmcsqn)*(cajsr + kmcsqn)) + 1));\n"
   "double Bcass = (1.0/(BSLmax*KmBSL/((KmBSL + cass)*(KmBSL + cass)) + BSRmax*KmBSR/((KmBSR + cass)*(KmBSR + cass)) + 1));\n"
   "double ki_diff = Acap*cm*(-IK1 - IKb - IKr - IKs + 2*INaK - Ito)/(F*vmyo) + JdiffK*vss/vmyo;\n"
   "double kss_diff = -Acap*ICaK*cm/(F*vss) - JdiffK;\n"
   "double nai_diff = Acap*cm*(-INa - 3*INaCa_i - 3*INaK - INaL - INab)/(F*vmyo) + JdiffNa*vss/vmyo;\n"
   "double nass_diff = Acap*cm*(-ICaNa - 3*INaCa_ss)/(F*vss) - JdiffNa;\n"
   "double Bcai = (1.0/(cmdnmax*kmcmdn/((cai + kmcmdn)*(cai + kmcmdn)) + kmtrpn*trpnmax/((cai + kmtrpn)*(cai + kmtrpn)) + 1));\n"
   "double cai_diff = Bcai*((1.0/2.0)*Acap*cm*(-ICab + 2*INaCa_i - IpCa)/(F*vmyo) + Jdiff*vss/vmyo - Jup*vnsr/vmyo);\n"
   "double fJrelp = (1.0/(1 + KmCaMK/CaMKa));\n"
   "double Jrel = Jrel_scaling_factor*(Jrelnp*(1 - fJrelp) + Jrelp*fJrelp);\n"
   "double __melodee_temp_046 = Jrel*JrelStiffConst > cajsr;\n"
   "double Jrel_001;\n"
   "if (__melodee_temp_046)\n"
   "{\n"
   "   double ryr_Jrel = cajsr/JrelStiffConst;\n"
   "   Jrel_001 = ryr_Jrel;\n"
   "}\n"
   "else\n"
   "{\n"
   "   Jrel_001 = Jrel;\n"
   "}\n"
   "double cajsr_diff = Bcajsr*(-Jrel_001 + Jtr);\n"
   "double cass_diff = Bcass*((1.0/2.0)*Acap*cm*(-ICaL + 2*INaCa_ss)/(F*vss) - Jdiff + Jrel_001*vjsr/vss);\n"
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
   "_state[_ii+C1_off*_nCells] += _dt*C1_diff;\n"
   "_state[_ii+C2_off*_nCells] += _dt*C2_diff;\n"
   "_state[_ii+CaMKt_off*_nCells] += _dt*CaMKt_diff;\n"
   "_state[_ii+Cbound_off*_nCells] += _dt*Cbound_diff;\n"
   "_state[_ii+D_off*_nCells] += _dt*D_diff;\n"
   "_state[_ii+IC1_off*_nCells] += _dt*IC1_diff;\n"
   "_state[_ii+IC2_off*_nCells] += _dt*IC2_diff;\n"
   "_state[_ii+IO_off*_nCells] += _dt*IO_diff;\n"
   "_state[_ii+IObound_off*_nCells] += _dt*IObound_diff;\n"
   "_state[_ii+O_off*_nCells] += _dt*O_diff;\n"
   "_state[_ii+Obound_off*_nCells] += _dt*Obound_diff;\n"
   "_state[_ii+cai_off*_nCells] += _dt*cai_diff;\n"
   "_state[_ii+cajsr_off*_nCells] += _dt*cajsr_diff;\n"
   "_state[_ii+cansr_off*_nCells] += _dt*cansr_diff;\n"
   "_state[_ii+cass_off*_nCells] += _dt*cass_diff;\n"
   "_state[_ii+d_off*_nCells] += _dt*d_diff;\n"
   "_state[_ii+fcaf_off*_nCells] += _dt*fcaf_diff;\n"
   "_state[_ii+fcafp_off*_nCells] += _dt*fcafp_diff;\n"
   "_state[_ii+fcas_off*_nCells] += _dt*fcas_diff;\n"
   "_state[_ii+ff_off*_nCells] += _dt*ff_diff;\n"
   "_state[_ii+ffp_off*_nCells] += _dt*ffp_diff;\n"
   "_state[_ii+fs_off*_nCells] += _dt*fs_diff;\n"
   "_state[_ii+jca_off*_nCells] += _dt*jca_diff;\n"
   "_state[_ii+ki_off*_nCells] += _dt*ki_diff;\n"
   "_state[_ii+kss_off*_nCells] += _dt*kss_diff;\n"
   "_state[_ii+nai_off*_nCells] += _dt*nai_diff;\n"
   "_state[_ii+nass_off*_nCells] += _dt*nass_diff;\n"
   "_state[_ii+nca_off*_nCells] += _dt*nca_diff;\n"
   "_state[_ii+Jrelnp_off*_nCells] += _Jrelnp_RLA*(Jrelnp+_Jrelnp_RLB);\n"
   "_state[_ii+Jrelp_off*_nCells] += _Jrelp_RLA*(Jrelp+_Jrelp_RLB);\n"
   "_state[_ii+a_off*_nCells] += _a_RLA*(a+_a_RLB);\n"
   "_state[_ii+ap_off*_nCells] += _ap_RLA*(ap+_ap_RLB);\n"
   "_state[_ii+hL_off*_nCells] += _hL_RLA*(hL+_hL_RLB);\n"
   "_state[_ii+hLp_off*_nCells] += _hLp_RLA*(hLp+_hLp_RLB);\n"
   "_state[_ii+hf_off*_nCells] += _hf_RLA*(hf+_hf_RLB);\n"
   "_state[_ii+hs_off*_nCells] += _hs_RLA*(hs+_hs_RLB);\n"
   "_state[_ii+hsp_off*_nCells] += _hsp_RLA*(hsp+_hsp_RLB);\n"
   "_state[_ii+iF_off*_nCells] += _iF_RLA*(iF+_iF_RLB);\n"
   "_state[_ii+iFp_off*_nCells] += _iFp_RLA*(iFp+_iFp_RLB);\n"
   "_state[_ii+iS_off*_nCells] += _iS_RLA*(iS+_iS_RLB);\n"
   "_state[_ii+iSp_off*_nCells] += _iSp_RLA*(iSp+_iSp_RLB);\n"
   "_state[_ii+j_off*_nCells] += _j_RLA*(j+_j_RLB);\n"
   "_state[_ii+jp_off*_nCells] += _jp_RLA*(jp+_jp_RLB);\n"
   "_state[_ii+m_off*_nCells] += _m_RLA*(m+_m_RLB);\n"
   "_state[_ii+mL_off*_nCells] += _mL_RLA*(mL+_mL_RLB);\n"
   "_state[_ii+xk1_off*_nCells] += _xk1_RLA*(xk1+_xk1_RLB);\n"
   "_state[_ii+xs1_off*_nCells] += _xs1_RLA*(xs1+_xs1_RLB);\n"
   "_state[_ii+xs2_off*_nCells] += _xs2_RLA*(xs2+_xs2_RLB);\n"
   "_dVm[_indexArray[_ii]] = -Iion_001;\n"
   "}\n";

   _program_code = ss.str();
   //cout << ss.str();
   nvrtcCreateProgram(&_program,
                      _program_code.c_str(),
                      "ohara_cipa_SCC_program",
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
   cuModuleGetFunction(&_kernel, _module, "ohara_cipa_SCC_kernel");
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
   _C1_off,
   _C2_off,
   _CaMKt_off,
   _Cbound_off,
   _D_off,
   _IC1_off,
   _IC2_off,
   _IO_off,
   _IObound_off,
   _Jrelnp_off,
   _Jrelp_off,
   _O_off,
   _Obound_off,
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
   _hL_off,
   _hLp_off,
   _hf_off,
   _hs_off,
   _hsp_off,
   _iF_off,
   _iFp_off,
   _iS_off,
   _iSp_off,
   _j_off,
   _jca_off,
   _jp_off,
   _ki_off,
   _kss_off,
   _m_off,
   _mL_off,
   _nai_off,
   _nass_off,
   _nca_off,
   _xk1_off,
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
   double cao = 1.8;
   double ko = 5.4000000000000004;
   double nao = 140;
   double L = 0.01;
   double rad = 0.0011000000000000001;
   double Ageo = 6.2800000000000002*L*rad + 6.2800000000000002*(rad*rad);
   double vcell = 3140.0*(rad*rad)*L;
   double Acap = 2*Ageo;
   double vjsr = 0.0047999999999999996*vcell;
   double vmyo = 0.68000000000000005*vcell;
   double vnsr = 0.055199999999999999*vcell;
   double vss = 0.02*vcell;
   double F = 96485;
   double R = 8314;
   double T = 310;
   double zca = 2;
   double zk = 1;
   double zna = 1;
   double PKNa = 0.018329999999999999;
   double cm = 1;
   double frt = F/(R*T);
   double ffrt = F*frt;
   double KmCaMK = 0.14999999999999999;
   double CaMKo = 0.050000000000000003;
   double KmCaM = 0.0015;
   double aCaMK = 0.050000000000000003;
   double bCaMK = 0.00068000000000000005;
   double Aff = 0.59999999999999998;
   double B_1 = 2*frt;
   double B_2 = frt;
   double B_3 = frt;
   double Kmn = 0.002;
   double PCa_b = 0.00010069999999999999;
   double k2n = 1000;
   double tjca = 75;
   double v0 = 0;
   double Afs = 1 - Aff;
   double __melodee_temp_030 = celltype == 2;
   double __melodee_temp_031;
   if (__melodee_temp_030)
   {
      __melodee_temp_031 = 2.5*PCa_b;
   }
   else
   {
      __melodee_temp_031 = PCa_b;
   }
   double __melodee_temp_032 = celltype == 1;
   double __melodee_temp_033;
   if (__melodee_temp_032)
   {
      __melodee_temp_033 = 1.2*PCa_b;
   }
   else
   {
      __melodee_temp_033 = __melodee_temp_031;
   }
   double PCa = __melodee_temp_033;
   double PCaK = 0.00035740000000000001*PCa;
   double PCaNa = 0.00125*PCa;
   double PCap = 1.1000000000000001*PCa;
   double PCaKp = 0.00035740000000000001*PCap;
   double PCaNap = 0.00125*PCap;
   double B = 2*frt;
   double PCab = 2.4999999999999999e-8;
   double ICab_v0 = 0;
   double GK1_b = 0.3239783999999998;
   double __melodee_temp_000 = celltype == 2;
   double __melodee_temp_001;
   if (__melodee_temp_000)
   {
      __melodee_temp_001 = 1.3*GK1_b;
   }
   else
   {
      __melodee_temp_001 = GK1_b;
   }
   double __melodee_temp_002 = celltype == 1;
   double __melodee_temp_003;
   if (__melodee_temp_002)
   {
      __melodee_temp_003 = 1.2*GK1_b;
   }
   else
   {
      __melodee_temp_003 = __melodee_temp_001;
   }
   double GK1 = __melodee_temp_003;
   double GKb_b = 0.0030000000000000001;
   double __melodee_temp_004 = celltype == 1;
   double __melodee_temp_005;
   if (__melodee_temp_004)
   {
      __melodee_temp_005 = 0.59999999999999998*GKb_b;
   }
   else
   {
      __melodee_temp_005 = GKb_b;
   }
   double GKb = __melodee_temp_005;
   double A1 = 0.0264;
   double A11 = 0.00078680000000000004;
   double A2 = 4.9860000000000002e-6;
   double A21 = 5.4550000000000003e-6;
   double A3 = 0.001214;
   double A31 = 0.005509;
   double A4 = 1.8539999999999999e-5;
   double A41 = 0.0014159999999999999;
   double A51 = 0.44919999999999999;
   double A52 = 0.31809999999999999;
   double A53 = 0.14899999999999999;
   double A61 = 0.012409999999999999;
   double A62 = 0.3226;
   double A63 = 0.0089779999999999999;
   double B1 = 4.6310000000000002e-5;
   double B11 = 1.5349999999999998e-8;
   double B2 = -0.0042259999999999997;
   double B21 = -0.16880000000000001;
   double B3 = 0.0085159999999999993;
   double B31 = 7.7710000000000006e-9;
   double B4 = -0.04641;
   double B41 = -0.02877;
   double B51 = 0.0085950000000000002;
   double B52 = 3.613e-8;
   double B53 = 0.0046680000000000003;
   double B61 = 0.17249999999999999;
   double B62 = -0.00065749999999999999;
   double B63 = -0.02215;
   double D_diff = 0;
   double GKr_b = 0.04658545454545456;
   double Kmax = 0;
   double Kt = 0;
   double Ku = 0;
   double Temp = 37;
   double Vhalf = 1;
   double halfmax = 1;
   double n = 1;
   double q1 = 4.843;
   double q11 = 4.9420000000000002;
   double q2 = 4.2300000000000004;
   double q21 = 4.1559999999999997;
   double q3 = 4.9619999999999997;
   double q31 = 4.2199999999999998;
   double q4 = 3.7690000000000001;
   double q41 = 1.4590000000000001;
   double q51 = 5;
   double q52 = 4.6630000000000003;
   double q53 = 2.4119999999999999;
   double q61 = 5.5679999999999996;
   double q62 = 5;
   double q63 = 5.6820000000000004;
   double __melodee_temp_006 = celltype == 2;
   double __melodee_temp_007;
   if (__melodee_temp_006)
   {
      __melodee_temp_007 = 0.80000000000000004*GKr_b;
   }
   else
   {
      __melodee_temp_007 = GKr_b;
   }
   double __melodee_temp_008 = celltype == 1;
   double __melodee_temp_009;
   if (__melodee_temp_008)
   {
      __melodee_temp_009 = 1.3*GKr_b;
   }
   else
   {
      __melodee_temp_009 = __melodee_temp_007;
   }
   double GKr = __melodee_temp_009;
   double _expensive_functions_025 = log(q63);
   double _expensive_functions_026 = exp(-1.0/10.0*_expensive_functions_025*(Temp - 20));
   double _expensive_functions_027 = log(q53);
   double _expensive_functions_028 = exp((1.0/10.0)*_expensive_functions_027*(Temp - 20));
   double GKs_b = 0.006358000000000001;
   double txs1_max = 817.29999999999995;
   double __melodee_temp_010 = celltype == 1;
   double __melodee_temp_011;
   if (__melodee_temp_010)
   {
      __melodee_temp_011 = 1.3999999999999999*GKs_b;
   }
   else
   {
      __melodee_temp_011 = GKs_b;
   }
   double GKs = __melodee_temp_011;
   double Ahf = 0.98999999999999999;
   double GNa = 75;
   double hssV1 = 82.900000000000006;
   double hssV2 = 6.0860000000000003;
   double mssV1 = 39.57;
   double mssV2 = 9.8710000000000004;
   double mtD1 = 6.7649999999999997;
   double mtD2 = 8.5519999999999996;
   double mtV1 = 11.640000000000001;
   double mtV2 = 34.770000000000003;
   double mtV3 = 77.420000000000002;
   double mtV4 = 5.9550000000000001;
   double shift_INa_inact = 0;
   double Ahs = 1 - Ahf;
   double Gncx_b = 0.00080000000000000004;
   double KmCaAct = 0.00014999999999999999;
   double kasymm = 12.5;
   double kcaoff = 5000.0;
   double kcaon = 1500000.0;
   double kna1 = 15;
   double kna2 = 5;
   double kna3 = 88.120000000000005;
   double qca = 0.16700000000000001;
   double qna = 0.52239999999999998;
   double wca = 60000.0;
   double wna = 60000.0;
   double wnaca = 5000.0;
   double __melodee_temp_012 = celltype == 2;
   double __melodee_temp_013;
   if (__melodee_temp_012)
   {
      __melodee_temp_013 = 1.3999999999999999*Gncx_b;
   }
   else
   {
      __melodee_temp_013 = Gncx_b;
   }
   double __melodee_temp_014 = celltype == 1;
   double __melodee_temp_015;
   if (__melodee_temp_014)
   {
      __melodee_temp_015 = 1.1000000000000001*Gncx_b;
   }
   else
   {
      __melodee_temp_015 = __melodee_temp_013;
   }
   double Gncx = __melodee_temp_015;
   double h10_i = kasymm + 1 + nao*(1 + nao/kna2)/kna1;
   double h10_ss = kasymm + 1 + nao*(1 + nao/kna2)/kna1;
   double k2_i = kcaoff;
   double k2_ss = kcaoff;
   double k5_i = kcaoff;
   double k5_ss = kcaoff;
   double h11_i = (nao*nao)/(h10_i*kna1*kna2);
   double h11_ss = (nao*nao)/(h10_ss*kna1*kna2);
   double h12_i = (1.0/h10_i);
   double h12_ss = (1.0/h10_ss);
   double k1_i = cao*h12_i*kcaon;
   double k1_ss = cao*h12_ss*kcaon;
   double H = 9.9999999999999995e-8;
   double Khp = 1.698e-7;
   double Kki = 0.5;
   double Kko = 0.35820000000000002;
   double Kmgatp = 1.698e-7;
   double Knai0 = 9.0730000000000004;
   double Knao0 = 27.780000000000001;
   double Knap = 224;
   double Kxkur = 292;
   double MgADP = 0.050000000000000003;
   double MgATP = 9.8000000000000007;
   double Pnak_b = 30;
   double delta = -0.155;
   double eP = 4.2000000000000002;
   double k1m = 182.40000000000001;
   double k1p = 949.5;
   double k2m = 39.399999999999999;
   double k2p = 687.20000000000005;
   double k3m = 79300;
   double k3p = 1899;
   double k4m = 40;
   double k4p = 639;
   double __melodee_temp_016 = celltype == 2;
   double __melodee_temp_017;
   if (__melodee_temp_016)
   {
      __melodee_temp_017 = 0.69999999999999996*Pnak_b;
   }
   else
   {
      __melodee_temp_017 = Pnak_b;
   }
   double __melodee_temp_018 = celltype == 1;
   double __melodee_temp_019;
   if (__melodee_temp_018)
   {
      __melodee_temp_019 = 0.90000000000000002*Pnak_b;
   }
   else
   {
      __melodee_temp_019 = __melodee_temp_017;
   }
   double Pnak = __melodee_temp_019;
   double a2 = k2p;
   double a4 = MgATP*k4p/(Kmgatp*(1 + MgATP/Kmgatp));
   double b1 = MgADP*k1m;
   double GNaL_b = 0.019957499999999975;
   double thL = 200;
   double __melodee_temp_028 = celltype == 1;
   double __melodee_temp_029;
   if (__melodee_temp_028)
   {
      __melodee_temp_029 = 0.59999999999999998*GNaL_b;
   }
   else
   {
      __melodee_temp_029 = GNaL_b;
   }
   double GNaL = __melodee_temp_029;
   double thLp = 3*thL;
   double INab_B = frt;
   double PNab = 3.75e-10;
   double INab_v0 = 0;
   double GpCa = 0.00050000000000000001;
   double KmCap = 0.00050000000000000001;
   double Gto_b = 0.02;
   double __melodee_temp_020 = celltype == 1;
   double __melodee_temp_021;
   if (__melodee_temp_020)
   {
   }
   else
   {
      __melodee_temp_021 = 1;
   }
   double __melodee_temp_022 = celltype == 2;
   double __melodee_temp_023;
   if (__melodee_temp_022)
   {
      __melodee_temp_023 = 4*Gto_b;
   }
   else
   {
      __melodee_temp_023 = Gto_b;
   }
   double __melodee_temp_024 = celltype == 1;
   double __melodee_temp_025;
   if (__melodee_temp_024)
   {
      __melodee_temp_025 = 4*Gto_b;
   }
   else
   {
      __melodee_temp_025 = __melodee_temp_023;
   }
   double Gto = __melodee_temp_025;
   double Jup_b = 1.0;
   double __melodee_temp_026 = celltype == 1;
   double __melodee_temp_027;
   if (__melodee_temp_026)
   {
      __melodee_temp_027 = 1.3;
   }
   else
   {
      __melodee_temp_027 = 1;
   }
   double upScale = __melodee_temp_027;
   double BSLmax = 1.1240000000000001;
   double BSRmax = 0.047;
   double KmBSL = 0.0086999999999999994;
   double KmBSR = 0.00087000000000000001;
   double cmdnmax_b = 0.050000000000000003;
   double csqnmax = 10;
   double kmcmdn = 0.0023800000000000002;
   double kmcsqn = 0.80000000000000004;
   double kmtrpn = 0.00050000000000000001;
   double trpnmax = 0.070000000000000007;
   double __melodee_temp_044 = celltype == 1;
   double __melodee_temp_045;
   if (__melodee_temp_044)
   {
      __melodee_temp_045 = 1.3*cmdnmax_b;
   }
   else
   {
      __melodee_temp_045 = cmdnmax_b;
   }
   double cmdnmax = __melodee_temp_045;
   double Jrel_scaling_factor = 1.0;
   double bt = 4.75;
   double a_rel = 0.5*bt;
   double btp = 1.25*bt;
   double a_relp = 0.5*btp;
   double __melodee_temp_049 = celltype == 2;
   double __melodee_temp_053 = celltype == 2;
   double _expensive_functions_075 = sqrt(ko);
   double _expensive_functions_077 = log(q2);
   double _expensive_functions_078 = exp((1.0/10.0)*_expensive_functions_077*(Temp - 20));
   double _expensive_functions_080 = log(q61);
   double _expensive_functions_081 = exp((1.0/10.0)*_expensive_functions_080*(Temp - 20));
   double _expensive_functions_083 = log(q1);
   double _expensive_functions_084 = exp((1.0/10.0)*_expensive_functions_083*(Temp - 20));
   double _expensive_functions_086 = log(q51);
   double _expensive_functions_087 = exp((1.0/10.0)*_expensive_functions_086*(Temp - 20));
   double _expensive_functions_089 = log(q1);
   double _expensive_functions_090 = exp((1.0/10.0)*_expensive_functions_089*(Temp - 20));
   double _expensive_functions_092 = log(q41);
   double _expensive_functions_093 = exp((1.0/10.0)*_expensive_functions_092*(Temp - 20));
   double _expensive_functions_095 = log(q62);
   double _expensive_functions_096 = exp((1.0/10.0)*_expensive_functions_095*(Temp - 20));
   double _expensive_functions_098 = log(q2);
   double _expensive_functions_099 = exp((1.0/10.0)*_expensive_functions_098*(Temp - 20));
   double _expensive_functions_101 = log(q31);
   double _expensive_functions_102 = exp((1.0/10.0)*_expensive_functions_101*(Temp - 20));
   double _expensive_functions_104 = log(q52);
   double _expensive_functions_105 = exp((1.0/10.0)*_expensive_functions_104*(Temp - 20));
   double _expensive_functions_108 = log(q21);
   double _expensive_functions_109 = exp((1.0/10.0)*_expensive_functions_108*(Temp - 20));
   double _expensive_functions_111 = log(q51);
   double _expensive_functions_112 = exp((1.0/10.0)*_expensive_functions_111*(Temp - 20));
   double _expensive_functions_114 = log(q11);
   double _expensive_functions_115 = exp((1.0/10.0)*_expensive_functions_114*(Temp - 20));
   double _expensive_functions_117 = log(q61);
   double _expensive_functions_118 = exp((1.0/10.0)*_expensive_functions_117*(Temp - 20));
   double _expensive_functions_120 = log(q11);
   double _expensive_functions_121 = exp((1.0/10.0)*_expensive_functions_120*(Temp - 20));
   double _expensive_functions_123 = log(q4);
   double _expensive_functions_124 = exp((1.0/10.0)*_expensive_functions_123*(Temp - 20));
   double _expensive_functions_126 = log(q52);
   double _expensive_functions_127 = exp((1.0/10.0)*_expensive_functions_126*(Temp - 20));
   double _expensive_functions_129 = log(q21);
   double _expensive_functions_130 = exp((1.0/10.0)*_expensive_functions_129*(Temp - 20));
   double _expensive_functions_132 = log(q3);
   double _expensive_functions_133 = exp((1.0/10.0)*_expensive_functions_132*(Temp - 20));
   double _expensive_functions_135 = log(q62);
   double _expensive_functions_136 = exp((1.0/10.0)*_expensive_functions_135*(Temp - 20));
   double _expensive_functions_138 = log(q3);
   double _expensive_functions_139 = exp((1.0/10.0)*_expensive_functions_138*(Temp - 20));
   double _expensive_functions_141 = log(q53);
   double _expensive_functions_142 = exp((1.0/10.0)*_expensive_functions_141*(Temp - 20));
   double _expensive_functions_144 = log(q4);
   double _expensive_functions_145 = exp((1.0/10.0)*_expensive_functions_144*(Temp - 20));
   double _expensive_functions_147 = log(q63);
   double _expensive_functions_148 = exp((1.0/10.0)*_expensive_functions_147*(Temp - 20));
   double _expensive_functions_155 = log(q63);
   double _expensive_functions_156 = exp(-1.0/10.0*_expensive_functions_155*(Temp - 20));
   double _expensive_functions_157 = log(q53);
   double _expensive_functions_158 = exp((1.0/10.0)*_expensive_functions_157*(Temp - 20));
   double _expensive_functions_160 = log(q31);
   double _expensive_functions_161 = exp((1.0/10.0)*_expensive_functions_160*(Temp - 20));
   double _expensive_functions_163 = log(q63);
   double _expensive_functions_164 = exp((1.0/10.0)*_expensive_functions_163*(Temp - 20));
   double _expensive_functions_166 = log(q41);
   double _expensive_functions_167 = exp((1.0/10.0)*_expensive_functions_166*(Temp - 20));
   double _expensive_functions_169 = log(q53);
   double _expensive_functions_170 = exp((1.0/10.0)*_expensive_functions_169*(Temp - 20));
   double _expensive_functions_175 = sqrt(ko);
   double _expensive_functions_182 = exp(-_dt/thL);
   double _hL_RLA = _expensive_functions_182 - 1;
   double _expensive_functions_183 = exp(-_dt/thLp);
   double _hLp_RLA = _expensive_functions_183 - 1;
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
      real C1=load(state_[__jj].C1);
      real C2=load(state_[__jj].C2);
      real CaMKt=load(state_[__jj].CaMKt);
      real Cbound=load(state_[__jj].Cbound);
      real D=load(state_[__jj].D);
      real IC1=load(state_[__jj].IC1);
      real IC2=load(state_[__jj].IC2);
      real IO=load(state_[__jj].IO);
      real IObound=load(state_[__jj].IObound);
      real Jrelnp=load(state_[__jj].Jrelnp);
      real Jrelp=load(state_[__jj].Jrelp);
      real O=load(state_[__jj].O);
      real Obound=load(state_[__jj].Obound);
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
      real hL=load(state_[__jj].hL);
      real hLp=load(state_[__jj].hLp);
      real hf=load(state_[__jj].hf);
      real hs=load(state_[__jj].hs);
      real hsp=load(state_[__jj].hsp);
      real iF=load(state_[__jj].iF);
      real iFp=load(state_[__jj].iFp);
      real iS=load(state_[__jj].iS);
      real iSp=load(state_[__jj].iSp);
      real j=load(state_[__jj].j);
      real jca=load(state_[__jj].jca);
      real jp=load(state_[__jj].jp);
      real ki=load(state_[__jj].ki);
      real kss=load(state_[__jj].kss);
      real m=load(state_[__jj].m);
      real mL=load(state_[__jj].mL);
      real nai=load(state_[__jj].nai);
      real nass=load(state_[__jj].nass);
      real nca=load(state_[__jj].nca);
      real xk1=load(state_[__jj].xk1);
      real xs1=load(state_[__jj].xs1);
      real xs2=load(state_[__jj].xs2);
      //get the gate updates (diagonalized exponential integrator)
      real v = V;
      real vfrt = frt*v;
      real _expensive_functions = exp((1.0/10.0)*v - 1);
      real Afcaf = 0.29999999999999999 + 0.59999999999999998/(_expensive_functions + 1);
      real Afcas = 1 - Afcaf;
      real U_1 = B_1*(v - v0);
      real __melodee_temp_036 = -9.9999999999999995e-8 >= U_1 && U_1 >= 9.9999999999999995e-8;
      real f = Aff*ff + Afs*fs;
      real fp = Aff*ffp + Afs*fs;
      real _expensive_functions_014 = exp(-0.049115913555992145*v);
      real _expensive_functions_015 = exp(0.014423770373575654*v);
      real txk1 = 122.2/(0.0019352007631390235*_expensive_functions_014 + 30.433647575249029*_expensive_functions_015);
      real _expensive_functions_016 = exp((-2.5537999999999998*ko - v - 144.59)/(1.5691999999999999*ko + 3.8115000000000001));
      real xk1ss = (1.0/(_expensive_functions_016 + 1));
      real _expensive_functions_034 = exp(-1.0/31.0*v);
      real _expensive_functions_035 = exp((1.0/20.0)*v - 5.0/2.0);
      real txs2 = (1.0/(0.0022561357010639103*_expensive_functions_034 + 0.01*_expensive_functions_035));
      real _expensive_functions_036 = exp(-0.11195700850873264*v);
      real xs1ss = (1.0/(0.27288596035656526*_expensive_functions_036 + 1));
      real _expensive_functions_037 = exp(0.056179775280898875*v);
      real _expensive_functions_038 = exp(-1.0/230.0*v - 21.0/23.0);
      real txs1 = txs1_max + (1.0/(0.0035040677630748581*_expensive_functions_037 + 0.001292*_expensive_functions_038));
      real xs2ss = xs1ss;
      real _expensive_functions_039 = exp((hssV1 - shift_INa_inact + v)/hssV2);
      real hss = (1.0/(_expensive_functions_039 + 1));
      real _expensive_functions_040 = exp(-0.16431153466973381*shift_INa_inact + 0.16431153466973381*v);
      real hssp = (1.0/(2281075.8166971938*_expensive_functions_040 + 1));
      real _expensive_functions_041 = exp((-mssV1 - v)/mssV2);
      real mss = (1.0/(_expensive_functions_041 + 1));
      real _expensive_functions_042 = exp(0.15910898965791567*shift_INa_inact - 0.15910898965791567*v);
      real _expensive_functions_043 = exp(-0.049333991119881598*shift_INa_inact + 0.049333991119881598*v);
      real thf = (1.0/(1.183856958289087e-5*_expensive_functions_042 + 6.3055491858172754*_expensive_functions_043));
      real _expensive_functions_044 = exp(0.035650623885918005*shift_INa_inact - 0.035650623885918005*v);
      real _expensive_functions_045 = exp(-0.017649135192375574*shift_INa_inact + 0.017649135192375574*v);
      real ths = (1.0/(0.0051646702353817919*_expensive_functions_044 + 0.36987619372096325*_expensive_functions_045));
      real _expensive_functions_046 = exp(-0.02600780234070221*shift_INa_inact + 0.02600780234070221*v);
      real _expensive_functions_047 = exp(0.12075836251660427*shift_INa_inact - 0.12075836251660427*v);
      real tj = 2.0379999999999998 + (1.0/(0.31319363947387729*_expensive_functions_046 + 1.1315282095590072e-7*_expensive_functions_047));
      real _expensive_functions_048 = exp((mtV1 + v)/mtV2);
      real _expensive_functions_049 = exp((-mtV3 - v)/mtV4);
      real tm = (1.0/(_expensive_functions_048*mtD1 + _expensive_functions_049*mtD2));
      real jss = hss;
      real thsp = 3*ths;
      real tjp = 1.46*tj;
      real _expensive_functions_052 = exp(0.13354700854700854*v);
      real hLss = (1.0/(120578.15595522427*_expensive_functions_052 + 1));
      real _expensive_functions_053 = exp(0.13354700854700854*v);
      real hLssp = (1.0/(275969.29038698709*_expensive_functions_053 + 1));
      real _expensive_functions_054 = exp(-0.18996960486322187*v);
      real mLss = (1.0/(0.00029157958563553099*_expensive_functions_054 + 1));
      real tmL = tm;
      real _expensive_functions_056 = exp(-0.067476383265856948*v);
      real ass = (1.0/(2.6316508161673635*_expensive_functions_056 + 1));
      real _expensive_functions_057 = exp(-0.067476383265856948*v);
      real assp = (1.0/(5.1674284622306663*_expensive_functions_057 + 1));
      if (__melodee_temp_020)
      {
         real _expensive_functions_058 = exp((1.0/5.0)*v + 14);
         __melodee_temp_021 = 1 - 0.94999999999999996/(_expensive_functions_058 + 1);
      }
      else
      {
      }
      real delta_epi = __melodee_temp_021;
      real _expensive_functions_058 = exp(0.062932662051604776*v);
      real _expensive_functions_059 = exp(-4.6425255338904359*v);
      real dti_develop = 1.3540000000000001 + 0.0001/(2.6591269045230603e-5*_expensive_functions_058 + 4.5541779737128264e+24*_expensive_functions_059);
      real _expensive_functions_060 = exp((1.0/20.0)*v + 7.0/2.0);
      real dti_recover = 1 - 0.5/(_expensive_functions_060 + 1);
      real _expensive_functions_061 = exp(0.17510068289266328*v);
      real iss = (1.0/(2194.970764538301*_expensive_functions_061 + 1));
      real _expensive_functions_062 = exp(-0.034035137876343539*v);
      real _expensive_functions_063 = exp(0.034035137876343539*v);
      real ta = 1.0515000000000001/(3.5/(30.069572727397507*_expensive_functions_063 + 1) + (1.0/(2.2621017070578837*_expensive_functions_062 + 1.2089000000000001)));
      real _expensive_functions_064 = exp(-1.0/100.0*v - 1);
      real _expensive_functions_065 = exp(0.060277275467148887*v);
      real tiF_b = 4.5620000000000003 + (1.0/(0.39329999999999998*_expensive_functions_064 + 1.6300896349780942*_expensive_functions_065));
      real _expensive_functions_066 = exp(-0.016934801016088061*v);
      real _expensive_functions_067 = exp(0.12377769525931426*v);
      real tiS_b = 23.620000000000001 + (1.0/(0.00027617763953377436*_expensive_functions_066 + 0.024208962804604526*_expensive_functions_067));
      real tiF = delta_epi*tiF_b;
      real tiS = delta_epi*tiS_b;
      real tiFp = dti_develop*dti_recover*tiF;
      real tiSp = dti_develop*dti_recover*tiS;
      real tau_rel_temp = bt/(1 + 0.0123/cajsr);
      real __melodee_temp_047 = tau_rel_temp < 0.001;
      real __melodee_temp_048;
      if (__melodee_temp_047)
      {
         __melodee_temp_048 = 0.001;
      }
      else
      {
         __melodee_temp_048 = tau_rel_temp;
      }
      real tau_rel = __melodee_temp_048;
      real tau_relp_temp = btp/(1 + 0.0123/cajsr);
      real __melodee_temp_051 = tau_relp_temp < 0.001;
      real __melodee_temp_052;
      if (__melodee_temp_051)
      {
         __melodee_temp_052 = 0.001;
      }
      else
      {
         __melodee_temp_052 = tau_relp_temp;
      }
      real tau_relp = __melodee_temp_052;
      real CaMKb = CaMKo*(1 - CaMKt)/(KmCaM/cass + 1);
      real CaMKa = CaMKb + CaMKt;
      real fICaLp = (1.0/(1 + KmCaMK/CaMKa));
      real _expensive_functions_071 = exp(2*vfrt);
      real A_1 = 4*ffrt*(_expensive_functions_071*cass - 0.34100000000000003*cao)/B_1;
      real __melodee_temp_037;
      if (__melodee_temp_036)
      {
         __melodee_temp_037 = A_1*(1 - 0.5*U_1);
      }
      else
      {
         real _expensive_functions_074 = exp(U_1);
         __melodee_temp_037 = A_1*U_1/(_expensive_functions_074 - 1);
      }
      real PhiCaL = __melodee_temp_037;
      real fca = Afcaf*fcaf + Afcas*fcas;
      real fcap = Afcaf*fcafp + Afcas*fcas;
      real ICaL = PCa*PhiCaL*d*(1 - fICaLp)*(f*(1 - nca) + fca*jca*nca) + PCap*PhiCaL*d*fICaLp*(fcap*jca*nca + fp*(1 - nca));
      real Jrel_inf_temp = -ICaL*a_rel/(1 + 25.62890625/(cajsr*cajsr*cajsr*cajsr*cajsr*cajsr*cajsr*cajsr));
      real __melodee_temp_050;
      if (__melodee_temp_049)
      {
         __melodee_temp_050 = 1.7*Jrel_inf_temp;
      }
      else
      {
         __melodee_temp_050 = Jrel_inf_temp;
      }
      real Jrel_inf = __melodee_temp_050;
      real Jrel_temp = -ICaL*a_relp/(1 + 25.62890625/(cajsr*cajsr*cajsr*cajsr*cajsr*cajsr*cajsr*cajsr));
      real __melodee_temp_054;
      if (__melodee_temp_053)
      {
         __melodee_temp_054 = 1.7*Jrel_temp;
      }
      else
      {
         __melodee_temp_054 = Jrel_temp;
      }
      real Jrel_infp = __melodee_temp_054;
      real _expensive_functions_178 = exp(-_dt/tau_rel);
      real _Jrelnp_RLA = _expensive_functions_178 - 1;
      real _Jrelnp_RLB = -Jrel_inf;
      real _expensive_functions_179 = exp(-_dt/tau_relp);
      real _Jrelp_RLA = _expensive_functions_179 - 1;
      real _Jrelp_RLB = -Jrel_infp;
      real _expensive_functions_180 = exp(-_dt/ta);
      real _a_RLA = _expensive_functions_180 - 1;
      real _a_RLB = -ass;
      real _expensive_functions_181 = exp(-_dt/ta);
      real _ap_RLA = _expensive_functions_181 - 1;
      real _ap_RLB = -assp;
      real _hL_RLB = -hLss;
      real _hLp_RLB = -hLssp;
      real _expensive_functions_184 = exp(-_dt/thf);
      real _hf_RLA = _expensive_functions_184 - 1;
      real _hf_RLB = -hss;
      real _expensive_functions_185 = exp(-_dt/ths);
      real _hs_RLA = _expensive_functions_185 - 1;
      real _hs_RLB = -hss;
      real _expensive_functions_186 = exp(-_dt/thsp);
      real _hsp_RLA = _expensive_functions_186 - 1;
      real _hsp_RLB = -hssp;
      real _expensive_functions_187 = exp(-_dt/tiF);
      real _iF_RLA = _expensive_functions_187 - 1;
      real _iF_RLB = -iss;
      real _expensive_functions_188 = exp(-_dt/tiFp);
      real _iFp_RLA = _expensive_functions_188 - 1;
      real _iFp_RLB = -iss;
      real _expensive_functions_189 = exp(-_dt/tiS);
      real _iS_RLA = _expensive_functions_189 - 1;
      real _iS_RLB = -iss;
      real _expensive_functions_190 = exp(-_dt/tiSp);
      real _iSp_RLA = _expensive_functions_190 - 1;
      real _iSp_RLB = -iss;
      real _expensive_functions_191 = exp(-_dt/tj);
      real _j_RLA = _expensive_functions_191 - 1;
      real _j_RLB = -jss;
      real _expensive_functions_192 = exp(-_dt/tjp);
      real _jp_RLA = _expensive_functions_192 - 1;
      real _jp_RLB = -jss;
      real _expensive_functions_193 = exp(-_dt/tm);
      real _m_RLA = _expensive_functions_193 - 1;
      real _m_RLB = -mss;
      real _expensive_functions_194 = exp(-_dt/tmL);
      real _mL_RLA = _expensive_functions_194 - 1;
      real _mL_RLB = -mLss;
      real _expensive_functions_195 = exp(-_dt/txk1);
      real _xk1_RLA = _expensive_functions_195 - 1;
      real _xk1_RLB = -xk1ss;
      real _expensive_functions_196 = exp(-_dt/txs1);
      real _xs1_RLA = _expensive_functions_196 - 1;
      real _xs1_RLB = -xs1ss;
      real _expensive_functions_197 = exp(-_dt/txs2);
      real _xs2_RLA = _expensive_functions_197 - 1;
      real _xs2_RLB = -xs2ss;
      //get the other differential updates
      real _expensive_functions_001 = exp(-0.23640661938534277*v);
      real dss = (1.0/(0.39398514226669484*_expensive_functions_001 + 1));
      real _expensive_functions_002 = exp(0.27056277056277056*v);
      real fss = (1.0/(199.86038496778565*_expensive_functions_002 + 1));
      real _expensive_functions_003 = exp(0.089999999999999997*v);
      real _expensive_functions_004 = exp(-0.050000000000000003*v);
      real td = 0.59999999999999998 + (1.0/(3.5254214873653824*_expensive_functions_003 + 0.74081822068171788*_expensive_functions_004));
      real _expensive_functions_005 = exp((1.0/7.0)*v - 4.0/7.0);
      real _expensive_functions_006 = exp(4.0/7.0 - 1.0/7.0*v);
      real tfcaf = 7 + (1.0/(0.040000000000000001*_expensive_functions_005 + 0.040000000000000001*_expensive_functions_006));
      real _expensive_functions_007 = exp(-1.0/3.0*v);
      real _expensive_functions_008 = exp((1.0/7.0)*v);
      real tfcas = 100 + (1.0/(0.00012*_expensive_functions_007 + 0.00012*_expensive_functions_008));
      real _expensive_functions_009 = exp(-1.0/10.0*v - 2);
      real _expensive_functions_010 = exp((1.0/10.0)*v + 2);
      real tff = 7 + (1.0/(0.0044999999999999997*_expensive_functions_009 + 0.0044999999999999997*_expensive_functions_010));
      real _expensive_functions_011 = exp(-1.0/4.0*v - 5.0/4.0);
      real _expensive_functions_012 = exp((1.0/6.0)*v + 5.0/6.0);
      real tfs = 1000 + (1.0/(3.4999999999999997e-5*_expensive_functions_011 + 3.4999999999999997e-5*_expensive_functions_012));
      real U_2 = B_2*(v - v0);
      real U_3 = B_3*(v - v0);
      real d_diff = (-d + dss)/td;
      real fcass = fss;
      real ff_diff = (-ff + fss)/tff;
      real fs_diff = (-fs + fss)/tfs;
      real tfcafp = 2.5*tfcaf;
      real tffp = 2.5*tff;
      real __melodee_temp_034 = -9.9999999999999995e-8 >= U_3 && U_3 >= 9.9999999999999995e-8;
      real __melodee_temp_038 = -9.9999999999999995e-8 >= U_2 && U_2 >= 9.9999999999999995e-8;
      real fcaf_diff = (-fcaf + fcass)/tfcaf;
      real fcafp_diff = (-fcafp + fcass)/tfcafp;
      real fcas_diff = (-fcas + fcass)/tfcas;
      real ffp_diff = (-ffp + fss)/tffp;
      real jca_diff = (fcass - jca)/tjca;
      real U = B*(-ICab_v0 + v);
      real __melodee_temp_040 = -9.9999999999999995e-8 >= U && U >= 9.9999999999999995e-8;
      real _expensive_functions_013 = exp(-0.27388602127883704*ko + 0.10534077741493732*v);
      real rk1 = (1.0/(69220.632210676704*_expensive_functions_013 + 1));
      real _expensive_functions_017 = exp(-0.054525627044711013*v);
      real xkb = (1.0/(2.2023634509492389*_expensive_functions_017 + 1));
      real _expensive_functions_018 = exp(0.14729709824716453*Vhalf - 0.14729709824716453*v);
      real _expensive_functions_019 = log(D);
      real _expensive_functions_020 = exp(_expensive_functions_019*n);
      real _expensive_functions_021 = log(D);
      real _expensive_functions_022 = exp(_expensive_functions_021*n);
      real _expensive_functions_023 = exp(B53*v);
      real _expensive_functions_024 = exp(-B63*v);
      real IObound_diff = -A53*IObound*Ku*_expensive_functions_023*_expensive_functions_024*_expensive_functions_026*_expensive_functions_028/A63 + Cbound*Kt/(_expensive_functions_018 + 1) + IO*Kmax*Ku*_expensive_functions_022/(_expensive_functions_020 + halfmax) - IObound*Kt;
      real _expensive_functions_029 = exp(0.14729709824716453*Vhalf - 0.14729709824716453*v);
      real _expensive_functions_030 = log(D);
      real _expensive_functions_031 = exp(_expensive_functions_030*n);
      real _expensive_functions_032 = log(D);
      real _expensive_functions_033 = exp(_expensive_functions_032*n);
      real Obound_diff = Cbound*Kt/(_expensive_functions_029 + 1) + Kmax*Ku*O*_expensive_functions_033/(_expensive_functions_031 + halfmax) - Kt*Obound - Ku*Obound;
      real hca = exp(F*qca*v/(R*T));
      real hna = exp(F*qna*v/(R*T));
      real h7_i = 1 + nao*(1 + (1.0/hna))/kna3;
      real h7_ss = 1 + nao*(1 + (1.0/hna))/kna3;
      real h8_i = nao/(h7_i*hna*kna3);
      real h8_ss = nao/(h7_ss*hna*kna3);
      real h9_i = (1.0/h7_i);
      real h9_ss = (1.0/h7_ss);
      real k3p_i = h9_i*wca;
      real k3p_ss = h9_ss*wca;
      real k3pp_i = h8_i*wnaca;
      real k3pp_ss = h8_ss*wnaca;
      real k8_i = h11_i*h8_i*wna;
      real k8_ss = h11_ss*h8_ss*wna;
      real k3_i = k3p_i + k3pp_i;
      real k3_ss = k3p_ss + k3pp_ss;
      real _expensive_functions_050 = exp((1.0/3.0)*F*delta*v/(R*T));
      real Knai = Knai0*_expensive_functions_050;
      real _expensive_functions_051 = exp((1.0/3.0)*F*v*(1 - delta)/(R*T));
      real Knao = Knao0*_expensive_functions_051;
      real a3 = (ko*ko)*k3p/((Kko*Kko)*(((1 + ko/Kko)*(1 + ko/Kko)) + ((1 + nao/Knao)*(1 + nao/Knao)*(1 + nao/Knao)) - 1));
      real b2 = (nao*nao*nao)*k2m/((Knao*Knao*Knao)*(((1 + ko/Kko)*(1 + ko/Kko)) + ((1 + nao/Knao)*(1 + nao/Knao)*(1 + nao/Knao)) - 1));
      real INab_U = INab_B*(-INab_v0 + v);
      real __melodee_temp_042 = -9.9999999999999995e-8 >= INab_U && INab_U >= 9.9999999999999995e-8;
      real _expensive_functions_055 = exp(0.0066137566137566143*v);
      real AiF = (1.0/(0.24348537187522867*_expensive_functions_055 + 1));
      real AiS = 1 - AiF;
      real Jtr = -1.0/100.0*cajsr + (1.0/100.0)*cansr;
      real _expensive_functions_068 = log(ko/ki);
      real EK = R*T*_expensive_functions_068/F;
      real _expensive_functions_069 = log(nao/nai);
      real ENa = R*T*_expensive_functions_069/F;
      real _expensive_functions_070 = log((PKNa*nao + ko)/(PKNa*nai + ki));
      real EKs = R*T*_expensive_functions_070/F;
      real Jdiff = -5.0*cai + 5.0*cass;
      real JdiffK = -1.0/2.0*ki + (1.0/2.0)*kss;
      real JdiffNa = -1.0/2.0*nai + (1.0/2.0)*nass;
      real CaMKt_diff = CaMKb*aCaMK*(CaMKb + CaMKt) - CaMKt*bCaMK;
      real km2n = jca;
      real _expensive_functions_072 = exp(vfrt);
      real A_2 = 0.75*ffrt*(_expensive_functions_072*nass - nao)/B_2;
      real _expensive_functions_073 = exp(vfrt);
      real A_3 = 0.75*ffrt*(_expensive_functions_073*kss - ko)/B_3;
      real anca = (1.0/(k2n/km2n + ((Kmn/cass + 1)*(Kmn/cass + 1)*(Kmn/cass + 1)*(Kmn/cass + 1))));
      real __melodee_temp_035;
      if (__melodee_temp_034)
      {
         __melodee_temp_035 = A_3*(1 - 0.5*U_3);
      }
      else
      {
         real _expensive_functions_074 = exp(U_3);
         __melodee_temp_035 = A_3*U_3/(_expensive_functions_074 - 1);
      }
      real PhiCaK = __melodee_temp_035;
      real __melodee_temp_039;
      if (__melodee_temp_038)
      {
         __melodee_temp_039 = A_2*(1 - 0.5*U_2);
      }
      else
      {
         real _expensive_functions_074 = exp(U_2);
         __melodee_temp_039 = A_2*U_2/(_expensive_functions_074 - 1);
      }
      real PhiCaNa = __melodee_temp_039;
      real nca_diff = anca*k2n - km2n*nca;
      real ICaK = PCaK*PhiCaK*d*(1 - fICaLp)*(f*(1 - nca) + fca*jca*nca) + PCaKp*PhiCaK*d*fICaLp*(fcap*jca*nca + fp*(1 - nca));
      real ICaNa = PCaNa*PhiCaNa*d*(1 - fICaLp)*(f*(1 - nca) + fca*jca*nca) + PCaNap*PhiCaNa*d*fICaLp*(fcap*jca*nca + fp*(1 - nca));
      real _expensive_functions_074 = exp(2*vfrt);
      real A = 4*PCab*ffrt*(_expensive_functions_074*cai - 0.34100000000000003*cao)/B;
      real __melodee_temp_041;
      if (__melodee_temp_040)
      {
         __melodee_temp_041 = A*(1 - 0.5*U);
      }
      else
      {
         _expensive_functions_075 = exp(U);
         __melodee_temp_041 = A*U/(_expensive_functions_075 - 1);
      }
      real ICab = __melodee_temp_041;
      real IK1 = GK1*_expensive_functions_075*rk1*xk1*(-EK + v);
      real IKb = GKb*xkb*(-EK + v);
      real _expensive_functions_076 = exp(B2*v);
      real _expensive_functions_079 = exp(B61*v);
      real _expensive_functions_082 = exp(B1*v);
      real _expensive_functions_085 = exp(B51*v);
      real C1_diff = -A1*C1*_expensive_functions_082*_expensive_functions_084 + A2*C2*_expensive_functions_076*_expensive_functions_078 - A51*C1*_expensive_functions_085*_expensive_functions_087 + A61*IC1*_expensive_functions_079*_expensive_functions_081;
      real _expensive_functions_088 = exp(B1*v);
      real _expensive_functions_091 = exp(B41*v);
      real _expensive_functions_094 = exp(B62*v);
      real _expensive_functions_097 = exp(B2*v);
      real _expensive_functions_100 = exp(B31*v);
      real _expensive_functions_103 = exp(B52*v);
      real C2_diff = A1*C1*_expensive_functions_088*_expensive_functions_090 - A2*C2*_expensive_functions_097*_expensive_functions_099 - A31*C2*_expensive_functions_100*_expensive_functions_102 + A41*O*_expensive_functions_091*_expensive_functions_093 - A52*C2*_expensive_functions_103*_expensive_functions_105 + A62*IC2*_expensive_functions_094*_expensive_functions_096;
      real _expensive_functions_106 = exp(0.14729709824716453*Vhalf - 0.14729709824716453*v);
      real Cbound_diff = -2*Cbound*Kt/(_expensive_functions_106 + 1) + IObound*Kt + Kt*Obound;
      real _expensive_functions_107 = exp(B21*v);
      real _expensive_functions_110 = exp(B51*v);
      real _expensive_functions_113 = exp(B11*v);
      real _expensive_functions_116 = exp(B61*v);
      real IC1_diff = -A11*IC1*_expensive_functions_113*_expensive_functions_115 + A21*IC2*_expensive_functions_107*_expensive_functions_109 + A51*C1*_expensive_functions_110*_expensive_functions_112 - A61*IC1*_expensive_functions_116*_expensive_functions_118;
      real _expensive_functions_119 = exp(B11*v);
      real _expensive_functions_122 = exp(B4*v);
      real _expensive_functions_125 = exp(B52*v);
      real _expensive_functions_128 = exp(B21*v);
      real _expensive_functions_131 = exp(B3*v);
      real _expensive_functions_134 = exp(B62*v);
      real IC2_diff = A11*IC1*_expensive_functions_119*_expensive_functions_121 - A21*IC2*_expensive_functions_128*_expensive_functions_130 - A3*IC2*_expensive_functions_131*_expensive_functions_133 + A4*IO*_expensive_functions_122*_expensive_functions_124 + A52*C2*_expensive_functions_125*_expensive_functions_127 - A62*IC2*_expensive_functions_134*_expensive_functions_136;
      real _expensive_functions_137 = exp(B3*v);
      real _expensive_functions_140 = exp(B53*v);
      real _expensive_functions_143 = exp(B4*v);
      real _expensive_functions_146 = exp(B63*v);
      real _expensive_functions_149 = log(D);
      real _expensive_functions_150 = exp(_expensive_functions_149*n);
      real _expensive_functions_151 = log(D);
      real _expensive_functions_152 = exp(_expensive_functions_151*n);
      real _expensive_functions_153 = exp(B53*v);
      real _expensive_functions_154 = exp(-B63*v);
      real IO_diff = A3*IC2*_expensive_functions_137*_expensive_functions_139 - A4*IO*_expensive_functions_143*_expensive_functions_145 + A53*O*_expensive_functions_140*_expensive_functions_142 + A53*IObound*Ku*_expensive_functions_153*_expensive_functions_154*_expensive_functions_156*_expensive_functions_158/A63 - A63*IO*_expensive_functions_146*_expensive_functions_148 - IO*Kmax*Ku*_expensive_functions_152/(_expensive_functions_150 + halfmax);
      real _expensive_functions_159 = exp(B31*v);
      real _expensive_functions_162 = exp(B63*v);
      real _expensive_functions_165 = exp(B41*v);
      real _expensive_functions_168 = exp(B53*v);
      real _expensive_functions_171 = log(D);
      real _expensive_functions_172 = exp(_expensive_functions_171*n);
      real _expensive_functions_173 = log(D);
      real _expensive_functions_174 = exp(_expensive_functions_173*n);
      real O_diff = A31*C2*_expensive_functions_159*_expensive_functions_161 - A41*O*_expensive_functions_165*_expensive_functions_167 - A53*O*_expensive_functions_168*_expensive_functions_170 + A63*IO*_expensive_functions_162*_expensive_functions_164 - Kmax*Ku*O*_expensive_functions_174/(_expensive_functions_172 + halfmax) + Ku*Obound;
      real IKr = 0.43033148291193518*GKr*O*_expensive_functions_175*(-EK + v);
      real _expensive_functions_176 = pow((1.0/cai), 1.3999999999999999);
      real KsCa = 1 + 0.59999999999999998/(6.4818210260626455e-7*_expensive_functions_176 + 1);
      real IKs = GKs*KsCa*xs1*xs2*(-EKs + v);
      real fINap = (1.0/(1 + KmCaMK/CaMKa));
      real h = Ahf*hf + Ahs*hs;
      real hp = Ahf*hf + Ahs*hsp;
      real INa = (m*m*m)*GNa*(-ENa + v)*(fINap*hp*jp + h*j*(1 - fINap));
      real allo_i = (1.0/((KmCaAct*KmCaAct)/(cai*cai) + 1));
      real allo_ss = (1.0/((KmCaAct*KmCaAct)/(cass*cass) + 1));
      real h4_i = 1 + nai*(1 + nai/kna2)/kna1;
      real h4_ss = 1 + nass*(1 + nass/kna2)/kna1;
      real h1_i = 1 + nai*(hna + 1)/kna3;
      real h1_ss = 1 + nass*(hna + 1)/kna3;
      real h5_i = (nai*nai)/(h4_i*kna1*kna2);
      real h5_ss = (nass*nass)/(h4_ss*kna1*kna2);
      real h6_i = (1.0/h4_i);
      real h6_ss = (1.0/h4_ss);
      real h2_i = hna*nai/(h1_i*kna3);
      real h2_ss = hna*nass/(h1_ss*kna3);
      real h3_i = (1.0/h1_i);
      real h3_ss = (1.0/h1_ss);
      real k6_i = cai*h6_i*kcaon;
      real k6_ss = cass*h6_ss*kcaon;
      real k4p_i = h3_i*wca/hca;
      real k4p_ss = h3_ss*wca/hca;
      real k4pp_i = h2_i*wnaca;
      real k4pp_ss = h2_ss*wnaca;
      real k7_i = h2_i*h5_i*wna;
      real k7_ss = h2_ss*h5_ss*wna;
      real k4_i = k4p_i + k4pp_i;
      real k4_ss = k4p_ss + k4pp_ss;
      real x1_i = k2_i*k4_i*(k6_i + k7_i) + k5_i*k7_i*(k2_i + k3_i);
      real x1_ss = k2_ss*k4_ss*(k6_ss + k7_ss) + k5_ss*k7_ss*(k2_ss + k3_ss);
      real x2_i = k1_i*k7_i*(k4_i + k5_i) + k4_i*k6_i*(k1_i + k8_i);
      real x2_ss = k1_ss*k7_ss*(k4_ss + k5_ss) + k4_ss*k6_ss*(k1_ss + k8_ss);
      real x3_i = k1_i*k3_i*(k6_i + k7_i) + k6_i*k8_i*(k2_i + k3_i);
      real x3_ss = k1_ss*k3_ss*(k6_ss + k7_ss) + k6_ss*k8_ss*(k2_ss + k3_ss);
      real x4_i = k2_i*k8_i*(k4_i + k5_i) + k3_i*k5_i*(k1_i + k8_i);
      real x4_ss = k2_ss*k8_ss*(k4_ss + k5_ss) + k3_ss*k5_ss*(k1_ss + k8_ss);
      real E1_i = x1_i/(x1_i + x2_i + x3_i + x4_i);
      real E1_ss = x1_ss/(x1_ss + x2_ss + x3_ss + x4_ss);
      real E2_i = x2_i/(x1_i + x2_i + x3_i + x4_i);
      real E2_ss = x2_ss/(x1_ss + x2_ss + x3_ss + x4_ss);
      real E3_i = x3_i/(x1_i + x2_i + x3_i + x4_i);
      real E3_ss = x3_ss/(x1_ss + x2_ss + x3_ss + x4_ss);
      real E4_i = x4_i/(x1_i + x2_i + x3_i + x4_i);
      real E4_ss = x4_ss/(x1_ss + x2_ss + x3_ss + x4_ss);
      real JncxCa_i = -E1_i*k1_i + E2_i*k2_i;
      real JncxCa_ss = -E1_ss*k1_ss + E2_ss*k2_ss;
      real JncxNa_i = -3*E1_i*k8_i - E2_i*k3pp_i + E3_i*k4pp_i + 3*E4_i*k7_i;
      real JncxNa_ss = -3*E1_ss*k8_ss - E2_ss*k3pp_ss + E3_ss*k4pp_ss + 3*E4_ss*k7_ss;
      real INaCa_i = 0.80000000000000004*Gncx*allo_i*(JncxCa_i*zca + JncxNa_i*zna);
      real INaCa_ss = 0.20000000000000001*Gncx*allo_ss*(JncxCa_ss*zca + JncxNa_ss*zna);
      real P = eP/(H/Khp + 1 + ki/Kxkur + nai/Knap);
      real a1 = (nai*nai*nai)*k1p/((Knai*Knai*Knai)*(((1 + ki/Kki)*(1 + ki/Kki)) + ((1 + nai/Knai)*(1 + nai/Knai)*(1 + nai/Knai)) - 1));
      real b3 = H*P*k3m/(1 + MgATP/Kmgatp);
      real b4 = (ki*ki)*k4m/((Kki*Kki)*(((1 + ki/Kki)*(1 + ki/Kki)) + ((1 + nai/Knai)*(1 + nai/Knai)*(1 + nai/Knai)) - 1));
      real x1 = a1*a2*a4 + a1*a2*b3 + a2*b3*b4 + b2*b3*b4;
      real x2 = a1*a2*a3 + a2*a3*b4 + a3*b1*b4 + b1*b2*b4;
      real x3 = a2*a3*a4 + a3*a4*b1 + a4*b1*b2 + b1*b2*b3;
      real x4 = a1*a3*a4 + a1*a4*b2 + a1*b2*b3 + b2*b3*b4;
      real E1 = x1/(x1 + x2 + x3 + x4);
      real E2 = x2/(x1 + x2 + x3 + x4);
      real E3 = x3/(x1 + x2 + x3 + x4);
      real E4 = x4/(x1 + x2 + x3 + x4);
      real JnakK = -2*E3*a1 + 2*E4*b1;
      real JnakNa = 3*E1*a3 - 3*E2*b3;
      real INaK = Pnak*(JnakK*zk + JnakNa*zna);
      real fINaLp = (1.0/(1 + KmCaMK/CaMKa));
      real INaL = GNaL*mL*(-ENa + v)*(fINaLp*hLp + hL*(1 - fINaLp));
      real _expensive_functions_177 = exp(vfrt);
      real INab_A = PNab*ffrt*(_expensive_functions_177*nai - nao)/INab_B;
      real __melodee_temp_043;
      if (__melodee_temp_042)
      {
         __melodee_temp_043 = INab_A*(1 - 0.5*INab_U);
      }
      else
      {
         _expensive_functions_178 = exp(INab_U);
         __melodee_temp_043 = INab_A*INab_U/(_expensive_functions_178 - 1);
      }
      real INab = __melodee_temp_043;
      real IpCa = GpCa*cai/(KmCap + cai);
      real fItop = (1.0/(1 + KmCaMK/CaMKa));
      real i = AiF*iF + AiS*iS;
      real ip = AiF*iFp + AiS*iSp;
      real Ito = Gto*(-EK + v)*(a*i*(1 - fItop) + ap*fItop*ip);
      real Jleak = 0.00026249999999999998*cansr;
      real fJupp = (1.0/(1 + KmCaMK/CaMKa));
      real Jupnp = 0.0043750000000000004*cai*upScale/(cai + 0.00092000000000000003);
      real Jupp = 0.01203125*cai*upScale/(cai + 0.00075000000000000002);
      real Jup = Jup_b*(-Jleak + Jupnp*(1 - fJupp) + Jupp*fJupp);
      real cansr_diff = -Jtr*vjsr/vnsr + Jup;
      real Bcajsr = (1.0/(csqnmax*kmcsqn/((cajsr + kmcsqn)*(cajsr + kmcsqn)) + 1));
      real Bcass = (1.0/(BSLmax*KmBSL/((KmBSL + cass)*(KmBSL + cass)) + BSRmax*KmBSR/((KmBSR + cass)*(KmBSR + cass)) + 1));
      real ki_diff = Acap*cm*(-IK1 - IKb - IKr - IKs + 2*INaK - Ito)/(F*vmyo) + JdiffK*vss/vmyo;
      real kss_diff = -Acap*ICaK*cm/(F*vss) - JdiffK;
      real nai_diff = Acap*cm*(-INa - 3*INaCa_i - 3*INaK - INaL - INab)/(F*vmyo) + JdiffNa*vss/vmyo;
      real nass_diff = Acap*cm*(-ICaNa - 3*INaCa_ss)/(F*vss) - JdiffNa;
      real Bcai = (1.0/(cmdnmax*kmcmdn/((cai + kmcmdn)*(cai + kmcmdn)) + kmtrpn*trpnmax/((cai + kmtrpn)*(cai + kmtrpn)) + 1));
      real cai_diff = Bcai*((1.0/2.0)*Acap*cm*(-ICab + 2*INaCa_i - IpCa)/(F*vmyo) + Jdiff*vss/vmyo - Jup*vnsr/vmyo);
      real fJrelp = (1.0/(1 + KmCaMK/CaMKa));
      real Jrel = Jrel_scaling_factor*(Jrelnp*(1 - fJrelp) + Jrelp*fJrelp);
      real __melodee_temp_046 = Jrel*JrelStiffConst > cajsr;
      real Jrel_001;
      if (__melodee_temp_046)
      {
         real ryr_Jrel = cajsr/JrelStiffConst;
         Jrel_001 = ryr_Jrel;
      }
      else
      {
         Jrel_001 = Jrel;
      }
      real cajsr_diff = Bcajsr*(-Jrel_001 + Jtr);
      real cass_diff = Bcass*((1.0/2.0)*Acap*cm*(-ICaL + 2*INaCa_ss)/(F*vss) - Jdiff + Jrel_001*vjsr/vss);
      //get Iion
      real Iion = ICaK + ICaL + ICaNa + ICab + IK1 + IKb + IKr + IKs + INa + INaCa_i + INaCa_ss + INaK + INaL + INab + IpCa + Ito;
      real Iion_001 = Iion;
      //Do the markov update (1 step rosenbrock with gauss siedel)
      //EDIT_STATE
      C1 += _dt*C1_diff;
      C2 += _dt*C2_diff;
      CaMKt += _dt*CaMKt_diff;
      Cbound += _dt*Cbound_diff;
      D += _dt*D_diff;
      IC1 += _dt*IC1_diff;
      IC2 += _dt*IC2_diff;
      IO += _dt*IO_diff;
      IObound += _dt*IObound_diff;
      O += _dt*O_diff;
      Obound += _dt*Obound_diff;
      cai += _dt*cai_diff;
      cajsr += _dt*cajsr_diff;
      cansr += _dt*cansr_diff;
      cass += _dt*cass_diff;
      d += _dt*d_diff;
      fcaf += _dt*fcaf_diff;
      fcafp += _dt*fcafp_diff;
      fcas += _dt*fcas_diff;
      ff += _dt*ff_diff;
      ffp += _dt*ffp_diff;
      fs += _dt*fs_diff;
      jca += _dt*jca_diff;
      ki += _dt*ki_diff;
      kss += _dt*kss_diff;
      nai += _dt*nai_diff;
      nass += _dt*nass_diff;
      nca += _dt*nca_diff;
      Jrelnp += _Jrelnp_RLA*(Jrelnp+_Jrelnp_RLB);
      Jrelp += _Jrelp_RLA*(Jrelp+_Jrelp_RLB);
      a += _a_RLA*(a+_a_RLB);
      ap += _ap_RLA*(ap+_ap_RLB);
      hL += _hL_RLA*(hL+_hL_RLB);
      hLp += _hLp_RLA*(hLp+_hLp_RLB);
      hf += _hf_RLA*(hf+_hf_RLB);
      hs += _hs_RLA*(hs+_hs_RLB);
      hsp += _hsp_RLA*(hsp+_hsp_RLB);
      iF += _iF_RLA*(iF+_iF_RLB);
      iFp += _iFp_RLA*(iFp+_iFp_RLB);
      iS += _iS_RLA*(iS+_iS_RLB);
      iSp += _iSp_RLA*(iSp+_iSp_RLB);
      j += _j_RLA*(j+_j_RLB);
      jp += _jp_RLA*(jp+_jp_RLB);
      m += _m_RLA*(m+_m_RLB);
      mL += _mL_RLA*(mL+_mL_RLB);
      xk1 += _xk1_RLA*(xk1+_xk1_RLB);
      xs1 += _xs1_RLA*(xs1+_xs1_RLB);
      xs2 += _xs2_RLA*(xs2+_xs2_RLB);
      store(state_[__jj].C1, C1);
      store(state_[__jj].C2, C2);
      store(state_[__jj].CaMKt, CaMKt);
      store(state_[__jj].Cbound, Cbound);
      store(state_[__jj].D, D);
      store(state_[__jj].IC1, IC1);
      store(state_[__jj].IC2, IC2);
      store(state_[__jj].IO, IO);
      store(state_[__jj].IObound, IObound);
      store(state_[__jj].Jrelnp, Jrelnp);
      store(state_[__jj].Jrelp, Jrelp);
      store(state_[__jj].O, O);
      store(state_[__jj].Obound, Obound);
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
      store(state_[__jj].hL, hL);
      store(state_[__jj].hLp, hLp);
      store(state_[__jj].hf, hf);
      store(state_[__jj].hs, hs);
      store(state_[__jj].hsp, hsp);
      store(state_[__jj].iF, iF);
      store(state_[__jj].iFp, iFp);
      store(state_[__jj].iS, iS);
      store(state_[__jj].iSp, iSp);
      store(state_[__jj].j, j);
      store(state_[__jj].jca, jca);
      store(state_[__jj].jp, jp);
      store(state_[__jj].ki, ki);
      store(state_[__jj].kss, kss);
      store(state_[__jj].m, m);
      store(state_[__jj].mL, mL);
      store(state_[__jj].nai, nai);
      store(state_[__jj].nass, nass);
      store(state_[__jj].nca, nca);
      store(state_[__jj].xk1, xk1);
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
   return "ohara_cipa_SCC";
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


   double V_init = -88.00190465;
   double V = V_init;
   double CaMKt_init = 0.0125840447;
   double CaMKt = CaMKt_init;
   double d_init = 2.3400000000000002e-9;
   double d = d_init;
   double ff_init = 0.99999999090000002;
   double ff = ff_init;
   double fs_init = 0.91024127769999996;
   double fs = fs_init;
   double fcaf_init = 0.99999999090000002;
   double fcaf = fcaf_init;
   double fcafp_init = 0.99999999090000002;
   double fcafp = fcafp_init;
   double fcas_init = 0.99980467770000003;
   double fcas = fcas_init;
   double ffp_init = 0.99999999090000002;
   double ffp = ffp_init;
   double jca_init = 0.99997383120000005;
   double jca = jca_init;
   double nca_init = 0.0027494140440000002;
   double nca = nca_init;
   double xk1_init = 0.99675975939999995;
   double xk1 = xk1_init;
   double D_init = 0;
   double D = D_init;
   double C1_init = 1.8014499999999999e-8;
   double C1 = C1_init;
   double C2_init = 8.2661899999999995e-5;
   double C2 = C2_init;
   double Cbound_init = 0;
   double Cbound = Cbound_init;
   double IC1_init = 0.999637;
   double IC1 = IC1_init;
   double IC2_init = 6.8320799999999998e-5;
   double IC2 = IC2_init;
   double IO_init = 5.6762299999999997e-5;
   double IO = IO_init;
   double IObound_init = 0;
   double IObound = IObound_init;
   double O_init = 0.00015551000000000001;
   double O = O_init;
   double Obound_init = 0;
   double Obound = Obound_init;
   double xs1_init = 0.2707758025;
   double xs1 = xs1_init;
   double xs2_init = 0.00019285034259999999;
   double xs2 = xs2_init;
   double hf_init = 0.69810719129999999;
   double hf = hf_init;
   double hs_init = 0.6980895801;
   double hs = hs_init;
   double m_init = 0.0073441211019999999;
   double m = m_init;
   double hsp_init = 0.45494855249999999;
   double hsp = hsp_init;
   double j_init = 0.69799084320000004;
   double j = j_init;
   double jp_init = 0.6979245865;
   double jp = jp_init;
   double hL_init = 0.50085488550000001;
   double hL = hL_init;
   double mL_init = 0.00018826172729999999;
   double mL = mL_init;
   double hLp_init = 0.26930653570000002;
   double hLp = hLp_init;
   double a_init = 0.0010010976869999999;
   double a = a_init;
   double ap_init = 0.00051008629340000002;
   double ap = ap_init;
   double iF_init = 0.99955417449999995;
   double iF = iF_init;
   double iS_init = 0.58650617360000001;
   double iS = iS_init;
   double iFp_init = 0.99955418230000004;
   double iFp = iFp_init;
   double iSp_init = 0.63933994819999995;
   double iSp = iSp_init;
   double cansr_init = 1.619574538;
   double cansr = cansr_init;
   double ki_init = 144.6555918;
   double ki = ki_init;
   double kss_init = 144.65556509999999;
   double kss = kss_init;
   double nai_init = 7.2680044979999998;
   double nai = nai_init;
   double nass_init = 7.2680899769999998;
   double nass = nass_init;
   double cajsr_init = 1.5712340140000001;
   double cajsr = cajsr_init;
   double cass_init = 8.4900000000000004e-5;
   double cass = cass_init;
   double cai_init = 8.6000000000000003e-5;
   double cai = cai_init;
   double Jrelnp_init = 2.4999999999999999e-7;
   double Jrelnp = Jrelnp_init;
   double Jrelp_init = 3.1199999999999999e-7;
   double Jrelp = Jrelp_init;
   for (int iCell=0; iCell<nCells_; iCell++)
   {
      READ_STATE(C1,iCell) = C1;
      READ_STATE(C2,iCell) = C2;
      READ_STATE(CaMKt,iCell) = CaMKt;
      READ_STATE(Cbound,iCell) = Cbound;
      READ_STATE(D,iCell) = D;
      READ_STATE(IC1,iCell) = IC1;
      READ_STATE(IC2,iCell) = IC2;
      READ_STATE(IO,iCell) = IO;
      READ_STATE(IObound,iCell) = IObound;
      READ_STATE(Jrelnp,iCell) = Jrelnp;
      READ_STATE(Jrelp,iCell) = Jrelp;
      READ_STATE(O,iCell) = O;
      READ_STATE(Obound,iCell) = Obound;
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
      READ_STATE(hL,iCell) = hL;
      READ_STATE(hLp,iCell) = hLp;
      READ_STATE(hf,iCell) = hf;
      READ_STATE(hs,iCell) = hs;
      READ_STATE(hsp,iCell) = hsp;
      READ_STATE(iF,iCell) = iF;
      READ_STATE(iFp,iCell) = iFp;
      READ_STATE(iS,iCell) = iS;
      READ_STATE(iSp,iCell) = iSp;
      READ_STATE(j,iCell) = j;
      READ_STATE(jca,iCell) = jca;
      READ_STATE(jp,iCell) = jp;
      READ_STATE(ki,iCell) = ki;
      READ_STATE(kss,iCell) = kss;
      READ_STATE(m,iCell) = m;
      READ_STATE(mL,iCell) = mL;
      READ_STATE(nai,iCell) = nai;
      READ_STATE(nass,iCell) = nass;
      READ_STATE(nca,iCell) = nca;
      READ_STATE(xk1,iCell) = xk1;
      READ_STATE(xs1,iCell) = xs1;
      READ_STATE(xs2,iCell) = xs2;
      __Vm[__indexArray[iCell]] = V_init;
   }
}

enum varHandles
{
   C1_handle,
   C2_handle,
   CaMKt_handle,
   Cbound_handle,
   D_handle,
   IC1_handle,
   IC2_handle,
   IO_handle,
   IObound_handle,
   Jrelnp_handle,
   Jrelp_handle,
   O_handle,
   Obound_handle,
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
   hL_handle,
   hLp_handle,
   hf_handle,
   hs_handle,
   hsp_handle,
   iF_handle,
   iFp_handle,
   iS_handle,
   iSp_handle,
   j_handle,
   jca_handle,
   jp_handle,
   ki_handle,
   kss_handle,
   m_handle,
   mL_handle,
   nai_handle,
   nass_handle,
   nca_handle,
   xk1_handle,
   xs1_handle,
   xs2_handle,
   NUMHANDLES
};

const string ThisReaction::getUnit(const std::string& varName) const
{
   if(0) {}
   else if (varName == "C1") { return "1"; }
   else if (varName == "C2") { return "1"; }
   else if (varName == "CaMKt") { return "mM"; }
   else if (varName == "Cbound") { return "1"; }
   else if (varName == "D") { return "1"; }
   else if (varName == "IC1") { return "1"; }
   else if (varName == "IC2") { return "1"; }
   else if (varName == "IO") { return "1"; }
   else if (varName == "IObound") { return "1"; }
   else if (varName == "Jrelnp") { return "1"; }
   else if (varName == "Jrelp") { return "1"; }
   else if (varName == "O") { return "1"; }
   else if (varName == "Obound") { return "1"; }
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
   else if (varName == "hL") { return "1"; }
   else if (varName == "hLp") { return "1"; }
   else if (varName == "hf") { return "1"; }
   else if (varName == "hs") { return "1"; }
   else if (varName == "hsp") { return "1"; }
   else if (varName == "iF") { return "1"; }
   else if (varName == "iFp") { return "1"; }
   else if (varName == "iS") { return "1"; }
   else if (varName == "iSp") { return "1"; }
   else if (varName == "j") { return "1"; }
   else if (varName == "jca") { return "1"; }
   else if (varName == "jp") { return "1"; }
   else if (varName == "ki") { return "mM"; }
   else if (varName == "kss") { return "mM"; }
   else if (varName == "m") { return "1"; }
   else if (varName == "mL") { return "1"; }
   else if (varName == "nai") { return "mM"; }
   else if (varName == "nass") { return "mM"; }
   else if (varName == "nca") { return "1"; }
   else if (varName == "xk1") { return "1"; }
   else if (varName == "xs1") { return "1"; }
   else if (varName == "xs2") { return "1"; }
   return "INVALID";
}

int ThisReaction::getVarHandle(const std::string& varName) const
{
   if (0) {}
   else if (varName == "C1") { return C1_handle; }
   else if (varName == "C2") { return C2_handle; }
   else if (varName == "CaMKt") { return CaMKt_handle; }
   else if (varName == "Cbound") { return Cbound_handle; }
   else if (varName == "D") { return D_handle; }
   else if (varName == "IC1") { return IC1_handle; }
   else if (varName == "IC2") { return IC2_handle; }
   else if (varName == "IO") { return IO_handle; }
   else if (varName == "IObound") { return IObound_handle; }
   else if (varName == "Jrelnp") { return Jrelnp_handle; }
   else if (varName == "Jrelp") { return Jrelp_handle; }
   else if (varName == "O") { return O_handle; }
   else if (varName == "Obound") { return Obound_handle; }
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
   else if (varName == "hL") { return hL_handle; }
   else if (varName == "hLp") { return hLp_handle; }
   else if (varName == "hf") { return hf_handle; }
   else if (varName == "hs") { return hs_handle; }
   else if (varName == "hsp") { return hsp_handle; }
   else if (varName == "iF") { return iF_handle; }
   else if (varName == "iFp") { return iFp_handle; }
   else if (varName == "iS") { return iS_handle; }
   else if (varName == "iSp") { return iSp_handle; }
   else if (varName == "j") { return j_handle; }
   else if (varName == "jca") { return jca_handle; }
   else if (varName == "jp") { return jp_handle; }
   else if (varName == "ki") { return ki_handle; }
   else if (varName == "kss") { return kss_handle; }
   else if (varName == "m") { return m_handle; }
   else if (varName == "mL") { return mL_handle; }
   else if (varName == "nai") { return nai_handle; }
   else if (varName == "nass") { return nass_handle; }
   else if (varName == "nca") { return nca_handle; }
   else if (varName == "xk1") { return xk1_handle; }
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
   else if (varHandle == C1_handle) { READ_STATE(C1,iCell) = value; }
   else if (varHandle == C2_handle) { READ_STATE(C2,iCell) = value; }
   else if (varHandle == CaMKt_handle) { READ_STATE(CaMKt,iCell) = value; }
   else if (varHandle == Cbound_handle) { READ_STATE(Cbound,iCell) = value; }
   else if (varHandle == D_handle) { READ_STATE(D,iCell) = value; }
   else if (varHandle == IC1_handle) { READ_STATE(IC1,iCell) = value; }
   else if (varHandle == IC2_handle) { READ_STATE(IC2,iCell) = value; }
   else if (varHandle == IO_handle) { READ_STATE(IO,iCell) = value; }
   else if (varHandle == IObound_handle) { READ_STATE(IObound,iCell) = value; }
   else if (varHandle == Jrelnp_handle) { READ_STATE(Jrelnp,iCell) = value; }
   else if (varHandle == Jrelp_handle) { READ_STATE(Jrelp,iCell) = value; }
   else if (varHandle == O_handle) { READ_STATE(O,iCell) = value; }
   else if (varHandle == Obound_handle) { READ_STATE(Obound,iCell) = value; }
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
   else if (varHandle == hL_handle) { READ_STATE(hL,iCell) = value; }
   else if (varHandle == hLp_handle) { READ_STATE(hLp,iCell) = value; }
   else if (varHandle == hf_handle) { READ_STATE(hf,iCell) = value; }
   else if (varHandle == hs_handle) { READ_STATE(hs,iCell) = value; }
   else if (varHandle == hsp_handle) { READ_STATE(hsp,iCell) = value; }
   else if (varHandle == iF_handle) { READ_STATE(iF,iCell) = value; }
   else if (varHandle == iFp_handle) { READ_STATE(iFp,iCell) = value; }
   else if (varHandle == iS_handle) { READ_STATE(iS,iCell) = value; }
   else if (varHandle == iSp_handle) { READ_STATE(iSp,iCell) = value; }
   else if (varHandle == j_handle) { READ_STATE(j,iCell) = value; }
   else if (varHandle == jca_handle) { READ_STATE(jca,iCell) = value; }
   else if (varHandle == jp_handle) { READ_STATE(jp,iCell) = value; }
   else if (varHandle == ki_handle) { READ_STATE(ki,iCell) = value; }
   else if (varHandle == kss_handle) { READ_STATE(kss,iCell) = value; }
   else if (varHandle == m_handle) { READ_STATE(m,iCell) = value; }
   else if (varHandle == mL_handle) { READ_STATE(mL,iCell) = value; }
   else if (varHandle == nai_handle) { READ_STATE(nai,iCell) = value; }
   else if (varHandle == nass_handle) { READ_STATE(nass,iCell) = value; }
   else if (varHandle == nca_handle) { READ_STATE(nca,iCell) = value; }
   else if (varHandle == xk1_handle) { READ_STATE(xk1,iCell) = value; }
   else if (varHandle == xs1_handle) { READ_STATE(xs1,iCell) = value; }
   else if (varHandle == xs2_handle) { READ_STATE(xs2,iCell) = value; }
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   if (0) {}
   else if (varHandle == C1_handle) { return READ_STATE(C1,iCell); }
   else if (varHandle == C2_handle) { return READ_STATE(C2,iCell); }
   else if (varHandle == CaMKt_handle) { return READ_STATE(CaMKt,iCell); }
   else if (varHandle == Cbound_handle) { return READ_STATE(Cbound,iCell); }
   else if (varHandle == D_handle) { return READ_STATE(D,iCell); }
   else if (varHandle == IC1_handle) { return READ_STATE(IC1,iCell); }
   else if (varHandle == IC2_handle) { return READ_STATE(IC2,iCell); }
   else if (varHandle == IO_handle) { return READ_STATE(IO,iCell); }
   else if (varHandle == IObound_handle) { return READ_STATE(IObound,iCell); }
   else if (varHandle == Jrelnp_handle) { return READ_STATE(Jrelnp,iCell); }
   else if (varHandle == Jrelp_handle) { return READ_STATE(Jrelp,iCell); }
   else if (varHandle == O_handle) { return READ_STATE(O,iCell); }
   else if (varHandle == Obound_handle) { return READ_STATE(Obound,iCell); }
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
   else if (varHandle == hL_handle) { return READ_STATE(hL,iCell); }
   else if (varHandle == hLp_handle) { return READ_STATE(hLp,iCell); }
   else if (varHandle == hf_handle) { return READ_STATE(hf,iCell); }
   else if (varHandle == hs_handle) { return READ_STATE(hs,iCell); }
   else if (varHandle == hsp_handle) { return READ_STATE(hsp,iCell); }
   else if (varHandle == iF_handle) { return READ_STATE(iF,iCell); }
   else if (varHandle == iFp_handle) { return READ_STATE(iFp,iCell); }
   else if (varHandle == iS_handle) { return READ_STATE(iS,iCell); }
   else if (varHandle == iSp_handle) { return READ_STATE(iSp,iCell); }
   else if (varHandle == j_handle) { return READ_STATE(j,iCell); }
   else if (varHandle == jca_handle) { return READ_STATE(jca,iCell); }
   else if (varHandle == jp_handle) { return READ_STATE(jp,iCell); }
   else if (varHandle == ki_handle) { return READ_STATE(ki,iCell); }
   else if (varHandle == kss_handle) { return READ_STATE(kss,iCell); }
   else if (varHandle == m_handle) { return READ_STATE(m,iCell); }
   else if (varHandle == mL_handle) { return READ_STATE(mL,iCell); }
   else if (varHandle == nai_handle) { return READ_STATE(nai,iCell); }
   else if (varHandle == nass_handle) { return READ_STATE(nass,iCell); }
   else if (varHandle == nca_handle) { return READ_STATE(nca,iCell); }
   else if (varHandle == xk1_handle) { return READ_STATE(xk1,iCell); }
   else if (varHandle == xs1_handle) { return READ_STATE(xs1,iCell); }
   else if (varHandle == xs2_handle) { return READ_STATE(xs2,iCell); }
   return NAN;
}

double ThisReaction::getValue(int iCell, int varHandle, double V) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   const double C1=READ_STATE(C1,iCell);
   const double C2=READ_STATE(C2,iCell);
   const double CaMKt=READ_STATE(CaMKt,iCell);
   const double Cbound=READ_STATE(Cbound,iCell);
   const double D=READ_STATE(D,iCell);
   const double IC1=READ_STATE(IC1,iCell);
   const double IC2=READ_STATE(IC2,iCell);
   const double IO=READ_STATE(IO,iCell);
   const double IObound=READ_STATE(IObound,iCell);
   const double Jrelnp=READ_STATE(Jrelnp,iCell);
   const double Jrelp=READ_STATE(Jrelp,iCell);
   const double O=READ_STATE(O,iCell);
   const double Obound=READ_STATE(Obound,iCell);
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
   const double hL=READ_STATE(hL,iCell);
   const double hLp=READ_STATE(hLp,iCell);
   const double hf=READ_STATE(hf,iCell);
   const double hs=READ_STATE(hs,iCell);
   const double hsp=READ_STATE(hsp,iCell);
   const double iF=READ_STATE(iF,iCell);
   const double iFp=READ_STATE(iFp,iCell);
   const double iS=READ_STATE(iS,iCell);
   const double iSp=READ_STATE(iSp,iCell);
   const double j=READ_STATE(j,iCell);
   const double jca=READ_STATE(jca,iCell);
   const double jp=READ_STATE(jp,iCell);
   const double ki=READ_STATE(ki,iCell);
   const double kss=READ_STATE(kss,iCell);
   const double m=READ_STATE(m,iCell);
   const double mL=READ_STATE(mL,iCell);
   const double nai=READ_STATE(nai,iCell);
   const double nass=READ_STATE(nass,iCell);
   const double nca=READ_STATE(nca,iCell);
   const double xk1=READ_STATE(xk1,iCell);
   const double xs1=READ_STATE(xs1,iCell);
   const double xs2=READ_STATE(xs2,iCell);
   if (0) {}
   else if (varHandle == C1_handle)
   {
      return C1;
   }
   else if (varHandle == C2_handle)
   {
      return C2;
   }
   else if (varHandle == CaMKt_handle)
   {
      return CaMKt;
   }
   else if (varHandle == Cbound_handle)
   {
      return Cbound;
   }
   else if (varHandle == D_handle)
   {
      return D;
   }
   else if (varHandle == IC1_handle)
   {
      return IC1;
   }
   else if (varHandle == IC2_handle)
   {
      return IC2;
   }
   else if (varHandle == IO_handle)
   {
      return IO;
   }
   else if (varHandle == IObound_handle)
   {
      return IObound;
   }
   else if (varHandle == Jrelnp_handle)
   {
      return Jrelnp;
   }
   else if (varHandle == Jrelp_handle)
   {
      return Jrelp;
   }
   else if (varHandle == O_handle)
   {
      return O;
   }
   else if (varHandle == Obound_handle)
   {
      return Obound;
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
   else if (varHandle == hL_handle)
   {
      return hL;
   }
   else if (varHandle == hLp_handle)
   {
      return hLp;
   }
   else if (varHandle == hf_handle)
   {
      return hf;
   }
   else if (varHandle == hs_handle)
   {
      return hs;
   }
   else if (varHandle == hsp_handle)
   {
      return hsp;
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
   else if (varHandle == j_handle)
   {
      return j;
   }
   else if (varHandle == jca_handle)
   {
      return jca;
   }
   else if (varHandle == jp_handle)
   {
      return jp;
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
   fieldNames.push_back("C1");
   fieldUnits.push_back(getUnit("C1"));
   fieldNames.push_back("C2");
   fieldUnits.push_back(getUnit("C2"));
   fieldNames.push_back("CaMKt");
   fieldUnits.push_back(getUnit("CaMKt"));
   fieldNames.push_back("Cbound");
   fieldUnits.push_back(getUnit("Cbound"));
   fieldNames.push_back("D");
   fieldUnits.push_back(getUnit("D"));
   fieldNames.push_back("IC1");
   fieldUnits.push_back(getUnit("IC1"));
   fieldNames.push_back("IC2");
   fieldUnits.push_back(getUnit("IC2"));
   fieldNames.push_back("IO");
   fieldUnits.push_back(getUnit("IO"));
   fieldNames.push_back("IObound");
   fieldUnits.push_back(getUnit("IObound"));
   fieldNames.push_back("Jrelnp");
   fieldUnits.push_back(getUnit("Jrelnp"));
   fieldNames.push_back("Jrelp");
   fieldUnits.push_back(getUnit("Jrelp"));
   fieldNames.push_back("O");
   fieldUnits.push_back(getUnit("O"));
   fieldNames.push_back("Obound");
   fieldUnits.push_back(getUnit("Obound"));
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
   fieldNames.push_back("hL");
   fieldUnits.push_back(getUnit("hL"));
   fieldNames.push_back("hLp");
   fieldUnits.push_back(getUnit("hLp"));
   fieldNames.push_back("hf");
   fieldUnits.push_back(getUnit("hf"));
   fieldNames.push_back("hs");
   fieldUnits.push_back(getUnit("hs"));
   fieldNames.push_back("hsp");
   fieldUnits.push_back(getUnit("hsp"));
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
   fieldNames.push_back("jp");
   fieldUnits.push_back(getUnit("jp"));
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
   fieldNames.push_back("xs1");
   fieldUnits.push_back(getUnit("xs1"));
   fieldNames.push_back("xs2");
   fieldUnits.push_back(getUnit("xs2"));
}

}
