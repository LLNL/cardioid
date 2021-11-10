/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "TT06_SCC.hh"
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
   "_fCass_RLA",
   "inward_rectifier_potassium_current_i_Kitot",
    NULL
};


   REACTION_FACTORY(TT06_SCC)(OBJECT* obj, const double _dt, const int numPoints, const ThreadTeam&)
   {
      TT06_SCC::ThisReaction* reaction = new TT06_SCC::ThisReaction(numPoints, _dt);

      //override the defaults
      //EDIT_PARAMETERS
      double celltype;
      double g_CaL;
      double g_K1;
      double g_Kr;
      double g_Ks;
      double g_Na;
      double g_bca;
      double g_bna;
      double g_pCa;
      double g_pK;
      double g_to;
      setDefault(celltype, 0);
      setDefault(g_CaL, 3.9799999999999998e-5);
      setDefault(g_bca, 0.00059199999999999997);
      setDefault(g_pCa, 0.12379999999999999);
      setDefault(g_Na, 14.837999999999999);
      setDefault(g_K1, 5.4050000000000002);
      setDefault(g_pK, 0.0146);
      setDefault(g_Kr, 0.153);
      double __melodee_temp_002 = celltype == 1;
      if (__melodee_temp_002)
      {
         g_Ks = 0.098000000000000004;
      }
      else
      {
         g_Ks = 0.39200000000000002;
      }
      setDefault(g_Ks, g_Ks);
      setDefault(g_bna, 0.00029);
      double __melodee_temp_004 = celltype == 0;
      if (__melodee_temp_004)
      {
         g_to = 0.072999999999999995;
      }
      else
      {
         g_to = 0.29399999999999998;
      }
      setDefault(g_to, g_to);
      reaction->celltype = celltype;
      reaction->g_CaL = g_CaL;
      reaction->g_K1 = g_K1;
      reaction->g_Kr = g_Kr;
      reaction->g_Ks = g_Ks;
      reaction->g_Na = g_Na;
      reaction->g_bca = g_bca;
      reaction->g_bna = g_bna;
      reaction->g_pCa = g_pCa;
      reaction->g_pK = g_pK;
      reaction->g_to = g_to;
      bool reusingInterpolants = false;
      string fitName;
      objectGet(obj, "fit", fitName, "");
      int funcCount = sizeof(reaction->_interpolant)/sizeof(reaction->_interpolant[0])-1; //BGQ_HACKFIX, compiler bug with zero length arrays
      if (fitName != "")
      {
         OBJECT* fitObj = objectFind(fitName, "FIT");
         double _fit_dt; objectGet(fitObj, "dt", _fit_dt, "nan");
         double _fit_g_K1; objectGet(fitObj, "g_K1", _fit_g_K1, "nan");
         if (1
            && _fit_dt == _dt
            && _fit_g_K1 == reaction->g_K1
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
            outfile << "   dt = " << _dt << ";\n";
            outfile << "   g_K1 = " << reaction->g_K1 << ";\n";
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

namespace TT06_SCC 
{

void ThisReaction::createInterpolants(const double _dt) {

   {
      int _numPoints = (1e-1 - 1e-7)/1e-6;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double Ca_ss = 1e-7 + (1e-1 - 1e-7)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = Ca_ss;
         double tau_fCass = 2 + 80/(400.0*(Ca_ss*Ca_ss) + 1);
         double _expensive_functions_049 = exp(-_dt/tau_fCass);
         double _fCass_RLA = _expensive_functions_049 - 1;
         _outputs[_ii] = _fCass_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[0].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _fCass_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double VEK = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = VEK;
         double K_o = 5.4000000000000004;
         double _expensive_functions_035 = exp(0.059999999999999998*VEK);
         double alpha_K1 = 0.10000000000000001/(6.1442123533282098e-6*_expensive_functions_035 + 1);
         double _expensive_functions_036 = exp(-0.5*VEK);
         double _expensive_functions_037 = exp(0.10000000000000001*VEK);
         double _expensive_functions_038 = exp(0.00020000000000000001*VEK);
         double beta_K1 = (0.36787944117144233*_expensive_functions_037 + 3.0606040200802673*_expensive_functions_038)/(_expensive_functions_036 + 1);
         double xK1_inf = alpha_K1/(alpha_K1 + beta_K1);
         double _expensive_functions_039 = sqrt(K_o);
         double i_K1 = 0.43033148291193518*VEK*_expensive_functions_039*g_K1*xK1_inf;
         double inward_rectifier_potassium_current_i_Kitot = i_K1;
         _outputs[_ii] = inward_rectifier_potassium_current_i_Kitot;
      }
      double relError = 1e-3;
      double actualTolerance = _interpolant[1].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for inward_rectifier_potassium_current_i_Kitot: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
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

   "   Ca_SR_off,\n"
   "   Ca_i_off,\n"
   "   Ca_ss_off,\n"
   "   K_i_off,\n"
   "   Na_i_off,\n"
   "   R_prime_off,\n"
   "   Xr1_off,\n"
   "   Xr2_off,\n"
   "   Xs_off,\n"
   "   d_off,\n"
   "   f_off,\n"
   "   f2_off,\n"
   "   fCass_off,\n"
   "   h_off,\n"
   "   j_off,\n"
   "   m_off,\n"
   "   r_off,\n"
   "   s_off,\n"
   "   NUMSTATES\n"
   "};\n"
   "extern \"C\"\n"
   "__global__ void TT06_SCC_kernel(const double* _Vm, const double* _iStim, double* _dVm, double* _state) {\n"
   "const double _dt = " << __cachedDt << ";\n"
   "const int _nCells = " << nCells_ << ";\n"

   "const double celltype = " << celltype << ";\n"
   "const double g_CaL = " << g_CaL << ";\n"
   "const double g_K1 = " << g_K1 << ";\n"
   "const double g_Kr = " << g_Kr << ";\n"
   "const double g_Ks = " << g_Ks << ";\n"
   "const double g_Na = " << g_Na << ";\n"
   "const double g_bca = " << g_bca << ";\n"
   "const double g_bna = " << g_bna << ";\n"
   "const double g_pCa = " << g_pCa << ";\n"
   "const double g_pK = " << g_pK << ";\n"
   "const double g_to = " << g_to << ";\n"
   "const int _ii = threadIdx.x + blockIdx.x*blockDim.x;\n"
   "if (_ii >= _nCells) { return; }\n"
   "const double V = _Vm[_ii];\n"
   "double _ratPoly;\n"

   "const double Ca_SR = _state[_ii+Ca_SR_off*_nCells];\n"
   "const double Ca_i = _state[_ii+Ca_i_off*_nCells];\n"
   "const double Ca_ss = _state[_ii+Ca_ss_off*_nCells];\n"
   "const double K_i = _state[_ii+K_i_off*_nCells];\n"
   "const double Na_i = _state[_ii+Na_i_off*_nCells];\n"
   "const double R_prime = _state[_ii+R_prime_off*_nCells];\n"
   "const double Xr1 = _state[_ii+Xr1_off*_nCells];\n"
   "const double Xr2 = _state[_ii+Xr2_off*_nCells];\n"
   "const double Xs = _state[_ii+Xs_off*_nCells];\n"
   "const double d = _state[_ii+d_off*_nCells];\n"
   "const double f = _state[_ii+f_off*_nCells];\n"
   "const double f2 = _state[_ii+f2_off*_nCells];\n"
   "const double fCass = _state[_ii+fCass_off*_nCells];\n"
   "const double h = _state[_ii+h_off*_nCells];\n"
   "const double j = _state[_ii+j_off*_nCells];\n"
   "const double m = _state[_ii+m_off*_nCells];\n"
   "const double r = _state[_ii+r_off*_nCells];\n"
   "const double s = _state[_ii+s_off*_nCells];\n"
   "//get the gate updates (diagonalized exponential integrator)\n"
   "double _expensive_functions = exp(-1.0/13.0*V - 35.0/13.0);\n"
   "double alpha_d = 0.25 + 1.3999999999999999/(_expensive_functions + 1);\n"
   "double _expensive_functions_001 = exp((1.0/5.0)*V + 1);\n"
   "double beta_d = 1.3999999999999999/(_expensive_functions_001 + 1);\n"
   "double _expensive_functions_002 = exp(-0.13333333333333333*V);\n"
   "double d_inf = (1.0/(0.34415378686541237*_expensive_functions_002 + 1));\n"
   "double _expensive_functions_003 = exp(5.0/2.0 - 1.0/20.0*V);\n"
   "double gamma_d = (1.0/(_expensive_functions_003 + 1));\n"
   "double tau_d = alpha_d*beta_d + gamma_d;\n"
   "double _expensive_functions_004 = exp((1.0/7.0)*V + 5);\n"
   "double f2_inf = 0.33000000000000002 + 0.67000000000000004/(_expensive_functions_004 + 1);\n"
   "double _expensive_functions_005 = exp(5.0/2.0 - 1.0/10.0*V);\n"
   "double _expensive_functions_006 = exp((1.0/10.0)*V + 3);\n"
   "double _expensive_functions_007 = exp(-1.0/240.0*((V + 27)*(V + 27)));\n"
   "double tau_f2 = 562*_expensive_functions_007 + 80/(_expensive_functions_006 + 1) + 31/(_expensive_functions_005 + 1);\n"
   "double fCass_inf = 0.40000000000000002 + 0.59999999999999998/(400.0*(Ca_ss*Ca_ss) + 1);\n"
   "double _expensive_functions_008 = exp((1.0/7.0)*V + 20.0/7.0);\n"
   "double f_inf = (1.0/(_expensive_functions_008 + 1));\n"
   "double _expensive_functions_009 = exp((1.0/10.0)*V + 3);\n"
   "double _expensive_functions_010 = exp(13.0/10.0 - 1.0/10.0*V);\n"
   "double _expensive_functions_011 = exp(-1.0/225.0*((V + 27)*(V + 27)));\n"
   "double tau_f = 1102.5*_expensive_functions_011 + 20 + 200/(_expensive_functions_010 + 1) + 180/(_expensive_functions_009 + 1);\n"
   "double __melodee_temp_000 = V < -40;\n"
   "double alpha_h;\n"
   "double beta_h;\n"
   "if (__melodee_temp_000)\n"
   "{\n"
   "   double _expensive_functions_013 = exp(-0.14705882352941177*V);\n"
   "   alpha_h = 4.4312679295805147e-7*_expensive_functions_013;\n"
   "   double _expensive_functions_014 = exp(0.34849999999999998*V);\n"
   "   double _expensive_functions_015 = exp(0.079000000000000001*V);\n"
   "   beta_h = 310000*_expensive_functions_014 + 2.7000000000000002*_expensive_functions_015;\n"
   "}\n"
   "else\n"
   "{\n"
   "   alpha_h = 0;\n"
   "   double _expensive_functions_013 = exp(-0.0900900900900901*V);\n"
   "   beta_h = 0.77000000000000002/(0.049758141083938695*_expensive_functions_013 + 0.13);\n"
   "}\n"
   "double _expensive_functions_013 = exp(0.13458950201884254*V);\n"
   "double h_inf = 4.3210917837689708e-9/((_expensive_functions_013 + 6.5735011856460265e-5)*(_expensive_functions_013 + 6.5735011856460265e-5));\n"
   "double tau_h = (1.0/(alpha_h + beta_h));\n"
   "double __melodee_temp_001 = V < -40;\n"
   "double alpha_j;\n"
   "double beta_j;\n"
   "if (__melodee_temp_001)\n"
   "{\n"
   "   double _expensive_functions_014 = exp(0.311*V);\n"
   "   double _expensive_functions_015 = exp(0.24440000000000001*V);\n"
   "   double _expensive_functions_016 = exp(-0.043909999999999998*V);\n"
   "   alpha_j = (-25428*_expensive_functions_015 - 6.9480000000000002e-6*_expensive_functions_016)*(V + 37.780000000000001)/(50262745825.953987*_expensive_functions_014 + 1);\n"
   "   double _expensive_functions_017 = exp(-0.13780000000000001*V);\n"
   "   double _expensive_functions_018 = exp(-0.01052*V);\n"
   "   beta_j = 0.024240000000000001*_expensive_functions_018/(0.003960868339904256*_expensive_functions_017 + 1);\n"
   "}\n"
   "else\n"
   "{\n"
   "   alpha_j = 0;\n"
   "   double _expensive_functions_014 = exp(-0.10000000000000001*V);\n"
   "   double _expensive_functions_015 = exp(0.057000000000000002*V);\n"
   "   beta_j = 0.59999999999999998*_expensive_functions_015/(0.040762203978366204*_expensive_functions_014 + 1);\n"
   "}\n"
   "double _expensive_functions_014 = exp(0.13458950201884254*V);\n"
   "double j_inf = 4.3210917837689708e-9/((_expensive_functions_014 + 6.5735011856460265e-5)*(_expensive_functions_014 + 6.5735011856460265e-5));\n"
   "double tau_j = (1.0/(alpha_j + beta_j));\n"
   "double _expensive_functions_015 = exp(-1.0/5.0*V - 12);\n"
   "double alpha_m = (1.0/(_expensive_functions_015 + 1));\n"
   "double _expensive_functions_016 = exp((1.0/5.0)*V + 7);\n"
   "double _expensive_functions_017 = exp((1.0/200.0)*V - 1.0/4.0);\n"
   "double beta_m = 0.10000000000000001/(_expensive_functions_017 + 1) + 0.10000000000000001/(_expensive_functions_016 + 1);\n"
   "double _expensive_functions_018 = exp(-0.11074197120708749*V);\n"
   "double m_inf = (1.0/(0.0018422115811651339*_expensive_functions_018 + 1)/(0.0018422115811651339*_expensive_functions_018 + 1));\n"
   "double tau_m = alpha_m*beta_m;\n"
   "double _expensive_functions_020 = exp(-1.0/10.0*V - 9.0/2.0);\n"
   "double alpha_xr1 = 450/(_expensive_functions_020 + 1);\n"
   "double _expensive_functions_021 = exp(0.086956521739130432*V);\n"
   "double beta_xr1 = 6/(13.581324522578193*_expensive_functions_021 + 1);\n"
   "double _expensive_functions_022 = exp(-1.0/7.0*V - 26.0/7.0);\n"
   "double xr1_inf = (1.0/(_expensive_functions_022 + 1));\n"
   "double tau_xr1 = alpha_xr1*beta_xr1;\n"
   "double _expensive_functions_023 = exp(-1.0/20.0*V - 3);\n"
   "double alpha_xr2 = 3/(_expensive_functions_023 + 1);\n"
   "double _expensive_functions_024 = exp((1.0/20.0)*V - 3);\n"
   "double beta_xr2 = 1.1200000000000001/(_expensive_functions_024 + 1);\n"
   "double _expensive_functions_025 = exp((1.0/24.0)*V + 11.0/3.0);\n"
   "double xr2_inf = (1.0/(_expensive_functions_025 + 1));\n"
   "double tau_xr2 = alpha_xr2*beta_xr2;\n"
   "double _expensive_functions_027 = exp(5.0/6.0 - 1.0/6.0*V);\n"
   "double _expensive_functions_028 = pow(_expensive_functions_027 + 1, -1.0/2.0);\n"
   "double alpha_xs = 1400*_expensive_functions_028;\n"
   "double _expensive_functions_029 = exp((1.0/15.0)*V - 7.0/3.0);\n"
   "double beta_xs = (1.0/(_expensive_functions_029 + 1));\n"
   "double _expensive_functions_030 = exp(-1.0/14.0*V - 5.0/14.0);\n"
   "double xs_inf = (1.0/(_expensive_functions_030 + 1));\n"
   "double tau_xs = alpha_xs*beta_xs + 80;\n"
   "double _expensive_functions_033 = exp(10.0/3.0 - 1.0/6.0*V);\n"
   "double r_inf = (1.0/(_expensive_functions_033 + 1));\n"
   "double _expensive_functions_034 = exp(-1.0/1800.0*((V + 40)*(V + 40)));\n"
   "double tau_r = 9.5*_expensive_functions_034 + 0.80000000000000004;\n"
   "double __melodee_temp_003 = celltype == 0;\n"
   "double s_inf;\n"
   "double tau_s;\n"
   "if (__melodee_temp_003)\n"
   "{\n"
   "   double _expensive_functions_035 = exp((1.0/5.0)*V + 28.0/5.0);\n"
   "   s_inf = (1.0/(_expensive_functions_035 + 1));\n"
   "   double _expensive_functions_036 = exp(-1.0/1000.0*((V + 67)*(V + 67)));\n"
   "   tau_s = 1000*_expensive_functions_036 + 8;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double _expensive_functions_035 = exp((1.0/5.0)*V + 4);\n"
   "   s_inf = (1.0/(_expensive_functions_035 + 1));\n"
   "   double _expensive_functions_036 = exp((1.0/5.0)*V - 4);\n"
   "   double _expensive_functions_037 = exp(-1.0/320.0*((V + 45)*(V + 45)));\n"
   "   tau_s = 85*_expensive_functions_037 + 3 + 5/(_expensive_functions_036 + 1);\n"
   "}\n"
   "double _expensive_functions_043 = exp(-_dt/tau_xr1);\n"
   "double _Xr1_RLA = _expensive_functions_043 - 1;\n"
   "double _Xr1_RLB = -xr1_inf;\n"
   "double _expensive_functions_044 = exp(-_dt/tau_xr2);\n"
   "double _Xr2_RLA = _expensive_functions_044 - 1;\n"
   "double _Xr2_RLB = -xr2_inf;\n"
   "double _expensive_functions_045 = exp(-_dt/tau_xs);\n"
   "double _Xs_RLA = _expensive_functions_045 - 1;\n"
   "double _Xs_RLB = -xs_inf;\n"
   "double _expensive_functions_046 = exp(-_dt/tau_d);\n"
   "double _d_RLA = _expensive_functions_046 - 1;\n"
   "double _d_RLB = -d_inf;\n"
   "double _expensive_functions_047 = exp(-_dt/tau_f);\n"
   "double _f_RLA = _expensive_functions_047 - 1;\n"
   "double _f_RLB = -f_inf;\n"
   "double _expensive_functions_048 = exp(-_dt/tau_f2);\n"
   "double _f2_RLA = _expensive_functions_048 - 1;\n"
   "double _f2_RLB = -f2_inf;\n"
   ""; generateInterpString(ss,_interpolant[0], "Ca_ss"); ss << "\n"
   "double _fCass_RLA = _ratPoly;\n"
   "double _fCass_RLB = -fCass_inf;\n"
   "double _expensive_functions_050 = exp(-_dt/tau_h);\n"
   "double _h_RLA = _expensive_functions_050 - 1;\n"
   "double _h_RLB = -h_inf;\n"
   "double _expensive_functions_051 = exp(-_dt/tau_j);\n"
   "double _j_RLA = _expensive_functions_051 - 1;\n"
   "double _j_RLB = -j_inf;\n"
   "double _expensive_functions_052 = exp(-_dt/tau_m);\n"
   "double _m_RLA = _expensive_functions_052 - 1;\n"
   "double _m_RLB = -m_inf;\n"
   "double _expensive_functions_053 = exp(-_dt/tau_r);\n"
   "double _r_RLA = _expensive_functions_053 - 1;\n"
   "double _r_RLB = -r_inf;\n"
   "double _expensive_functions_054 = exp(-_dt/tau_s);\n"
   "double _s_RLA = _expensive_functions_054 - 1;\n"
   "double _s_RLB = -s_inf;\n"
   "//get the other differential updates\n"
   "double Cm = 0.185;\n"
   "double F = 96485.341499999995;\n"
   "double R = 8314.4719999999998;\n"
   "double T = 310;\n"
   "double V_c = 0.016403999999999998;\n"
   "double factor_fix = 1;\n"
   "double Ca_o = 2;\n"
   "double Buf_c = 0.20000000000000001;\n"
   "double Buf_sr = 10;\n"
   "double Buf_ss = 0.40000000000000002;\n"
   "double EC = 1.5;\n"
   "double K_buf_c = 0.001;\n"
   "double K_buf_sr = 0.29999999999999999;\n"
   "double K_buf_ss = 0.00025000000000000001;\n"
   "double K_up = 0.00025000000000000001;\n"
   "double V_leak = 0.00036000000000000002;\n"
   "double V_rel = 0.10199999999999999;\n"
   "double V_sr = 0.0010939999999999999;\n"
   "double V_ss = 5.4679999999999998e-5;\n"
   "double V_xfer = 0.0038;\n"
   "double Vmax_up = 0.0063749999999999996;\n"
   "double k1_prime = 0.14999999999999999;\n"
   "double k2_prime = 0.044999999999999998;\n"
   "double k3 = 0.059999999999999998;\n"
   "double k4 = 0.0050000000000000001;\n"
   "double max_sr = 2.5;\n"
   "double min_sr = 1;\n"
   "double i_CalTerm1 = (F*F)*(4*V - 60)/(R*T);\n"
   "double i_CalTerm2 = exp(F*(2*V - 30)/(R*T));\n"
   "double __melodee_temp_005 = V == 15;\n"
   "double i_CalTerm3;\n"
   "if (__melodee_temp_005)\n"
   "{\n"
   "   i_CalTerm3 = 2*F;\n"
   "}\n"
   "else\n"
   "{\n"
   "   i_CalTerm3 = i_CalTerm1/(i_CalTerm2 - 1);\n"
   "}\n"
   "double i_CalTerm4 = i_CalTerm2*i_CalTerm3;\n"
   "double _expensive_functions_012 = log(Ca_o/Ca_i);\n"
   "double E_Ca = 0.5*R*T*_expensive_functions_012/F;\n"
   "double i_b_Ca = g_bca*(-E_Ca + V);\n"
   "double K_pCa = 0.00050000000000000001;\n"
   "double i_p_Ca = Ca_i*g_pCa/(Ca_i + K_pCa);\n"
   "double K_NaCa = 1000;\n"
   "double K_sat = 0.10000000000000001;\n"
   "double Km_Ca = 1.3799999999999999;\n"
   "double Km_Nai = 87.5;\n"
   "double alpha = 2.5;\n"
   "double gamma = 0.34999999999999998;\n"
   "double exp_gamma_VFRT = exp(F*V*gamma/(R*T));\n"
   "double exp_gamma_m1_VFRT = exp(F*V*(gamma - 1)/(R*T));\n"
   "double K_o = 5.4000000000000004;\n"
   "double _expensive_functions_019 = exp(-0.16722408026755853*V);\n"
   "double i_p_K_term = (1.0/(65.405215741938321*_expensive_functions_019 + 1));\n"
   "double _expensive_functions_026 = log(K_o/K_i);\n"
   "double E_K = R*T*_expensive_functions_026/F;\n"
   "double P_kna = 0.029999999999999999;\n"
   "double Na_o = 140;\n"
   "double K_mNa = 40;\n"
   "double K_mk = 1;\n"
   "double P_NaK = 2.7240000000000002;\n"
   "double _expensive_functions_031 = exp(-F*V/(R*T));\n"
   "double _expensive_functions_032 = exp(-0.10000000000000001*F*V/(R*T));\n"
   "double i_NaK_term = K_o*P_NaK/((K_mk + K_o)*(0.035299999999999998*_expensive_functions_031 + 0.1245*_expensive_functions_032 + 1));\n"
   "double i_NaK = Na_i*i_NaK_term/(K_mNa + Na_i);\n"
   "double i_Naitot = 3*i_NaK;\n"
   "double i_Kitot = -2*i_NaK;\n"
   "double VEK = -E_K + V;\n"
   "double Ca_i_bufc = (1.0/(Buf_c*K_buf_c/((Ca_i + K_buf_c)*(Ca_i + K_buf_c)) + 1));\n"
   "double Ca_sr_bufsr = (1.0/(Buf_sr*K_buf_sr/((Ca_SR + K_buf_sr)*(Ca_SR + K_buf_sr)) + 1));\n"
   "double Ca_ss_bufss = (1.0/(Buf_ss*K_buf_ss/((Ca_ss + K_buf_ss)*(Ca_ss + K_buf_ss)) + 1));\n"
   "double i_leak = V_leak*(Ca_SR - Ca_i);\n"
   "double i_up = Vmax_up/(1 + (K_up*K_up)/(Ca_i*Ca_i));\n"
   "double i_xfer = V_xfer*(-Ca_i + Ca_ss);\n"
   "double kcasr = max_sr - (max_sr - min_sr)/(1 + (EC*EC)/(Ca_SR*Ca_SR));\n"
   "double k1 = k1_prime/kcasr;\n"
   "double k2 = k2_prime*kcasr;\n"
   "double O = (Ca_ss*Ca_ss)*R_prime*k1/((Ca_ss*Ca_ss)*k1 + k3);\n"
   "double R_prime_diff = -Ca_ss*R_prime*k2 + k4*(1 - R_prime);\n"
   "double i_rel = O*V_rel*(Ca_SR - Ca_ss);\n"
   "double Ca_SR_diff = Ca_sr_bufsr*(-i_leak - i_rel + i_up);\n"
   "double i_CaL = d*f*f2*fCass*g_CaL*(-Ca_o*i_CalTerm3 + 0.25*Ca_ss*i_CalTerm4);\n"
   "double i_NaCa = K_NaCa*(-(Na_o*Na_o*Na_o)*Ca_i*alpha*exp_gamma_m1_VFRT + (Na_i*Na_i*Na_i)*Ca_o*exp_gamma_VFRT)/(((Km_Nai*Km_Nai*Km_Nai) + (Na_o*Na_o*Na_o))*(Ca_o + Km_Ca)*(K_sat*exp_gamma_m1_VFRT + 1));\n"
   "double sodium_calcium_exchanger_current_i_Naitot = 3*i_NaCa;\n"
   ""; generateInterpString(ss,_interpolant[1], "VEK"); ss << "\n"
   "double inward_rectifier_potassium_current_i_Kitot = _ratPoly;\n"
   "double i_p_K = VEK*g_pK*i_p_K_term;\n"
   "double potassium_pump_current_i_Kitot = i_p_K;\n"
   "double _expensive_functions_040 = sqrt(K_o);\n"
   "double i_Kr = 0.43033148291193518*VEK*Xr1*Xr2*_expensive_functions_040*g_Kr;\n"
   "double rapid_time_dependent_potassium_current_i_Kitot = i_Kr;\n"
   "double _expensive_functions_041 = log(Na_o/Na_i);\n"
   "double E_Na = R*T*_expensive_functions_041/F;\n"
   "double _expensive_functions_042 = log((K_o + Na_o*P_kna)/(K_i + Na_i*P_kna));\n"
   "double E_Ks = R*T*_expensive_functions_042/F;\n"
   "double i_Ks = (Xs*Xs)*g_Ks*(-E_Ks + V);\n"
   "double slow_time_dependent_potassium_current_i_Kitot = i_Ks;\n"
   "double i_b_Na = g_bna*(-E_Na + V);\n"
   "double sodium_background_current_i_Naitot = i_b_Na;\n"
   "double i_to = VEK*g_to*r*s;\n"
   "double transient_outward_current_i_Kitot = i_to;\n"
   "double i_Kitot_001 = i_Kitot + inward_rectifier_potassium_current_i_Kitot + potassium_pump_current_i_Kitot + rapid_time_dependent_potassium_current_i_Kitot + slow_time_dependent_potassium_current_i_Kitot + transient_outward_current_i_Kitot;\n"
   "double Ca_i_diff = Ca_i_bufc*(-1.0/2.0*Cm*(-2*i_NaCa + i_b_Ca + i_p_Ca)/(F*V_c) + i_xfer + V_sr*(i_leak - i_up)/V_c);\n"
   "double Ca_ss_diff = Ca_ss_bufss*(-1.0/2.0*Cm*i_CaL/(F*V_ss) - V_c*i_xfer/V_ss + V_sr*i_rel/V_ss);\n"
   "double i_Na = (m*m*m)*g_Na*h*j*(-E_Na + V);\n"
   "double fast_sodium_current_i_Naitot = i_Na;\n"
   "double K_i_diff = -Cm*factor_fix*i_Kitot_001/(F*V_c);\n"
   "double i_Naitot_001 = fast_sodium_current_i_Naitot + i_Naitot + sodium_background_current_i_Naitot + sodium_calcium_exchanger_current_i_Naitot;\n"
   "double Na_i_diff = -Cm*factor_fix*i_Naitot_001/(F*V_c);\n"
   "//get Iion\n"
   "double i_Caitot = i_b_Ca;\n"
   "double calcium_pump_current_i_Caitot = i_p_Ca;\n"
   "double L_type_Ca_current_i_Caitot = i_CaL;\n"
   "double sodium_calcium_exchanger_current_i_Caitot = -2*i_NaCa;\n"
   "double i_Caitot_001 = L_type_Ca_current_i_Caitot + calcium_pump_current_i_Caitot + i_Caitot + sodium_calcium_exchanger_current_i_Caitot;\n"
   "double Iion = i_Caitot_001 + i_Kitot_001 + i_Naitot_001;\n"
   "double Iion_001 = Iion;\n"
   "//Do the markov update (1 step rosenbrock with gauss siedel)\n"
   "int _count=0;\n"
   "do\n"
   "{\n"
   "   _count++;\n"
   "} while (_count<50);\n"
   "//EDIT_STATE\n"
   "_state[_ii+Ca_SR_off*_nCells] += _dt*Ca_SR_diff;\n"
   "_state[_ii+Ca_i_off*_nCells] += _dt*Ca_i_diff;\n"
   "_state[_ii+Ca_ss_off*_nCells] += _dt*Ca_ss_diff;\n"
   "_state[_ii+K_i_off*_nCells] += _dt*K_i_diff;\n"
   "_state[_ii+Na_i_off*_nCells] += _dt*Na_i_diff;\n"
   "_state[_ii+R_prime_off*_nCells] += _dt*R_prime_diff;\n"
   "_state[_ii+Xr1_off*_nCells] += _Xr1_RLA*(Xr1+_Xr1_RLB);\n"
   "_state[_ii+Xr2_off*_nCells] += _Xr2_RLA*(Xr2+_Xr2_RLB);\n"
   "_state[_ii+Xs_off*_nCells] += _Xs_RLA*(Xs+_Xs_RLB);\n"
   "_state[_ii+d_off*_nCells] += _d_RLA*(d+_d_RLB);\n"
   "_state[_ii+f_off*_nCells] += _f_RLA*(f+_f_RLB);\n"
   "_state[_ii+f2_off*_nCells] += _f2_RLA*(f2+_f2_RLB);\n"
   "_state[_ii+fCass_off*_nCells] += _fCass_RLA*(fCass+_fCass_RLB);\n"
   "_state[_ii+h_off*_nCells] += _h_RLA*(h+_h_RLB);\n"
   "_state[_ii+j_off*_nCells] += _j_RLA*(j+_j_RLB);\n"
   "_state[_ii+m_off*_nCells] += _m_RLA*(m+_m_RLB);\n"
   "_state[_ii+r_off*_nCells] += _r_RLA*(r+_r_RLB);\n"
   "_state[_ii+s_off*_nCells] += _s_RLA*(s+_s_RLB);\n"
   "_dVm[_ii] = -Iion_001;\n"
   "}\n";

   _program_code = ss.str();
   //cout << ss.str();
   nvrtcCreateProgram(&_program,
                      _program_code.c_str(),
                      "TT06_SCC_program",
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
   cuModuleGetFunction(&_kernel, _module, "TT06_SCC_kernel");
}

void ThisReaction::calc(double dt,
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
         const double* VmRaw = Vm_m.useOn(GPU).raw();
         const double* iStimRaw = iStim_m.useOn(GPU).raw();
         double* dVmRaw = dVm_m.useOn(GPU).raw();
         double* stateRaw= stateTransport_.readwrite(GPU).raw();
         void* args[] = { &VmRaw,
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
   _Ca_SR_off,
   _Ca_i_off,
   _Ca_ss_off,
   _K_i_off,
   _Na_i_off,
   _R_prime_off,
   _Xr1_off,
   _Xr2_off,
   _Xs_off,
   _d_off,
   _f_off,
   _f2_off,
   _fCass_off,
   _h_off,
   _j_off,
   _m_off,
   _r_off,
   _s_off,
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
                ro_mgarray_ptr<double> ___Vm,
                ro_mgarray_ptr<double> ___iStim,
                wo_mgarray_ptr<double> ___dVm)
{
   ro_array_ptr<double> __Vm = ___Vm.useOn(CPU);
   ro_array_ptr<double> __iStim = ___iStim.useOn(CPU);
   wo_array_ptr<double> __dVm = ___dVm.useOn(CPU);

   //define the constants
   double Cm = 0.185;
   double F = 96485.341499999995;
   double R = 8314.4719999999998;
   double T = 310;
   double V_c = 0.016403999999999998;
   double factor_fix = 1;
   double Ca_o = 2;
   double Buf_c = 0.20000000000000001;
   double Buf_sr = 10;
   double Buf_ss = 0.40000000000000002;
   double EC = 1.5;
   double K_buf_c = 0.001;
   double K_buf_sr = 0.29999999999999999;
   double K_buf_ss = 0.00025000000000000001;
   double K_up = 0.00025000000000000001;
   double V_leak = 0.00036000000000000002;
   double V_rel = 0.10199999999999999;
   double V_sr = 0.0010939999999999999;
   double V_ss = 5.4679999999999998e-5;
   double V_xfer = 0.0038;
   double Vmax_up = 0.0063749999999999996;
   double k1_prime = 0.14999999999999999;
   double k2_prime = 0.044999999999999998;
   double k3 = 0.059999999999999998;
   double k4 = 0.0050000000000000001;
   double max_sr = 2.5;
   double min_sr = 1;
   double K_pCa = 0.00050000000000000001;
   double K_NaCa = 1000;
   double K_sat = 0.10000000000000001;
   double Km_Ca = 1.3799999999999999;
   double Km_Nai = 87.5;
   double alpha = 2.5;
   double gamma = 0.34999999999999998;
   double K_o = 5.4000000000000004;
   double P_kna = 0.029999999999999999;
   double Na_o = 140;
   double K_mNa = 40;
   double K_mk = 1;
   double P_NaK = 2.7240000000000002;
   double __melodee_temp_003 = celltype == 0;
   double _expensive_functions_040 = sqrt(K_o);
   for (unsigned __jj=0; __jj<(nCells_+width-1)/width; __jj++)
   {
      const int __ii = __jj*width;
      //set Vm
      const real V = load(&__Vm[__ii]);
      const real iStim = load(&__iStim[__ii]);

      //set all state variables
      real Ca_SR=load(state_[__jj].Ca_SR);
      real Ca_i=load(state_[__jj].Ca_i);
      real Ca_ss=load(state_[__jj].Ca_ss);
      real K_i=load(state_[__jj].K_i);
      real Na_i=load(state_[__jj].Na_i);
      real R_prime=load(state_[__jj].R_prime);
      real Xr1=load(state_[__jj].Xr1);
      real Xr2=load(state_[__jj].Xr2);
      real Xs=load(state_[__jj].Xs);
      real d=load(state_[__jj].d);
      real f=load(state_[__jj].f);
      real f2=load(state_[__jj].f2);
      real fCass=load(state_[__jj].fCass);
      real h=load(state_[__jj].h);
      real j=load(state_[__jj].j);
      real m=load(state_[__jj].m);
      real r=load(state_[__jj].r);
      real s=load(state_[__jj].s);
      //get the gate updates (diagonalized exponential integrator)
      real _expensive_functions = exp(-1.0/13.0*V - 35.0/13.0);
      real alpha_d = 0.25 + 1.3999999999999999/(_expensive_functions + 1);
      real _expensive_functions_001 = exp((1.0/5.0)*V + 1);
      real beta_d = 1.3999999999999999/(_expensive_functions_001 + 1);
      real _expensive_functions_002 = exp(-0.13333333333333333*V);
      real d_inf = (1.0/(0.34415378686541237*_expensive_functions_002 + 1));
      real _expensive_functions_003 = exp(5.0/2.0 - 1.0/20.0*V);
      real gamma_d = (1.0/(_expensive_functions_003 + 1));
      real tau_d = alpha_d*beta_d + gamma_d;
      real _expensive_functions_004 = exp((1.0/7.0)*V + 5);
      real f2_inf = 0.33000000000000002 + 0.67000000000000004/(_expensive_functions_004 + 1);
      real _expensive_functions_005 = exp(5.0/2.0 - 1.0/10.0*V);
      real _expensive_functions_006 = exp((1.0/10.0)*V + 3);
      real _expensive_functions_007 = exp(-1.0/240.0*((V + 27)*(V + 27)));
      real tau_f2 = 562*_expensive_functions_007 + 80/(_expensive_functions_006 + 1) + 31/(_expensive_functions_005 + 1);
      real fCass_inf = 0.40000000000000002 + 0.59999999999999998/(400.0*(Ca_ss*Ca_ss) + 1);
      real _expensive_functions_008 = exp((1.0/7.0)*V + 20.0/7.0);
      real f_inf = (1.0/(_expensive_functions_008 + 1));
      real _expensive_functions_009 = exp((1.0/10.0)*V + 3);
      real _expensive_functions_010 = exp(13.0/10.0 - 1.0/10.0*V);
      real _expensive_functions_011 = exp(-1.0/225.0*((V + 27)*(V + 27)));
      real tau_f = 1102.5*_expensive_functions_011 + 20 + 200/(_expensive_functions_010 + 1) + 180/(_expensive_functions_009 + 1);
      real __melodee_temp_000 = V < -40;
      real alpha_h;
      real beta_h;
      if (__melodee_temp_000)
      {
         real _expensive_functions_013 = exp(-0.14705882352941177*V);
         alpha_h = 4.4312679295805147e-7*_expensive_functions_013;
         real _expensive_functions_014 = exp(0.34849999999999998*V);
         real _expensive_functions_015 = exp(0.079000000000000001*V);
         beta_h = 310000*_expensive_functions_014 + 2.7000000000000002*_expensive_functions_015;
      }
      else
      {
         alpha_h = 0;
         real _expensive_functions_013 = exp(-0.0900900900900901*V);
         beta_h = 0.77000000000000002/(0.049758141083938695*_expensive_functions_013 + 0.13);
      }
      real _expensive_functions_013 = exp(0.13458950201884254*V);
      real h_inf = 4.3210917837689708e-9/((_expensive_functions_013 + 6.5735011856460265e-5)*(_expensive_functions_013 + 6.5735011856460265e-5));
      real tau_h = (1.0/(alpha_h + beta_h));
      real __melodee_temp_001 = V < -40;
      real alpha_j;
      real beta_j;
      if (__melodee_temp_001)
      {
         real _expensive_functions_014 = exp(0.311*V);
         real _expensive_functions_015 = exp(0.24440000000000001*V);
         real _expensive_functions_016 = exp(-0.043909999999999998*V);
         alpha_j = (-25428*_expensive_functions_015 - 6.9480000000000002e-6*_expensive_functions_016)*(V + 37.780000000000001)/(50262745825.953987*_expensive_functions_014 + 1);
         real _expensive_functions_017 = exp(-0.13780000000000001*V);
         real _expensive_functions_018 = exp(-0.01052*V);
         beta_j = 0.024240000000000001*_expensive_functions_018/(0.003960868339904256*_expensive_functions_017 + 1);
      }
      else
      {
         alpha_j = 0;
         real _expensive_functions_014 = exp(-0.10000000000000001*V);
         real _expensive_functions_015 = exp(0.057000000000000002*V);
         beta_j = 0.59999999999999998*_expensive_functions_015/(0.040762203978366204*_expensive_functions_014 + 1);
      }
      real _expensive_functions_014 = exp(0.13458950201884254*V);
      real j_inf = 4.3210917837689708e-9/((_expensive_functions_014 + 6.5735011856460265e-5)*(_expensive_functions_014 + 6.5735011856460265e-5));
      real tau_j = (1.0/(alpha_j + beta_j));
      real _expensive_functions_015 = exp(-1.0/5.0*V - 12);
      real alpha_m = (1.0/(_expensive_functions_015 + 1));
      real _expensive_functions_016 = exp((1.0/5.0)*V + 7);
      real _expensive_functions_017 = exp((1.0/200.0)*V - 1.0/4.0);
      real beta_m = 0.10000000000000001/(_expensive_functions_017 + 1) + 0.10000000000000001/(_expensive_functions_016 + 1);
      real _expensive_functions_018 = exp(-0.11074197120708749*V);
      real m_inf = (1.0/(0.0018422115811651339*_expensive_functions_018 + 1)/(0.0018422115811651339*_expensive_functions_018 + 1));
      real tau_m = alpha_m*beta_m;
      real _expensive_functions_020 = exp(-1.0/10.0*V - 9.0/2.0);
      real alpha_xr1 = 450/(_expensive_functions_020 + 1);
      real _expensive_functions_021 = exp(0.086956521739130432*V);
      real beta_xr1 = 6/(13.581324522578193*_expensive_functions_021 + 1);
      real _expensive_functions_022 = exp(-1.0/7.0*V - 26.0/7.0);
      real xr1_inf = (1.0/(_expensive_functions_022 + 1));
      real tau_xr1 = alpha_xr1*beta_xr1;
      real _expensive_functions_023 = exp(-1.0/20.0*V - 3);
      real alpha_xr2 = 3/(_expensive_functions_023 + 1);
      real _expensive_functions_024 = exp((1.0/20.0)*V - 3);
      real beta_xr2 = 1.1200000000000001/(_expensive_functions_024 + 1);
      real _expensive_functions_025 = exp((1.0/24.0)*V + 11.0/3.0);
      real xr2_inf = (1.0/(_expensive_functions_025 + 1));
      real tau_xr2 = alpha_xr2*beta_xr2;
      real _expensive_functions_027 = exp(5.0/6.0 - 1.0/6.0*V);
      real _expensive_functions_028 = pow(_expensive_functions_027 + 1, -1.0/2.0);
      real alpha_xs = 1400*_expensive_functions_028;
      real _expensive_functions_029 = exp((1.0/15.0)*V - 7.0/3.0);
      real beta_xs = (1.0/(_expensive_functions_029 + 1));
      real _expensive_functions_030 = exp(-1.0/14.0*V - 5.0/14.0);
      real xs_inf = (1.0/(_expensive_functions_030 + 1));
      real tau_xs = alpha_xs*beta_xs + 80;
      real _expensive_functions_033 = exp(10.0/3.0 - 1.0/6.0*V);
      real r_inf = (1.0/(_expensive_functions_033 + 1));
      real _expensive_functions_034 = exp(-1.0/1800.0*((V + 40)*(V + 40)));
      real tau_r = 9.5*_expensive_functions_034 + 0.80000000000000004;
      real s_inf;
      real tau_s;
      if (__melodee_temp_003)
      {
         real _expensive_functions_035 = exp((1.0/5.0)*V + 28.0/5.0);
         s_inf = (1.0/(_expensive_functions_035 + 1));
         real _expensive_functions_036 = exp(-1.0/1000.0*((V + 67)*(V + 67)));
         tau_s = 1000*_expensive_functions_036 + 8;
      }
      else
      {
         real _expensive_functions_035 = exp((1.0/5.0)*V + 4);
         s_inf = (1.0/(_expensive_functions_035 + 1));
         real _expensive_functions_036 = exp((1.0/5.0)*V - 4);
         real _expensive_functions_037 = exp(-1.0/320.0*((V + 45)*(V + 45)));
         tau_s = 85*_expensive_functions_037 + 3 + 5/(_expensive_functions_036 + 1);
      }
      real _expensive_functions_043 = exp(-_dt/tau_xr1);
      real _Xr1_RLA = _expensive_functions_043 - 1;
      real _Xr1_RLB = -xr1_inf;
      real _expensive_functions_044 = exp(-_dt/tau_xr2);
      real _Xr2_RLA = _expensive_functions_044 - 1;
      real _Xr2_RLB = -xr2_inf;
      real _expensive_functions_045 = exp(-_dt/tau_xs);
      real _Xs_RLA = _expensive_functions_045 - 1;
      real _Xs_RLB = -xs_inf;
      real _expensive_functions_046 = exp(-_dt/tau_d);
      real _d_RLA = _expensive_functions_046 - 1;
      real _d_RLB = -d_inf;
      real _expensive_functions_047 = exp(-_dt/tau_f);
      real _f_RLA = _expensive_functions_047 - 1;
      real _f_RLB = -f_inf;
      real _expensive_functions_048 = exp(-_dt/tau_f2);
      real _f2_RLA = _expensive_functions_048 - 1;
      real _f2_RLB = -f2_inf;
      real _fCass_RLA = _interpolant[0].eval(Ca_ss);
      real _fCass_RLB = -fCass_inf;
      real _expensive_functions_050 = exp(-_dt/tau_h);
      real _h_RLA = _expensive_functions_050 - 1;
      real _h_RLB = -h_inf;
      real _expensive_functions_051 = exp(-_dt/tau_j);
      real _j_RLA = _expensive_functions_051 - 1;
      real _j_RLB = -j_inf;
      real _expensive_functions_052 = exp(-_dt/tau_m);
      real _m_RLA = _expensive_functions_052 - 1;
      real _m_RLB = -m_inf;
      real _expensive_functions_053 = exp(-_dt/tau_r);
      real _r_RLA = _expensive_functions_053 - 1;
      real _r_RLB = -r_inf;
      real _expensive_functions_054 = exp(-_dt/tau_s);
      real _s_RLA = _expensive_functions_054 - 1;
      real _s_RLB = -s_inf;
      //get the other differential updates
      real i_CalTerm1 = (F*F)*(4*V - 60)/(R*T);
      real i_CalTerm2 = exp(F*(2*V - 30)/(R*T));
      real __melodee_temp_005 = V == 15;
      real i_CalTerm3;
      if (__melodee_temp_005)
      {
         i_CalTerm3 = 2*F;
      }
      else
      {
         i_CalTerm3 = i_CalTerm1/(i_CalTerm2 - 1);
      }
      real i_CalTerm4 = i_CalTerm2*i_CalTerm3;
      real _expensive_functions_012 = log(Ca_o/Ca_i);
      real E_Ca = 0.5*R*T*_expensive_functions_012/F;
      real i_b_Ca = g_bca*(-E_Ca + V);
      real i_p_Ca = Ca_i*g_pCa/(Ca_i + K_pCa);
      real exp_gamma_VFRT = exp(F*V*gamma/(R*T));
      real exp_gamma_m1_VFRT = exp(F*V*(gamma - 1)/(R*T));
      real _expensive_functions_019 = exp(-0.16722408026755853*V);
      real i_p_K_term = (1.0/(65.405215741938321*_expensive_functions_019 + 1));
      real _expensive_functions_026 = log(K_o/K_i);
      real E_K = R*T*_expensive_functions_026/F;
      real _expensive_functions_031 = exp(-F*V/(R*T));
      real _expensive_functions_032 = exp(-0.10000000000000001*F*V/(R*T));
      real i_NaK_term = K_o*P_NaK/((K_mk + K_o)*(0.035299999999999998*_expensive_functions_031 + 0.1245*_expensive_functions_032 + 1));
      real i_NaK = Na_i*i_NaK_term/(K_mNa + Na_i);
      real i_Naitot = 3*i_NaK;
      real i_Kitot = -2*i_NaK;
      real VEK = -E_K + V;
      real Ca_i_bufc = (1.0/(Buf_c*K_buf_c/((Ca_i + K_buf_c)*(Ca_i + K_buf_c)) + 1));
      real Ca_sr_bufsr = (1.0/(Buf_sr*K_buf_sr/((Ca_SR + K_buf_sr)*(Ca_SR + K_buf_sr)) + 1));
      real Ca_ss_bufss = (1.0/(Buf_ss*K_buf_ss/((Ca_ss + K_buf_ss)*(Ca_ss + K_buf_ss)) + 1));
      real i_leak = V_leak*(Ca_SR - Ca_i);
      real i_up = Vmax_up/(1 + (K_up*K_up)/(Ca_i*Ca_i));
      real i_xfer = V_xfer*(-Ca_i + Ca_ss);
      real kcasr = max_sr - (max_sr - min_sr)/(1 + (EC*EC)/(Ca_SR*Ca_SR));
      real k1 = k1_prime/kcasr;
      real k2 = k2_prime*kcasr;
      real O = (Ca_ss*Ca_ss)*R_prime*k1/((Ca_ss*Ca_ss)*k1 + k3);
      real R_prime_diff = -Ca_ss*R_prime*k2 + k4*(1 - R_prime);
      real i_rel = O*V_rel*(Ca_SR - Ca_ss);
      real Ca_SR_diff = Ca_sr_bufsr*(-i_leak - i_rel + i_up);
      real i_CaL = d*f*f2*fCass*g_CaL*(-Ca_o*i_CalTerm3 + 0.25*Ca_ss*i_CalTerm4);
      real i_NaCa = K_NaCa*(-(Na_o*Na_o*Na_o)*Ca_i*alpha*exp_gamma_m1_VFRT + (Na_i*Na_i*Na_i)*Ca_o*exp_gamma_VFRT)/(((Km_Nai*Km_Nai*Km_Nai) + (Na_o*Na_o*Na_o))*(Ca_o + Km_Ca)*(K_sat*exp_gamma_m1_VFRT + 1));
      real sodium_calcium_exchanger_current_i_Naitot = 3*i_NaCa;
      real inward_rectifier_potassium_current_i_Kitot = _interpolant[1].eval(VEK);
      real i_p_K = VEK*g_pK*i_p_K_term;
      real potassium_pump_current_i_Kitot = i_p_K;
      real i_Kr = 0.43033148291193518*VEK*Xr1*Xr2*_expensive_functions_040*g_Kr;
      real rapid_time_dependent_potassium_current_i_Kitot = i_Kr;
      real _expensive_functions_041 = log(Na_o/Na_i);
      real E_Na = R*T*_expensive_functions_041/F;
      real _expensive_functions_042 = log((K_o + Na_o*P_kna)/(K_i + Na_i*P_kna));
      real E_Ks = R*T*_expensive_functions_042/F;
      real i_Ks = (Xs*Xs)*g_Ks*(-E_Ks + V);
      real slow_time_dependent_potassium_current_i_Kitot = i_Ks;
      real i_b_Na = g_bna*(-E_Na + V);
      real sodium_background_current_i_Naitot = i_b_Na;
      real i_to = VEK*g_to*r*s;
      real transient_outward_current_i_Kitot = i_to;
      real i_Kitot_001 = i_Kitot + inward_rectifier_potassium_current_i_Kitot + potassium_pump_current_i_Kitot + rapid_time_dependent_potassium_current_i_Kitot + slow_time_dependent_potassium_current_i_Kitot + transient_outward_current_i_Kitot;
      real Ca_i_diff = Ca_i_bufc*(-1.0/2.0*Cm*(-2*i_NaCa + i_b_Ca + i_p_Ca)/(F*V_c) + i_xfer + V_sr*(i_leak - i_up)/V_c);
      real Ca_ss_diff = Ca_ss_bufss*(-1.0/2.0*Cm*i_CaL/(F*V_ss) - V_c*i_xfer/V_ss + V_sr*i_rel/V_ss);
      real i_Na = (m*m*m)*g_Na*h*j*(-E_Na + V);
      real fast_sodium_current_i_Naitot = i_Na;
      real K_i_diff = -Cm*factor_fix*i_Kitot_001/(F*V_c);
      real i_Naitot_001 = fast_sodium_current_i_Naitot + i_Naitot + sodium_background_current_i_Naitot + sodium_calcium_exchanger_current_i_Naitot;
      real Na_i_diff = -Cm*factor_fix*i_Naitot_001/(F*V_c);
      //get Iion
      real i_Caitot = i_b_Ca;
      real calcium_pump_current_i_Caitot = i_p_Ca;
      real L_type_Ca_current_i_Caitot = i_CaL;
      real sodium_calcium_exchanger_current_i_Caitot = -2*i_NaCa;
      real i_Caitot_001 = L_type_Ca_current_i_Caitot + calcium_pump_current_i_Caitot + i_Caitot + sodium_calcium_exchanger_current_i_Caitot;
      real Iion = i_Caitot_001 + i_Kitot_001 + i_Naitot_001;
      real Iion_001 = Iion;
      //Do the markov update (1 step rosenbrock with gauss siedel)
      //EDIT_STATE
      Ca_SR += _dt*Ca_SR_diff;
      Ca_i += _dt*Ca_i_diff;
      Ca_ss += _dt*Ca_ss_diff;
      K_i += _dt*K_i_diff;
      Na_i += _dt*Na_i_diff;
      R_prime += _dt*R_prime_diff;
      Xr1 += _Xr1_RLA*(Xr1+_Xr1_RLB);
      Xr2 += _Xr2_RLA*(Xr2+_Xr2_RLB);
      Xs += _Xs_RLA*(Xs+_Xs_RLB);
      d += _d_RLA*(d+_d_RLB);
      f += _f_RLA*(f+_f_RLB);
      f2 += _f2_RLA*(f2+_f2_RLB);
      fCass += _fCass_RLA*(fCass+_fCass_RLB);
      h += _h_RLA*(h+_h_RLB);
      j += _j_RLA*(j+_j_RLB);
      m += _m_RLA*(m+_m_RLB);
      r += _r_RLA*(r+_r_RLB);
      s += _s_RLA*(s+_s_RLB);
      store(state_[__jj].Ca_SR, Ca_SR);
      store(state_[__jj].Ca_i, Ca_i);
      store(state_[__jj].Ca_ss, Ca_ss);
      store(state_[__jj].K_i, K_i);
      store(state_[__jj].Na_i, Na_i);
      store(state_[__jj].R_prime, R_prime);
      store(state_[__jj].Xr1, Xr1);
      store(state_[__jj].Xr2, Xr2);
      store(state_[__jj].Xs, Xs);
      store(state_[__jj].d, d);
      store(state_[__jj].f, f);
      store(state_[__jj].f2, f2);
      store(state_[__jj].fCass, fCass);
      store(state_[__jj].h, h);
      store(state_[__jj].j, j);
      store(state_[__jj].m, m);
      store(state_[__jj].r, r);
      store(state_[__jj].s, s);
      simdops::store(&__dVm.raw()[__ii],-Iion_001);
   }
}
#endif //USE_CUDA
   
string ThisReaction::methodName() const
{
   return "TT06_SCC";
}

void ThisReaction::initializeMembraneVoltage(wo_mgarray_ptr<double> __Vm_m)
{
   assert(__Vm_m.size() >= nCells_);

   wo_array_ptr<double> __Vm = __Vm_m.useOn(CPU);
#ifdef USE_CUDA
#define READ_STATE(state,index) (stateData[_##state##_off*nCells_+index])
   wo_array_ptr<double> stateData = stateTransport_.useOn(CPU);
#else //USE_CUDA
#define READ_STATE(state,index) (state_[index/width].state[index % width])
   state_.resize((nCells_+width-1)/width);
#endif //USE_CUDA


   double V_init = -83;
   double V = V_init;
   double Ca_i_init = 2.0000000000000002e-5;
   double Ca_i = Ca_i_init;
   double R_prime_init = 0.98680000000000001;
   double R_prime = R_prime_init;
   double Ca_SR_init = 3.1549999999999998;
   double Ca_SR = Ca_SR_init;
   double Ca_ss_init = 0.00017000000000000001;
   double Ca_ss = Ca_ss_init;
   double d_init = 3.1640000000000002e-5;
   double d = d_init;
   double f2_init = 0.9778;
   double f2 = f2_init;
   double fCass_init = 0.99529999999999996;
   double fCass = fCass_init;
   double f_init = 0.96089999999999998;
   double f = f_init;
   double h_init = 0.55000000000000004;
   double h = h_init;
   double j_init = 0.66000000000000003;
   double j = j_init;
   double m_init = 0.0015499999999999999;
   double m = m_init;
   double K_i_init = 138.40000000000001;
   double K_i = K_i_init;
   double Xr1_init = 0.0044799999999999996;
   double Xr1 = Xr1_init;
   double Xr2_init = 0.47599999999999998;
   double Xr2 = Xr2_init;
   double Xs_init = 0.0086999999999999994;
   double Xs = Xs_init;
   double Na_i_init = 10.355;
   double Na_i = Na_i_init;
   double r_init = 2.2350000000000002e-8;
   double r = r_init;
   double s_init = 0.60119999999999996;
   double s = s_init;
   for (int iCell=0; iCell<nCells_; iCell++)
   {
      READ_STATE(Ca_SR,iCell) = Ca_SR;
      READ_STATE(Ca_i,iCell) = Ca_i;
      READ_STATE(Ca_ss,iCell) = Ca_ss;
      READ_STATE(K_i,iCell) = K_i;
      READ_STATE(Na_i,iCell) = Na_i;
      READ_STATE(R_prime,iCell) = R_prime;
      READ_STATE(Xr1,iCell) = Xr1;
      READ_STATE(Xr2,iCell) = Xr2;
      READ_STATE(Xs,iCell) = Xs;
      READ_STATE(d,iCell) = d;
      READ_STATE(f,iCell) = f;
      READ_STATE(f2,iCell) = f2;
      READ_STATE(fCass,iCell) = fCass;
      READ_STATE(h,iCell) = h;
      READ_STATE(j,iCell) = j;
      READ_STATE(m,iCell) = m;
      READ_STATE(r,iCell) = r;
      READ_STATE(s,iCell) = s;
   }

   __Vm.assign(__Vm.size(), V_init);
}

enum varHandles
{
   Ca_SR_handle,
   Ca_i_handle,
   Ca_ss_handle,
   K_i_handle,
   Na_i_handle,
   R_prime_handle,
   Xr1_handle,
   Xr2_handle,
   Xs_handle,
   d_handle,
   f_handle,
   f2_handle,
   fCass_handle,
   h_handle,
   j_handle,
   m_handle,
   r_handle,
   s_handle,
   i_CaL_handle,
   i_K1_handle,
   i_Kr_handle,
   i_Ks_handle,
   i_Na_handle,
   i_NaCa_handle,
   i_NaK_handle,
   i_b_Ca_handle,
   i_b_Na_handle,
   i_leak_handle,
   i_p_Ca_handle,
   i_p_K_handle,
   i_rel_handle,
   i_to_handle,
   i_up_handle,
   i_xfer_handle,
   NUMHANDLES
};

const string ThisReaction::getUnit(const std::string& varName) const
{
   if(0) {}
   else if (varName == "Ca_SR") { return "uM"; }
   else if (varName == "Ca_i") { return "uM"; }
   else if (varName == "Ca_ss") { return "uM"; }
   else if (varName == "K_i") { return "mM"; }
   else if (varName == "Na_i") { return "mM"; }
   else if (varName == "R_prime") { return "1"; }
   else if (varName == "Xr1") { return "1"; }
   else if (varName == "Xr2") { return "1"; }
   else if (varName == "Xs") { return "1"; }
   else if (varName == "d") { return "1"; }
   else if (varName == "f") { return "1"; }
   else if (varName == "f2") { return "1"; }
   else if (varName == "fCass") { return "1"; }
   else if (varName == "h") { return "1"; }
   else if (varName == "i_CaL") { return "V/s"; }
   else if (varName == "i_K1") { return "INVALID"; }
   else if (varName == "i_Kr") { return "INVALID"; }
   else if (varName == "i_Ks") { return "INVALID"; }
   else if (varName == "i_Na") { return "INVALID"; }
   else if (varName == "i_NaCa") { return "V/s"; }
   else if (varName == "i_NaK") { return "INVALID"; }
   else if (varName == "i_b_Ca") { return "V/s"; }
   else if (varName == "i_b_Na") { return "INVALID"; }
   else if (varName == "i_leak") { return "INVALID"; }
   else if (varName == "i_p_Ca") { return "V/s"; }
   else if (varName == "i_p_K") { return "INVALID"; }
   else if (varName == "i_rel") { return "INVALID"; }
   else if (varName == "i_to") { return "INVALID"; }
   else if (varName == "i_up") { return "INVALID"; }
   else if (varName == "i_xfer") { return "INVALID"; }
   else if (varName == "j") { return "1"; }
   else if (varName == "m") { return "1"; }
   else if (varName == "r") { return "1"; }
   else if (varName == "s") { return "1"; }
   return "INVALID";
}

int ThisReaction::getVarHandle(const std::string& varName) const
{
   if (0) {}
   else if (varName == "Ca_SR") { return Ca_SR_handle; }
   else if (varName == "Ca_i") { return Ca_i_handle; }
   else if (varName == "Ca_ss") { return Ca_ss_handle; }
   else if (varName == "K_i") { return K_i_handle; }
   else if (varName == "Na_i") { return Na_i_handle; }
   else if (varName == "R_prime") { return R_prime_handle; }
   else if (varName == "Xr1") { return Xr1_handle; }
   else if (varName == "Xr2") { return Xr2_handle; }
   else if (varName == "Xs") { return Xs_handle; }
   else if (varName == "d") { return d_handle; }
   else if (varName == "f") { return f_handle; }
   else if (varName == "f2") { return f2_handle; }
   else if (varName == "fCass") { return fCass_handle; }
   else if (varName == "h") { return h_handle; }
   else if (varName == "i_CaL") { return i_CaL_handle; }
   else if (varName == "i_K1") { return i_K1_handle; }
   else if (varName == "i_Kr") { return i_Kr_handle; }
   else if (varName == "i_Ks") { return i_Ks_handle; }
   else if (varName == "i_Na") { return i_Na_handle; }
   else if (varName == "i_NaCa") { return i_NaCa_handle; }
   else if (varName == "i_NaK") { return i_NaK_handle; }
   else if (varName == "i_b_Ca") { return i_b_Ca_handle; }
   else if (varName == "i_b_Na") { return i_b_Na_handle; }
   else if (varName == "i_leak") { return i_leak_handle; }
   else if (varName == "i_p_Ca") { return i_p_Ca_handle; }
   else if (varName == "i_p_K") { return i_p_K_handle; }
   else if (varName == "i_rel") { return i_rel_handle; }
   else if (varName == "i_to") { return i_to_handle; }
   else if (varName == "i_up") { return i_up_handle; }
   else if (varName == "i_xfer") { return i_xfer_handle; }
   else if (varName == "j") { return j_handle; }
   else if (varName == "m") { return m_handle; }
   else if (varName == "r") { return r_handle; }
   else if (varName == "s") { return s_handle; }
   return -1;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readwrite(CPU);
#endif //USE_CUDA



   if (0) {}
   else if (varHandle == Ca_SR_handle) { READ_STATE(Ca_SR,iCell) = value; }
   else if (varHandle == Ca_i_handle) { READ_STATE(Ca_i,iCell) = value; }
   else if (varHandle == Ca_ss_handle) { READ_STATE(Ca_ss,iCell) = value; }
   else if (varHandle == K_i_handle) { READ_STATE(K_i,iCell) = value; }
   else if (varHandle == Na_i_handle) { READ_STATE(Na_i,iCell) = value; }
   else if (varHandle == R_prime_handle) { READ_STATE(R_prime,iCell) = value; }
   else if (varHandle == Xr1_handle) { READ_STATE(Xr1,iCell) = value; }
   else if (varHandle == Xr2_handle) { READ_STATE(Xr2,iCell) = value; }
   else if (varHandle == Xs_handle) { READ_STATE(Xs,iCell) = value; }
   else if (varHandle == d_handle) { READ_STATE(d,iCell) = value; }
   else if (varHandle == f_handle) { READ_STATE(f,iCell) = value; }
   else if (varHandle == f2_handle) { READ_STATE(f2,iCell) = value; }
   else if (varHandle == fCass_handle) { READ_STATE(fCass,iCell) = value; }
   else if (varHandle == h_handle) { READ_STATE(h,iCell) = value; }
   else if (varHandle == j_handle) { READ_STATE(j,iCell) = value; }
   else if (varHandle == m_handle) { READ_STATE(m,iCell) = value; }
   else if (varHandle == r_handle) { READ_STATE(r,iCell) = value; }
   else if (varHandle == s_handle) { READ_STATE(s,iCell) = value; }
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   if (0) {}
   else if (varHandle == Ca_SR_handle) { return READ_STATE(Ca_SR,iCell); }
   else if (varHandle == Ca_i_handle) { return READ_STATE(Ca_i,iCell); }
   else if (varHandle == Ca_ss_handle) { return READ_STATE(Ca_ss,iCell); }
   else if (varHandle == K_i_handle) { return READ_STATE(K_i,iCell); }
   else if (varHandle == Na_i_handle) { return READ_STATE(Na_i,iCell); }
   else if (varHandle == R_prime_handle) { return READ_STATE(R_prime,iCell); }
   else if (varHandle == Xr1_handle) { return READ_STATE(Xr1,iCell); }
   else if (varHandle == Xr2_handle) { return READ_STATE(Xr2,iCell); }
   else if (varHandle == Xs_handle) { return READ_STATE(Xs,iCell); }
   else if (varHandle == d_handle) { return READ_STATE(d,iCell); }
   else if (varHandle == f_handle) { return READ_STATE(f,iCell); }
   else if (varHandle == f2_handle) { return READ_STATE(f2,iCell); }
   else if (varHandle == fCass_handle) { return READ_STATE(fCass,iCell); }
   else if (varHandle == h_handle) { return READ_STATE(h,iCell); }
   else if (varHandle == j_handle) { return READ_STATE(j,iCell); }
   else if (varHandle == m_handle) { return READ_STATE(m,iCell); }
   else if (varHandle == r_handle) { return READ_STATE(r,iCell); }
   else if (varHandle == s_handle) { return READ_STATE(s,iCell); }
   return NAN;
}

double ThisReaction::getValue(int iCell, int varHandle, double V) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   const double Ca_SR=READ_STATE(Ca_SR,iCell);
   const double Ca_i=READ_STATE(Ca_i,iCell);
   const double Ca_ss=READ_STATE(Ca_ss,iCell);
   const double K_i=READ_STATE(K_i,iCell);
   const double Na_i=READ_STATE(Na_i,iCell);
   const double R_prime=READ_STATE(R_prime,iCell);
   const double Xr1=READ_STATE(Xr1,iCell);
   const double Xr2=READ_STATE(Xr2,iCell);
   const double Xs=READ_STATE(Xs,iCell);
   const double d=READ_STATE(d,iCell);
   const double f=READ_STATE(f,iCell);
   const double f2=READ_STATE(f2,iCell);
   const double fCass=READ_STATE(fCass,iCell);
   const double h=READ_STATE(h,iCell);
   const double j=READ_STATE(j,iCell);
   const double m=READ_STATE(m,iCell);
   const double r=READ_STATE(r,iCell);
   const double s=READ_STATE(s,iCell);
   if (0) {}
   else if (varHandle == Ca_SR_handle)
   {
      return Ca_SR;
   }
   else if (varHandle == Ca_i_handle)
   {
      return Ca_i;
   }
   else if (varHandle == Ca_ss_handle)
   {
      return Ca_ss;
   }
   else if (varHandle == K_i_handle)
   {
      return K_i;
   }
   else if (varHandle == Na_i_handle)
   {
      return Na_i;
   }
   else if (varHandle == R_prime_handle)
   {
      return R_prime;
   }
   else if (varHandle == Xr1_handle)
   {
      return Xr1;
   }
   else if (varHandle == Xr2_handle)
   {
      return Xr2;
   }
   else if (varHandle == Xs_handle)
   {
      return Xs;
   }
   else if (varHandle == d_handle)
   {
      return d;
   }
   else if (varHandle == f_handle)
   {
      return f;
   }
   else if (varHandle == f2_handle)
   {
      return f2;
   }
   else if (varHandle == fCass_handle)
   {
      return fCass;
   }
   else if (varHandle == h_handle)
   {
      return h;
   }
   else if (varHandle == i_CaL_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double Ca_o = 2;
      double i_CalTerm1 = (F*F)*(4*V - 60)/(R*T);
      double i_CalTerm2 = exp(F*(2*V - 30)/(R*T));
      double __melodee_temp_005 = V == 15;
      double i_CalTerm3;
      if (__melodee_temp_005)
      {
         i_CalTerm3 = 2*F;
      }
      else
      {
         i_CalTerm3 = i_CalTerm1/(i_CalTerm2 - 1);
      }
      double i_CalTerm4 = i_CalTerm2*i_CalTerm3;
      double i_CaL = d*f*f2*fCass*g_CaL*(-Ca_o*i_CalTerm3 + 0.25*Ca_ss*i_CalTerm4);
      return i_CaL;
   }
   else if (varHandle == i_K1_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double K_o = 5.4000000000000004;
      double _expensive_functions_026 = log(K_o/K_i);
      double E_K = R*T*_expensive_functions_026/F;
      double VEK = -E_K + V;
      double _expensive_functions_035 = exp(0.059999999999999998*VEK);
      double alpha_K1 = 0.10000000000000001/(6.1442123533282098e-6*_expensive_functions_035 + 1);
      double _expensive_functions_036 = exp(-0.5*VEK);
      double _expensive_functions_037 = exp(0.10000000000000001*VEK);
      double _expensive_functions_038 = exp(0.00020000000000000001*VEK);
      double beta_K1 = (0.36787944117144233*_expensive_functions_037 + 3.0606040200802673*_expensive_functions_038)/(_expensive_functions_036 + 1);
      double xK1_inf = alpha_K1/(alpha_K1 + beta_K1);
      double _expensive_functions_039 = sqrt(K_o);
      double i_K1 = 0.43033148291193518*VEK*_expensive_functions_039*g_K1*xK1_inf;
      return i_K1;
   }
   else if (varHandle == i_Kr_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double K_o = 5.4000000000000004;
      double _expensive_functions_026 = log(K_o/K_i);
      double E_K = R*T*_expensive_functions_026/F;
      double VEK = -E_K + V;
      double _expensive_functions_040 = sqrt(K_o);
      double i_Kr = 0.43033148291193518*VEK*Xr1*Xr2*_expensive_functions_040*g_Kr;
      return i_Kr;
   }
   else if (varHandle == i_Ks_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double K_o = 5.4000000000000004;
      double P_kna = 0.029999999999999999;
      double Na_o = 140;
      double _expensive_functions_042 = log((K_o + Na_o*P_kna)/(K_i + Na_i*P_kna));
      double E_Ks = R*T*_expensive_functions_042/F;
      double i_Ks = (Xs*Xs)*g_Ks*(-E_Ks + V);
      return i_Ks;
   }
   else if (varHandle == i_Na_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double Na_o = 140;
      double _expensive_functions_041 = log(Na_o/Na_i);
      double E_Na = R*T*_expensive_functions_041/F;
      double i_Na = (m*m*m)*g_Na*h*j*(-E_Na + V);
      return i_Na;
   }
   else if (varHandle == i_NaCa_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double Ca_o = 2;
      double K_NaCa = 1000;
      double K_sat = 0.10000000000000001;
      double Km_Ca = 1.3799999999999999;
      double Km_Nai = 87.5;
      double alpha = 2.5;
      double gamma = 0.34999999999999998;
      double exp_gamma_VFRT = exp(F*V*gamma/(R*T));
      double exp_gamma_m1_VFRT = exp(F*V*(gamma - 1)/(R*T));
      double Na_o = 140;
      double i_NaCa = K_NaCa*(-(Na_o*Na_o*Na_o)*Ca_i*alpha*exp_gamma_m1_VFRT + (Na_i*Na_i*Na_i)*Ca_o*exp_gamma_VFRT)/(((Km_Nai*Km_Nai*Km_Nai) + (Na_o*Na_o*Na_o))*(Ca_o + Km_Ca)*(K_sat*exp_gamma_m1_VFRT + 1));
      return i_NaCa;
   }
   else if (varHandle == i_NaK_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double K_o = 5.4000000000000004;
      double K_mNa = 40;
      double K_mk = 1;
      double P_NaK = 2.7240000000000002;
      double _expensive_functions_031 = exp(-F*V/(R*T));
      double _expensive_functions_032 = exp(-0.10000000000000001*F*V/(R*T));
      double i_NaK_term = K_o*P_NaK/((K_mk + K_o)*(0.035299999999999998*_expensive_functions_031 + 0.1245*_expensive_functions_032 + 1));
      double i_NaK = Na_i*i_NaK_term/(K_mNa + Na_i);
      return i_NaK;
   }
   else if (varHandle == i_b_Ca_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double Ca_o = 2;
      double _expensive_functions_012 = log(Ca_o/Ca_i);
      double E_Ca = 0.5*R*T*_expensive_functions_012/F;
      double i_b_Ca = g_bca*(-E_Ca + V);
      return i_b_Ca;
   }
   else if (varHandle == i_b_Na_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double Na_o = 140;
      double _expensive_functions_041 = log(Na_o/Na_i);
      double E_Na = R*T*_expensive_functions_041/F;
      double i_b_Na = g_bna*(-E_Na + V);
      return i_b_Na;
   }
   else if (varHandle == i_leak_handle)
   {
      double V_leak = 0.00036000000000000002;
      double i_leak = V_leak*(Ca_SR - Ca_i);
      return i_leak;
   }
   else if (varHandle == i_p_Ca_handle)
   {
      double K_pCa = 0.00050000000000000001;
      double i_p_Ca = Ca_i*g_pCa/(Ca_i + K_pCa);
      return i_p_Ca;
   }
   else if (varHandle == i_p_K_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double K_o = 5.4000000000000004;
      double _expensive_functions_019 = exp(-0.16722408026755853*V);
      double i_p_K_term = (1.0/(65.405215741938321*_expensive_functions_019 + 1));
      double _expensive_functions_026 = log(K_o/K_i);
      double E_K = R*T*_expensive_functions_026/F;
      double VEK = -E_K + V;
      double i_p_K = VEK*g_pK*i_p_K_term;
      return i_p_K;
   }
   else if (varHandle == i_rel_handle)
   {
      double EC = 1.5;
      double V_rel = 0.10199999999999999;
      double k1_prime = 0.14999999999999999;
      double k3 = 0.059999999999999998;
      double max_sr = 2.5;
      double min_sr = 1;
      double kcasr = max_sr - (max_sr - min_sr)/(1 + (EC*EC)/(Ca_SR*Ca_SR));
      double k1 = k1_prime/kcasr;
      double O = (Ca_ss*Ca_ss)*R_prime*k1/((Ca_ss*Ca_ss)*k1 + k3);
      double i_rel = O*V_rel*(Ca_SR - Ca_ss);
      return i_rel;
   }
   else if (varHandle == i_to_handle)
   {
      double F = 96485.341499999995;
      double R = 8314.4719999999998;
      double T = 310;
      double K_o = 5.4000000000000004;
      double _expensive_functions_026 = log(K_o/K_i);
      double E_K = R*T*_expensive_functions_026/F;
      double VEK = -E_K + V;
      double i_to = VEK*g_to*r*s;
      return i_to;
   }
   else if (varHandle == i_up_handle)
   {
      double K_up = 0.00025000000000000001;
      double Vmax_up = 0.0063749999999999996;
      double i_up = Vmax_up/(1 + (K_up*K_up)/(Ca_i*Ca_i));
      return i_up;
   }
   else if (varHandle == i_xfer_handle)
   {
      double V_xfer = 0.0038;
      double i_xfer = V_xfer*(-Ca_i + Ca_ss);
      return i_xfer;
   }
   else if (varHandle == j_handle)
   {
      return j;
   }
   else if (varHandle == m_handle)
   {
      return m;
   }
   else if (varHandle == r_handle)
   {
      return r;
   }
   else if (varHandle == s_handle)
   {
      return s;
   }
   return NAN;
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.clear();
   fieldUnits.clear();
   fieldNames.push_back("Ca_SR");
   fieldUnits.push_back(getUnit("Ca_SR"));
   fieldNames.push_back("Ca_i");
   fieldUnits.push_back(getUnit("Ca_i"));
   fieldNames.push_back("Ca_ss");
   fieldUnits.push_back(getUnit("Ca_ss"));
   fieldNames.push_back("K_i");
   fieldUnits.push_back(getUnit("K_i"));
   fieldNames.push_back("Na_i");
   fieldUnits.push_back(getUnit("Na_i"));
   fieldNames.push_back("R_prime");
   fieldUnits.push_back(getUnit("R_prime"));
   fieldNames.push_back("Xr1");
   fieldUnits.push_back(getUnit("Xr1"));
   fieldNames.push_back("Xr2");
   fieldUnits.push_back(getUnit("Xr2"));
   fieldNames.push_back("Xs");
   fieldUnits.push_back(getUnit("Xs"));
   fieldNames.push_back("d");
   fieldUnits.push_back(getUnit("d"));
   fieldNames.push_back("f");
   fieldUnits.push_back(getUnit("f"));
   fieldNames.push_back("f2");
   fieldUnits.push_back(getUnit("f2"));
   fieldNames.push_back("fCass");
   fieldUnits.push_back(getUnit("fCass"));
   fieldNames.push_back("h");
   fieldUnits.push_back(getUnit("h"));
   fieldNames.push_back("j");
   fieldUnits.push_back(getUnit("j"));
   fieldNames.push_back("m");
   fieldUnits.push_back(getUnit("m"));
   fieldNames.push_back("r");
   fieldUnits.push_back(getUnit("r"));
   fieldNames.push_back("s");
   fieldUnits.push_back(getUnit("s"));
}

}
