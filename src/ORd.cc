/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "ORd.hh"
#include "object_cc.hh"
#include "mpiUtils.h"
#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>

using namespace std;

namespace scanReaction 
{


   
   Reaction* scanORd(OBJECT* obj, const int numPoints, const double _dt)
   {
      ORd::ThisReaction* reaction = new ORd::ThisReaction(numPoints, _dt);

      //override the defaults
#define setDefault(name, value) objectGet(obj, #name, reaction->name, #value)
      //EDIT_PARAMETERS
      setDefault(celltype, 0);
      setDefault(JrelStiffConst, 0.00500000000000000);
#undef setDefault

      bool reusingInterpolants = false;
      string fitName;
      objectGet(obj, "fit", fitName, "");
      int funcCount = sizeof(reaction->_interpolant)/sizeof(reaction->_interpolant[0]);
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
               outfile << obj->name << "_interpFunc" << _ii << " ";
            }
            outfile << ";\n";
            outfile << "}\n";

            for (int _ii=0; _ii<funcCount; _ii++)
            {
               outfile << obj->name << "_interpFunc" << _ii << " FUNCTION { "
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

      reaction->constructKernel();
      return reaction;
   }

}

namespace ORd 
{
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
    for (int _ii=interp.numDenom_+interp.numNumer_-2; _ii>=0; _ii--)
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
   "ORd_program\n"
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
   "   xrf_off,\n"
   "   xrs_off,\n"
   "   xs1_off,\n"
   "   xs2_off,\n"
   "   NUMSTATES\n"
   "};\n"
   "__global__ void ORd_kernel(const double* _Vm, const double* _iStim, double* _dVm, double* _state) {\n"
   "const double _dt = " << __cachedDt << ";\n"
   "const int _nCells = " << nCells_ << ";\n"

   "const double JrelStiffConst = " << JrelStiffConst << ";\n"
   "const double celltype = " << celltype << ";\n"
   "const int _ii = threadIdx.x + blockIdx.x*blockDim.x;\n"
   "if (_ii >= _nCells) { return; }\n"
   "const double V = _Vm[_ii];\n"
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
   "const double xrf = _state[_ii+xrf_off*_nCells];\n"
   "const double xrs = _state[_ii+xrs_off*_nCells];\n"
   "const double xs1 = _state[_ii+xs1_off*_nCells];\n"
   "const double xs2 = _state[_ii+xs2_off*_nCells];\n"
   "double Ahf = 0.990000000000000;\n"
   "double Ahs = -Ahf + 1.0;\n"
   "double h = Ahf*hf + Ahs*hs;\n"
   "double hp = Ahf*hf + Ahs*hsp;\n"
   "double GNa = 75;\n"
   "double v = V;\n"
   "double nao = 140.000000000000;\n"
   "double cao = 1.80000000000000;\n"
   "double ko = 5.40000000000000;\n"
   "double R = 8314.00000000000;\n"
   "double T = 310.000000000000;\n"
   "double F = 96485.0000000000;\n"
   "double L = 0.0100000000000000;\n"
   "double rad = 0.00110000000000000;\n"
   "double vcell = 3140.0*rad*rad*L;\n"
   "double Ageo = 6.28*L*rad + 6.28*rad*rad;\n"
   "double Acap = 2*Ageo;\n"
   "double vmyo = 0.68*vcell;\n"
   "double vnsr = 0.0552*vcell;\n"
   "double vjsr = 0.0048*vcell;\n"
   "double vss = 0.02*vcell;\n"
   "double _expensive_functions = log(nao/nai);\n"
   "double ENa = R*T*_expensive_functions/F;\n"
   "double _expensive_functions_001 = log(ko/ki);\n"
   "double EK = R*T*_expensive_functions_001/F;\n"
   "double PKNa = 0.0183300000000000;\n"
   "double _expensive_functions_002 = log((PKNa*nao + ko)/(PKNa*nai + ki));\n"
   "double EKs = R*T*_expensive_functions_002/F;\n"
   "double vffrt = F*F*v/(R*T);\n"
   "double vfrt = F*v/(R*T);\n"
   "double KmCaMK = 0.150000000000000;\n"
   "double aCaMK = 0.0500000000000000;\n"
   "double bCaMK = 0.000680000000000000;\n"
   "double CaMKo = 0.0500000000000000;\n"
   "double KmCaM = 0.00150000000000000;\n"
   "double CaMKb = CaMKo*(-CaMKt + 1.0)/(KmCaM/cass + 1.0);\n"
   "double CaMKa = CaMKb + CaMKt;\n"
   "double CaMKt_diff = CaMKb*aCaMK*(CaMKb + CaMKt) - CaMKt*bCaMK;\n"
   "double _expensive_functions_003 = exp(-0.189969604863222*v);\n"
   "double mLss = 1.0/(0.000291579585635531*_expensive_functions_003 + 1.0);\n"
   "double _expensive_functions_004 = exp(0.133547008547009*v);\n"
   "double hLss = 1.0/(120578.155955224*_expensive_functions_004 + 1.0);\n"
   "double thL = 200.000000000000;\n"
   "double _expensive_functions_005 = exp(0.133547008547009*v);\n"
   "double hLpss = 1.0/(275969.290386987*_expensive_functions_005 + 1.0);\n"
   "double thLp = 3.0*thL;\n"
   "double GNaL = 0.00750000000000000;\n"
   "double __melodee_temp_000 = celltype == 1;\n"
   "double GNaL_001;\n"
   "if (__melodee_temp_000)\n"
   "{\n"
   "   GNaL_001 = 0.6*GNaL;\n"
   "}\n"
   "else\n"
   "{\n"
   "   GNaL_001 = GNaL;\n"
   "}\n"
   "double fINaLp = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double INaL = GNaL_001*mL*(-ENa + v)*(fINaLp*hLp + hL*(-fINaLp + 1.0));\n"
   "double _expensive_functions_006 = exp(-0.0674763832658569*v);\n"
   "double ass = 1.0/(2.63165081616736*_expensive_functions_006 + 1.0);\n"
   "double _expensive_functions_007 = exp(-0.0340351378763435*v);\n"
   "double _expensive_functions_008 = exp(0.0340351378763435*v);\n"
   "double ta = 1.0515/(3.5/(30.0695727273975*_expensive_functions_008 + 1.0) + 1.0/(2.26210170705788*_expensive_functions_007 + 1.2089));\n"
   "double _expensive_functions_009 = exp(0.175100682892663*v);\n"
   "double iss = 1.0/(2194.9707645383*_expensive_functions_009 + 1.0);\n"
   "double __melodee_temp_001 = celltype == 1;\n"
   "double delta_epi;\n"
   "if (__melodee_temp_001)\n"
   "{\n"
   "   double _expensive_functions_010 = exp(0.2*v);\n"
   "   delta_epi = 1.0 - 0.95/(1202604.28416478*_expensive_functions_010 + 1.0);\n"
   "}\n"
   "else\n"
   "{\n"
   "   delta_epi = 1.00000000000000;\n"
   "}\n"
   "double _expensive_functions_010 = exp(-0.01*v);\n"
   "double _expensive_functions_011 = exp(0.0602772754671489*v);\n"
   "double tiF = 4.562 + 1.0/(0.144686984212728*_expensive_functions_010 + 1.63008963497809*_expensive_functions_011);\n"
   "double _expensive_functions_012 = exp(-0.0169348010160881*v);\n"
   "double _expensive_functions_013 = exp(0.123777695259314*v);\n"
   "double tiS = 23.62 + 1.0/(0.000276177639533774*_expensive_functions_012 + 0.0242089628046045*_expensive_functions_013);\n"
   "double tiF_001 = delta_epi*tiF;\n"
   "double tiS_001 = delta_epi*tiS;\n"
   "double _expensive_functions_014 = exp(0.00661375661375661*v);\n"
   "double AiF = 1.0/(0.243485371875229*_expensive_functions_014 + 1.0);\n"
   "double AiS = -AiF + 1.0;\n"
   "double i = AiF*iF + AiS*iS;\n"
   "double _expensive_functions_015 = exp(-0.0674763832658569*v);\n"
   "double apss = 1.0/(5.16742846223067*_expensive_functions_015 + 1.0);\n"
   "double _expensive_functions_016 = exp(0.0629326620516048*v);\n"
   "double _expensive_functions_017 = exp(-4.64252553389044*v);\n"
   "double dti_develop = 1.354 + 0.0001/(2.65912690452306e-5*_expensive_functions_016 + 4.55417797371283e+24*_expensive_functions_017);\n"
   "double _expensive_functions_018 = exp(0.05*v);\n"
   "double dti_recover = 1.0 - 0.5/(33.1154519586923*_expensive_functions_018 + 1.0);\n"
   "double tiFp = dti_develop*dti_recover*tiF_001;\n"
   "double tiSp = dti_develop*dti_recover*tiS_001;\n"
   "double ip = AiF*iFp + AiS*iSp;\n"
   "double Gto = 0.0200000000000000;\n"
   "double __melodee_temp_003 = celltype == 1;\n"
   "double Gto_001;\n"
   "if (__melodee_temp_003)\n"
   "{\n"
   "   Gto_001 = 4.0*Gto;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_002 = celltype == 2;\n"
   "   if (__melodee_temp_002)\n"
   "   {\n"
   "      Gto_001 = 4.0*Gto;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      Gto_001 = Gto;\n"
   "   }\n"
   "}\n"
   "double fItop = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double Ito = Gto_001*(-EK + v)*(a*i*(-fItop + 1.0) + ap*fItop*ip);\n"
   "double _expensive_functions_019 = exp(-0.236406619385343*v);\n"
   "double dss = 1.0/(0.393985142266695*_expensive_functions_019 + 1.0);\n"
   "double _expensive_functions_020 = exp(0.09*v);\n"
   "double _expensive_functions_021 = exp(-0.05*v);\n"
   "double td = 0.6 + 1.0/(3.52542148736538*_expensive_functions_020 + 0.740818220681718*_expensive_functions_021);\n"
   "double _expensive_functions_022 = exp(0.270562770562771*v);\n"
   "double fss = 1.0/(199.860384967786*_expensive_functions_022 + 1.0);\n"
   "double _expensive_functions_023 = exp(0.1*v);\n"
   "double _expensive_functions_024 = exp(-0.1*v);\n"
   "double tff = 7.0 + 1.0/(0.0332507524451879*_expensive_functions_023 + 0.000609008774564757*_expensive_functions_024);\n"
   "double _expensive_functions_025 = exp(-0.25*v);\n"
   "double _expensive_functions_026 = exp(0.166666666666667*v);\n"
   "double tfs = 1000.0 + 1.0/(1.00276678901067e-5*_expensive_functions_025 + 8.05341561812489e-5*_expensive_functions_026);\n"
   "double Aff = 0.600000000000000;\n"
   "double Afs = -Aff + 1.0;\n"
   "double f = Aff*ff + Afs*fs;\n"
   "double fcass = fss;\n"
   "double _expensive_functions_027 = exp(-0.142857142857143*v);\n"
   "double _expensive_functions_028 = exp(0.142857142857143*v);\n"
   "double tfcaf = 7.0 + 1.0/(0.0708317980974062*_expensive_functions_027 + 0.0225887248803104*_expensive_functions_028);\n"
   "double _expensive_functions_029 = exp(0.142857142857143*v);\n"
   "double _expensive_functions_030 = exp(-0.333333333333333*v);\n"
   "double tfcas = 100.0 + 1.0/(0.00012*_expensive_functions_029 + 0.00012*_expensive_functions_030);\n"
   "double _expensive_functions_031 = exp(0.1*v);\n"
   "double Afcaf = 0.3 + 0.6/(0.367879441171442*_expensive_functions_031 + 1.0);\n"
   "double Afcas = -Afcaf + 1.0;\n"
   "double fca = Afcaf*fcaf + Afcas*fcas;\n"
   "double tjca = 75.0000000000000;\n"
   "double tffp = 2.5*tff;\n"
   "double fp = Aff*ffp + Afs*fs;\n"
   "double tfcafp = 2.5*tfcaf;\n"
   "double fcap = Afcaf*fcafp + Afcas*fcas;\n"
   "double Kmn = 0.00200000000000000;\n"
   "double k2n = 1000.00000000000;\n"
   "double km2n = 1.0*jca;\n"
   "double anca = 1.0/(k2n/km2n + (Kmn/cass + 1.0)*(Kmn/cass + 1.0)*(Kmn/cass + 1.0)*(Kmn/cass + 1.0));\n"
   "double _expensive_functions_032 = exp(2.0*vfrt);\n"
   "double _expensive_functions_033 = exp(2.0*vfrt);\n"
   "double PhiCaL = 4.0*vffrt*(_expensive_functions_033*cass - 0.341*cao)/(_expensive_functions_032 - 1.0);\n"
   "double _expensive_functions_034 = exp(1.0*vfrt);\n"
   "double _expensive_functions_035 = exp(1.0*vfrt);\n"
   "double PhiCaNa = 1.0*vffrt*(0.75*_expensive_functions_035*nass - 0.75*nao)/(_expensive_functions_034 - 1.0);\n"
   "double _expensive_functions_036 = exp(1.0*vfrt);\n"
   "double _expensive_functions_037 = exp(1.0*vfrt);\n"
   "double PhiCaK = 1.0*vffrt*(0.75*_expensive_functions_037*kss - 0.75*ko)/(_expensive_functions_036 - 1.0);\n"
   "double zca = 2.00000000000000;\n"
   "double PCa = 0.000100000000000000;\n"
   "double __melodee_temp_005 = celltype == 1;\n"
   "double PCa_001;\n"
   "if (__melodee_temp_005)\n"
   "{\n"
   "   PCa_001 = 1.2*PCa;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_004 = celltype == 2;\n"
   "   if (__melodee_temp_004)\n"
   "   {\n"
   "      PCa_001 = 2.5*PCa;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      PCa_001 = PCa;\n"
   "   }\n"
   "}\n"
   "double PCap = 1.1*PCa_001;\n"
   "double PCaNa = 0.00125*PCa_001;\n"
   "double PCaK = 0.0003574*PCa_001;\n"
   "double PCaNap = 0.00125*PCap;\n"
   "double PCaKp = 0.0003574*PCap;\n"
   "double fICaLp = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double ICaL = PCa_001*PhiCaL*d*(-fICaLp + 1.0)*(f*(-nca + 1.0) + fca*jca*nca) + PCap*PhiCaL*d*fICaLp*(fcap*jca*nca + fp*(-nca + 1.0));\n"
   "double ICaNa = PCaNa*PhiCaNa*d*(-fICaLp + 1.0)*(f*(-nca + 1.0) + fca*jca*nca) + PCaNap*PhiCaNa*d*fICaLp*(fcap*jca*nca + fp*(-nca + 1.0));\n"
   "double ICaK = PCaK*PhiCaK*d*(-fICaLp + 1.0)*(f*(-nca + 1.0) + fca*jca*nca) + PCaKp*PhiCaK*d*fICaLp*(fcap*jca*nca + fp*(-nca + 1.0));\n"
   "double _expensive_functions_038 = exp(-0.147297098247165*v);\n"
   "double xrss = 1.0/(0.292873088723775*_expensive_functions_038 + 1.0);\n"
   "double _expensive_functions_039 = exp(0.258464719565779*v);\n"
   "double _expensive_functions_040 = exp(-0.0490677134445535*v);\n"
   "double txrf = 12.98 + 1.0/(0.000102023931289489*_expensive_functions_039 + 0.000429929608919291*_expensive_functions_040);\n"
   "double _expensive_functions_041 = exp(0.135961930659415*v);\n"
   "double _expensive_functions_042 = exp(-0.038550501156515*v);\n"
   "double txrs = 1.865 + 1.0/(0.000592242003680939*_expensive_functions_041 + 3.54996611180246e-5*_expensive_functions_042);\n"
   "double _expensive_functions_043 = exp(0.0261711593823606*v);\n"
   "double Axrf = 1.0/(4.19729909473472*_expensive_functions_043 + 1.0);\n"
   "double Axrs = -Axrf + 1.0;\n"
   "double xr = Axrf*xrf + Axrs*xrs;\n"
   "double _expensive_functions_044 = exp(0.0133333333333333*v);\n"
   "double _expensive_functions_045 = exp(0.0333333333333333*v);\n"
   "double rkr = 1.0/((2.08200908407846*_expensive_functions_044 + 1.0)*(0.716531310573789*_expensive_functions_045 + 1.0));\n"
   "double GKr = 0.0460000000000000;\n"
   "double __melodee_temp_007 = celltype == 1;\n"
   "double GKr_001;\n"
   "if (__melodee_temp_007)\n"
   "{\n"
   "   GKr_001 = 1.3*GKr;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_006 = celltype == 2;\n"
   "   if (__melodee_temp_006)\n"
   "   {\n"
   "      GKr_001 = 0.8*GKr;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      GKr_001 = GKr;\n"
   "   }\n"
   "}\n"
   "double _expensive_functions_046 = sqrt(ko);\n"
   "double IKr = 0.430331482911935*GKr_001*_expensive_functions_046*rkr*xr*(-EK + v);\n"
   "double _expensive_functions_047 = exp(-0.111957008508733*v);\n"
   "double xs1ss = 1.0/(0.272885960356565*_expensive_functions_047 + 1.0);\n"
   "double _expensive_functions_048 = exp(0.0561797752808989*v);\n"
   "double _expensive_functions_049 = exp(-0.00434782608695652*v);\n"
   "double txs1 = 817.3 + 1.0/(0.00350406776307486*_expensive_functions_048 + 0.000518480908358166*_expensive_functions_049);\n"
   "double xs2ss = xs1ss;\n"
   "double _expensive_functions_050 = exp(-0.032258064516129*v);\n"
   "double _expensive_functions_051 = exp(0.05*v);\n"
   "double txs2 = 1.0/(0.00225613570106391*_expensive_functions_050 + 0.000820849986238988*_expensive_functions_051);\n"
   "double _expensive_functions_052 = pow(1.0/cai, 1.4);\n"
   "double KsCa = 1.0 + 0.6/(6.48182102606265e-7*_expensive_functions_052 + 1.0);\n"
   "double GKs = 0.00340000000000000;\n"
   "double __melodee_temp_008 = celltype == 1;\n"
   "double GKs_001;\n"
   "if (__melodee_temp_008)\n"
   "{\n"
   "   GKs_001 = 1.4*GKs;\n"
   "}\n"
   "else\n"
   "{\n"
   "   GKs_001 = GKs;\n"
   "}\n"
   "double IKs = GKs_001*KsCa*xs1*xs2*(-EKs + v);\n"
   "double _expensive_functions_053 = exp((-2.5538*ko - v - 144.59)/(1.5692*ko + 3.8115));\n"
   "double xk1ss = 1.0/(_expensive_functions_053 + 1.0);\n"
   "double _expensive_functions_054 = exp(-0.0491159135559921*v);\n"
   "double _expensive_functions_055 = exp(0.0144237703735757*v);\n"
   "double txk1 = 122.2/(0.00193520076313902*_expensive_functions_054 + 30.433647575249*_expensive_functions_055);\n"
   "double _expensive_functions_056 = exp(-0.273886021278837*ko + 0.105340777414937*v);\n"
   "double rk1 = 1.0/(69220.6322106767*_expensive_functions_056 + 1.0);\n"
   "double GK1 = 0.190800000000000;\n"
   "double __melodee_temp_010 = celltype == 1;\n"
   "double GK1_001;\n"
   "if (__melodee_temp_010)\n"
   "{\n"
   "   GK1_001 = 1.2*GK1;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_009 = celltype == 2;\n"
   "   if (__melodee_temp_009)\n"
   "   {\n"
   "      GK1_001 = 1.3*GK1;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      GK1_001 = GK1;\n"
   "   }\n"
   "}\n"
   "double _expensive_functions_057 = sqrt(ko);\n"
   "double IK1 = GK1_001*_expensive_functions_057*rk1*xk1*(-EK + v);\n"
   "double kna1 = 15.0000000000000;\n"
   "double kna2 = 5.00000000000000;\n"
   "double kna3 = 88.1200000000000;\n"
   "double kasymm = 12.5000000000000;\n"
   "double wna = 60000.0000000000;\n"
   "double wca = 60000.0000000000;\n"
   "double wnaca = 5000.00000000000;\n"
   "double kcaon = 1500000.00000000;\n"
   "double kcaoff = 5000.00000000000;\n"
   "double qna = 0.522400000000000;\n"
   "double qca = 0.167000000000000;\n"
   "double hca = exp(F*qca*v/(R*T));\n"
   "double hna = exp(F*qna*v/(R*T));\n"
   "double h1 = 1 + nai*(hna + 1)/kna3;\n"
   "double h2 = hna*nai/(h1*kna3);\n"
   "double h3 = 1.0/h1;\n"
   "double h4 = 1.0 + nai*(1 + nai/kna2)/kna1;\n"
   "double h5 = nai*nai/(h4*kna1*kna2);\n"
   "double h6 = 1.0/h4;\n"
   "double h7 = 1.0 + nao*(1.0 + 1.0/hna)/kna3;\n"
   "double h8 = nao/(h7*hna*kna3);\n"
   "double h9 = 1.0/h7;\n"
   "double h10 = kasymm + 1.0 + nao*(1.0 + nao/kna2)/kna1;\n"
   "double h11 = nao*nao/(h10*kna1*kna2);\n"
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
   "double KmCaAct = 0.000150000000000000;\n"
   "double allo = 1.0/(KmCaAct*KmCaAct/cai*cai + 1.0);\n"
   "double zna = 1.00000000000000;\n"
   "double JncxNa = -3.0*E1*k8 - E2*k3pp + E3*k4pp + 3.0*E4*k7;\n"
   "double JncxCa = -E1*k1 + E2*k2;\n"
   "double Gncx = 0.000800000000000000;\n"
   "double __melodee_temp_012 = celltype == 1;\n"
   "double Gncx_001;\n"
   "if (__melodee_temp_012)\n"
   "{\n"
   "   Gncx_001 = 1.1*Gncx;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_011 = celltype == 2;\n"
   "   if (__melodee_temp_011)\n"
   "   {\n"
   "      Gncx_001 = 1.4*Gncx;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      Gncx_001 = Gncx;\n"
   "   }\n"
   "}\n"
   "double INaCa_i = 0.8*Gncx_001*allo*(JncxCa*zca + JncxNa*zna);\n"
   "double h1_001 = 1 + nass*(hna + 1)/kna3;\n"
   "double h2_001 = hna*nass/(h1_001*kna3);\n"
   "double h3_001 = 1.0/h1_001;\n"
   "double h4_001 = 1.0 + nass*(1 + nass/kna2)/kna1;\n"
   "double h5_001 = nass*nass/(h4_001*kna1*kna2);\n"
   "double h6_001 = 1.0/h4_001;\n"
   "double h7_001 = 1.0 + nao*(1.0 + 1.0/hna)/kna3;\n"
   "double h8_001 = nao/(h7_001*hna*kna3);\n"
   "double h9_001 = 1.0/h7_001;\n"
   "double h10_001 = kasymm + 1.0 + nao*(1 + nao/kna2)/kna1;\n"
   "double h11_001 = nao*nao/(h10_001*kna1*kna2);\n"
   "double h12_001 = 1.0/h10_001;\n"
   "double k1_001 = cao*h12_001*kcaon;\n"
   "double k2_001 = kcaoff;\n"
   "double k3p_001 = h9_001*wca;\n"
   "double k3pp_001 = h8_001*wnaca;\n"
   "double k3_001 = k3p_001 + k3pp_001;\n"
   "double k4p_001 = h3_001*wca/hca;\n"
   "double k4pp_001 = h2_001*wnaca;\n"
   "double k4_001 = k4p_001 + k4pp_001;\n"
   "double k5_001 = kcaoff;\n"
   "double k6_001 = cass*h6_001*kcaon;\n"
   "double k7_001 = h2_001*h5_001*wna;\n"
   "double k8_001 = h11_001*h8_001*wna;\n"
   "double x1_001 = k2_001*k4_001*(k6_001 + k7_001) + k5_001*k7_001*(k2_001 + k3_001);\n"
   "double x2_001 = k1_001*k7_001*(k4_001 + k5_001) + k4_001*k6_001*(k1_001 + k8_001);\n"
   "double x3_001 = k1_001*k3_001*(k6_001 + k7_001) + k6_001*k8_001*(k2_001 + k3_001);\n"
   "double x4_001 = k2_001*k8_001*(k4_001 + k5_001) + k3_001*k5_001*(k1_001 + k8_001);\n"
   "double E1_001 = x1_001/(x1_001 + x2_001 + x3_001 + x4_001);\n"
   "double E2_001 = x2_001/(x1_001 + x2_001 + x3_001 + x4_001);\n"
   "double E3_001 = x3_001/(x1_001 + x2_001 + x3_001 + x4_001);\n"
   "double E4_001 = x4_001/(x1_001 + x2_001 + x3_001 + x4_001);\n"
   "double KmCaAct_001 = 0.000150000000000000;\n"
   "double allo_001 = 1.0/(KmCaAct_001*KmCaAct_001/cass*cass + 1.0);\n"
   "double JncxNa_001 = -3.0*E1_001*k8_001 - E2_001*k3pp_001 + E3_001*k4pp_001 + 3.0*E4_001*k7_001;\n"
   "double JncxCa_001 = -E1_001*k1_001 + E2_001*k2_001;\n"
   "double INaCa_ss = 0.2*Gncx_001*allo_001*(JncxCa_001*zca + JncxNa_001*zna);\n"
   "double k1p = 949.500000000000;\n"
   "double k1m = 182.400000000000;\n"
   "double k2p = 687.200000000000;\n"
   "double k2m = 39.4000000000000;\n"
   "double k3p_002 = 1899.00000000000;\n"
   "double k3m = 79300.0000000000;\n"
   "double k4p_002 = 639.000000000000;\n"
   "double k4m = 40.0000000000000;\n"
   "double Knai0 = 9.07300000000000;\n"
   "double Knao0 = 27.7800000000000;\n"
   "double delta = -0.155000000000000;\n"
   "double _expensive_functions_058 = exp(0.333333333333333*F*delta*v/(R*T));\n"
   "double Knai = Knai0*_expensive_functions_058;\n"
   "double _expensive_functions_059 = exp(0.333333333333333*F*v*(-delta + 1.0)/(R*T));\n"
   "double Knao = Knao0*_expensive_functions_059;\n"
   "double Kki = 0.500000000000000;\n"
   "double Kko = 0.358200000000000;\n"
   "double MgADP = 0.0500000000000000;\n"
   "double MgATP = 9.80000000000000;\n"
   "double Kmgatp = 1.69800000000000e-7;\n"
   "double H = 1.00000000000000e-7;\n"
   "double eP = 4.20000000000000;\n"
   "double Khp = 1.69800000000000e-7;\n"
   "double Knap = 224.000000000000;\n"
   "double Kxkur = 292.000000000000;\n"
   "double P = eP/(H/Khp + 1.0 + ki/Kxkur + nai/Knap);\n"
   "double a1 = k1p*(nai/Knai)*(nai/Knai)*(nai/Knai)/((1.0 + ki/Kki)*(1.0 + ki/Kki) + (1.0 + nai/Knai)*(1.0 + nai/Knai)*(1.0 + nai/Knai) - 1.0);\n"
   "double b1 = MgADP*k1m;\n"
   "double a2 = k2p;\n"
   "double b2 = k2m*(nao/Knao)*(nao/Knao)*(nao/Knao)/((1.0 + ko/Kko)*(1.0 + ko/Kko) + (1.0 + nao/Knao)*(1.0 + nao/Knao)*(1.0 + nao/Knao) - 1.0);\n"
   "double a3 = k3p_002*(ko/Kko)*(ko/Kko)/((1.0 + ko/Kko)*(1.0 + ko/Kko) + (1.0 + nao/Knao)*(1.0 + nao/Knao)*(1.0 + nao/Knao) - 1.0);\n"
   "double b3 = H*P*k3m/(1.0 + MgATP/Kmgatp);\n"
   "double a4 = MgATP*k4p_002/(Kmgatp*(1.0 + MgATP/Kmgatp));\n"
   "double b4 = k4m*(ki/Kki)*(ki/Kki)/((1.0 + ki/Kki)*(1.0 + ki/Kki) + (1.0 + nai/Knai)*(1.0 + nai/Knai)*(1.0 + nai/Knai) - 1.0);\n"
   "double x1_002 = a1*a2*a4 + a1*a2*b3 + a2*b3*b4 + b2*b3*b4;\n"
   "double x2_002 = a1*a2*a3 + a2*a3*b4 + a3*b1*b4 + b1*b2*b4;\n"
   "double x3_002 = a2*a3*a4 + a3*a4*b1 + a4*b1*b2 + b1*b2*b3;\n"
   "double x4_002 = a1*a3*a4 + a1*a4*b2 + a1*b2*b3 + b2*b3*b4;\n"
   "double E1_002 = x1_002/(x1_002 + x2_002 + x3_002 + x4_002);\n"
   "double E2_002 = x2_002/(x1_002 + x2_002 + x3_002 + x4_002);\n"
   "double E3_002 = x3_002/(x1_002 + x2_002 + x3_002 + x4_002);\n"
   "double E4_002 = x4_002/(x1_002 + x2_002 + x3_002 + x4_002);\n"
   "double zk = 1.00000000000000;\n"
   "double JnakNa = 3.0*E1_002*a3 - 3.0*E2_002*b3;\n"
   "double JnakK = -2.0*E3_002*a1 + 2.0*E4_002*b1;\n"
   "double Pnak = 30;\n"
   "double __melodee_temp_014 = celltype == 1;\n"
   "double Pnak_001;\n"
   "if (__melodee_temp_014)\n"
   "{\n"
   "   Pnak_001 = 0.9*Pnak;\n"
   "}\n"
   "else\n"
   "{\n"
   "   double __melodee_temp_013 = celltype == 2;\n"
   "   if (__melodee_temp_013)\n"
   "   {\n"
   "      Pnak_001 = 0.7*Pnak;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      Pnak_001 = Pnak;\n"
   "   }\n"
   "}\n"
   "double INaK = Pnak_001*(JnakK*zk + JnakNa*zna);\n"
   "double _expensive_functions_060 = exp(-0.054525627044711*v);\n"
   "double xkb = 1.0/(2.20236345094924*_expensive_functions_060 + 1.0);\n"
   "double GKb = 0.00300000000000000;\n"
   "double __melodee_temp_015 = celltype == 1;\n"
   "double GKb_001;\n"
   "if (__melodee_temp_015)\n"
   "{\n"
   "   GKb_001 = 0.6*GKb;\n"
   "}\n"
   "else\n"
   "{\n"
   "   GKb_001 = GKb;\n"
   "}\n"
   "double IKb = GKb_001*xkb*(-EK + v);\n"
   "double PNab = 3.75000000000000e-10;\n"
   "double _expensive_functions_061 = exp(vfrt);\n"
   "double _expensive_functions_062 = exp(vfrt);\n"
   "double INab = PNab*vffrt*(_expensive_functions_062*nai - nao)/(_expensive_functions_061 - 1.0);\n"
   "double PCab = 2.50000000000000e-8;\n"
   "double _expensive_functions_063 = exp(2.0*vfrt);\n"
   "double _expensive_functions_064 = exp(2.0*vfrt);\n"
   "double ICab = 4.0*PCab*vffrt*(_expensive_functions_064*cai - 0.341*cao)/(_expensive_functions_063 - 1.0);\n"
   "double GpCa = 0.000500000000000000;\n"
   "double IpCa = GpCa*cai/(cai + 0.0005);\n"
   "double JdiffNa = -0.5*nai + 0.5*nass;\n"
   "double JdiffK = -0.5*ki + 0.5*kss;\n"
   "double Jdiff = -5.0*cai + 5.0*cass;\n"
   "double bt = 4.75000000000000;\n"
   "double a_rel = 0.5*bt;\n"
   "double Jrel_inf = -ICaL*a_rel/(25.62890625*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr) + 1.0);\n"
   "double __melodee_temp_016 = celltype == 2;\n"
   "double Jrel_inf_001;\n"
   "if (__melodee_temp_016)\n"
   "{\n"
   "   Jrel_inf_001 = 1.7*Jrel_inf;\n"
   "}\n"
   "else\n"
   "{\n"
   "   Jrel_inf_001 = Jrel_inf;\n"
   "}\n"
   "double tau_rel = bt/(1.0 + 0.0123/cajsr);\n"
   "double __melodee_temp_017 = tau_rel < 0.001;\n"
   "double tau_rel_001;\n"
   "if (__melodee_temp_017)\n"
   "{\n"
   "   tau_rel_001 = 0.00100000000000000;\n"
   "}\n"
   "else\n"
   "{\n"
   "   tau_rel_001 = tau_rel;\n"
   "}\n"
   "double Jrelnp_diff = (Jrel_inf_001 - Jrelnp)/tau_rel_001;\n"
   "double btp = 1.25*bt;\n"
   "double a_relp = 0.5*btp;\n"
   "double Jrelp_inf = -ICaL*a_relp/(25.62890625*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr)*(1.0/cajsr) + 1.0);\n"
   "double __melodee_temp_018 = celltype == 2;\n"
   "double Jrelp_inf_001;\n"
   "if (__melodee_temp_018)\n"
   "{\n"
   "   Jrelp_inf_001 = 1.7*Jrelp_inf;\n"
   "}\n"
   "else\n"
   "{\n"
   "   Jrelp_inf_001 = Jrelp_inf;\n"
   "}\n"
   "double tau_relp = btp/(1.0 + 0.0123/cajsr);\n"
   "double __melodee_temp_019 = tau_relp < 0.001;\n"
   "double tau_relp_001;\n"
   "if (__melodee_temp_019)\n"
   "{\n"
   "   tau_relp_001 = 0.00100000000000000;\n"
   "}\n"
   "else\n"
   "{\n"
   "   tau_relp_001 = tau_relp;\n"
   "}\n"
   "double Jrelp_diff = (-Jrelp + Jrelp_inf_001)/tau_relp_001;\n"
   "double fJrelp = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double Jrel = Jrelnp*(-fJrelp + 1.0) + Jrelp*fJrelp;\n"
   "double __melodee_temp_020 = Jrel*JrelStiffConst > cajsr;\n"
   "double Jrel_001;\n"
   "if (__melodee_temp_020)\n"
   "{\n"
   "   Jrel_001 = cajsr/JrelStiffConst;\n"
   "}\n"
   "else\n"
   "{\n"
   "   Jrel_001 = Jrel;\n"
   "}\n"
   "double Jupnp = 0.004375*cai/(cai + 0.00092);\n"
   "double Jupp = 0.01203125*cai/(cai + 0.00075);\n"
   "double __melodee_temp_021 = celltype == 1;\n"
   "double Jupnp_001;\n"
   "double Jupp_001;\n"
   "if (__melodee_temp_021)\n"
   "{\n"
   "   Jupnp_001 = 1.3*Jupnp;\n"
   "   Jupp_001 = 1.3*Jupp;\n"
   "}\n"
   "else\n"
   "{\n"
   "   Jupnp_001 = Jupnp;\n"
   "   Jupp_001 = Jupp;\n"
   "}\n"
   "double fJupp = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double Jleak = 0.0002625*cansr;\n"
   "double Jup = -Jleak + Jupnp_001*(-fJupp + 1.0) + Jupp_001*fJupp;\n"
   "double Jtr = -0.01*cajsr + 0.01*cansr;\n"
   "double cmdnmax = 0.0500000000000000;\n"
   "double __melodee_temp_022 = celltype == 1;\n"
   "double cmdnmax_001;\n"
   "if (__melodee_temp_022)\n"
   "{\n"
   "   cmdnmax_001 = 1.3*cmdnmax;\n"
   "}\n"
   "else\n"
   "{\n"
   "   cmdnmax_001 = cmdnmax;\n"
   "}\n"
   "double kmcmdn = 0.00238000000000000;\n"
   "double trpnmax = 0.0700000000000000;\n"
   "double kmtrpn = 0.000500000000000000;\n"
   "double BSRmax = 0.0470000000000000;\n"
   "double KmBSR = 0.000870000000000000;\n"
   "double BSLmax = 1.12400000000000;\n"
   "double KmBSL = 0.00870000000000000;\n"
   "double csqnmax = 10.0000000000000;\n"
   "double kmcsqn = 0.800000000000000;\n"
   "double nass_diff = Acap*(-ICaNa - 3.0*INaCa_ss)/(F*vss) - JdiffNa;\n"
   "double ki_diff = Acap*(-IK1 - IKb - IKr - IKs + 2.0*INaK - Ito)/(F*vmyo) + JdiffK*vss/vmyo;\n"
   "double kss_diff = -Acap*ICaK/(F*vss) - JdiffK;\n"
   "double Bcai = 1.0/(cmdnmax_001*kmcmdn*1.0/(cai + kmcmdn)/(cai + kmcmdn) + kmtrpn*trpnmax*1.0/(cai + kmtrpn)/(cai + kmtrpn) + 1.0);\n"
   "double cai_diff = Bcai*(0.5*Acap*(-ICab + 2.0*INaCa_i - IpCa)/(F*vmyo) + Jdiff*vss/vmyo - Jup*vnsr/vmyo);\n"
   "double Bcass = 1.0/(BSLmax*KmBSL*1.0/(KmBSL + cass)/(KmBSL + cass) + BSRmax*KmBSR*1.0/(KmBSR + cass)/(KmBSR + cass) + 1.0);\n"
   "double cass_diff = Bcass*(0.5*Acap*(-ICaL + 2.0*INaCa_ss)/(F*vss) - Jdiff + Jrel_001*vjsr/vss);\n"
   "double cansr_diff = -Jtr*vjsr/vnsr + Jup;\n"
   "double Bcajsr = 1.0/(csqnmax*kmcsqn*1.0/(cajsr + kmcsqn)/(cajsr + kmcsqn) + 1.0);\n"
   "double cajsr_diff = Bcajsr*(-Jrel_001 + Jtr);\n"
   "double _expensive_functions_065 = exp(-0.101306858474319*v);\n"
   "double mss = 1.0/(0.0181567590194406*_expensive_functions_065 + 1.0);\n"
   "double _expensive_functions_066 = exp(0.0287604256542997*v);\n"
   "double _expensive_functions_067 = exp(-0.167926112510495*v);\n"
   "double tm = 1.0/(9.45490463856472*_expensive_functions_066 + 1.93141135585369e-5*_expensive_functions_067);\n"
   "double _expensive_functions_068 = exp(0.164311534669734*v);\n"
   "double hss = 1.0/(823588.444960838*_expensive_functions_068 + 1);\n"
   "double _expensive_functions_069 = exp(-0.159108989657916*v);\n"
   "double _expensive_functions_070 = exp(0.0493339911198816*v);\n"
   "double thf = 1.0/(1.18385695828909e-5*_expensive_functions_069 + 6.30554918581728*_expensive_functions_070);\n"
   "double jss = hss;\n"
   "double _expensive_functions_071 = exp(0.0260078023407022*v);\n"
   "double _expensive_functions_072 = exp(-0.120758362516604*v);\n"
   "double tj = 2.038 + 1.0/(0.313193639473877*_expensive_functions_071 + 1.13152820955901e-7*_expensive_functions_072);\n"
   "double _expensive_functions_073 = exp(-0.035650623885918*v);\n"
   "double _expensive_functions_074 = exp(0.0176491351923756*v);\n"
   "double ths = 1.0/(0.00516467023538179*_expensive_functions_073 + 0.369876193720963*_expensive_functions_074);\n"
   "double _expensive_functions_075 = exp(0.164311534669734*v);\n"
   "double hspss = 1.0/(2281075.81669719*_expensive_functions_075 + 1);\n"
   "double thsp = 3.0*ths;\n"
   "double tjp = 1.46*tj;\n"
   "double fINap = 1.0/(1.0 + KmCaMK/CaMKa);\n"
   "double INa = m*m*m*GNa*(-ENa + v)*(fINap*hp*jp + h*j*(-fINap + 1.0));\n"
   "double tmL = tm;\n"
   "double nai_diff = Acap*(-INa - 3.0*INaCa_i - 3.0*INaK - INaL - INab)/(F*vmyo) + JdiffNa*vss/vmyo;\n"
   "double Iion = ICaK + ICaL + ICaNa + ICab + IK1 + IKb + IKr + IKs + INa + INaCa_i + INaCa_ss + INaK + INaL + INab + IpCa + Ito;\n"
   "double Iion_001 = Iion;\n"
   "double _a_RLA = exp(-_dt/ta);\n"
   "double _a_RLB = -ass*(_a_RLA - 1);\n"
   "double _ap_RLA = exp(-_dt/ta);\n"
   "double _ap_RLB = -apss*(_ap_RLA - 1);\n"
   "double _d_RLA = exp(-_dt/td);\n"
   "double _d_RLB = -dss*(_d_RLA - 1);\n"
   "double _fcaf_RLA = exp(-_dt/tfcaf);\n"
   "double _fcaf_RLB = -fcass*(_fcaf_RLA - 1);\n"
   "double _fcafp_RLA = exp(-_dt/tfcafp);\n"
   "double _fcafp_RLB = -fcass*(_fcafp_RLA - 1);\n"
   "double _fcas_RLA = exp(-_dt/tfcas);\n"
   "double _fcas_RLB = -fcass*(_fcas_RLA - 1);\n"
   "double _ff_RLA = exp(-_dt/tff);\n"
   "double _ff_RLB = -fss*(_ff_RLA - 1);\n"
   "double _ffp_RLA = exp(-_dt/tffp);\n"
   "double _ffp_RLB = -fss*(_ffp_RLA - 1);\n"
   "double _fs_RLA = exp(-_dt/tfs);\n"
   "double _fs_RLB = -fss*(_fs_RLA - 1);\n"
   "double _hL_RLA = exp(-_dt/thL);\n"
   "double _hL_RLB = -hLss*(_hL_RLA - 1);\n"
   "double _hLp_RLA = exp(-_dt/thLp);\n"
   "double _hLp_RLB = -hLpss*(_hLp_RLA - 1);\n"
   "double _hf_RLA = exp(-_dt/thf);\n"
   "double _hf_RLB = -hss*(_hf_RLA - 1);\n"
   "double _hs_RLA = exp(-_dt/ths);\n"
   "double _hs_RLB = -hss*(_hs_RLA - 1);\n"
   "double _hsp_RLA = exp(-_dt/thsp);\n"
   "double _hsp_RLB = -hspss*(_hsp_RLA - 1);\n"
   "double _iF_RLA = exp(-_dt/tiF_001);\n"
   "double _iF_RLB = -iss*(_iF_RLA - 1);\n"
   "double _iFp_RLA = exp(-_dt/tiFp);\n"
   "double _iFp_RLB = -iss*(_iFp_RLA - 1);\n"
   "double _iS_RLA = exp(-_dt/tiS_001);\n"
   "double _iS_RLB = -iss*(_iS_RLA - 1);\n"
   "double _iSp_RLA = exp(-_dt/tiSp);\n"
   "double _iSp_RLB = -iss*(_iSp_RLA - 1);\n"
   "double _j_RLA = exp(-_dt/tj);\n"
   "double _j_RLB = -jss*(_j_RLA - 1);\n"
   "double _jca_RLA = exp(-_dt/tjca);\n"
   "double _jca_RLB = -fcass*(_jca_RLA - 1);\n"
   "double _jp_RLA = exp(-_dt/tjp);\n"
   "double _jp_RLB = -jss*(_jp_RLA - 1);\n"
   "double _m_RLA = exp(-_dt/tm);\n"
   "double _m_RLB = -mss*(_m_RLA - 1);\n"
   "double _mL_RLA = exp(-_dt/tmL);\n"
   "double _mL_RLB = -mLss*(_mL_RLA - 1);\n"
   "double _nca_RLA = exp(-_dt*km2n);\n"
   "double _nca_RLB = -anca*k2n*(_nca_RLA - 1)/km2n;\n"
   "double _xk1_RLA = exp(-_dt/txk1);\n"
   "double _xk1_RLB = -xk1ss*(_xk1_RLA - 1);\n"
   "double _xrf_RLA = exp(-_dt/txrf);\n"
   "double _xrf_RLB = -xrss*(_xrf_RLA - 1);\n"
   "double _xrs_RLA = exp(-_dt/txrs);\n"
   "double _xrs_RLB = -xrss*(_xrs_RLA - 1);\n"
   "double _xs1_RLA = exp(-_dt/txs1);\n"
   "double _xs1_RLB = -xs1ss*(_xs1_RLA - 1);\n"
   "double _xs2_RLA = exp(-_dt/txs2);\n"
   "double _xs2_RLB = -xs2ss*(_xs2_RLA - 1);\n"
   "\n\n//EDIT STATE\n"
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
   "_state[_ii+a_off*_nCells] = _a_RLA*a + _a_RLB;\n"
   "_state[_ii+ap_off*_nCells] = _ap_RLA*ap + _ap_RLB;\n"
   "_state[_ii+d_off*_nCells] = _d_RLA*d + _d_RLB;\n"
   "_state[_ii+fcaf_off*_nCells] = _fcaf_RLA*fcaf + _fcaf_RLB;\n"
   "_state[_ii+fcafp_off*_nCells] = _fcafp_RLA*fcafp + _fcafp_RLB;\n"
   "_state[_ii+fcas_off*_nCells] = _fcas_RLA*fcas + _fcas_RLB;\n"
   "_state[_ii+ff_off*_nCells] = _ff_RLA*ff + _ff_RLB;\n"
   "_state[_ii+ffp_off*_nCells] = _ffp_RLA*ffp + _ffp_RLB;\n"
   "_state[_ii+fs_off*_nCells] = _fs_RLA*fs + _fs_RLB;\n"
   "_state[_ii+hL_off*_nCells] = _hL_RLA*hL + _hL_RLB;\n"
   "_state[_ii+hLp_off*_nCells] = _hLp_RLA*hLp + _hLp_RLB;\n"
   "_state[_ii+hf_off*_nCells] = _hf_RLA*hf + _hf_RLB;\n"
   "_state[_ii+hs_off*_nCells] = _hs_RLA*hs + _hs_RLB;\n"
   "_state[_ii+hsp_off*_nCells] = _hsp_RLA*hsp + _hsp_RLB;\n"
   "_state[_ii+iF_off*_nCells] = _iF_RLA*iF + _iF_RLB;\n"
   "_state[_ii+iFp_off*_nCells] = _iFp_RLA*iFp + _iFp_RLB;\n"
   "_state[_ii+iS_off*_nCells] = _iS_RLA*iS + _iS_RLB;\n"
   "_state[_ii+iSp_off*_nCells] = _iSp_RLA*iSp + _iSp_RLB;\n"
   "_state[_ii+j_off*_nCells] = _j_RLA*j + _j_RLB;\n"
   "_state[_ii+jca_off*_nCells] = _jca_RLA*jca + _jca_RLB;\n"
   "_state[_ii+jp_off*_nCells] = _jp_RLA*jp + _jp_RLB;\n"
   "_state[_ii+m_off*_nCells] = _m_RLA*m + _m_RLB;\n"
   "_state[_ii+mL_off*_nCells] = _mL_RLA*mL + _mL_RLB;\n"
   "_state[_ii+nca_off*_nCells] = _nca_RLA*nca + _nca_RLB;\n"
   "_state[_ii+xk1_off*_nCells] = _xk1_RLA*xk1 + _xk1_RLB;\n"
   "_state[_ii+xrf_off*_nCells] = _xrf_RLA*xrf + _xrf_RLB;\n"
   "_state[_ii+xrs_off*_nCells] = _xrs_RLA*xrs + _xrs_RLB;\n"
   "_state[_ii+xs1_off*_nCells] = _xs1_RLA*xs1 + _xs1_RLB;\n"
   "_state[_ii+xs2_off*_nCells] = _xs2_RLA*xs2 + _xs2_RLB;\n"
   "_dVm[_ii] = -Iion_001;\n"
   "}\n";

   _kernel_program = ss.str();
   cout << ss.str();
}

string ThisReaction::methodName() const
{
   return "ORd";
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
   _xrf_off,
   _xrs_off,
   _xs1_off,
   _xs2_off,
   NUMSTATES
};
   
ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   stateTransport_.setup(std::vector<double>(nCells_*NUMSTATES));
   __cachedDt = __dt;
   blockSize_ = -1;
}

void ThisReaction::createInterpolants(const double _dt) {

}

void ThisReaction::calc(double _dt, const VectorDouble32& __Vm,
                       const vector<double>& __iStim , VectorDouble32& __dVm)
{
   vector<double>& state(stateTransport_.modifyOnDevice());

   const double* VmRaw=&__Vm[0];
   double* dVmRaw=&__dVm[0];
   const double* iStimRaw=&__iStim[0];
   double* stateDataRaw=&state[0];

   #if _OPENMP >= 201511
   //#pragma omp target data use_device_ptr(VmRaw, dVmRaw, iStimRaw, stateDataRaw)
   #endif
   {
      jitify::Program _program = _kernel_cache.program(_kernel_program, 0);

      int errorCode=-1;
      if (blockSize_ == -1) { blockSize_ = 1024; }
      while(1)
      {
         int errorCode = _program
            .kernel("ORd_kernel")
            .instantiate()
            .configure(dim3((nCells_+blockSize_-1)/blockSize_,1,1),dim3(blockSize_,1,1))
            .launch(VmRaw,iStimRaw,dVmRaw,stateDataRaw);
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
      }
   }
}

void ThisReaction::initializeMembraneVoltage(VectorDouble32& __Vm)
{
   assert(__Vm.size() >= nCells_);


   double V_init = -85;
   double V = V_init;
   double m_init = 0;
   double m = m_init;
   double hf_init = 1;
   double hf = hf_init;
   double hs_init = 1;
   double hs = hs_init;
   double j_init = 1;
   double j = j_init;
   double hsp_init = 1;
   double hsp = hsp_init;
   double jp_init = 1;
   double jp = jp_init;
   double nai_init = 7;
   double nai = nai_init;
   double nass_init = 7;
   double nass = nass_init;
   double ki_init = 145;
   double ki = ki_init;
   double kss_init = 145;
   double kss = kss_init;
   double cai_init = 0.000100000000000000;
   double cai = cai_init;
   double cass_init = 0.000100000000000000;
   double cass = cass_init;
   double cansr_init = 1.20000000000000;
   double cansr = cansr_init;
   double cajsr_init = 1.20000000000000;
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
      
   
   vector<double>& stateData(stateTransport_.modifyOnHost());
   
   for (int iCell=0; iCell<nCells_; iCell++)
   {
      stateData[_CaMKt_off*nCells_+iCell] = CaMKt;
      stateData[_Jrelnp_off*nCells_+iCell] = Jrelnp;
      stateData[_Jrelp_off*nCells_+iCell] = Jrelp;
      stateData[_a_off*nCells_+iCell] = a;
      stateData[_ap_off*nCells_+iCell] = ap;
      stateData[_cai_off*nCells_+iCell] = cai;
      stateData[_cajsr_off*nCells_+iCell] = cajsr;
      stateData[_cansr_off*nCells_+iCell] = cansr;
      stateData[_cass_off*nCells_+iCell] = cass;
      stateData[_d_off*nCells_+iCell] = d;
      stateData[_fcaf_off*nCells_+iCell] = fcaf;
      stateData[_fcafp_off*nCells_+iCell] = fcafp;
      stateData[_fcas_off*nCells_+iCell] = fcas;
      stateData[_ff_off*nCells_+iCell] = ff;
      stateData[_ffp_off*nCells_+iCell] = ffp;
      stateData[_fs_off*nCells_+iCell] = fs;
      stateData[_hL_off*nCells_+iCell] = hL;
      stateData[_hLp_off*nCells_+iCell] = hLp;
      stateData[_hf_off*nCells_+iCell] = hf;
      stateData[_hs_off*nCells_+iCell] = hs;
      stateData[_hsp_off*nCells_+iCell] = hsp;
      stateData[_iF_off*nCells_+iCell] = iF;
      stateData[_iFp_off*nCells_+iCell] = iFp;
      stateData[_iS_off*nCells_+iCell] = iS;
      stateData[_iSp_off*nCells_+iCell] = iSp;
      stateData[_j_off*nCells_+iCell] = j;
      stateData[_jca_off*nCells_+iCell] = jca;
      stateData[_jp_off*nCells_+iCell] = jp;
      stateData[_ki_off*nCells_+iCell] = ki;
      stateData[_kss_off*nCells_+iCell] = kss;
      stateData[_m_off*nCells_+iCell] = m;
      stateData[_mL_off*nCells_+iCell] = mL;
      stateData[_nai_off*nCells_+iCell] = nai;
      stateData[_nass_off*nCells_+iCell] = nass;
      stateData[_nca_off*nCells_+iCell] = nca;
      stateData[_xk1_off*nCells_+iCell] = xk1;
      stateData[_xrf_off*nCells_+iCell] = xrf;
      stateData[_xrs_off*nCells_+iCell] = xrs;
      stateData[_xs1_off*nCells_+iCell] = xs1;
      stateData[_xs2_off*nCells_+iCell] = xs2;
   }

   __Vm.assign(__Vm.size(), V_init);
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
   xrf_handle,
   xrs_handle,
   xs1_handle,
   xs2_handle,
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
   else if (varName == "xrf") { return xrf_handle; }
   else if (varName == "xrs") { return xrs_handle; }
   else if (varName == "xs1") { return xs1_handle; }
   else if (varName == "xs2") { return xs2_handle; }
   return -1;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
   vector<double>& stateData(stateTransport_.modifyOnHost());
   
   if (0) {}
   else if (varHandle == CaMKt_handle) { stateData[_CaMKt_off*nCells_+iCell] = value; }
   else if (varHandle == Jrelnp_handle) { stateData[_Jrelnp_off*nCells_+iCell] = value; }
   else if (varHandle == Jrelp_handle) { stateData[_Jrelp_off*nCells_+iCell] = value; }
   else if (varHandle == a_handle) { stateData[_a_off*nCells_+iCell] = value; }
   else if (varHandle == ap_handle) { stateData[_ap_off*nCells_+iCell] = value; }
   else if (varHandle == cai_handle) { stateData[_cai_off*nCells_+iCell] = value; }
   else if (varHandle == cajsr_handle) { stateData[_cajsr_off*nCells_+iCell] = value; }
   else if (varHandle == cansr_handle) { stateData[_cansr_off*nCells_+iCell] = value; }
   else if (varHandle == cass_handle) { stateData[_cass_off*nCells_+iCell] = value; }
   else if (varHandle == d_handle) { stateData[_d_off*nCells_+iCell] = value; }
   else if (varHandle == fcaf_handle) { stateData[_fcaf_off*nCells_+iCell] = value; }
   else if (varHandle == fcafp_handle) { stateData[_fcafp_off*nCells_+iCell] = value; }
   else if (varHandle == fcas_handle) { stateData[_fcas_off*nCells_+iCell] = value; }
   else if (varHandle == ff_handle) { stateData[_ff_off*nCells_+iCell] = value; }
   else if (varHandle == ffp_handle) { stateData[_ffp_off*nCells_+iCell] = value; }
   else if (varHandle == fs_handle) { stateData[_fs_off*nCells_+iCell] = value; }
   else if (varHandle == hL_handle) { stateData[_hL_off*nCells_+iCell] = value; }
   else if (varHandle == hLp_handle) { stateData[_hLp_off*nCells_+iCell] = value; }
   else if (varHandle == hf_handle) { stateData[_hf_off*nCells_+iCell] = value; }
   else if (varHandle == hs_handle) { stateData[_hs_off*nCells_+iCell] = value; }
   else if (varHandle == hsp_handle) { stateData[_hsp_off*nCells_+iCell] = value; }
   else if (varHandle == iF_handle) { stateData[_iF_off*nCells_+iCell] = value; }
   else if (varHandle == iFp_handle) { stateData[_iFp_off*nCells_+iCell] = value; }
   else if (varHandle == iS_handle) { stateData[_iS_off*nCells_+iCell] = value; }
   else if (varHandle == iSp_handle) { stateData[_iSp_off*nCells_+iCell] = value; }
   else if (varHandle == j_handle) { stateData[_j_off*nCells_+iCell] = value; }
   else if (varHandle == jca_handle) { stateData[_jca_off*nCells_+iCell] = value; }
   else if (varHandle == jp_handle) { stateData[_jp_off*nCells_+iCell] = value; }
   else if (varHandle == ki_handle) { stateData[_ki_off*nCells_+iCell] = value; }
   else if (varHandle == kss_handle) { stateData[_kss_off*nCells_+iCell] = value; }
   else if (varHandle == m_handle) { stateData[_m_off*nCells_+iCell] = value; }
   else if (varHandle == mL_handle) { stateData[_mL_off*nCells_+iCell] = value; }
   else if (varHandle == nai_handle) { stateData[_nai_off*nCells_+iCell] = value; }
   else if (varHandle == nass_handle) { stateData[_nass_off*nCells_+iCell] = value; }
   else if (varHandle == nca_handle) { stateData[_nca_off*nCells_+iCell] = value; }
   else if (varHandle == xk1_handle) { stateData[_xk1_off*nCells_+iCell] = value; }
   else if (varHandle == xrf_handle) { stateData[_xrf_off*nCells_+iCell] = value; }
   else if (varHandle == xrs_handle) { stateData[_xrs_off*nCells_+iCell] = value; }
   else if (varHandle == xs1_handle) { stateData[_xs1_off*nCells_+iCell] = value; }
   else if (varHandle == xs2_handle) { stateData[_xs2_off*nCells_+iCell] = value; }
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
   const vector<double>& stateData(stateTransport_.readOnHost());

   if (0) {}
   else if (varHandle == CaMKt_handle) { return stateData[_CaMKt_off*nCells_+iCell]; }
   else if (varHandle == Jrelnp_handle) { return stateData[_Jrelnp_off*nCells_+iCell]; }
   else if (varHandle == Jrelp_handle) { return stateData[_Jrelp_off*nCells_+iCell]; }
   else if (varHandle == a_handle) { return stateData[_a_off*nCells_+iCell]; }
   else if (varHandle == ap_handle) { return stateData[_ap_off*nCells_+iCell]; }
   else if (varHandle == cai_handle) { return stateData[_cai_off*nCells_+iCell]; }
   else if (varHandle == cajsr_handle) { return stateData[_cajsr_off*nCells_+iCell]; }
   else if (varHandle == cansr_handle) { return stateData[_cansr_off*nCells_+iCell]; }
   else if (varHandle == cass_handle) { return stateData[_cass_off*nCells_+iCell]; }
   else if (varHandle == d_handle) { return stateData[_d_off*nCells_+iCell]; }
   else if (varHandle == fcaf_handle) { return stateData[_fcaf_off*nCells_+iCell]; }
   else if (varHandle == fcafp_handle) { return stateData[_fcafp_off*nCells_+iCell]; }
   else if (varHandle == fcas_handle) { return stateData[_fcas_off*nCells_+iCell]; }
   else if (varHandle == ff_handle) { return stateData[_ff_off*nCells_+iCell]; }
   else if (varHandle == ffp_handle) { return stateData[_ffp_off*nCells_+iCell]; }
   else if (varHandle == fs_handle) { return stateData[_fs_off*nCells_+iCell]; }
   else if (varHandle == hL_handle) { return stateData[_hL_off*nCells_+iCell]; }
   else if (varHandle == hLp_handle) { return stateData[_hLp_off*nCells_+iCell]; }
   else if (varHandle == hf_handle) { return stateData[_hf_off*nCells_+iCell]; }
   else if (varHandle == hs_handle) { return stateData[_hs_off*nCells_+iCell]; }
   else if (varHandle == hsp_handle) { return stateData[_hsp_off*nCells_+iCell]; }
   else if (varHandle == iF_handle) { return stateData[_iF_off*nCells_+iCell]; }
   else if (varHandle == iFp_handle) { return stateData[_iFp_off*nCells_+iCell]; }
   else if (varHandle == iS_handle) { return stateData[_iS_off*nCells_+iCell]; }
   else if (varHandle == iSp_handle) { return stateData[_iSp_off*nCells_+iCell]; }
   else if (varHandle == j_handle) { return stateData[_j_off*nCells_+iCell]; }
   else if (varHandle == jca_handle) { return stateData[_jca_off*nCells_+iCell]; }
   else if (varHandle == jp_handle) { return stateData[_jp_off*nCells_+iCell]; }
   else if (varHandle == ki_handle) { return stateData[_ki_off*nCells_+iCell]; }
   else if (varHandle == kss_handle) { return stateData[_kss_off*nCells_+iCell]; }
   else if (varHandle == m_handle) { return stateData[_m_off*nCells_+iCell]; }
   else if (varHandle == mL_handle) { return stateData[_mL_off*nCells_+iCell]; }
   else if (varHandle == nai_handle) { return stateData[_nai_off*nCells_+iCell]; }
   else if (varHandle == nass_handle) { return stateData[_nass_off*nCells_+iCell]; }
   else if (varHandle == nca_handle) { return stateData[_nca_off*nCells_+iCell]; }
   else if (varHandle == xk1_handle) { return stateData[_xk1_off*nCells_+iCell]; }
   else if (varHandle == xrf_handle) { return stateData[_xrf_off*nCells_+iCell]; }
   else if (varHandle == xrs_handle) { return stateData[_xrs_off*nCells_+iCell]; }
   else if (varHandle == xs1_handle) { return stateData[_xs1_off*nCells_+iCell]; }
   else if (varHandle == xs2_handle) { return stateData[_xs2_off*nCells_+iCell]; }
   return NAN;
}

double ThisReaction::getValue(int iCell, int varHandle, double V) const
{
   const vector<double>& stateData(stateTransport_.readOnHost());

   const double CaMKt=stateData[_CaMKt_off*nCells_+iCell];
   const double Jrelnp=stateData[_Jrelnp_off*nCells_+iCell];
   const double Jrelp=stateData[_Jrelp_off*nCells_+iCell];
   const double a=stateData[_a_off*nCells_+iCell];
   const double ap=stateData[_ap_off*nCells_+iCell];
   const double cai=stateData[_cai_off*nCells_+iCell];
   const double cajsr=stateData[_cajsr_off*nCells_+iCell];
   const double cansr=stateData[_cansr_off*nCells_+iCell];
   const double cass=stateData[_cass_off*nCells_+iCell];
   const double d=stateData[_d_off*nCells_+iCell];
   const double fcaf=stateData[_fcaf_off*nCells_+iCell];
   const double fcafp=stateData[_fcafp_off*nCells_+iCell];
   const double fcas=stateData[_fcas_off*nCells_+iCell];
   const double ff=stateData[_ff_off*nCells_+iCell];
   const double ffp=stateData[_ffp_off*nCells_+iCell];
   const double fs=stateData[_fs_off*nCells_+iCell];
   const double hL=stateData[_hL_off*nCells_+iCell];
   const double hLp=stateData[_hLp_off*nCells_+iCell];
   const double hf=stateData[_hf_off*nCells_+iCell];
   const double hs=stateData[_hs_off*nCells_+iCell];
   const double hsp=stateData[_hsp_off*nCells_+iCell];
   const double iF=stateData[_iF_off*nCells_+iCell];
   const double iFp=stateData[_iFp_off*nCells_+iCell];
   const double iS=stateData[_iS_off*nCells_+iCell];
   const double iSp=stateData[_iSp_off*nCells_+iCell];
   const double j=stateData[_j_off*nCells_+iCell];
   const double jca=stateData[_jca_off*nCells_+iCell];
   const double jp=stateData[_jp_off*nCells_+iCell];
   const double ki=stateData[_ki_off*nCells_+iCell];
   const double kss=stateData[_kss_off*nCells_+iCell];
   const double m=stateData[_m_off*nCells_+iCell];
   const double mL=stateData[_mL_off*nCells_+iCell];
   const double nai=stateData[_nai_off*nCells_+iCell];
   const double nass=stateData[_nass_off*nCells_+iCell];
   const double nca=stateData[_nca_off*nCells_+iCell];
   const double xk1=stateData[_xk1_off*nCells_+iCell];
   const double xrf=stateData[_xrf_off*nCells_+iCell];
   const double xrs=stateData[_xrs_off*nCells_+iCell];
   const double xs1=stateData[_xs1_off*nCells_+iCell];
   const double xs2=stateData[_xs2_off*nCells_+iCell];
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
