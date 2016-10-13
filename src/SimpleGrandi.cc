/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "SimpleGrandi.hh"
#include "object_cc.hh"
#include "Anatomy.hh"
#include <cmath>

#define pi 3.141592653589793238462643383279502884197169399375105820974944592307816406286

using namespace std;

namespace scanReaction 
{

#define setDefault(name, value) objectGet(obj, #name, reaction->name, #value)
   
   Reaction* scanSimpleGrandi(OBJECT* obj, const Anatomy& anatomy)
   {
      SimpleGrandi::ThisReaction* reaction = new SimpleGrandi::ThisReaction(anatomy);

      //override the defaults
      //EDIT_FLAGS
      setDefault(AF, 0);
      setDefault(ISO, 0);
      setDefault(RA, 1);

      //EDIT_PARAMETERS
      /*setDefault(ks, 25); // [1/ms]      
      setDefault(Vmax_SRCaP, 5.3114e-3);  // [mM/msec] (286 umol/L cytosol/sec)
      setDefault(GCaB, 6.0643e-4);    // [uA/uF] 3
      setDefault(GClCa, 0.0548);   // [mS/uF]
      setDefault(GClB, 9e-3);        // [mS/uF]
      setDefault(gkp, 0.002);
      setDefault(IbarNaK, 1.26);     // [uA/uF]
      setDefault(IbarSLCaP, 0.0471); // IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
      */

      return reaction;
   }
#undef setDefault

}

namespace SimpleGrandi 
{

inline double pow(const double x, const int p)
{
   double ret=1;
   if (p > 0) 
   {
      for (int ii=0; ii<p; ii++) 
      {
         ret *= x;
      }
   }
   else
   {
      for (int ii=0; ii<-p; ii++) 
      {
         ret /= x;
      }
   }
   return ret;
}
   
string ThisReaction::methodName() const
{
   return "SimpleGrandi";
}

const char* varNames[] = 
{
   //EDIT_STATE
   "m",
   "h",
   "j",
   "mL",
   "hL",
   "xtf",
   "ytf",
   "xkr",
   "xks",
   "xkur",
   "ykur",
   "d",
   "f",
   "fcaBj",
   "fcaBsl",
   "RyRr",
   "RyRo",
   "RyRi",
   "NaBj",
   "NaBsl",
   "TnCL",
   "TnCHc",
   "TnCHm",
   "CaM",
   "Myc",
   "Mym",
   "SRB",
   "SLLj",
   "SLLsl",
   "SLHj",
   "SLHsl",
   "Csqnb",
   "Naj",
   "Nasl",
   "Nai",
   "Ki",
   "Casr",
   "Caj",
   "Casl",
   "Cai"
};
#define NUMVARS (sizeof(varNames)/sizeof(char*))

int getVarOffset(const std::string& varName)
{
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      if (varNames[ivar] == varName) 
      {
         return ivar;
      }
   }
   assert(0 && "Control should never get here.");
   return -1;
}

void assertStateOrderAndVarNamesAgree(void)
{
   State s;
#define checkVarOrder(x) assert(reinterpret_cast<double*>(&s)+getVarOffset(#x) == &s . x)

   int STATIC_ASSERT_checkAllDouble[(NUMVARS == sizeof(s)/sizeof(double))? 1: 0];

   //EDIT_STATE
   checkVarOrder(m);
   checkVarOrder(h);
   checkVarOrder(j);
   checkVarOrder(mL);
   checkVarOrder(hL);
   checkVarOrder(xtf);
   checkVarOrder(ytf);
   checkVarOrder(xkr);
   checkVarOrder(xks);
   checkVarOrder(xkur);
   checkVarOrder(ykur);
   checkVarOrder(d);
   checkVarOrder(f);
   checkVarOrder(fcaBj);
   checkVarOrder(fcaBsl);
   checkVarOrder(RyRr);
   checkVarOrder(RyRo);
   checkVarOrder(RyRi);
   checkVarOrder(NaBj);
   checkVarOrder(NaBsl);
   checkVarOrder(TnCL);
   checkVarOrder(TnCHc);
   checkVarOrder(TnCHm);
   checkVarOrder(CaM);
   checkVarOrder(Myc);
   checkVarOrder(Mym);
   checkVarOrder(SRB);
   checkVarOrder(SLLj);
   checkVarOrder(SLLsl);
   checkVarOrder(SLHj);
   checkVarOrder(SLHsl);
   checkVarOrder(Csqnb);
   checkVarOrder(Naj);
   checkVarOrder(Nasl);
   checkVarOrder(Nai);
   checkVarOrder(Ki);
   checkVarOrder(Casr);
   checkVarOrder(Caj);
   checkVarOrder(Casl);
   checkVarOrder(Cai);
}
   
ThisReaction::ThisReaction(const Anatomy& anatomy)
: nCells_(anatomy.nLocal())
{
   assertStateOrderAndVarNamesAgree();
   state_.resize(nCells_);
   perCellFlags_.resize(nCells_);
   perCellParameters_.resize(nCells_);
}

void ThisReaction::calc(double dt, const VectorDouble32& Vm,
                       const vector<double>& iStim , VectorDouble32& dVm)
{
   for (unsigned ii=0; ii<nCells_; ++ii)
   {

      //set Vm
      const double v = Vm[ii];
      const double istim = iStim[ii];

      //set all state variables
      //EDIT_STATE
      const double m= state_[ii].m;
      const double h= state_[ii].h;
      const double j= state_[ii].j;
      const double mL= state_[ii].mL;
      const double hL= state_[ii].hL;
      const double xtf= state_[ii].xtf;
      const double ytf= state_[ii].ytf;
      const double xkr= state_[ii].xkr;
      const double xks= state_[ii].xks;
      const double xkur= state_[ii].xkur;
      const double ykur= state_[ii].ykur;
      const double d= state_[ii].d;
      const double f= state_[ii].f;
      const double fcaBj= state_[ii].fcaBj;
      const double fcaBsl= state_[ii].fcaBsl;
      const double RyRr= state_[ii].RyRr;
      const double RyRo= state_[ii].RyRo;
      const double RyRi= state_[ii].RyRi;
      const double NaBj= state_[ii].NaBj;
      const double NaBsl= state_[ii].NaBsl;
      const double TnCL= state_[ii].TnCL;
      const double TnCHc= state_[ii].TnCHc;
      const double TnCHm= state_[ii].TnCHm;
      const double CaM= state_[ii].CaM;
      const double Myc= state_[ii].Myc;
      const double Mym= state_[ii].Mym;
      const double SRB= state_[ii].SRB;
      const double SLLj= state_[ii].SLLj;
      const double SLLsl= state_[ii].SLLsl;
      const double SLHj= state_[ii].SLHj;
      const double SLHsl= state_[ii].SLHsl;
      const double Csqnb= state_[ii].Csqnb;
      const double Naj= state_[ii].Naj;
      const double Nasl= state_[ii].Nasl;
      const double Nai= state_[ii].Nai;
      const double Ki= state_[ii].Ki;
      const double Casr= state_[ii].Casr;
      const double Caj= state_[ii].Caj;
      const double Casl= state_[ii].Casl;
      const double Cai= state_[ii].Cai;

      //set per-cell flags
      //EDIT_PERCELL_FLAGS
      
      //set per-cell parameters
      //EDIT_PERCELL_PARAMETERS

      // Constants
      const double R = 8314.0;       // [J/kmol*K]  
      const double Frdy = 96485.0;   // [C/mol]  
      const double Temp = 310.0;     // [K]
      const double FoRT = Frdy/R/Temp;
      const double Cmem = 1.1e-10;   // [F] membrane capacitance 1.3810e-10;//
      const double Qpow = (Temp-310.0)/10.0;

      // Cell geometry
      const double cellLength = 100.0;     // cell length [um]113;//100
      const double cellRadius = 10.25;   // cell radius [um]12;//10.25
      const double Vcell = pi*cellRadius*cellRadius*cellLength*1.0e-15;    // [L]
      const double Vmyo = 0.65*Vcell;
      const double Vsr = 0.035*Vcell;
      const double Vsl = 0.02*Vcell;
      const double Vjunc = 0.0539*.01*Vcell; 
      const double Jca_juncsl =1.0/1.2134e12; // [L/msec] = 8.2413e-13
      const double Jca_slmyo = 1.0/2.68510e11; // [L/msec] = 3.2743e-12
      const double Jna_juncsl = 1.0/(1.6382e12/3.0*100.0); // [L/msec] = 6.1043e-13
      const double Jna_slmyo = 1.0/(1.8308e10/3.0*100.0);  // [L/msec] = 5.4621e-11

      // Fractional currents in compartments
      const double Fjunc = 0.11;
      const double Fsl = 1-Fjunc;
      const double Fjunc_CaL = 0.9;
      const double Fsl_CaL = 1-Fjunc_CaL;

      // Fixed ion concentrations     
      const double Cli = 15.0;   // Intracellular Cl  [mM]
      const double Clo = 150.0;  // Extracellular Cl  [mM]
      const double Ko = 5.4;   // Extracellular K   [mM]
      const double Nao = 140.0;  // Extracellular Na  [mM]
      const double Cao = 1.8;  // Extracellular Ca  [mM]
      const double Mgi = 1.0;    // Intracellular Mg  [mM]

      // Nernst Potentials
      const double ENa_junc = (1.0/FoRT)*log(Nao/Naj);     // [mV]
      const double ENa_sl = (1.0/FoRT)*log(Nao/Nasl);       // [mV]
      const double EK = (1.0/FoRT)*log(Ko/Ki);	        // [mV]
      const double pNaK=0.01833;
      const double EKs=(1.0/FoRT)*log((Ko+pNaK*Nao)/(Ki+pNaK*Nai));
      const double ECa_junc = (1.0/FoRT/2.0)*log(Cao/Caj);   // [mV]
      const double ECa_sl = (1.0/FoRT/2.0)*log(Cao/Casl);     // [mV]
      const double ECl = (1.0/FoRT)*log(Cli/Clo);            // [mV]

      //// Membrane Currents
      const double mss=1.0/(pow(1.0+exp(-(56.86+v)/9.03),2.0));
      const double taum=0.1292*exp(-pow((v+45.79)/15.54,2.0))+0.06487*exp(-pow((v-4.823)/51.12,2.0));
      double ah,bh,aj,bj;
      if (v >= -40)
	{
	  ah=0.0;
	  bh=0.77/(0.13*(1.0+exp(-(v+10.66)/11.1)));
	  aj=0.0;
	  bj=((0.6*exp(0.057*v))/(1.0+exp(-0.1*(v+32.0))));
	}
      else
	{
	  ah=0.057*exp(-(v+80.0)/6.8);
	  bh=2.7*exp(0.079*v)+3.1e5*exp(0.3485*v);
	  aj=((-2.5428e4*exp(0.2444*v)-6.948e-6*exp(-0.04391*v))*(v+37.78))/(1.0+exp(0.311*(v+79.23)));
	  bj=(0.02424*exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14)));  
	}
      const double tauh=1.0/(ah+bh);
      const double hss=1.0/(pow(1.0+exp((v+71.55)/7.43),2.0));
      const double tauj=1.0/(aj+bj);
      const double jss=1.0/(pow(1.0+exp((v+71.55)/7.43),2.0));
      const double diff_m=(mss-m)/taum;
      const double diff_h=(hss-h)/tauh;
      const double diff_j=(jss-j)/tauj;
      const double GNa=23.0*(1.0-0.1*AF);
      const double INa_junc=Fjunc*GNa*m*m*m*h*j*(v-ENa_junc);
      const double INa_sl=Fsl*GNa*m*m*m*h*j*(v-ENa_sl);

      const double aml=0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
      const double bml=0.08*exp(-v/11.0);
      const double taumL=1.0/(aml+bml);
      const double mLss=aml*taumL;
      const double hlinf=1.0/(1.0+exp((v+91.0)/6.1));
      const double tauhl=600.0;
      const double diff_mL=aml*(1.0-mL)-bml*mL;
      const double diff_hL=(hlinf-hL)/tauhl;
      const double GNaL=0.0025*(1.0-AF);
      const double INaL_junc=Fjunc*GNaL*mL*mL*mL*hL*(v-ENa_junc);
      const double INaL_sl=Fsl*GNaL*mL*mL*mL*hL*(v-ENa_sl);

      const double GNaB=0.597e-3;
      const double INaBk_junc=Fjunc*GNaB*(v-ENa_junc);
      const double INaBk_sl=Fsl*GNaB*(v-ENa_sl);

      const double KmNaip=11.0*(1.0-0.25*ISO);
      const double KmKo=1.5;
      const double sigma=(exp(Nao/67.3)-1.0)/7.0;
      const double fnak=1.0/(1.0+0.1245*exp(-0.1*v*FoRT)+0.0365*sigma*exp(-v*FoRT));
      const double IbarNaK=1.26;
      const double INAK_junc=Fjunc*IbarNaK*fnak*Ko/(1.0+pow(KmNaip/Naj,4.0))/(Ko+KmKo);
      const double INAK_sl=Fsl*IbarNaK*fnak*Ko/(1.0+pow(KmNaip/Nasl,4.0))/(Ko+KmKo);
      const double INAK=INAK_junc+INAK_sl;

      const double gkr =0.035*sqrt(Ko/5.4);
      const double xrss=1.0/(1.0+exp(-(v+10.0)/5.0));
      const double tauxr=550.0/(1.0+exp((-22.0-v)/9.0))*6.0/(1.0+exp((v-(-11.0))/9.0))+230.0/(1.0+exp((v-(-40.0))/20.0));
      const double diff_xkr=(xrss-xkr)/tauxr;
      const double rkr=1.0/(1.0+exp((v+74.0)/24.0));
      const double IKr=gkr*xkr*rkr*(v-EK);

      const double gks_junc=(1.0+1.0*AF+2.0*ISO)*0.0035;
      const double gks_sl=(1.0+1.0*AF+2.0*ISO)*0.0035;
      const double xsss=1.0/(1.0+exp(-(v+40.0*ISO+3.8)/14.25));
      const double tauxs=990.1/(1.0+exp(-(v+40.0*ISO+2.436)/14.12));
      const double diff_xks=(xsss-xks)/tauxs;
      const double IKs_junc=Fjunc*gks_junc*xks*xks*(v-EKs);
      const double IKs_sl=Fsl*gks_sl*xks*xks*(v-EKs);
      const double IKs=IKs_junc+IKs_sl;

      const double kp_kp=1.0/(1.0+exp(7.488-v/5.98));
      const double gkp=0.002;
      const double IKp_junc=Fjunc*gkp*kp_kp*(v-EK);
      const double IKp_sl=Fsl*gkp*kp_kp*(v-EK);
      const double IKp=IKp_junc+IKp_sl;

      const double xtss=1.0/(1.0+exp(-(v+1.0)/11.0));
      const double tauxtf=3.5*exp(-(pow(v/30.0,2.0)))+1.5;
      const double ytss=1.0/(1.0+exp((v+40.5)/11.5));
      const double tauytf=25.635*exp(-(pow((v+52.45)/15.8827,2.0)))+24.14;
      const double diff_xtf=(xtss-xtf)/tauxtf;
      const double diff_ytf=(ytss-ytf)/tauytf;
      const double GtoFast=(1.0-0.7*AF)*0.165*1.0;
      const double Ito=GtoFast*xtf*ytf*(v-EK);

      const double xkurss=1.0/(1.0+exp((v+6.0)/-8.6));
      const double tauxkur=9.0/(1.0+exp((v+5.0)/12.0))+0.5;
      const double ykurss=1.0/(1.0+exp((v+7.5)/10.0));
      const double tauykur=590.0/(1.0+exp((v+60.0)/10.0))+3050.0;
      const double diff_xkur=(xkurss-xkur)/tauxkur;
      const double diff_ykur=(ykurss-ykur)/tauykur;
      const double Gkur=(1.0-0.5*AF)*(1.0+2.0*ISO)*0.045*(1.0+0.2*RA);
      const double IKur=Gkur*xkur*ykur*(v-EK);

      const double aki=1.02/(1+exp(0.2385*(v-EK-59.215)));
      const double bki=(0.49124*exp(0.08032*(v+5.476-EK))+exp(0.06175*(v-EK-594.31)))/(1.0+exp(-0.5143*(v-EK+4.753)));
      const double kiss=aki/(aki+bki);
      const double IK1 =(1.0+1.0*AF)*0.0525*sqrt(Ko/5.4)*kiss*(v-EK);

      const double KdClCa=100.0e-3;
      const double GClCa =0.0548;
      const double IClCa_junc=Fjunc*GClCa/(1.0+KdClCa/Caj)*(v-ECl);
      const double IClCa_sl=Fsl*GClCa/(1.0+KdClCa/Casl)*(v-ECl);
      const double IClCa=IClCa_junc+IClCa_sl;

      const double GClB=9.0e-3;
      const double IClbk=GClB*(v-ECl);

      const double dss=1.0/(1.0+exp(-(v+3.0*ISO+9.0)/6.0));
      const double taud=1.0*dss*(1.0-exp(-(v+3.0*ISO+9.0)/6.0))/(0.035*(v+3.0*ISO+9.0));
      const double fss=1.0/(1.0+exp((v+3.0*ISO+30.0)/7.0))+0.2/(1+exp((50.0-v-3.0*ISO)/20.0));
      const double tauf=1.0/(0.0197*exp(-pow(0.0337*(v+3.0*ISO+25.0),2.0))+0.02);
      const double diff_d=(dss-d)/taud;
      const double diff_f=(fss-f)/tauf;
      const double diff_fcaBj=1.7*Caj*(1.0-fcaBj)-11.9e-3*fcaBj;
      const double diff_fcaBsl=1.7*Casl*(1-fcaBsl)-11.9e-3*fcaBsl;
      const double pNa=(1.0+0.5*ISO)*(1.0-0.5*AF)*7.5e-9;
      const double pCa=(1.0+0.5*ISO)*(1.0-0.5*AF)*2.7e-4;
      const double pK=(1.0+0.5*ISO)*(1.0-0.5*AF)*1.35e-7;
      const double Q10CaL=1.8;
      const double ibarca_j=pCa*4.0*(v*Frdy*FoRT)*(0.341*Caj*exp(2.0*v*FoRT)-0.341*Cao)/(exp(2.0*v*FoRT)-1.0);
      const double ibarca_sl=pCa*4.0*(v*Frdy*FoRT)*(0.341*Casl*exp(2.0*v*FoRT)-0.341*Cao)/(exp(2.0*v*FoRT)-1.0);
      const double ibark=pK*(v*Frdy*FoRT)*(0.75*Ki*exp(v*FoRT)-0.75*Ko)/(exp(v*FoRT)-1.0);
      const double ibarna_j=pNa*(v*Frdy*FoRT) *(0.75*Naj*exp(v*FoRT)-0.75*Nao)/(exp(v*FoRT)-1.0);
      const double ibarna_sl=pNa*(v*Frdy*FoRT) *(0.75*Nasl*exp(v*FoRT)-0.75*Nao)/(exp(v*FoRT)-1.0);
      const double ICa_junc=(Fjunc_CaL*ibarca_j*d*f*(1-fcaBj)*pow(Q10CaL,Qpow))*0.45;
      const double ICa_sl=(Fsl_CaL*ibarca_sl*d*f*(1-fcaBsl)*pow(Q10CaL,Qpow))*0.45;
      const double ICaNa_junc=(Fjunc_CaL*ibarna_j*d*f*(1-fcaBj)*pow(Q10CaL,Qpow))*0.45;
      const double ICaNa_sl=(Fsl_CaL*ibarna_sl*d*f*(1-fcaBsl)*pow(Q10CaL,Qpow))*0.45;
      const double ICaK=(ibark*d*f*(Fjunc_CaL*(1-fcaBj)+Fsl_CaL*(1-fcaBsl))*pow(Q10CaL,Qpow))*0.45;

      const double KmCai=3.59e-3;
      const double KmCao=1.3;
      const double KmNai=12.29;
      const double KmNao=87.5;
      const double ksat=0.27;
      const double nu=0.35;
      const double Kdact =0.384e-3;
      const double Q10NCX=1.57;
      const double Ka_junc=1.0/(1.0+pow(Kdact/Caj,2.0));
      const double Ka_sl=1.0/(1.0+pow(Kdact/Casl,2.0));
      const double s1_junc=exp(nu*v*FoRT)*Naj*Naj*Naj*Cao;
      const double s1_sl=exp(nu*v*FoRT)*Nasl*Nasl*Nasl*Cao;
      const double s2_junc=exp((nu-1.0)*v*FoRT)*Nao*Nao*Nao*Caj;
      const double s3_junc=KmCai*Nao*Nao*Nao*(1+pow(Naj/KmNai,3.0))+KmNao*KmNao*KmNao*Caj*(1.0+Caj/KmCai)+KmCao*Naj*Naj*Naj+Naj*Naj*Naj*Cao+Nao*Nao*Nao*Caj;
      const double s2_sl=exp((nu-1.0)*v*FoRT)*Nao*Nao*Nao*Casl;
      const double s3_sl=KmCai*Nao*Nao*Nao*(1.0+pow(Nasl/KmNai,3.0))+KmNao*KmNao*KmNao*Casl*(1.0+Casl/KmCai)+KmCao*Nasl*Nasl*Nasl+Nasl*Nasl*Nasl*Cao+Nao*Nao*Nao*Casl;
      const double IbarNCX=(1.0+0.4*AF)*3.15;
      const double Incx_junc=Fjunc*IbarNCX*pow(Q10NCX,Qpow)*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1.0+ksat*exp((nu-1.0)*v*FoRT));
      const double Incx_sl=Fsl*IbarNCX*pow(Q10NCX,Qpow)*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1.0+ksat*exp((nu-1.0)*v*FoRT));

      const double Q10SLCaP=2.35;
      const double KmPCa =0.5e-3;
      const double IbarSLCaP= 0.0471;
      const double Ipca_junc=Fjunc*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(Caj,1.6)/(pow(KmPCa,1.6)+pow(Caj,1.6));
      const double Ipca_sl=Fsl*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(Casl,1.6)/(pow(KmPCa,1.6)+pow(Casl,1.6));

      const double GCaB=6.0643e-4;
      const double Icabk_junc=Fjunc*GCaB*(v-ECa_junc);
      const double Icabk_sl=Fsl*GCaB*(v-ECa_sl);

      const double Q10SRCaP=2.6;
      const double Vmax_SRCaP=5.3114e-3;
      const double Kmf=(2.5-1.25*ISO)*0.246e-3;
      const double Kmr=1.7;
      const double hillSRCaP=1.787;
      const double ks=25.0;
      const double koCa=10.0+20.0*AF+10.0*ISO*(1.0-AF);
      const double kom=0.06;
      const double kiCa=0.5;
      const double kim=0.005;
      const double ec50SR=0.45;
      const double MaxSR=15.0;
      const double MinSR=1.0;
      const double kCaSR=MaxSR-(MaxSR-MinSR)/(1.0+pow(ec50SR/Casr,2.5));
      const double koSRCa=koCa/kCaSR;
      const double kiSRCa=kiCa*kCaSR;
      const double RI=1.0-RyRr-RyRo-RyRi;
      const double diff_RyRr=(kim*RI-kiSRCa*Caj*RyRr)-(koSRCa*Caj*Caj*RyRr-kom*RyRo);
      const double diff_RyRo=(koSRCa*Caj*Caj*RyRr-kom*RyRo)-(kiSRCa*Caj*RyRo-kim*RyRi);
      const double diff_RyRi=(kiSRCa*Caj*RyRo-kim*RyRi)-(kom*RyRi-koSRCa*Caj*Caj*RI);
      const double JSRCarel=ks*RyRo*(Casr-Caj);
      const double Jserca=1.0*pow(Q10SRCaP,Qpow)*Vmax_SRCaP*(pow(Cai/Kmf,hillSRCaP)-pow(Casr/Kmr,hillSRCaP))/(1+pow(Cai/Kmf,hillSRCaP)+pow(Casr/Kmr,hillSRCaP));
      const double JSRleak=(1.0+0.25*AF)*5.348e-6*(Casr-Caj);

      const double Bmax_Naj=7.561;
      const double Bmax_Nasl=1.65;
      const double koff_na=1.0e-3;
      const double kon_na=0.1e-3;
      const double diff_NaBj=kon_na*Naj*(Bmax_Naj-NaBj)-koff_na*NaBj;
      const double diff_NaBsl=kon_na*Nasl*(Bmax_Nasl-NaBsl)-koff_na*NaBsl;

      const double Bmax_TnClow=70.0e-3;
      const double koff_tncl=(1.0+0.5*ISO)*19.6e-3;
      const double kon_tncl=32.7;
      const double Bmax_TnChigh=140e-3;
      const double koff_tnchca=0.032e-3;
      const double kon_tnchca=2.37;
      const double koff_tnchmg=3.33e-3;
      const double kon_tnchmg=3.0e-3;
      const double Bmax_CaM=24.0e-3;
      const double koff_cam=238.0e-3;
      const double kon_cam=34.0;
      const double Bmax_myosin=140e-3;
      const double koff_myoca=0.46e-3;
      const double kon_myoca=13.8;
      const double koff_myomg=0.057e-3;
      const double kon_myomg=0.0157;
      const double Bmax_SR=19.0*.9e-3;
      const double koff_sr=60.0e-3;
      const double kon_sr=100.0;
      const double diff_TnCL=kon_tncl*Cai*(Bmax_TnClow-TnCL)-koff_tncl*TnCL;
      const double diff_TnCHc=kon_tnchca*Cai*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchca*TnCHc;
      const double diff_TnCHm=kon_tnchmg*Mgi*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchmg*TnCHm;
      const double diff_CaM=kon_cam*Cai*(Bmax_CaM-CaM)-koff_cam*CaM;
      const double diff_Myc=kon_myoca*Cai*(Bmax_myosin-Myc-Mym)-koff_myoca*Myc;
      const double diff_Mym=kon_myomg*Mgi*(Bmax_myosin-Myc-Mym)-koff_myomg*Mym;
      const double diff_SRB=kon_sr*Cai*(Bmax_SR-SRB)-koff_sr*SRB;
      const double JCaB_cytsol=diff_TnCL+diff_TnCHc+diff_TnCHm+diff_CaM+diff_Myc+diff_Mym+diff_SRB;

      const double Bmax_SLlowsl=37.4e-3*Vmyo/Vsl;
      const double Bmax_SLlowj=4.6e-3*Vmyo/Vjunc*0.1;
      const double koff_sll=1300e-3;
      const double kon_sll=100.0;
      const double Bmax_SLhighsl=13.4e-3*Vmyo/Vsl;
      const double Bmax_SLhighj=1.65e-3*Vmyo/Vjunc*0.1;
      const double koff_slh=30.0e-3;
      const double kon_slh=100.0;
      const double diff_SLLj=kon_sll*Caj*(Bmax_SLlowj-SLLj)-koff_sll*SLLj;
      const double diff_SLLsl=kon_sll*Casl*(Bmax_SLlowsl-SLLsl)-koff_sll*SLLsl;
      const double diff_SLHj=kon_slh*Caj*(Bmax_SLhighj-SLHj)-koff_slh*SLHj;
      const double diff_SLHsl=kon_slh*Casl*(Bmax_SLhighsl-SLHsl)-koff_slh*SLHsl;
      const double JCaB_junc=diff_SLLj+diff_SLHj;
      const double JCaB_sl=diff_SLLsl+diff_SLHsl;

      const double Bmax_Csqn=140.0e-3*Vmyo/Vsr;
      const double koff_csqn=65.0;
      const double kon_csqn=100.0;
      const double KmCsqnb=koff_csqn/kon_csqn;
      const double Buff_Csqnb=1.0/(1.0+Bmax_Csqn*KmCsqnb/pow(KmCsqnb+Casr,2.0));
      const double diff_Csqnb=0;
      const double diff_Casr=Buff_Csqnb*(Jserca-(JSRleak*Vmyo/Vsr+JSRCarel));

      const double INa_tot_junc=INa_junc+INaBk_junc+3.0*Incx_junc+3.0*INAK_junc+ICaNa_junc+INaL_junc;
      const double INa_tot_sl=INa_sl+INaBk_sl+3.0*Incx_sl+3.0*INAK_sl+ICaNa_sl+INaL_sl;
      const double diff_Naj=-INa_tot_junc*Cmem/(Vjunc*Frdy)+Jna_juncsl/Vjunc*(Nasl-Naj)-diff_NaBj;
      const double diff_Nasl=-INa_tot_sl*Cmem/(Vsl*Frdy)+Jna_juncsl/Vsl*(Naj-Nasl)+Jna_slmyo/Vsl*(Nai-Nasl)-diff_NaBsl;
      const double diff_Nai=Jna_slmyo/Vmyo*(Nasl-Nai);

      const double IK_tot=Ito+IKr+IKs+IK1-2.0*INAK+ICaK+IKp+IKur;
      const double diff_Ki=-IK_tot*Cmem/(Vmyo*Frdy);

      const double ICa_tot_junc=ICa_junc+Icabk_junc+Ipca_junc-2.0*Incx_junc;
      const double ICa_tot_sl=ICa_sl+Icabk_sl+Ipca_sl-2.0*Incx_sl;
      const double diff_Caj=-ICa_tot_junc*Cmem/(Vjunc*2.0*Frdy)+Jca_juncsl/Vjunc*(Casl-Caj)-JCaB_junc+(JSRCarel)*Vsr/Vjunc+JSRleak*Vmyo/Vjunc;
      const double diff_Casl=-ICa_tot_sl*Cmem/(Vsl*2.0*Frdy)+Jca_juncsl/Vsl*(Caj-Casl)+Jca_slmyo/Vsl*(Cai-Casl)-JCaB_sl;
      const double diff_Cai=-Jserca*Vsr/Vmyo-JCaB_cytsol+Jca_slmyo/Vmyo*(Casl-Cai);

      const double INa_tot=INa_tot_junc+INa_tot_sl;
      const double ICl_tot=IClCa+IClbk;
      const double ICa_tot=ICa_tot_junc+ICa_tot_sl;
      const double Itot=INa_tot+ICl_tot+ICa_tot+IK_tot;

      const double diff_v = -Itot;
      dVm[ii] = diff_v;

      if (1) 
      {
         bool foundError=false;
#define CHECK_BLOWUP(x) do { if (!isfinite(x)) { fprintf(stderr, "Error in node %d, variable " #x " = %g\n", ii, (x)); foundError=true; } } while(0)
         CHECK_BLOWUP(diff_v);
            
         //EDIT_STATE
         CHECK_BLOWUP(diff_m);
         CHECK_BLOWUP(diff_h);
         CHECK_BLOWUP(diff_j);
         CHECK_BLOWUP(diff_mL);
         CHECK_BLOWUP(diff_hL);
         CHECK_BLOWUP(diff_xtf);
         CHECK_BLOWUP(diff_ytf);
         CHECK_BLOWUP(diff_xkr);
         CHECK_BLOWUP(diff_xks);
         CHECK_BLOWUP(diff_xkur);
         CHECK_BLOWUP(diff_ykur);
         CHECK_BLOWUP(diff_d);
         CHECK_BLOWUP(diff_f);
         CHECK_BLOWUP(diff_fcaBj);
         CHECK_BLOWUP(diff_fcaBsl);
         CHECK_BLOWUP(diff_RyRr);
         CHECK_BLOWUP(diff_RyRo);
         CHECK_BLOWUP(diff_RyRi);
         CHECK_BLOWUP(diff_NaBj);
         CHECK_BLOWUP(diff_NaBsl);
         CHECK_BLOWUP(diff_TnCL);
         CHECK_BLOWUP(diff_TnCHc);
         CHECK_BLOWUP(diff_TnCHm);
         CHECK_BLOWUP(diff_CaM);
         CHECK_BLOWUP(diff_Myc);
         CHECK_BLOWUP(diff_Mym);
         CHECK_BLOWUP(diff_SRB);
         CHECK_BLOWUP(diff_SLLj);
         CHECK_BLOWUP(diff_SLLsl);
         CHECK_BLOWUP(diff_SLHj);
         CHECK_BLOWUP(diff_SLHsl);
         CHECK_BLOWUP(diff_Csqnb);
         CHECK_BLOWUP(diff_Naj);
         CHECK_BLOWUP(diff_Nasl);
         CHECK_BLOWUP(diff_Nai);
         CHECK_BLOWUP(diff_Ki);
         CHECK_BLOWUP(diff_Casr);
         CHECK_BLOWUP(diff_Caj);
         CHECK_BLOWUP(diff_Casl);
         CHECK_BLOWUP(diff_Cai);


#undef CHECK_BLOWUP
         
         if (foundError) 
         {
#define PRINT_STATE(x) do { fprintf(stderr, "node %d: " #x " = %g\n", ii, (x)); } while(0)
            //EDIT_STATE
            PRINT_STATE(m);
            PRINT_STATE(h);
            PRINT_STATE(j);
            PRINT_STATE(mL);
            PRINT_STATE(hL);
            PRINT_STATE(xtf);
            PRINT_STATE(ytf);
            PRINT_STATE(xkr);
            PRINT_STATE(xks);
            PRINT_STATE(xkur);
            PRINT_STATE(ykur);
            PRINT_STATE(d);
            PRINT_STATE(f);
            PRINT_STATE(fcaBj);
            PRINT_STATE(fcaBsl);
            PRINT_STATE(RyRr);
            PRINT_STATE(RyRo);
            PRINT_STATE(RyRi);
            PRINT_STATE(NaBj);
            PRINT_STATE(NaBsl);
            PRINT_STATE(TnCL);
            PRINT_STATE(TnCHc);
            PRINT_STATE(TnCHm);
            PRINT_STATE(CaM);
            PRINT_STATE(Myc);
            PRINT_STATE(Mym);
            PRINT_STATE(SRB);
            PRINT_STATE(SLLj);
            PRINT_STATE(SLLsl);
            PRINT_STATE(SLHj);
            PRINT_STATE(SLHsl);
            PRINT_STATE(Csqnb);
            PRINT_STATE(Naj);
            PRINT_STATE(Nasl);
            PRINT_STATE(Nai);
            PRINT_STATE(Ki);
            PRINT_STATE(Casr);
            PRINT_STATE(Caj);
            PRINT_STATE(Casl);
            PRINT_STATE(Cai);
            
#undef PRINT_STATE
            
            exit(255);
         }
      }
      
      
      //EDIT_STATE
      state_[ii].m = mss-(mss-m)*exp(-dt/taum);
      state_[ii].h = hss-(hss-h)*exp(-dt/tauh);
      state_[ii].j = jss-(jss-j)*exp(-dt/tauj);
      state_[ii].mL = mLss-(mLss-mL)*exp(-dt/taumL);
      state_[ii].hL = hlinf-(hlinf-hL)*exp(-dt/tauhl);
      state_[ii].xtf = xtss-(xtss-xtf)*exp(-dt/tauxtf);
      state_[ii].ytf = ytss-(ytss-ytf)*exp(-dt/tauytf);
      state_[ii].xkr = xrss-(xrss-xkr)*exp(-dt/tauxr);
      state_[ii].xks = xsss-(xsss-xks)*exp(-dt/tauxs);
      state_[ii].xkur = xkurss-(xkurss-xkur)*exp(-dt/tauxkur);
      state_[ii].ykur = ykurss-(ykurss-ykur)*exp(-dt/tauykur);
      state_[ii].d = dss-(dss-d)*exp(-dt/taud);
      state_[ii].f = fss-(fss-f)*exp(-dt/tauf);
      state_[ii].fcaBj += diff_fcaBj*dt;
      state_[ii].fcaBsl += diff_fcaBsl*dt;
      state_[ii].RyRr += diff_RyRr*dt;
      state_[ii].RyRo += diff_RyRo*dt;
      state_[ii].RyRi += diff_RyRi*dt;
      state_[ii].NaBj += diff_NaBj*dt;
      state_[ii].NaBsl += diff_NaBsl*dt;
      state_[ii].TnCL += diff_TnCL*dt;
      state_[ii].TnCHc += diff_TnCHc*dt;
      state_[ii].TnCHm += diff_TnCHm*dt;
      state_[ii].CaM += diff_CaM*dt;
      state_[ii].Myc += diff_Myc*dt;
      state_[ii].Mym += diff_Mym*dt;
      state_[ii].SRB += diff_SRB*dt;
      state_[ii].SLLj += diff_SLLj*dt;
      state_[ii].SLLsl += diff_SLLsl*dt;
      state_[ii].SLHj += diff_SLHj*dt;
      state_[ii].SLHsl += diff_SLHsl*dt;
      state_[ii].Csqnb += diff_Csqnb*dt;
      state_[ii].Naj += diff_Naj*dt;
      state_[ii].Nasl += diff_Nasl*dt;
      state_[ii].Nai += diff_Nai*dt;
      state_[ii].Ki += diff_Ki*dt;
      state_[ii].Casr += diff_Casr*dt;
      state_[ii].Caj += diff_Caj*dt;
      state_[ii].Casl += diff_Casl*dt;
      state_[ii].Cai += diff_Cai*dt;
      
   }
}

void ThisReaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= nCells_);
   Vm.assign(Vm.size(), -87.84);
   State initState;
   //EDIT_STATE
   initState.m=0.0;
   initState.h=1.0;
   initState.j=1.0;
   initState.mL=0.0;
   initState.hL=1.0;
   initState.xtf=0.0;
   initState.ytf=1.0;
   initState.xkr=0.0;
   initState.xks= 0.0;
   initState.xkur=0.0;
   initState.ykur=1.0;
   initState.d=0.0;
   initState.f=1.0;
   initState.fcaBj=0.025;
   initState.fcaBsl=0.015;
   initState.RyRr=1.0;
   initState.RyRo=0.0;
   initState.RyRi=0.0;
   initState.NaBj=3.5;
   initState.NaBsl=0.8;
   initState.TnCL=0.01;
   initState.TnCHc=0.1;
   initState.TnCHm=0.01;
   initState.CaM=3.0e-4;
   initState.Myc=1.3e-3;
   initState.Mym=0.14;
   initState.SRB=2.0e-3;
   initState.SLLj=0.01;
   initState.SLLsl=0.1;
   initState.SLHj=7.3e-3;
   initState.SLHsl=7.3e-2;
   initState.Csqnb=1.25;
   initState.Naj=9.136;
   initState.Nasl=9.136;
   initState.Nai=9.136;
   initState.Ki=120.0;
   initState.Casr=0.01;
   initState.Caj=1.7e-4;
   initState.Casl=1.0e-4;
   initState.Cai=1.0e-4;

   state_.resize(nCells_);
   state_.assign(state_.size(), initState);

}

const string ThisReaction::getUnit(const std::string& varName) const
{
   //deliberatly broken for now, if this code still is being used past 2016-11-01 something has gone wrong.
   return "1";
}

   
#define HANDLE_OFFSET 1000
int ThisReaction::getVarHandle(const std::string& varName) const
{
   return getVarOffset(varName)+HANDLE_OFFSET;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
   reinterpret_cast<double*>(&state_[iCell])[varHandle-HANDLE_OFFSET] = value;
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
   return reinterpret_cast<const double*>(&state_[iCell])[varHandle-HANDLE_OFFSET];
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.resize(NUMVARS);
   fieldUnits.resize(NUMVARS);
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      fieldNames[ivar] = varNames[ivar];
      fieldUnits[ivar] = getUnit(fieldNames[ivar]);
   }
}

}
