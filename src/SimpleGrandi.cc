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
      setDefault(RA, 0);

      //EDIT_PARAMETERS
      setDefault(ks, 25); // [1/ms]      
      setDefault(Vmax_SRCaP, 5.3114e-3);  // [mM/msec] (286 umol/L cytosol/sec)
      setDefault(GCaB, 6.0643e-4);    // [uA/uF] 3
      setDefault(GClCa, 0.0548);   // [mS/uF]
      setDefault(GClB, 9e-3);        // [mS/uF]
      setDefault(gkp, 0.002);
      setDefault(IbarNaK, 1.26);     // [uA/uF]
      setDefault(IbarSLCaP, 0.0471); // IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
      
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

ThisReaction::ThisReaction(const Anatomy& anatomy)
: nCells_(anatomy.nLocal())
{
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
      const double d= state_[ii].d;
      const double f= state_[ii].f;
      const double fcaBj= state_[ii].fcaBj;
      const double fcaBsl= state_[ii].fcaBsl;
      const double xtof= state_[ii].xtof;
      const double ytof= state_[ii].ytof;
      const double xkr= state_[ii].xkr;
      const double xks= state_[ii].xks;
      const double xkur= state_[ii].xkur;
      const double ykur= state_[ii].ykur;
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
      const double R = 8314;       // [J/kmol*K]  
      const double Frdy = 96485;   // [C/mol]  
      const double Temp = 310;     // [K]
      const double FoRT = Frdy/R/Temp;
      const double Cmem = 1.1e-10;   // [F] membrane capacitance 1.3810e-10;//
      const double Qpow = (Temp-310)/10;

      // Cell geometry
      const double cellLength = 100;     // cell length [um]113;//100
      const double cellRadius = 10.25;   // cell radius [um]12;//10.25
      const double Vcell = pi*cellRadius*cellRadius*cellLength*1e-15;    // [L]
      const double Vmyo = 0.65*Vcell;
      const double Vsr = 0.035*Vcell;
      const double Vsl = 0.02*Vcell;
      const double Vjunc = 1*0.0539*.01*Vcell; 
      const double J_ca_juncsl =1/1.2134e12; // [L/msec] = 8.2413e-13
      const double J_ca_slmyo = 1/2.68510e11; // [L/msec] = 3.2743e-12
      const double J_na_juncsl = 1/(1.6382e12/3*100); // [L/msec] = 6.1043e-13
      const double J_na_slmyo = 1/(1.8308e10/3*100);  // [L/msec] = 5.4621e-11

      // Fractional currents in compartments
      const double Fjunc = 0.11;
      const double Fsl = 1-Fjunc;
      const double Fjunc_CaL = 0.9;
      const double Fsl_CaL = 1-Fjunc_CaL;

      // Fixed ion concentrations     
      const double Cli = 15;   // Intracellular Cl  [mM]
      const double Clo = 150;  // Extracellular Cl  [mM]
      const double Ko = 5.4;   // Extracellular K   [mM]
      const double Nao = 140;  // Extracellular Na  [mM]
      const double Cao = 1.8;  // Extracellular Ca  [mM]
      const double Mgi = 1;    // Intracellular Mg  [mM]

      // Nernst Potentials
      const double ena_junc = (1/FoRT)*log(Nao/Naj);     // [mV]
      const double ena_sl = (1/FoRT)*log(Nao/Nasl);       // [mV]
      const double ek = (1/FoRT)*log(Ko/Ki);	        // [mV]
      const double eca_junc = (1/FoRT/2)*log(Cao/Caj);   // [mV]
      const double eca_sl = (1/FoRT/2)*log(Cao/Casl);     // [mV]
      const double ecl = (1/FoRT)*log(Cli/Clo);            // [mV]

      // Na transport parameters
      const double GNa=23*(1-0.1*AF);  // [mS/uF]
      const double GNaB = 0.597e-3;    // [mS/uF] 
      const double KmNaip = 11*(1-0.25*ISO);         // [mM]11
      const double KmKo =1.5;         // [mM]1.5

      //// K current parameters
      const double pNaK = 0.01833;      

      // Cl current parameters
      const double KdClCa = 100e-3;    // [mM]

      // I_Ca parameters
      const double pNa = (1+0.5*ISO)*(1-0.5*AF)*0.75e-8;       // [cm/sec]
      const double pCa = (1+0.5*ISO)*(1-0.5*AF)*2.7e-4;       // [cm/sec]
      const double pK = (1+0.5*ISO)*(1-0.5*AF)*1.35e-7;        // [cm/sec]
      const double Q10CaL = 1.8;       

      //// Ca transport parameters
      const double IbarNCX = (1+0.4*AF)*3.15;      // [uA/uF]5.5 before - 9 in rabbit
      const double KmCai = 3.59e-3;    // [mM]
      const double KmCao = 1.3;        // [mM]
      const double KmNai = 12.29;      // [mM]
      const double KmNao = 87.5;       // [mM]
      const double ksat = 0.27;        // [none]  
      const double nu = 0.35;          // [none]
      const double Kdact =0.384e-3;   // [mM] 0.256 rabbit
      const double Q10NCX = 1.57;      // [none]
      const double KmPCa =0.5e-3;     // [mM] 
      const double Q10SLCaP = 2.35;    // [none]

      // SR flux parameters
      const double Q10SRCaP = 2.6;          // [none]
      const double Kmf = (2.5-1.25*ISO)*0.246e-3;          // [mM] default
      const double Kmr = 1.7;               // [mM]L cytosol
      const double hillSRCaP = 1.787;       // [mM]
      const double koCa = 10+20*AF+10*ISO*(1-AF);               // [mM^-2 1/ms]   //default 10   modified 20
      const double kom = 0.06;              // [1/ms]     
      const double kiCa = 0.5;              // [1/mM/ms]
      const double kim = 0.005;             // [1/ms]
      const double ec50SR = 0.45;           // [mM]

      // Buffering parameters
      // koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
      const double Bmax_Naj = 7.561;       // [mM] // Na buffering
      const double Bmax_Nasl = 1.65;       // [mM]
      const double koff_na = 1e-3;         // [1/ms]
      const double kon_na = 0.1e-3;        // [1/mM/ms]
      const double Bmax_TnClow = 70e-3;    // [mM]                      // TnC low affinity
      const double koff_tncl = (1+0.5*ISO)*19.6e-3;    // [1/ms] 
      const double kon_tncl = 32.7;        // [1/mM/ms]
      const double Bmax_TnChigh = 140e-3;  // [mM]                      // TnC high affinity 
      const double koff_tnchca = 0.032e-3; // [1/ms] 
      const double kon_tnchca = 2.37;      // [1/mM/ms]
      const double koff_tnchmg = 3.33e-3;  // [1/ms] 
      const double kon_tnchmg = 3e-3;      // [1/mM/ms]
      const double Bmax_CaM = 24e-3;       // [mM] **? about setting to 0 in c-code**   // CaM buffering
      const double koff_cam = 238e-3;      // [1/ms] 
      const double kon_cam = 34;           // [1/mM/ms]
      const double Bmax_myosin = 140e-3;   // [mM]                      // Myosin buffering
      const double koff_myoca = 0.46e-3;   // [1/ms]
      const double kon_myoca = 13.8;       // [1/mM/ms]
      const double koff_myomg = 0.057e-3;  // [1/ms]
      const double kon_myomg = 0.0157;     // [1/mM/ms]
      const double Bmax_SR = 19*.9e-3;     // [mM] (Bers text says 47e-3) 19e-3
      const double koff_sr = 60e-3;        // [1/ms]
      const double kon_sr = 100;           // [1/mM/ms]
      const double Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        // [mM]    // SL buffering
      const double Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    // [mM]    //Fei *0.1!!! junction reduction factor
      const double koff_sll = 1300e-3;     // [1/ms]
      const double kon_sll = 100;          // [1/mM/ms]
      const double Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       // [mM] 
      const double Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  // [mM] //Fei *0.1!!! junction reduction factor
      const double koff_slh = 30e-3;       // [1/ms]
      const double kon_slh = 100;          // [1/mM/ms]
      const double Bmax_Csqn = 140e-3*Vmyo/Vsr;            // [mM] // Bmax_Csqn = 2.6;      // Csqn buffering
      const double koff_csqn = 65;         // [1/ms] 
      const double kon_csqn = 100;         // [1/mM/ms] 


      //// Membrane Currents
      const double mss = 1 / pow(1 + exp( -(56.86 + v) / 9.03 ),2);
      const double taum = 0.1292 * exp(-pow((v+45.79)/15.54,2)) + 0.06487 * exp(-pow((v-4.823)/51.12,2)); 
      const double ah = (v >= -40) * (0)
         + (v < -40) * (0.057 * exp( -(v + 80) / 6.8 )); 
      const double bh = (v >= -40) * (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 )))) 
         + (v < -40) * ((2.7 * exp( 0.079 * v) + 3.1*1e5 * exp(0.3485 * v))); 
      const double tauh = 1 / (ah + bh); 
      const double hss = 1 / (pow(1 + exp( (v + 71.55)/7.43 ),2));
 
      const double aj = (v >= -40) * (0) 
         +(v < -40) * (((-2.5428 * 1e4*exp(0.2444*v) - 6.948*1e-6 * exp(-0.04391*v)) * (v + 37.78)) / 
                       (1 + exp( 0.311 * (v + 79.23) )));
      const double bj = (v >= -40) * ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) ))) 
         + (v < -40) * ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) ))); 
      const double tauj = 1 / (aj + bj);
      const double jss = 1 / (pow(1 + exp( (v + 71.55)/7.43 ),2));         
 
      const double diff_m = (mss - m) / taum;
      const double diff_h = (hss - h) / tauh;
      const double diff_j = (jss - j) / tauj;
    
      const double I_Na_junc = Fjunc*GNa*pow(m,3)*h*j*(v-ena_junc);
      const double I_Na_sl = Fsl*GNa*pow(m,3)*h*j*(v-ena_sl);
      const double I_Na = I_Na_junc+I_Na_sl;


      // Late I_Na
      const double GNaL=0.0025*AF;
      const double aml = 0.32*(v+47.13)/(1-exp(-0.1*(v+47.13)));
      const double bml = 0.08*exp(-v/11);
      const double hlinf = 1/(1+exp((v+91)/6.1));
      const double tauhl=600;
      const double diff_mL = aml*(1-mL)-bml*mL;
      const double diff_hL = (hlinf-hL)/tauhl;      
      const double hLss = hlinf;
      const double tauhL = tauhl;
      const double taumL=1.0/(aml+bml);
      const double mLss=aml*taumL;
      
      const double I_NaL_junc = Fjunc*GNaL*pow(mL,3)*hL*(v-ena_junc);
      const double I_NaL_sl = Fsl*GNaL*pow(mL,3)*hL*(v-ena_sl);
      const double I_NaL = I_NaL_junc + I_NaL_sl;

      // I_nabk: Na Background Current
      const double I_nabk_junc = Fjunc*GNaB*(v-ena_junc);
      const double I_nabk_sl = Fsl*GNaB*(v-ena_sl);
      const double I_nabk = I_nabk_junc+I_nabk_sl;

      // I_nak: Na/K Pump Current
      const double sigma = (exp(Nao/67.3)-1)/7;
      const double fnak = 1/(1+0.1245*exp(-0.1*v*FoRT)+0.0365*sigma*exp(-v*FoRT));
      const double I_nak_junc = 1*Fjunc*IbarNaK*fnak*Ko /(1+pow(KmNaip/Naj,4)) /(Ko+KmKo);
      const double I_nak_sl = 1*Fsl*IbarNaK*fnak*Ko /(1+pow(KmNaip/Nasl,4)) /(Ko+KmKo);
      const double I_nak = I_nak_junc+I_nak_sl;

      //// I_kr: Rapidly Activating K Current
      const double gkr =0.035*sqrt(Ko/5.4);
      const double xrss = 1/(1+exp(-(v+10)/5));
      const double tauxr = 550/(1+exp((-22-v)/9))*6/(1+exp((v-(-11))/9))+230/(1+exp((v-(-40))/20));
      const double diff_xkr = (xrss-xkr)/tauxr;
      const double rkr = 1/(1+exp((v+74)/24));
      const double I_kr = gkr*xkr*rkr*(v-ek);
      const double xkrss = xrss;
      const double tauxkr = tauxr;

      //// I_ks: Slowly Activating K Current
      const double eks = (1/FoRT)*log((Ko+pNaK*Nao)/(Ki+pNaK*Nai));
      const double gks_junc=1*(1+1*AF+2*ISO)*0.0035*1;
      const double gks_sl=1*(1+1*AF+2*ISO)*0.0035*1; //FRA
      const double xsss = 1 / (1+exp(-(v+40*ISO + 3.8)/14.25)); // fitting Fra
      const double tauxs=990.1/(1+exp(-(v+40*ISO+2.436)/14.12));
      const double diff_xks = (xsss-xks)/tauxs;
      const double I_ks_junc = Fjunc*gks_junc*pow(xks,2)*(v-eks);
      const double I_ks_sl = Fsl*gks_sl*pow(xks,2)*(v-eks);
      const double I_ks = I_ks_junc+I_ks_sl;
      const double xksss = xsss;
      const double tauxks = tauxs;
      
      //I_kp: Plateau K current
      const double kp_kp = 1/(1+exp(7.488-v/5.98));
      const double I_kp_junc = Fjunc*gkp*kp_kp*(v-ek);
      const double I_kp_sl = Fsl*gkp*kp_kp*(v-ek);
      const double I_kp = I_kp_junc+I_kp_sl;

      //// I_to: Transient Outward K Current (slow and fast components)
      // modified for human myocytes

      const double GtoFast=(1.0-0.7*AF)*0.165*1.0; //nS/pF maleckar; //human atrium

      //11/12/09; changed Itof to that from maleckar/giles/2009; removed I_tos
      //atrium
      //equations for activation; 
      const double xtoss = ( (1)/ ( 1 + exp( -(v+1.0)/11.0 ) ) );
      const double tauxtof = 3.5*exp(-pow(v/30.0,2))+1.5;
      const double xtofss = xtoss;
      
      //equations for inactivation;
      const double ytoss = ( (1.0)/ ( 1 + exp( (v+40.5)/11.5) ) ) ;
      const double tauytof =25.635*exp(-(pow((v+52.45)/15.8827,2)))+24.14;//14.14
      const double ytofss = ytoss;
      
      const double diff_xtof = (xtoss-xtof)/tauxtof;
      const double diff_ytof = (ytoss-ytof)/tauytof;
      const double I_tof = 1.0*GtoFast*xtof*ytof*(v-ek);
      const double I_to = 1*I_tof;

      //// I_kur: Ultra rapid delayed rectifier Outward K Current
      //Equation for IKur; from Maleckar et al. 2009 - EG
      //atrium
      //equations for activation;
      const double Gkur = 1*(1.0-0.5*AF)*(1+2*ISO)* 0.045*(1+0.2*RA); //nS/pF maleckar 0.045
      const double xkurss = ( (1)/ ( 1 + exp( (v+6)/-8.6 ) ) );
      const double tauxkur = 9/(1+exp((v+5)/12.0))+0.5;

      //equations for inactivation;
      const double ykurss = ( (1)/ ( 1 + exp( (v+7.5)/10 ) ) );
      const double tauykur = 590/(1+exp((v+60)/10.0))+3050;
      const double diff_xkur = (xkurss-xkur)/tauxkur;
      const double diff_ykur = (ykurss-ykur)/tauykur;
      const double I_kur = 1*Gkur*xkur*ykur*(v-ek);

      //// I_ki: Time-Independent K Current
      const double aki = 1.02/(1+exp(0.2385*(v-ek-59.215)));
      const double bki =(0.49124*exp(0.08032*(v+5.476-ek)) + exp(0.06175*(v-ek-594.31))) /(1 + exp(-0.5143*(v-ek+4.753)));
      const double kiss = aki/(aki+bki);

      //I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*(v-ek);
      //SVP 11/11/09
      //multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
      //resting potential
      const double I_ki =(1+1*AF)*0.0525*sqrt(Ko/5.4)*kiss*(v-ek);

      // I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
      const double I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/Caj)*(v-ecl);
      const double I_ClCa_sl = Fsl*GClCa/(1+KdClCa/Casl)*(v-ecl);
      const double I_ClCa = I_ClCa_junc+I_ClCa_sl;
      const double I_Clbk = GClB*(v-ecl);

      const double GClCFTR=0;//4.9e-3*ISO;     // [mS/uF]
      const double I_ClCFTR = GClCFTR*(v-ecl);

      //// I_Ca: L-type Calcium Current
      const double dss = 1/(1+exp(-(v+3*ISO+9)/6)); //in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
      const double taud = 1*dss*(1-exp(-(v+3*ISO+9)/6))/(0.035*(v+3*ISO+9)); 
      const double fss = 1/(1+exp((v+3*ISO+30)/7))+0.2/(1+exp((50-v-3*ISO)/20)); // in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
      const double tauf = 1/(0.0197*exp( -pow(0.0337*(v+3*ISO+25),2) )+0.02);
      const double diff_d = (dss-d)/taud;
      const double diff_f = (fss-f)/tauf;
      const double diff_fcaBj = 1.7*Caj*(1-fcaBj)-1*11.9e-3*fcaBj; // fCa_junc   koff!!!!!!!!
      const double diff_fcaBsl = 1.7*Casl*(1-fcaBsl)-1*11.9e-3*fcaBsl; // fCa_sl
      // fcaCaMSL= 0.1/(1+(0.01/Casl));
      // fcaCaj= 0.1/(1+(0.01/Caj));
      const double fcaCaMSL=0;
      const double fcaCaj= 0;
      const double ibarca_j = pCa*4*(v*Frdy*FoRT) * (0.341*Caj*exp(2*v*FoRT)-0.341*Cao) /(exp(2*v*FoRT)-1);
      const double ibarca_sl = pCa*4*(v*Frdy*FoRT) * (0.341*Casl*exp(2*v*FoRT)-0.341*Cao) /(exp(2*v*FoRT)-1);
      const double ibark = pK*(v*Frdy*FoRT)*(0.75*Ki*exp(v*FoRT)-0.75*Ko) /(exp(v*FoRT)-1);
      const double ibarna_j = pNa*(v*Frdy*FoRT) *(0.75*Naj*exp(v*FoRT)-0.75*Nao)  /(exp(v*FoRT)-1);
      const double ibarna_sl = pNa*(v*Frdy*FoRT) *(0.75*Nasl*exp(v*FoRT)-0.75*Nao)  /(exp(v*FoRT)-1);
      const double I_Ca_junc = (Fjunc_CaL*ibarca_j*d*f*((1-fcaBj)+fcaCaj)*pow(Q10CaL,Qpow))*0.45;
      const double I_Ca_sl = (Fsl_CaL*ibarca_sl*d*f*((1-fcaBsl)+fcaCaMSL)*pow(Q10CaL,Qpow))*0.45;
      const double I_Ca = I_Ca_junc+I_Ca_sl;
      const double I_CaK = (ibark*d*f*(Fjunc_CaL*(fcaCaj+(1-fcaBj))+Fsl_CaL*(fcaCaMSL+(1-fcaBsl)))*pow(Q10CaL,Qpow))*0.45;
      const double I_CaNa_junc = (Fjunc_CaL*ibarna_j*d*f*((1-fcaBj)+fcaCaj)*pow(Q10CaL,Qpow))*0.45;
      const double I_CaNa_sl = (Fsl_CaL*ibarna_sl*d*f*((1-fcaBsl)+fcaCaMSL)*pow(Q10CaL,Qpow))*.45;
      const double I_CaNa = I_CaNa_junc+I_CaNa_sl;
      const double I_Catot = I_Ca+I_CaK+I_CaNa;

      // I_ncx: Na/Ca Exchanger flux
      const double Ka_junc = 1/(1+pow(Kdact/Caj,2));
      const double Ka_sl = 1/(1+pow(Kdact/Casl,2));
      const double s1_junc = exp(nu*v*FoRT)*pow(Naj,3)*Cao;
      const double s1_sl = exp(nu*v*FoRT)*pow(Nasl,3)*Cao;
      const double s2_junc = exp((nu-1)*v*FoRT)*pow(Nao,3)*Caj;
      const double s3_junc = KmCai*pow(Nao,3)*(1+pow(Naj/KmNai,3)) + pow(KmNao,3)*Caj*(1+Caj/KmCai)+KmCao*pow(Naj,3)+pow(Naj,3)*Cao+pow(Nao,3)*Caj;
      const double s2_sl = exp((nu-1)*v*FoRT)*pow(Nao,3)*Casl;
      const double s3_sl = KmCai*pow(Nao,3)*(1+pow(Nasl/KmNai,3)) + pow(KmNao,3)*Casl*(1+Casl/KmCai)+KmCao*pow(Nasl,3)+pow(Nasl,3)*Cao+pow(Nao,3)*Casl;


      const double I_ncx_junc = Fjunc*IbarNCX*pow(Q10NCX,Qpow)*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*v*FoRT));
      const double I_ncx_sl = Fsl*IbarNCX*pow(Q10NCX,Qpow)*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*v*FoRT));
      const double I_ncx = I_ncx_junc+I_ncx_sl;

      // I_pca: Sarcolemmal Ca Pump Current
      const double I_pca_junc = Fjunc*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(Caj,1.6)/(pow(KmPCa,1.6)+pow(Caj,1.6));
      const double I_pca_sl = Fsl*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(Casl,1.6)/(pow(KmPCa,1.6)+pow(Casl,1.6));
      const double I_pca = I_pca_junc+I_pca_sl;

      // I_cabk: Ca Background Current
      const double I_cabk_junc = Fjunc*GCaB*(v-eca_junc);
      const double I_cabk_sl = Fsl*GCaB*(v-eca_sl);
      const double I_cabk = I_cabk_junc+I_cabk_sl;

      //// SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
      const double MaxSR = 15;
      const double MinSR = 1;
      const double kCaSR = MaxSR - (MaxSR-MinSR)/(1+pow(ec50SR/Casr,2.5));
      const double koSRCa = (1)*koCa/kCaSR;//
      const double kiSRCa = kiCa*kCaSR;
      const double RI = 1-RyRr-RyRo-RyRi;
      const double diff_RyRr = (kim*RI-kiSRCa*Caj*RyRr)-(koSRCa*pow(Caj,2)*RyRr-kom*RyRo);   // R
      const double diff_RyRo = (koSRCa*pow(Caj,2)*RyRr-kom*RyRo)-(kiSRCa*Caj*RyRo-kim*RyRi);// O
      const double diff_RyRi = (kiSRCa*Caj*RyRo-kim*RyRi)-(kom*RyRi-koSRCa*pow(Caj,2)*RI);   // I
      const double J_SRCarel = ks*RyRo*(Casr-Caj);          // [mM/ms]

      const double J_serca = 1.0*pow(Q10SRCaP,Qpow)*Vmax_SRCaP*(pow(Cai/Kmf,hillSRCaP)-pow(Casr/Kmr,hillSRCaP))
         /(1+pow(Cai/Kmf,hillSRCaP)+pow(Casr/Kmr,hillSRCaP));
      const double J_SRleak = (1)*(1.0+0.25*AF)*5.348e-6*(Casr-Caj);           //   [mM/ms]


      //// Sodium and Calcium Buffering
      const double diff_NaBj = kon_na*Naj*(Bmax_Naj-NaBj)-koff_na*NaBj;        // NaBj      [mM/ms]
      const double diff_NaBsl = kon_na*Nasl*(Bmax_Nasl-NaBsl)-koff_na*NaBsl;       // NaBsl     [mM/ms]

      // Cytosolic Ca Buffers
      const double diff_TnCL = kon_tncl*Cai*(Bmax_TnClow-TnCL)-koff_tncl*TnCL;            // TnCL      [mM/ms]
      const double diff_TnCHc = kon_tnchca*Cai*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchca*TnCHc; // TnCHc     [mM/ms]
      const double diff_TnCHm = kon_tnchmg*Mgi*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchmg*TnCHm;   // TnCHm     [mM/ms]
      const double diff_CaM = kon_cam*Cai*(Bmax_CaM-CaM)-koff_cam*CaM;                 // CaM       [mM/ms]
      const double diff_Myc = kon_myoca*Cai*(Bmax_myosin-Myc-Mym)-koff_myoca*Myc;    // Myosin_ca [mM/ms]
      const double diff_Mym = kon_myomg*Mgi*(Bmax_myosin-Myc-Mym)-koff_myomg*Mym;      // Myosin_mg [mM/ms]
      const double diff_SRB = kon_sr*Cai*(Bmax_SR-SRB)-koff_sr*SRB;                    // SRB       [mM/ms]
      const double J_CaB_cytosol=diff_TnCL+diff_TnCHc+diff_TnCHm+diff_CaM+diff_Myc+diff_Mym+diff_SRB;

      // Junctional and SL Ca Buffers
      const double diff_SLLj = kon_sll*Caj*(Bmax_SLlowj-SLLj)-koff_sll*SLLj;       // SLLj      [mM/ms]
      const double diff_SLLsl = kon_sll*Casl*(Bmax_SLlowsl-SLLsl)-koff_sll*SLLsl;      // SLLsl     [mM/ms]
      const double diff_SLHj = kon_slh*Caj*(Bmax_SLhighj-SLHj)-koff_slh*SLHj;      // SLHj      [mM/ms]
      const double diff_SLHsl = kon_slh*Casl*(Bmax_SLhighsl-SLHsl)-koff_slh*SLHsl;     // SLHsl     [mM/ms]
      const double J_CaB_junction = diff_SLLj+diff_SLHj;
      const double J_CaB_sl = diff_SLLsl+diff_SLHsl;

      //// Ion concentrations
      // SR Ca Concentrations
      const double diff_Csqnb = kon_csqn*Casr*(Bmax_Csqn-Csqnb)-koff_csqn*Csqnb;       // Csqn      [mM/ms]
      const double diff_Casr = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-diff_Csqnb;         // Ca_sr     [mM/ms] //Ratio 3 leak current
      // diff_Casr=0;

      // Sodium Concentrations
      const double I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc+I_NaL_junc;   // [uA/uF]
      const double I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl+I_NaL_sl;   // [uA/uF]
      const double I_Na_tot_sl2 = 3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   // [uA/uF]
      const double I_Na_tot_junc2 = 3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   // [uA/uF]

      const double diff_Naj = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(Nasl-Naj)-diff_NaBj;
      const double diff_Nasl = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(Naj-Nasl)
         +J_na_slmyo/Vsl*(Nai-Nasl)-diff_NaBsl;
      //FluxNaSL=diff_Nasl;
      // diff_Naj = 0;
      // diff_Nasl = 0;
      const double diff_Nai = J_na_slmyo/Vmyo*(Nasl-Nai);             // [mM/msec] 
      // diff_Nai=0;
      
      // Potassium Concentration
      const double I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur;     // [uA/uF] //SVP: added IKur
      // diff_Ki = 0; //-I_K_tot*Cmem/(Vmyo*Frdy);           // [mM/msec]
      const double diff_Ki =0; // -I_K_tot*Cmem/(Vmyo*Frdy);

      // Calcium Concentrations
      const double I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   // [uA/uF]
      const double I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            // [uA/uF]
      const double diff_Caj = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(Casl-Caj)
         -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  // Ca_j
      const double diff_Casl = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(Caj-Casl)
         + J_ca_slmyo/Vsl*(Cai-Casl)-J_CaB_sl;   // Ca_sl
      // diff_Caj=0;
      // diff_Casl=0;
      // diff_Cai = -J_serca*Vsr/Vmyo-J_CaB_cytosol;//+J_ca_slmyo/Vmyo*(Casl-Cai);    // [mM/msec]
      const double diff_Cai = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(Casl-Cai);
      // diff_Cai=0;

      //calculate the stimulus current, Istim

      //// Membrane Potential
      ////
      const double I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          // [uA/uF]
      const double I_Cl_tot = I_ClCa+I_Clbk+I_ClCFTR;                        // [uA/uF]
      const double I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
      const double I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
      //diff_v = -(I_Ca_tot+I_K_tot+I_Na_tot-I_app);

      const double diff_v = -I_tot;
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
         CHECK_BLOWUP(diff_d);
         CHECK_BLOWUP(diff_f);
         CHECK_BLOWUP(diff_xtof);
         CHECK_BLOWUP(diff_ytof);
         CHECK_BLOWUP(diff_xkr);
         CHECK_BLOWUP(diff_xks);
         CHECK_BLOWUP(diff_xkur);
         CHECK_BLOWUP(diff_ykur);
         CHECK_BLOWUP(diff_fcaBj );
         CHECK_BLOWUP(diff_fcaBsl );
         CHECK_BLOWUP(diff_RyRr );
         CHECK_BLOWUP(diff_RyRo );
         CHECK_BLOWUP(diff_RyRi );
         CHECK_BLOWUP(diff_NaBj );
         CHECK_BLOWUP(diff_NaBsl );
         CHECK_BLOWUP(diff_TnCL );
         CHECK_BLOWUP(diff_TnCHc );
         CHECK_BLOWUP(diff_TnCHm );
         CHECK_BLOWUP(diff_CaM );
         CHECK_BLOWUP(diff_Myc );
         CHECK_BLOWUP(diff_Mym );
         CHECK_BLOWUP(diff_SRB );
         CHECK_BLOWUP(diff_SLLj );
         CHECK_BLOWUP(diff_SLLsl );
         CHECK_BLOWUP(diff_SLHj );
         CHECK_BLOWUP(diff_SLHsl );
         CHECK_BLOWUP(diff_Csqnb );
         CHECK_BLOWUP(diff_Naj );
         CHECK_BLOWUP(diff_Nasl );
         CHECK_BLOWUP(diff_Nai );
         CHECK_BLOWUP(diff_Ki );
         CHECK_BLOWUP(diff_Casr );
         CHECK_BLOWUP(diff_Caj );
         CHECK_BLOWUP(diff_Casl );
         CHECK_BLOWUP(diff_Cai );


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
            PRINT_STATE(d);
            PRINT_STATE(f);
            PRINT_STATE(xtof);
            PRINT_STATE(ytof);
            PRINT_STATE(xkr);
            PRINT_STATE(xks);
            PRINT_STATE(xkur);
            PRINT_STATE(ykur);
            PRINT_STATE(fcaBj );
            PRINT_STATE(fcaBsl );
            PRINT_STATE(RyRr );
            PRINT_STATE(RyRo );
            PRINT_STATE(RyRi );
            PRINT_STATE(NaBj );
            PRINT_STATE(NaBsl );
            PRINT_STATE(TnCL );
            PRINT_STATE(TnCHc );
            PRINT_STATE(TnCHm );
            PRINT_STATE(CaM );
            PRINT_STATE(Myc );
            PRINT_STATE(Mym );
            PRINT_STATE(SRB );
            PRINT_STATE(SLLj );
            PRINT_STATE(SLLsl );
            PRINT_STATE(SLHj );
            PRINT_STATE(SLHsl );
            PRINT_STATE(Csqnb );
            PRINT_STATE(Naj );
            PRINT_STATE(Nasl );
            PRINT_STATE(Nai );
            PRINT_STATE(Ki );
            PRINT_STATE(Casr );
            PRINT_STATE(Caj );
            PRINT_STATE(Casl );
            PRINT_STATE(Cai );
            
#undef PRINT_STATE
            
            exit(255);
         }
      }
      
      
      //EDIT_STATE
      state_[ii].m = mss-(mss-m)*exp(-dt/taum);
      state_[ii].h = hss-(hss-m)*exp(-dt/tauh);
      state_[ii].j = jss-(jss-m)*exp(-dt/tauj);
      state_[ii].mL = mLss-(mLss-m)*exp(-dt/taumL);
      state_[ii].hL = hLss-(hLss-m)*exp(-dt/tauhL);
      state_[ii].d = dss-(dss-m)*exp(-dt/taud);
      state_[ii].f = fss-(fss-m)*exp(-dt/tauf);
      state_[ii].xtof = xtofss-(xtofss-m)*exp(-dt/tauxtof);
      state_[ii].ytof = ytofss-(ytofss-m)*exp(-dt/tauytof);
      state_[ii].xkr = xkrss-(xkrss-m)*exp(-dt/tauxkr);
      state_[ii].xks = xksss-(xksss-m)*exp(-dt/tauxks);
      state_[ii].xkur = xkurss-(xkurss-m)*exp(-dt/tauxkur);
      state_[ii].ykur = ykurss-(ykurss-m)*exp(-dt/tauykur);
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
   Vm.assign(Vm.size(), -85.0);
   State initState;
   //EDIT_STATE
   initState.m=0.0;
   initState.h=1.0;
   initState.j=1.0;
   initState.mL=0.0;
   initState.hL=1.0;
   initState.d=0.0;
   initState.f=1.0;
   initState.fcaBj=0.025;
   initState.fcaBsl=0.015;
   initState.xtof=0.0;
   initState.ytof=1.0;
   initState.xkr=0.0;
   initState.xks= 0.0;
   initState.xkur=0.0;
   initState.ykur=1.0;
   initState.RyRr=1.0;
   initState.RyRo=0.0;
   initState.RyRi=0.0;
   initState.NaBj=3.5;
   initState.NaBsl=0.8;
   initState.TnCL=0.01;
   initState.TnCHc=0.1;
   initState.TnCHm=0.01;
   initState.CaM=3e-4;
   initState.Myc=1.3e-3;
   initState.Mym=0.14;
   initState.SRB=2.0e-3;
   initState.SLLj=0.01;
   initState.SLLsl=0.1;
   initState.SLHj=7.3e-3;
   initState.SLHsl=7.3e-2;
   initState.Csqnb= 1.25;
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

const char* varNames[] = 
{
   //EDIT_STATE
   "m",
   "h",
   "j",
   "mL",
   "hL",
   "d",
   "f",
   "fcaBj",
   "fcaBsl",
   "xtof",
   "ytof",
   "xkr",
   "xks",
   "xkur",
   "ykur",
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

#define HANDLE_OFFSET 1000
int ThisReaction::getVarHandle(const std::string& varName) const
{
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      if (varNames[ivar] == varName) 
      {
         return HANDLE_OFFSET+ivar;
      }
   }
   assert(0 && "Control should never get here.");
   return -1;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
   //assure all states are doubles.
   int STATIC_ASSERT_checkAllDouble[(NUMVARS == sizeof(state_[0])/sizeof(double))? 1: 0];
   reinterpret_cast<double*>(&state_[iCell])[varHandle-HANDLE_OFFSET] = value;
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
   //assure all states are doubles.
   int STATIC_ASSERT_checkAllDouble[(NUMVARS == sizeof(state_[0])/sizeof(double))? 1: 0];
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
