#include "OHaraRudy.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "OHaraRudyNameArray.h"
typedef struct derived_st 
{ 
  double ENa, EK, EKs, phiCaMK,dCaMKtrap;
} DERIVED; 
typedef struct current_st 
{ 
   double NaFast, NaL, to, CaL, CaNa, CaK, Kr, Ks, K1, 
   NaCai, NaCass, NaK, Nab, Cab, Kb, pCa;
} CURRENTS;
typedef struct flux_st    { double diffNa, diffK, diffCa, rel, upNP, upCaMK, up, tr;} FLUXES; 
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define sigm(x)   (1.0/(1.0 + (x)))
#define sigm2(x)  (1.0/SQ(1.0 + (x)))
#define sige(x)   (1.0/(1.0 + exp((x))))
#define SQ(x)   ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

//Physical Constants 
static double R   = 8314.0;  //J/kmol/K
static double F   = 96485; //Coulomb/mol

static double T   = 310; //K
static double Cm     = 1.0    ; // uF; 

//Cell Geometry Parameters 
static   double AFV     ;//= Acap/(F*Vmyo);
static   double AFVss   ;//= Acap/(F*Vss);
static   double VssVmyo ;// = Vss/Vmyo; 
static   double VnsrVmyo;//= Vnsr/Vmyo; 
static   double VjsrVss ;//= Vjsr/Vss; 
static   double VjsrVnsr;//= Vjsr/Vnsr; 
void OHaraRudyCellular()
{
   double L = 0.01; // cm 
   double r = 0.0011; // cm
   double RCG = 2.0; 
   double Vcell ; // uL   Vcell = pi*r^2*L  Volume of cell; 
   double Ageo  ;// cm^2  Ageo =  2*pi*(r^2 + r*L)  Surface Area of cell 
   double Acap  ;// cm^2  Acap = RCG*Ageo 
   double Vmyo  ;// uL   Vmyo = 0.6800*Vcell 
   double Vnsr  ;// uL   Vnsr = 0.0552*Vcell 
   double Vjsr  ;// uL   Vjsr = 0.0048*Vcell 
   double Vss   ;// uL   Vss  = 0.0200*Vcell 
   double pi = 3.14; 
   Vcell = pi *r*r*L * 1e3; 
   Ageo  = 2*pi*(r*r+r*L);
   Acap  = RCG*Ageo; 
   Vmyo = 0.68*Vcell; 
   Vnsr = 0.0552*Vcell; 
   Vjsr = 0.0048*Vcell; 
   Vss = 0.02*Vcell; 
//
   AFV = Acap/(F*Vmyo);
   AFVss = Acap/(F*Vss);
   VssVmyo  = Vss/Vmyo; 
   VnsrVmyo = Vnsr/Vmyo; 
   VjsrVss  = Vjsr/Vss; 
   VjsrVnsr  = Vjsr/Vnsr; 
}
// External Concentrations 

static double PRNaK = 0.01833; 
static double Nao =  140; //mM;
static double Cao =  1.8; // mM
static double Ko  =  5.4; // mM

static double RTF ;
static double FRT ;
static double RTFlogNao; 
static double RTFlogKo; 
static double RTFlogKoNao; 
void OHaraRudyPerComputedConstants()
{
   RTF = (R*T)/F;
   FRT = F/(R*T);
   RTFlogNao= RTF*log(Nao); 
   RTFlogKo= RTF*log(Ko); 
   RTFlogKoNao=RTF*log(Ko+PRNaK*Nao);
}
void OHaraRudySetup()
{
   OHaraRudyCellular();
   OHaraRudyPerComputedConstants();
}

void integrateState(STATE *S, STATE *D, double dt)
{
   double *s = (double *)S; 
   double *d = (double *)D; 
   for (int i=37;i<41;i++) 
   {
      s[i] += dt*d[i]; 
   }
}
void writeState(FILE *file,char *stateVarNames[], STATE *S, STATE *D, double t)
{
   double *s = (double *)S; 
   double *d = (double *)D; 
   for (int i=0;i<41;i++) 
   {
      fprintf(file,"%2d %-11s %8.5f %24.12e %28.16e\n",i,stateVarNames[i],t,s[i],d[i]); 
   }
}
void reversalPotentials(CELLPARMS *cP, STATE *state, DERIVED *derived)
{
   derived->ENa = RTFlogNao-RTF*log(state->Nai);
   derived->EK  = RTFlogKo -RTF*log(state->Ki);
   derived->EKs = RTFlogKoNao-RTF*log((state->Ki+PRNaK*state->Nai));
}

void CaMK(CELLPARMS *cP, STATE *state, DERIVED *derived)
{
   double aCaMK = 0.05; //1/ms
   double bCaMK = 0.00068; // 1/ms
   double CaMK0 = 0.05; 
   double KmCaMK = 0.15;  // mM
   double KmCaM  = 0.0015;  // mM
   double CaMKBound = CaMK0 * (1-state->CaMKtrap)*sigm(KmCaM/state->Cass); 
   double CaMKActive = CaMKBound + state->CaMKtrap; 
   derived->phiCaMK = sigm(KmCaMK/CaMKActive); 
   derived->dCaMKtrap = aCaMK * CaMKBound * CaMKActive - bCaMK*state->CaMKtrap; 
}
void  initialState(STATE *state)
{
   state->Vm  = -87.84;
   state->Nai = 7.23; 
   state->Nass =7.23;//*0.9; 
   state->Ki = 143.79; 
   state->Kss= 143.79;//*0.9; 
   state->Cai= 8.54e-5 ;
   state->Cass= 8.43e-5 ;
   state->Cansr=1.61; 
   state->Cajsr=1.56; 
   state->m = 0.0074621;
   state->hFast = 0.692591; 
   state->hSlow = 0.692574; 
   state->j = 0.692477; 
   state->hCaMKSlow = 0.448501;
   state->jCaMK = 0.692413; 
   state->mL = 0.000194015; 
   state->hL = 0.496116; 
   state->hLCaMK = 0.265885; 
   state->a = 0.00101185; 
   state->iFast = 0.999542; 
   state->iSlow = 0.589579; 
   state->aCaMK = 0.000515567; 
   state->iCaMKFast = 0.999542;
   state->iCaMKSlow = 0.641861;
   state->d = 2.43015e-9;
   state->fFast = 1.0; 
   state->fSlow = 0.910671; 
   state->fCaFast=1.0; 
   state->fCaSlow=0.99982; 
   state->jCa = 0.999977; 
   state->n = 0.00267171; 
   state->fCaMKFast = 1.0; 
   state->fCaCaMKFast = 1.0; 
   state->XrFast = 8.26608e-6; 
   state->XrSlow = 0.453268;   
   state->Xs1 = 0.270492; 
   state->Xs2 = 0.0001963;
   state->XK1 = 0.996801;
   state->JrelNP = 2.53943e-5;  //mM/ms
   state->JrelCaMK = 3.17262e-7;  //mM/ms
   state->CaMKtrap = 0.0124065; 
}
void  initialConditions1(STATE *state)
{
   state->Vm = -87.5;
   state->Nai = 7.00; 
   state->Nass =7.00; 
   state->Ki = 145.00; 
   state->Kss= 145.00; 
   state->Cai= 1.00e-4 ;
   state->Cass= 1.00e-4 ;
   state->Cansr=1.20; 
   state->Cajsr=1.20; 
   state->m = 0.00;
   state->hFast = 1.00; 
   state->hSlow = 1.00; 
   state->j = 1.00; 
   state->hCaMKSlow = 1.00;
   state->jCaMK = 1.00; 
   state->mL = 0.00; 
   state->hL = 1.00; 
   state->hLCaMK = 1.00; 
   state->a = 0.00; 
   state->iFast = 1.00; 
   state->iSlow = 1.00; 
   state->aCaMK = 0.00; 
   state->iCaMKFast = 1.00;
   state->iCaMKSlow = 1.00;
   state->d = 0.00;
   state->fFast = 1.0; 
   state->fSlow = 1.00; 
   state->fCaFast=1.00; 
   state->fCaSlow=1.00; 
   state->jCa = 1.00; 
   state->n = 0.00 ;
   state->fCaMKFast = 1.0; 
   state->fCaCaMKFast = 1.0; 
   state->XrFast = 0.00; 
   state->XrSlow = 0.00;   
   state->Xs1 = 0.00; 
   state->Xs2 = 0.00;
   state->XK1 = 1.00;
   state->JrelNP = 0.00;  //mM/ms
   state->JrelCaMK = 0.0;  //mM/ms
   state->CaMKtrap = 0.00; 
}


void INaFunc(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I, STATE *D, double dt)
{
#define hSlowTauC   0.009794   //code
   //#define hSlowTauC   0.009764   //paper

   double V = state->Vm; 

   // input
   double ENa=derived->ENa;
   double phiCaMK=derived->phiCaMK;

   double mTauR; 
   // INaFast calc
   {
      //Gates needed to calculate INaFast; 
      double m=state->m; 
      double hFast=state->hFast; 
      double hSlow=state->hSlow; 
      double j=state->j; 
      double hCaMKSlow=state->hCaMKSlow;
      double jCaMK = state->jCaMK;   

      double AhFast=0.99;
      double AhSlow = 1-AhFast; 
      double h     = AhFast*hFast     + AhSlow*hSlow; 
      double AhCaMKFast = AhFast; 
      double AhCaMKSlow = AhSlow; 
      double hCaMKFast = hFast; 
      double hCaMK = AhCaMKFast*hCaMKFast + AhCaMKSlow*hCaMKSlow; 

      I->NaFast = cP->GNaFast*(V-ENa)*m*m*m*((1-phiCaMK)*h*j+phiCaMK*hCaMK*jCaMK); 

      //double mMhu =  sige(-(V + 39.57)/9.871);  //OHR orginal  mMhu
      double mMhu = SQ(sige(-(V + 56.86)/9.030)); //TT06  mMhu
      mTauR = 6.765*exp((V + 11.64)/34.77)+8.552*exp(-(V + 77.42)/5.955); 
      double  dm = (mMhu-m)*mTauR;  // gate

      double hMhu = sige((V + 82.90)/6.086); 
      double hFastTauR = 1.432e-5 * exp(-(V + 1.196)/6.285) + 6.149 * exp((V + 0.5096)/20.27);
      double hSlowTauR = hSlowTauC * exp (-(V + 17.95)/28.05) + 0.3343 * exp((V+5.730)/56.66);
      double dhFast = (hMhu - hFast)* hFastTauR;  // gate
      double dhSlow = (hMhu - hSlow)* hSlowTauR;  // gate

      double jMhu = hMhu; 
      double jTau  = 2.038 + 1.0/(0.02136*exp(-(V + 100.6)/8.281) + 0.3052*exp((V + 0.9941)/38.45));
      double jTauR = 1.0/jTau; 
      double dj = (jMhu-j)*jTauR;  // gate

      double jCaMKMhu = jMhu; 
      double jCaMKTauR = jTauR/1.46; 
      double djCaMK = (jCaMKMhu-jCaMK)*jCaMKTauR;  // gate

      double hCaMKMhu = sige((V+89.1)/6.086);
      double hCaMKSlowTauR = hSlowTauR/3.0; 
      double dhCaMKSlow = (hCaMKMhu-hCaMKSlow)*hCaMKSlowTauR;  // gate
      D->m = dm; 
      D->hFast = dhFast; 
      D->hSlow = dhSlow; 
      D->j = dj;
      D->jCaMK = djCaMK; 
      D->hCaMKSlow=dhCaMKSlow; 
      //state->m += dt*dm;    //Forward Euler for original OR mGate 
      double tauRdt = 1-exp(-dt*mTauR); //Rush Larsen needed with TT06 mMhu; 
      state->m += (mMhu-m)*tauRdt;     // Rush Larsen needed with TT06 mMhu; 
      state->hFast += dt*dhFast; 
      state->hSlow += dt*dhSlow; 
      state->j     += dt*dj; 
      state->hCaMKSlow += dt*dhCaMKSlow; 
      state->jCaMK += dt*djCaMK;
   }

   // INaL Calculation 
   {
      //   Gates needed to calculate INaL; 
      double mL=state->mL;
      double hL=state->hL; 
      double hLCaMK=state->hLCaMK;     

      I->NaL= cP->GNaL*(V-ENa)*mL*((1-phiCaMK)*hL+phiCaMK*hLCaMK); 

      double mLMhu =  sige(-(V+42.85)/5.264);
      double mLTauR = mTauR; 
      double dmL = (mLMhu-mL)*mLTauR;  // gate

      double hLMhu = sige((V+87.61)/7.488);
      double hLTauR = 0.005; 
      double dhL = (hLMhu-hL)*hLTauR;  // gate

      double hLCaMKMhu = sige((V+93.81)/7.488);
      double hLCaMKTauR = hLTauR/3.0; 
      double dhLCaMK = (hLCaMKMhu-hLCaMK)*hLCaMKTauR;  // gate
      D->mL = dmL; 
      D->hL = dhL; 
      D->hLCaMK=dhLCaMK; 
      state->mL += dt*dmL; 
      state->hL += dt*dhL; 
      state->hLCaMK += dt*dhLCaMK; 
   }
}
void ItoFunc(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I, STATE *D, double dt )
{
#define iSlowTauC   1.780e-8 // code
#define aTauC0     18.4099   // code
#define aTauC1     29.3814   // code

   //#define iSlowTauC 1.7808e-8 // paper
   //#define aTauC0   18.41      // paper 
   //#define aTauC1   29.38      // paper 

   double V = state->Vm; 
   double EK = derived->EK; 

   double a     =state->a;
   double iFast = state->iFast; 
   double iSlow = state->iSlow; 
   double aCaMK = state->aCaMK;
   double iCaMKFast = state->iCaMKFast; 
   double iCaMKSlow = state->iCaMKSlow; 
   double phiCaMK   = derived->phiCaMK; 

   double AiFast = sige((V-213.6)/151.2); 
   double AiSlow = 1 - AiFast; 
   double AiCaMKFast = AiFast; 
   double AiCaMKSlow = AiSlow; 

   double i = AiFast * iFast + AiSlow * iSlow; 
   double iCaMK = AiCaMKFast * iCaMKFast + AiCaMKSlow * iCaMKSlow; 

   I->to=cP->Gto*(V-EK) * ((1-phiCaMK)*a*i + phiCaMK *aCaMK*iCaMK); 

   double aMhu = sige(-(V-14.34)/14.82); 
   //double aTauR = ((1./1.2089)*sige(-(V-18.4099)/29.3814) + 3.5*sige((V+100)/29.3814))/1.0515;
   double aTauR = ((1./1.2089)*sige(-(V-aTauC0)/29.3814) + 3.5*sige((V+100)/aTauC1))/1.0515;
   double da = (aMhu-a)*aTauR;  // gate

   double iMhu = sige((V+43.94)/5.711); 
   double delta = 1.0-cP->aDelta*sige((V+70)/5.0); 
   double iFastTau = 4.562 + 1/(0.3933*exp(-(V+100)/100) + 0.08004*exp((V+50)/16.59));
   iFastTau = iFastTau*delta; 
   //double iSlowTau = 23.62 + 1/(0.001416*exp(-(V+96.52)/59.05) + 1.7808e-8*exp((V+114.1)/8.079));
   double iSlowTau = 23.62 + 1/(0.001416*exp(-(V+96.52)/59.05) + iSlowTauC*exp((V+114.1)/8.079));
   iSlowTau = iSlowTau*delta; 
   double iFastTauR = 1/iFastTau; 
   double iSlowTauR = 1/iSlowTau; 

   double diSlow = (iMhu-iSlow)*iSlowTauR;  // gate
   double diFast = (iMhu-iFast)*iFastTauR;  // gate

   double aCaMKMhu = sige(-(V-24.34)/14.82); 
   double aCaMKTauR = aTauR; 
   double daCaMK = (aCaMKMhu-aCaMK)*aCaMKTauR;  // gate

   double iCaMKMhu = iMhu; 
   double deltaCaMKdevelop = 1.354 + 1e-4/(exp((V-167.4)/15.89)+exp(-(V-12.23)/0.2154));
   double deltaCaMKrecover = 1.0 - 0.5*sige((V+70)/20);
   double deltaCaMKR =  1/(deltaCaMKdevelop*deltaCaMKrecover); 
   double iCaMKFastTauR = iFastTauR * deltaCaMKR;
   double iCaMKSlowTauR = iSlowTauR * deltaCaMKR;
   double diCaMKSlow = (iCaMKMhu-iCaMKSlow)*iCaMKSlowTauR;  // gate
   double diCaMKFast = (iCaMKMhu-iCaMKFast)*iCaMKFastTauR;  // gate
   D->a  = da; 
   D->iSlow = diSlow; 
   D->iFast = diFast; 
   D->aCaMK = daCaMK; 
   D->iCaMKSlow = diCaMKSlow; 
   D->iCaMKFast = diCaMKFast; 
   state->a += dt*da; 
   state->iSlow += dt*diSlow; 
   state->iFast += dt*diFast; 
   state->aCaMK += dt*daCaMK; 
   state->iCaMKSlow += dt*diCaMKSlow; 
   state->iCaMKFast += dt*diCaMKFast; 
}
void ICaL_ICaNa_ICaK_Func(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I, STATE *D, double dt )
{
   double V = state->Vm; 
   double fFast     =state->fFast;
   double fSlow     =state->fSlow;
   double fCaFast     =state->fCaFast;
   double fCaSlow     =state->fCaSlow;
   double fCaMKFast     =state->fCaMKFast;
   double fCaCaMKFast     =state->fCaCaMKFast;
   double n     =state->n;
   double jCa     =state->jCa;
   double d     =state->d;
   double Cass     =state->Cass;
   double Nass     =state->Nass;
   double Kss     =state->Kss;
   double phiCaMK   = derived->phiCaMK; 
   double AfFast=0.6; 
   double AfSlow=1-AfFast;  
   double f = AfFast*fFast + AfSlow*fSlow; 

   double AfCaFast=0.3+0.6*sige((V-10)/10); 
   double AfCaSlow=1-AfCaFast;  
   double fCa = AfCaFast*fCaFast + AfCaSlow*fCaSlow; 

   double AfCaMKFast=AfFast; 
   double AfCaMKSlow=1-AfCaMKFast;  
   double fCaMKSlow = fSlow; 
   double fCaMK = AfCaMKFast*fCaMKFast + AfCaMKSlow*fCaMKSlow; 

   double AfCaCaMKFast = AfCaFast; 
   double AfCaCaMKSlow = AfCaSlow; 
   double fCaCaMKSlow = fCaSlow; 
   double fCaCaMK = AfCaCaMKFast*fCaCaMKFast + AfCaCaMKSlow*fCaCaMKSlow; 


   double chi1 = d*(1-phiCaMK)*(f    *(1-n)+fCa    *n*jCa); 
   double chi2 = d*   phiCaMK *(fCaMK*(1-n)+fCaCaMK*n*jCa); 
   double chi  = chi1 + 1.1 * chi2; 
   double gammaCai= 1.0 ;
   double gammaCao= 0.341; 
   double zCa = 2;
   double psiCa   = SQ(zCa)*V*F*FRT*(gammaCai*Cass*exp(zCa*V*FRT)-gammaCao*Cao)/(exp(zCa*V*FRT)-1.0); 
   double gammaNai= 0.75;
   double gammaNao= 0.75; 
   double zNa = 1;
   double psiCaNa = SQ(zNa)*V*F*FRT*(gammaNai*Nass*exp(zNa*V*FRT)-gammaNao*Nao)/(exp(zNa*V*FRT)-1.0); 
   double gammaKi= 0.75;
   double gammaKo= 0.75; 
   double zK = 1;
   double psiCaK  = SQ(zK )*V*F*FRT*(gammaKi * Kss*exp(zK*V*FRT)-gammaKo*Ko)/(exp(zK*V*FRT)-1.0); 
   /*
      double PCaL      =          PCa;
      double PCaNa     = 1.250e-3*PCa;
      double PCaK      = 3.574e-4*PCa;
      double PCaCaMK   = 1.10e0  *PCa;
      double PCaNaCaMK = 1.250e-3*PCaCaMK;
      double PCaKCaMK  = 3.574e-4*PCaCaMK;
      double  ICaLBar = PCa * psiCa; 
      double  ICaLCaMKBar = PCaCaMK * psiCa; 
      double ICaL = ICaLBar*chi1 + ICaLCaMKBar*chi2;
    */

   //I->CaL = cP->PCaL*psiCa*chi1 + PCaCaMK*psiCa*chi2; 
   I->CaL   = cP->PCaL*psiCa*chi; 

   //I->CaNa = cP->PCaNa*psiCaNa*chi1 + PCaNaCaMK*psiCaNa * chi2; 
   I->CaNa   = cP->PCaNa*psiCaNa*chi; 

   //I->CaK  = (cP->PCaK * chi1 + cP->PCaKCaMK * chi2)*psiCaK; 
   I->CaK    = cP->PCaK*psiCaK*chi; 


   double dMhu = sige(-(V+3.940)/4.230); 
   double dTau = 0.6 + 1/(exp(-0.05*(V+6.0)) + exp(0.09*(V+14.0))); 
   double dTauR = 1.0/dTau; 
   double dd = (dMhu-d)*dTauR;  // gate

   double fMhu = sige((V+19.58)/3.696); 
   double fFastTau = 7.0 + 1/(0.0045*exp(-(V+20)/10) + 0.0045*exp((V+20)/10)); 
   double fFastTauR = 1/fFastTau; 

   double fSlowTau  = 1000 + 1/(3.5e-5*exp(-(V+5)/4) + 3.5e-5*exp((V+5)/6)); 
   double fSlowTauR = 1/fSlowTau; 
   double dfSlow = (fMhu-fSlow)*fSlowTauR;  // gate
   double dfFast = (fMhu-fFast)*fFastTauR;  // gate

   double fCaMhu  =  fMhu; 
   double fCaFastTau = 7.0 + 1/(0.04*exp(-(V-4)/7) + 0.04*exp((V-4)/7)); 
   double fCaFastTauR = 1/fCaFastTau; 
   double fCaSlowTau = 100 + 1/(1.2e-4*exp(-(V)/3) + 1.2e-4*exp((V)/7)); 
   double fCaSlowTauR = 1/fCaSlowTau; 
   double dfCaSlow = (fCaMhu-fCaSlow)*fCaSlowTauR;  // gate
   double dfCaFast = (fCaMhu-fCaFast)*fCaFastTauR;  // gate

   double jCaMhu = fMhu; 
   double jCaTauR = 1/75.0; 
   double djCa = (jCaMhu-jCa)*jCaTauR;  // gate

   double fCaMKMhu = fMhu;
   double fCaMKFastTauR = fFastTauR/2.5; 
   double dfCaMKFast = (fCaMKMhu-fCaMKFast)*fCaMKFastTauR;  // gate

   double fCaCaMKFastMhu = fMhu; 
   double fCaCaMKFastTauR = fCaFastTauR/2.5; 
   double dfCaCaMKFast = (fCaCaMKFastMhu-fCaCaMKFast)*fCaCaMKFastTauR;  // gate

   double Kmn = 0.002; 
   double Kp2n=1000.0; 
   double Km2n=jCa*1.0; 
   double an = 1.0/(Kp2n/Km2n + SQ(SQ(1+Kmn/Cass)));
   double dn =  an*Kp2n -n * Km2n;  //gate
   D->d = dd; 
   D->fSlow = dfSlow; 
   D->fFast = dfFast; 
   D->fCaSlow = dfCaSlow; 
   D->fCaFast = dfCaFast; 
   D->jCa     = djCa; 
   D->fCaMKFast = dfCaMKFast; 
   D->fCaCaMKFast = dfCaCaMKFast; 
   D->n  = dn; 
   state->d += dt*dd; 
   state->fSlow += dt*dfSlow; 
   state->fFast += dt*dfFast; 
   state->fCaSlow += dt*dfCaSlow; 
   state->fCaFast += dt*dfCaFast; 
   state->jCa     += dt*djCa; 
   state->fCaMKFast += dt*dfCaMKFast; 
   state->fCaCaMKFast += dt*dfCaCaMKFast; 
   state->n  += dt*dn; 

}

void IKrFunc(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I, STATE *D, double dt)
{

   // EK, XrFast, XrSlow
   double V = state->Vm; 
   double EK = derived->EK; 
   double XrFast = state->XrFast; 
   double XrSlow = state->XrSlow; 
   double AXrFast = sige((V+54.81)/38.21); 
   double AXrSlow = 1.0 - AXrFast; 
   double Xr = AXrFast*XrFast + AXrSlow*XrSlow; 
   double RKr = sige((V+55)/75)*sige((V-10)/30); 
   I->Kr = cP->GKr * (V-EK) *sqrt(Ko/5.4) * Xr * RKr; 

   double XrMhu = sige(-(V+8.337)/6.789); 
   double XrFastTau = 12.98  + 1/(0.3652*exp((V-31.66)/3.869) + 4.123e-5*exp(-(V-47.78)/20.38));
   double XrSlowTau =  1.865 + 1/(0.06629*exp((V-34.70)/7.355) + 1.128e-5*exp(-(V-29.74)/25.94));
   double XrFastTauR = 1/XrFastTau; 
   double XrSlowTauR = 1/XrSlowTau; 
   double dXrSlow = (XrMhu - XrSlow)*XrSlowTauR;  // gate
   double dXrFast = (XrMhu - XrFast)*XrFastTauR;  // gate
   D->XrSlow = dXrSlow; 
   D->XrFast = dXrFast; 
   state->XrSlow += dt*dXrSlow; 
   state->XrFast += dt*dXrFast; 
}
void IKsFunc(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I, STATE *D, double dt)
{
   double V = state->Vm; 
   double EKs = derived->EKs; 
   double Xs1 = state->Xs1; 
   double Xs2 = state->Xs2; 
   double Cai = state->Cai; 

   double phi =  1 + 0.6*sigm(pow(3.8e-5/Cai,1.4)); 
   I->Ks = cP->GKs *(V-EKs) * phi * Xs1*Xs2;

   double XsMhu = sige(-(V+11.60)/8.932); 
   double Xs1Tau =  817.3 + 1/(2.326e-4*exp((V+48.28)/17.80) + 1.292e-3*exp(-(V+210.0)/230));
   double Xs1TauR =  1/Xs1Tau;
   double Xs2TauR = 1e-2*exp((V-50.0)/20.0) + 1.93e-2*exp(-(V+66.54)/31);
   double dXs1 = (XsMhu-Xs1)*Xs1TauR;  // gate
   double dXs2 = (XsMhu-Xs2)*Xs2TauR;  // gate
   D->Xs1 = dXs1; 
   D->Xs2 = dXs2; 
   state->Xs1 += dt*dXs1; 
   state->Xs2 += dt*dXs2; 
}
void IK1Func(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I, STATE *D, double dt)
{
   double V = state->Vm; 
   double EK = derived->EK; 
   double XK1 = state->XK1; 

   double RK1 = sige((V+105.8-2.6*Ko)/9.493); 
   double phi = sqrt(Ko)*RK1; 
   I->K1 = cP->GK1 * (V-EK) * XK1 * phi;

   double XK1Mhu = sige(-(V + 2.5538*Ko + 144.59)/(1.5692*Ko + 3.8115)); 
   double XK1TauR = (exp(-(V+127.2)/20.36) + exp((V+236.8)/69.33))/122.2;
   double dXK1 = (XK1Mhu-XK1)*XK1TauR; // gate
   D->XK1 = dXK1; 
   state->XK1 += dt*dXK1; 
}
//  All currents function have no dependence on state of gates. 
static inline double INaCaYFunc( double V,double CaY, double NaY  )
{
   double kNa1 = 15; //mM
   double kNa2 =  5; //mM
   double kNa3 =  88.12; //mM
   double kasymm =  12.5; //mM
   double wNa = 6e4 ; // Hz
   double wCa = 6e4 ; // Hz
   double wNaCa = 5e3; //Hz
   double kCaOn = 1.5e6 ; //mM/ms
   double kCaOff=5e3; // Hz
   double qNa = 0.5224; 
   double qCa = 0.1670; 
   double hCa = exp(qCa*V*FRT);
   double hNa = exp(qNa*V*FRT);
   double h1 = 1 + NaY*(1 + hNa)/kNa3; 
   double h2 = NaY *hNa/(kNa3*h1); 
   double h3 = 1/h1;
   double h4 = 1 +NaY*(1+NaY/kNa2)/kNa1;
   double h5 = SQ(NaY)/(h4*kNa1*kNa2); 
   double h6 = 1/h4; 
   double h7 = 1+Nao*(1+1/hNa)/kNa3; 
   double h8 = Nao/(kNa3*hNa*h7); 
   double h9 = 1/h7; 
   double h10 = kasymm + 1 + Nao*(1+Nao/kNa2)/kNa1; 
   double h11 = SQ(Nao)/(h10*kNa1*kNa2); 
   double h12 = 1/h10; 
   double k1 = h12*Cao*kCaOn; 
   double k2 = kCaOff; 
   double k3p = h9*wCa; 
   double k3pp=h8*wNaCa; 
   double k3 = k3p + k3pp; 
   double k4p = h3*wCa/hCa; 
   double k4pp= h2*wNaCa; 
   double k4 = k4p + k4pp; 
   double k5 = kCaOff; 
   double k6=h6*CaY*kCaOn;
   double k7 = h5*h2*wNa;
   double k8 = h8*h11*wNa; 
   double x1 = k2*k4*(k7+k6) + k5*k7*(k2+k3); 
   double x2 = k1*k7*(k4+k5) + k4*k6*(k1+k8); 
   double x3 = k1*k3*(k7+k6) + k8*k6*(k2+k3); 
   double x4 = k2*k8*(k4+k5) + k3*k5*(k1+k8); 
   double E1 = x1/(x1+x2+x3+x4); 
   double E2 = x2/(x1+x2+x3+x4); 
   double E3 = x3/(x1+x2+x3+x4); 
   double E4 = x4/(x1+x2+x3+x4); 
   double KmCaAct = 150e-6;//mM
   double alloY = sigm(SQ(KmCaAct/CaY)); 
   double JNaCaNaY = 3*(E4*k7-E1*k8) + E3*k4pp - E2*k3pp; 
   double JNaCaCaY = E2*k2 - E1*k1; 
   double zNa = 1; 
   double zCa = 2; 
   double INaCaY = alloY*(zNa*JNaCaNaY+ zCa*JNaCaCaY); 
   return INaCaY; 
}
void INaCaFunc(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I)
{
   double V = state->Vm; 
   //double GNaCai  = 0.8*GNaCa ;
   //double GNaCass = 0.2*GNaCa ;
   I->NaCai =  cP->GNaCai *INaCaYFunc(state->Vm,state->Cai, state->Nai   );
   I->NaCass = cP->GNaCass*INaCaYFunc(state->Vm,state->Cass, state->Nass  );
}
void INaKFunc(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I)
{
   double V = state->Vm; 
   double Nai = state->Nai; 
   double Ki = state->Ki; 

   double k1Plus = 949.5;//Hz
   double k1Minus = 182.4;//mM^-1
   double k2Plus = 687.2;//Hz
   double k2Minus = 39.4 ;//Hz
   double k3Plus = 1899.;//Hz
   double k3Minus = 79300 ;//Hz*mM^-2
   double k4Plus = 639.0;//Hz
   double k4Minus = 40.0 ;//Hz
   double KoNai = 9.073 ; //mM
   double KoNao = 27.78 ; //mM
   double delta = -0.1550;  
   double KNai = KoNai*exp(delta*V*FRT/3); 
   double KNao = KoNao*exp((1-delta)*V*FRT/3); 
   double KKi = 0.5;//mM
   double KKo = 0.3582; //mM
   double MgADP = 0.05; 
   double MgATP = 9.8; 
   double KMgATP = 1.698e-7 ; //mM
   double HPlus  = 1e-7 ; //mM
   double sigmaP = 4.2; //mM
   double KHP = 1.698e-7; //mM
   double KNaP = 224;// mM
   double KKP  = 292;// mM
   double P = sigmaP/(1+HPlus/KHP + Nai/KNaP + Ki/KKP); 
   double a1 = k1Plus*CUBE(Nai/KNai)/(CUBE(1+Nai/KNai)+SQ(1+Ki/KKi) -1.0);
   double b1  = k1Minus*MgADP; 
   double a2 = k2Plus;
   double factor = 1.0/(CUBE(1+Nao/KNao)+SQ(1+Ko/KKo) -1.0);
   double b2  = k2Minus*CUBE(Nao/KNao)*factor;
   double a3 = k3Plus*SQ(Ko/KKo)*factor;
   double b3  = k3Minus*P*HPlus        /(1+MgATP/KMgATP);
   double a4  = k4Plus*(MgATP/KMgATP)  /(1+MgATP/KMgATP);
   double b4 = k4Minus*SQ(Ki/KKi)/(CUBE(1+Nai/KNai)+SQ(1+Ki/KKi) -1.0);
   double x1 = a4*a1*a2 + b2*b4*b3 + a2*b4*b3 + b3*a1*a2; 
   double x2 = b2*b1*b4 + a1*a2*a3 + a3*b1*b4 + a2*a3*b4; 
   double x3 = a2*a3*a4 + b3*b2*b1 + b2*b1*a4 + a3*a4*b1; 
   double x4 = b4*b3*b2 + a3*a4*a1 + b2*a4*a1 + b3*b2*a1; 
   double norm = 1.0/(x1+x2+x3+x4); 
   double E1 = x1*norm;
   double E2 = x2*norm;
   double E3 = x3*norm;
   double E4 = x4*norm;
   double zNa=1; 
   double zK=1; 
   double JNaKNa=3*(E1*a3-E2*b3);
   double JNaKK = 2*(E4*b1-E3*a1); 
   I->NaK = cP->PNaK*(zNa*JNaKNa + zK*JNaKK); 
}
void INab_ICab_IKb_IpCaFunc(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I )
{
   double V = state->Vm; 
   double Nai = state->Nai; 
   double Cai = state->Cai; 
   double EK = derived->EK;

   double zNa=1.0;//
   I->Nab = cP->PNab*SQ(zNa) *V*FRT*F*(Nai*exp(zNa*V*FRT)-Nao)/(exp(zNa*V*FRT)-1); 
   double gammaCai= 1.0;
   double gammaCao=0.341;
   double zCa = 2; 
   I->Cab = cP->PCab*SQ(zCa)*V*FRT*F*(gammaCai*Cai*exp(zCa*V*FRT)-gammaCao*Cao)/(exp(zCa*V*FRT)-1); 
   double xKb = 1/(1+exp(-(V-14.48)/18.34)) ; 
   I->Kb = cP->GKb*xKb*(V-EK) ;
   I->pCa = cP->GpCa * Cai/(0.0005 + Cai); 
}
void integrateNonGates(CELLPARMS *cP, STATE *state, DERIVED *derived, CURRENTS *I, STATE *D, double Istim, double dt)
{

   double V = state->Vm; 
   double Nass = state->Nass; 
   double Nai  = state->Nai; 
   double Cass = state->Cass; 
   double Cajsr = state->Cajsr; 
   double Cai  = state->Cai; 
   double Kss  = state->Kss; 
   double Ki   = state->Ki; 

   double tauDiffNa=2.0; // ms
   double tauDiffK=2.0; // ms
   double tauDiffCa=0.2;// ms
   double JdiffNa = (Nass-Nai)/tauDiffNa;
   double JdiffCa = (Cass-Cai)/tauDiffCa;
   double JdiffK  = (Kss-Ki)/tauDiffK;

   double alpha0   = (-I->CaL) /(1.0+SQ(SQ(SQ(1.5/state->Cajsr))));
   double beta0    = 1.0/(1+0.0123/state->Cajsr); 
   double mhuJrelNP  = cP->alphaJrelNP * alpha0 ;
   double mhuJrelCaMK= cP->alphaJrelCaMK*alpha0;
   double tauJrelNP  = cP->betaJrelNP * beta0;
   double tauJrelCaMK= cP->betaJrelCaMK*beta0; 
   tauJrelNP = MAX(tauJrelNP,0.001); 
   tauJrelCaMK = MAX(tauJrelCaMK,0.001); 

   //cP->relABRatio = 0.5; 
   //double bRelRatio = 1.25; 
   //double aRelRatio = bRelRatio; 
   //double bTau = 4.75; // ms
   //double aRel = cP->JrelABRatio*bTau;
   //double aRel = 0.5*bTau;
   //double bTauCaMK = bRelRatio*bTau; 
   //double aRelCaMK = relRatio*bTauCaMK;
   //double aRelCaMK = 0.5*bTauCaMK;
   //double mhuJrelCaMK = (aRelCaMK/aRel)*mhuJrelNP;
   //double  tauRelCaMK = (bTauCaMK/bTau)*tauRelNP; 

   //double mhuJrelNP=aRel * (-I->CaL) /(1.0+SQ(SQ(SQ(1.5/state->Cajsr))));
   //double mhuJrelCaMK= aRelRatio*mhuJrelNP;
   //double tauRelNP   = bTau/(1+0.0123/state->Cajsr); 
   //double tauRelCaMK = bRelRatio*tauRelNP; 
   //tauRelNP = MAX(tauRelNP,0.001); 
   //tauRelCaMK = MAX(tauRelCaMK,0.001); 

   double phiRelCaMK=derived->phiCaMK;
   double Jrel = (1-phiRelCaMK)*state->JrelNP + phiRelCaMK*state->JrelCaMK; 

   double JupNP = cP->cJup*Cai/(0.00092 + Cai) ; 
   double deltaKmPLB = 0.00017 ; // mM
   double deltaJupCaMK=1.75; 
   double JupCaMK = cP->cJup*(1 + deltaJupCaMK)*Cai/(0.00092 -deltaKmPLB + Cai) ;

   double phiUpCaMK=phiRelCaMK;
   double Jleak = 0.0039375*state->Cansr/15.0;
   double Jup = (1-phiUpCaMK)*JupNP + phiUpCaMK*JupCaMK - Jleak; 
   double tautr = 100 ;// ms
   double Jtr = (state->Cansr-state->Cajsr)/tautr; 

   //double CMDN =  0.05 ;  // mM
   double TRPN =  0.07 ;  // mM
   double BSR  =  0.047;  // mM
   double BSL  =  1.124;  // mM
   double CSQN = 10.0  ;  // mM
   double KmCMDN = 0.00238; // mM
   double KmTRPN = 0.0005 ; // mM
   double KmBSR  = 0.00087; // mM
   double KmBSL  = 0.0087 ; // mM
   double KmCSQN = 0.8    ; // mM

   double dNai = -(I->NaFast + I->NaL + 3*I->NaCai + 3*I->NaK + I->Nab)*AFV + JdiffNa*VssVmyo;
   double dNass= -(I->CaNa + 3*I->NaCass)*AFVss - JdiffNa; 
   double dKi =  -(I->to + I->Kr + I->Ks + I->K1 + I->Kb + Istim - 2*I->NaK )*AFV + JdiffK*VssVmyo;
   double dKss = -I->CaK * AFVss - JdiffK; 
   double bCai = 1/(1 + (cP->CMDN/KmCMDN)*sigm2(Cai/KmCMDN)+(TRPN/KmTRPN)*sigm2(Cai/KmTRPN));
   double dCai= bCai*(-0.5*(I->pCa + I->Cab - 2*I->NaCai)*AFV - Jup*VnsrVmyo + JdiffCa*VssVmyo); 
   double bCass= 1/(1 + (BSR /KmBSR )*sigm2(Cass/KmBSR)+(BSL/KmBSL)*sigm2(Cass/KmBSL));
   double dCass=bCass*(-0.5*(I->CaL - 2*I->NaCass)*AFVss + Jrel*VjsrVss - JdiffCa); 
   double dCansr = Jup-Jtr*VjsrVnsr; 
   double bCajsr = 1/(1+(CSQN/KmCSQN)*sigm2(Cajsr/KmCSQN)); 
   double dCajsr =bCajsr*(Jtr-Jrel); 
   double dJrelNP = (mhuJrelNP-state->JrelNP)/tauJrelNP; 
   double dJrelCaMK = (mhuJrelCaMK-state->JrelCaMK)/tauJrelCaMK; 
   double dV = -(I->NaFast + I->NaL + I->to + I->CaL + I->CaNa + I->CaK + I->Kr + I->Ks + I->K1 + I->NaCai + I->NaCass + I->NaK + I->Nab + I->Cab + I->Kb + I->pCa )/Cm;
   D->Nai  = dNai; 
   D->Nass = dNass; 
   D->Ki   = dKi; 
   D->Kss  = dKss; 
   D->Cai  = dCai; 
   D->Cass = dCass; 
   D->Cansr= dCansr; 
   D->Cajsr= dCajsr; 
   D->JrelNP= dJrelNP; 
   D->JrelCaMK= dJrelCaMK; 
   D->CaMKtrap=derived->dCaMKtrap; 
   D->Vm = dV; 
   state->Nai  += dt*dNai; 
   state->Nass += dt*dNass; 
   state->Ki   += dt*dKi; 
   state->Kss  += dt*dKss; 
   state->Cai  += dt*dCai; 
   state->Cass += dt*dCass; 
   state->Cansr+= dt*dCansr; 
   state->Cajsr+= dt*dCajsr; 

   state->JrelNP+= dt*dJrelNP; 
   state->JrelCaMK+= dt*dJrelCaMK; 
   state->CaMKtrap+=dt*derived->dCaMKtrap; 

}
void OHaraRudyInitialState(int cellType, STATE *state)
{
   initialState(state); 
}
CELLPARMS *OHaraRudyCellParms(int cellType)
{
static   CELLPARMS endo = { 75., 0.0075, 0.0200, 1e-4, 1.250e-7, 3.574e-8, 0.0460, 0.1908, 0.0034, 6.4e-4, 1.6e-4,  30.0, 3.75e-10, 2.5e-8, 0.003, 0.0005, 2.375, 4.75, 2.96875, 5.9375, 0.004375,  0.05, 0.0};
static   CELLPARMS epi ; 
static   CELLPARMS M ; 
   if  (cellType == ENDO_CELL) return &endo; 
   if  (cellType == EPI_CELL)
   {
      epi=endo; 
      epi.GNaL *= 0.6; 
      epi.Gto *= 4.0; 
      epi.PCaL  *= 1.2;
      epi.PCaNa *= 1.2;
      epi.PCaK  *= 1.2;
      epi.GKr   *= 1.3; 
      epi.GKs   *= 1.4; 
      epi.GK1   *= 1.2; 
      epi.GNaCai*= 1.1; 
      epi.GNaCass*= 1.1; 
      epi.PNaK   *= 0.9; 
      epi.GKb    *= 0.6; 
      epi.alphaJrelNP*=1.0;
      epi.alphaJrelCaMK*=1.0;
      epi.cJup *= 1.3; 
      epi.CMDN *= 1.3 ; 
      epi.aDelta = 0.95; 
      return &epi; 
   }
   if  (cellType == M_CELL)
   {
      M=endo; 
      M.GNaL *= 1.0; 
      M.Gto *= 4.0; 
      M.PCaL  *= 2.5;
      M.PCaNa *= 2.5;
      M.PCaK  *= 2.5;
      M.GKr   *= 0.8; 
      M.GKs   *= 1.0; 
      M.GK1   *= 1.3; 
      M.GNaCai*= 1.4; 
      M.GNaCass*= 1.4; 
      M.PNaK   *= 0.7; 
      M.GKb    *= 1.0; 
      M.alphaJrelNP*=1.7;
      M.alphaJrelCaMK*=1.7;
      M.cJup *= 1.0; 
      M.CMDN *= 1.0 ; 
      M.aDelta = 0.00  ;
      return &M; 
   }

}
double OHaraRudyIntegrate(double dt, double Istim, STATE *state, CELLPARMS *cP, STATE *D)
{
   CURRENTS I; 
   //STATE Dlocal;  STATE *D = Dlocal; 

   DERIVED derived; 
   static double t=0; 


   reversalPotentials(cP,state,&derived); 
   CaMK(cP,state,&derived); 
   INaFunc(cP,state, &derived, &I, D, dt);
   ItoFunc(cP,state, &derived, &I, D, dt);
   ICaL_ICaNa_ICaK_Func(cP, state, &derived, &I, D, dt);
   IKrFunc(cP, state, &derived, &I, D, dt );
   IKsFunc(cP, state, &derived, &I, D, dt);
   IK1Func(cP, state, &derived, &I, D, dt);

   INaCaFunc(cP, state, &derived, &I );
   INaKFunc(cP, state, &derived, &I );
   INab_ICab_IKb_IpCaFunc(cP, state, &derived, &I );
   integrateNonGates(cP, state, &derived, &I, D, Istim, dt );
   t += dt; 
    //printf("%f %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",t,I.NaFast, I.NaL, I.to, I.CaL, I.CaNa, I.CaK, I.Kr, I.Ks, I.K1, I.NaCai, I.NaCass, I.NaK, I.Nab, I.Cab, I.Kb, I.pCa); 
   return D->Vm; 
} 

#ifdef SA
void main()
{
   STATE state; 
   DERIVED derived; 
   CURRENTS I; 
   STATE D; 
   CELLPARMS *cP; 
   OHaraRudySetup(); 
   OHaraRudyInitialState(ENDO_CELL,&state); 
   CELLPARMS *endo = OHaraRudyCellParms(ENDO_CELL); 
   CELLPARMS *epi = OHaraRudyCellParms(EPI_CELL); 
   double dt=1.0e-6;
   double t=0.0; 
   int loop =0; 
   FILE *file0 = fopen("voltage.dat","w"); 
   FILE *file1 = fopen("state.dat","w"); 
   FILE *file2 = fopen("current.dat","w"); 
   double Istim = 0.0; 
   double dVm = 0.0; 
   if (loop%1==0) fprintf(file0,"%13.8f %20.15f %20.15f %20.15f\n",t,state.Vm, dVm,Istim); 
   while (t<.10) 
   {

      double rt = fmod(t,1000.0); 
      Istim = 0.0; 
      if (rt<0.5) Istim = -80.0 ; // uA/uF; 
      //STATE stateLast = state; 
      dVm=OHaraRudyIntegrate(dt,Istim,&state,epi,&D); 
      state.Vm += dt*(dVm-Istim/Cm);
      D.Vm = dVm; 
      //integrateState(&state,&D,dt); 
      //if (loop%100==0) writeState(file1,stateVarNames, &stateLast , &D, t);
      //if (loop%100==0) fprintf(file2,"%f %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",t,I.NaFast, I.NaL, I.to, I.CaL, I.CaNa, I.CaK, I.Kr, I.Ks, I.K1, I.NaCai, I.NaCass, I.NaK, I.Nab, I.Cab, I.Kb, I.pCa); 
      t += dt; 
      if (loop%1==0) fprintf(file0,"%13.8f %20.15f %20.15f %20.15f\n",t,state.Vm, dVm,Istim); 
      loop++; 
   } 
}
#endif


