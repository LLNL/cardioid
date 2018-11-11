#include "TT06Func.hh"
#include <stdio.h>
#include <math.h>
#include <cassert>
#include "TT06Tau.hh"
#include "pade.hh"
#include "mpiUtils.h"
#include "TT06Gates.h"
#include "TT06NonGates.h"
#include "ThreadServer.hh"

using namespace std;


namespace TT06Func
{
double Xr1Mhu(double Vm, void *parms) ;
double Xr2Mhu(double Vm, void *parms) ;
double XsMhu(double Vm, void *parms) ;
double mMhu(double Vm, void *parms) ;
double hjMhu(double Vm, void *parms) ;
double rMhu(double Vm, void *parms) ;
double dMhu(double Vm, void *parms) ;
double fMhu(double Vm, void *parms) ;
double f2Mhu(double Vm, void *parms) ;
double jLMhu(double Vm, void *parms) ;
double sMhu0(double Vm, void *parms) ;
double sMhu1(double Vm, void *parms) ;

double Xr1TauR(double Vm, void *parms) ;
double Xr2TauR(double Vm, void *parms) ;
double XsTauR(double Vm, void *parms) ;
double mTauR(double Vm, void *parms) ;
double hTauR(double Vm, void *parms) ;
double jTauR(double Vm, void *parms) ;
double rTauR(double Vm, void *parms) ;
double dTauR(double Vm, void *parms) ;
double fTauR(double Vm, void *parms) ;
double f2TauR(double Vm, void *parms) ;
double jLTauR(double Vm, void *parms) ;
double sTauR0(double Vm, void *parms) ;
double sTauR1(double Vm, void *parms) ;

double hTauRMod(double Vm, void *parms) ;
double jTauRMod(double Vm, void *parms) ;



static TauRecipParms *jParms;
static TauRecipParms *hParms;
}


namespace TT06Func
{
void initCnst()
{
   jParms =makeTauRecipParms(-48.85,-17.6,jTauRecip); 
   hParms =makeTauRecipParms(-64.20,-23.3,hTauRecip); 
   initNonGateCnst(); 
}
// Update Gates; 

void updateGate(double dt, int nCellsTotal, int s_switch, double *VM, int offset, double **gate, PADE *fit, WORK &work)
{
   
  int offsetCell = work.offsetCell; 
  int nCell      = work.nCell; 
  int offsetEq =  work.offsetEq; 
  int nEq =  work.nEq; 

  PADE *gatefit = fit + gateFitOffset; 

  void *mhuParms ; 
  void *tauRParms ; 
  double (*mhuFunc)(double V,void *parms); 
  double (*tauRFunc)(double V,void *parms);
  for (int eq=offsetEq;eq<offsetEq+nEq;eq++) 
  {
   double *g = gate[eq];
   switch (eq)
   {
     case  0:
     {
        mhuParms  = (gatefit[2*eq+0].aparms); 
        tauRParms = (gatefit[2*eq+1].aparms); 
        mhuFunc =gatefit[2*eq+0].afunc; 
        tauRFunc=gatefit[2*eq+1].afunc; 

        for (int ii=offsetCell;ii<offsetCell+nCell;ii++) 
        {
            double Vm = VM[ii]; 
            double mhu =mhuFunc(Vm,mhuParms); 
            double tauR=tauRFunc(Vm,tauRParms); 
            g[ii] +=   (mhu - g[ii])*(1-exp(-dt*tauR));                   //mtauR can be very large   ~1140.0 use Rush_Larson to integrate. 
        }
     }
     break; 
     case  1:
     case  2:
     case  3:
     case  4:
     case  5:
     case  6:
     case  7:
     case  8:
     case  9:
     case 10:
     {
        mhuParms  = (gatefit[2*eq].aparms); 
        tauRParms = (gatefit[2*eq+1].aparms); 
        mhuFunc   = (gatefit[2*eq].afunc); 
        tauRFunc  = (gatefit[2*eq+1].afunc); 
        for (int ii=offsetCell;ii<offsetCell+nCell;ii++) 
        {
          double Vm = VM[ii]; 
          double mhu=mhuFunc(Vm,mhuParms); 
          double tauR=tauRFunc(Vm,tauRParms); 
          g[ii] +=  dt*(mhu - g[ii])*tauR; 
        }
     }
     break; 
     case 11:
     {
        if (s_switch == 0)
        {
           mhuParms  = (gatefit[2*eq+0].aparms); 
           tauRParms = (gatefit[2*eq+1].aparms); 
           mhuFunc   = (gatefit[2*eq+0].afunc); 
           tauRFunc  = (gatefit[2*eq+1].afunc); 
           int ii=offsetCell; 
           for (;ii<offsetCell+nCell;ii++) 
           {
              double Vm = VM[ii]; 
              double mhu =mhuFunc (Vm,mhuParms); 
              double tauR=tauRFunc(Vm,tauRParms); 
              g[ii] +=  dt*(mhu - g[ii])*tauR;     //sGate
           }
        }
        else
        {
           mhuParms  = (gatefit[2*eq+2].aparms); 
           tauRParms = (gatefit[2*eq+3].aparms); 
           mhuFunc   = (gatefit[2*eq+2].afunc); 
           tauRFunc  = (gatefit[2*eq+3].afunc); 
           int ii=offsetCell; 
           for (;ii<offsetCell+nCell;ii++) 
           {
              double Vm = VM[ii]; 
              double mhu =mhuFunc (Vm,mhuParms); 
              double tauR=tauRFunc(Vm,tauRParms); 
              g[ii] +=  dt*(mhu - g[ii])*tauR;     //sGate
           }
        }
      }
      break; 
      default:
        assert(0); 
   }
  }
}
typedef void  (*UPDATEGATE)(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a) ; 
/*
UPDATEGATE updateGateFuncs0[]={ update_mGate, update_hGate, update_jGate, update_Xr1Gate, update_Xr2Gate, update_XsGate, update_rGate, update_dGate, 
                        update_fGate, update_f2Gate, update_jLGate, update_s0Gate, update_s1Gate} ;
UPDATEGATE updateGateFuncs1[]={ update_mGate_v1, update_hGate_v1, update_jGate_v1, update_Xr1Gate_v1, update_Xr2Gate_v1, update_XsGate_v1, update_rGate_v1, update_dGate_v1, 
                          update_fGate_v1, update_f2Gate_v1, update_jLGate_v1, update_s0Gate_v1, update_s1Gate_v1} ;
*/
void updateGateFast(double dt, int nCellsTotal, int s_switch, double *Vm, int offset, double **gate, PADE *fit,  WORK &work)
{
#if 0
     UPDATEGATE updateGateFuncs[]={ update_mGate_v1, update_hGate_v1, update_jGate_v1, update_Xr1Gate_v1, update_Xr2Gate_v1, update_XsGate_v1, 
                                    update_rGate_v1, update_dGate_v1, update_fGate_v1, update_f2Gate_v1,  update_jLGate_v1, 
                                    update_s0Gate_v1,update_s1Gate_v1} ;
#else
     UPDATEGATE updateGateFuncs[]={ update_mGate_v2, update_hGate_v2, update_jGate_v2, update_Xr1Gate_v2, update_Xr2Gate_v2, update_XsGate_v2, 
                                    update_rGate_v2, update_dGate_v2, update_fGate_v2, update_f2Gate_v2,  update_jLGate_v2, 
                                    update_s0Gate_v2,update_s1Gate_v2} ;
#endif
    int nCell=work.nCell; 
    if (nCell ==0) return; 
    int offsetCell=work.offsetCell; 
    int offsetEq = work.offsetEq; 
    int nEq = work.nEq; 

   PADE *gatefit = fit + gateFitOffset; 

   for (int eq=offsetEq;eq<offsetEq+nEq;eq++)
   {
       if (eq < 11) 
       {
            double *mhu  = gatefit[2*eq+0].coef; 
            double *tauR = gatefit[2*eq+1].coef; 
            updateGateFuncs[eq](dt, nCell , &Vm[offsetCell], gate[eq]+offsetCell, mhu, tauR);
       }
       else 
       {
          if (s_switch == 0)
          {
             double *mhu  = gatefit[2*eq+0].coef;
             double *tauR = gatefit[2*eq+1].coef;
             updateGateFuncs[11](dt, nCell , &Vm[offsetCell], gate[eq]+offsetCell, mhu, tauR);
          }
          else
          {
             double *mhu  = gatefit[2*eq+2].coef;
             double *tauR = gatefit[2*eq+3].coef;
             updateGateFuncs[12](dt, nCell , &Vm[offsetCell], gate[eq]+offsetCell, mhu, tauR);
          }
       }
    }
}
void updateGateFast0(double dt, int nCellsTotal, int *cellTypeVector, double *Vm, int offset, double **gate, PADE *fit,  WORK &work)
{
#if 0
     UPDATEGATE updateGateFuncs[]={ update_mGate_v1, update_hGate_v1, update_jGate_v1, update_Xr1Gate_v1, update_Xr2Gate_v1, update_XsGate_v1, 
                                    update_rGate_v1, update_dGate_v1, update_fGate_v1, update_f2Gate_v1,  update_jLGate_v1, update_s0Gate_v1} ;
#else
     UPDATEGATE updateGateFuncs[]={ update_mGate_v2, update_hGate_v2, update_jGate_v2, update_Xr1Gate_v2, update_Xr2Gate_v2, update_XsGate_v2, 
                                    update_rGate_v2, update_dGate_v2, update_fGate_v2, update_f2Gate_v2,  update_jLGate_v2, update_s0Gate_v2} ;
#endif
    int nCell=work.nCell; 
    if ( nCell ==0)  return; 
    int offsetCell=work.offsetCell; 
    int offsetEq = work.offsetEq; 
    int nEq = work.nEq; 

   PADE *gatefit = fit + gateFitOffset; 

   for (int eq=offsetEq;eq<offsetEq+nEq;eq++)
   {
         double *mhu  = gatefit[2*eq+0].coef; 
         double *tauR = gatefit[2*eq+1].coef; 
         updateGateFuncs[eq](dt, nCell , &Vm[offsetCell], gate[eq]+offsetCell, mhu, tauR);
   }
}
void updateGateFast1(double dt, int nCellsTotal, int *cellTypeVector, double *Vm, int offset, double **gate, PADE *fit,  WORK &work)
{
#if 0
     UPDATEGATE updateGateFuncs[]={ update_mGate_v1, update_hGate_v1, update_jGate_v1, update_Xr1Gate_v1, update_Xr2Gate_v1, update_XsGate_v1, 
                                    update_rGate_v1, update_dGate_v1, update_fGate_v1, update_f2Gate_v1,  update_jLGate_v1, update_s1Gate_v1} ;
#else
     UPDATEGATE updateGateFuncs[]={ update_mGate_v2, update_hGate_v2, update_jGate_v2, update_Xr1Gate_v2, update_Xr2Gate_v2, update_XsGate_v2, 
                                    update_rGate_v2, update_dGate_v2, update_fGate_v2, update_f2Gate_v2,  update_jLGate_v2, update_s1Gate_v2} ;
#endif
    int nCell=work.nCell; 
   if ( nCell ==0)  return; 
    int offsetCell=work.offsetCell; 
    int offsetEq = work.offsetEq; 
    int nEq = work.nEq; 

   PADE *gatefit = fit + gateFitOffset; 

   for (int eq=offsetEq;eq<offsetEq+nEq;eq++)
   {
         double *mhu  = gatefit[2*eq+0].coef; 
         double *tauR = gatefit[2*eq+1].coef; 
         updateGateFuncs[eq](dt, nCell , &Vm[offsetCell], gate[eq]+offsetCell, mhu, tauR);
   }
}
double Xr1Mhu(double Vm, void *) 
{ 
   double mhu=1.0/(1.0+(exp(((-26.0 - Vm)/7.0))));
   return mhu ; 
}
double Xr1TauR(double Vm, void *) 
{ 
   double t1 = (1.0+(exp(((-45.0 - Vm)/10.0))))/450;
   double t2 = (1.0+(exp(((Vm+30.0)/11.5000))))/6.0;
   double tauR =  t1*t2;
   return tauR;
}
double Xr2Mhu(double Vm, void *) 
{ 
   double mhu=1.0/(1.0+(exp(((Vm+88.0)/24.0))));
   return mhu;
}
double Xr2TauR(double Vm, void *) 
{ 
   double t1 = (1.0+(exp(((-60.0 - Vm)/20.0))))/3.0;
   double t2 = (1.0+(exp(((Vm - 60.0)/20.0))))/1.120;
   double tauR =  t1*t2;
   return tauR;
}
double XsMhu(double Vm, void *) 
{ 
   double mhu = 1.0/(1.0+(exp(((-5.0 - Vm)/14.0))));
   return mhu;
}
double XsTauR(double Vm, void *) 
{ 
   double t1 = 1400.00/sqrt(1.0+exp((5.0 - Vm)/6.0));
   double t2 = 1.0/(1.0+(exp(((Vm - 35.0)/15.0))));
   double tauR  =  1.0/(t1*t2+80.0);
   return tauR;
}
double mMhu(double Vm, void *) 
{ 
   double mhu = 1.0/SQ(1.0+exp((-56.8600 - Vm)/9.03000));
   return mhu;
}
double mTauR(double Vm, void *) 
{ 
   double t1  = 1.0/(1.0+(exp(((- 60.0 - Vm)/5.0))));
   double t2  =  0.10000/(1.0+(exp(((Vm+35.0)/5.0))))+0.100000/(1.0+(exp(((Vm - 50.0)/200.0))));
   double tauR =  1.0/(t1*t2);
   return tauR;
}
double hjMhu(double Vm, void *) 
{ 
   double mhu = 1.0/SQ((1.0+(exp(((Vm+71.5500)/7.43000)))));
   return mhu;
}
double hTauR(double Vm, void *) 
{ 
   double t1  = (Vm<- 40.0 ?  0.0570000*(exp((- (Vm+80.0)/6.80000))) : 0.0);
   double t2  = (Vm<- 40.0 ?  2.70000*(exp(( 0.0790000*Vm)))+ 310000.*(exp(( 0.348500*Vm))) : 0.770000/( 0.130000*(1.0+(exp(((Vm+10.6600)/- 11.1000))))));
   double tauR = (t1+t2);
   return tauR;
}
double hTauRMod(double Vm, void *) 
{ 
   double dtauR,ddtauR; 
   double tauR = TauRecipMod(Vm,hParms,&dtauR,&ddtauR); 
   return tauR;
}
double jTauR(double Vm, void *) 
{ 
   double t1  = (Vm < -40.0 ? (( ( - 25428.0*(exp(( 0.244400*Vm))) -  6.94800e-06*(exp(( - 0.0439100*Vm))))*(Vm+37.7800))/1.0)/(1.0+(exp(( 0.311000*(Vm+79.2300))))) : 0.0);
   double t2 = (Vm < -40.0 ? ( 0.0242400*(exp(( - 0.0105200*Vm))))/(1.0+(exp(( - 0.137800*(Vm+40.1400))))) : ( 0.600000*(exp(( 0.0570000*Vm))))/(1.0+(exp(( - 0.100000*(Vm+32.0))))));
   double tauR  = (t1+t2);
   return tauR;
}
double jTauRMod(double Vm, void *) 
{ 
   double dtauR,ddtauR; 
   double tauR = TauRecipMod(Vm,jParms,&dtauR,&ddtauR); 
   return tauR;
}
double rMhu(double Vm, void *) 
{ 
   double mhu = 1.0/(1.0+(exp(((20.0 - Vm)/6.0))));
   return mhu;
}
double rTauR(double Vm, void *) 
{ 
   double tau =  9.50000*(exp((- SQ((Vm+40.0)))/1800.00))+0.800000;
   double tauR = 1.0/tau; 
   return tauR;
}
double dMhu(double Vm, void *) 
{ 
   double mhu = 1.0/(1.0+(exp(((- 8.0 - Vm)/7.50000))));
   return mhu;
}
double dTauR(double Vm, void *) 
{ 
   double t1  = 1.40000/(1.0+(exp(((- 35.0 - Vm)/13.0))))+0.250000;
   double t2 = 1.40000/(1.0+(exp(((Vm+5.0)/5.0))));
   double t3 = 1.0/(1.0+(exp(((50.0 - Vm)/20.0))));
   double tauR =  1/(t1*t2+t3);
   return tauR;
}
double fMhu(double Vm, void *) 
{ 
   double mhu = 1.0/(1.0+(exp(((Vm+20.0)/7.0))));
   return mhu;
}
double fTauR(double Vm, void *) 
{ 
   double tau =  1102.50*(exp((- SQ(Vm+27.0)/225.0)))+200.0/(1.0+(exp(((13.0 - Vm)/10.0))))+180.0/(1.0+(exp(((Vm+30.0)/10.0))))+20.0;
   double tauR = 1/tau; 
   return tauR;
}
double f2Mhu(double Vm, void *) 
{ 
   double mhu = 0.670000/(1.0+(exp(((Vm+35.0)/7.0))))+0.330000;
   return mhu;
}
double f2TauR(double Vm, void *) 
{ 
   double tau =  562.0*exp(-SQ((Vm+27.0))/240.0)+31.0/(1.0+(exp(((25.0 - Vm)/10.0))))+80.0/(1.0+(exp(((Vm+30.0)/10.0))));
   double tauR = 1/tau; 
   return tauR;
}
double jLMhu(double Vm, void *) 
{ 
   double mhu =exp(-(Vm+91.0)*2/6.1);  // jL_gate
   return mhu;
}
double jLTauR(double Vm, void *) 
{ 
   double tauR = 1.0/670.0 ; 
   return tauR;
}
double sMhu0(double Vm, void *) 
{ 
   double mhu = 1.0/(1.0+(exp(((Vm+28.0)/5.0))));
   return mhu;
}
double sMhu1(double Vm, void *) 
{ 
   double  mhu=1.00000/(1.00000+(exp(((Vm+20.0000)/5.00000))));
   return mhu;
}
double sTauR0(double Vm, void *) 
{ 
   double   tau =  1000.0*(exp((-SQ(Vm+67.0)/1000.0)))+8.0;
   double tauR = 1/tau; 
   return tauR;
}
double sTauR1(double Vm, void *) 
{ 
   double tau =  85.0000*(exp((- (pow((Vm+45.0000), 2.00000))/320.000)))+5.00000/(1.00000+(exp(((Vm - 20.0000)/5.00000))))+3.00000;
   double tauR = 1/tau; 
   return tauR;
}
}
namespace TT06Func
{
map<string,CellTypeParmsFull>getStandardCellTypes()
{
   map<string,CellTypeParmsFull> cellTypes; 
   string name; 
   {
      int indices[] = {30,31,100} ;
      name ="endoCellML"; 
      cellTypes[name].name=name;
      cellTypes[name].s_switch=0; 
      cellTypes[name].P_NaK=2.724;
      cellTypes[name].g_Ks = 0.392; 
      cellTypes[name].g_Kr = 0.153; 
      cellTypes[name].g_to = 0.073; 
      cellTypes[name].g_NaL=0.0; 
      cellTypes[name].Vm = -86.709;
      cellTypes[name].state["K_i"]      = STATE(nonGateVar, K_i     , 138.4   );
      cellTypes[name].state["Na_i"]     = STATE(nonGateVar, Na_i    , 10.355  );
      cellTypes[name].state["Ca_i"]     = STATE(nonGateVar, Ca_i    , 0.00013 );
      cellTypes[name].state["Ca_ss"]    = STATE(nonGateVar, Ca_ss   , 0.00036 );
      cellTypes[name].state["Ca_SR"]    = STATE(nonGateVar, Ca_SR   , 3.715   );
      cellTypes[name].state["R_prime"]  = STATE(nonGateVar, R_prime , 0.9068  );
      cellTypes[name].state["fCass"]    = STATE(nonGateVar, fCass   , 0.9953  );
      cellTypes[name].state["Xr1_gate"] = STATE(   GateVar, Xr1_gate, 0.00448 );
      cellTypes[name].state["Xr2_gate"] = STATE(   GateVar, Xr2_gate, 0.476   );
      cellTypes[name].state["Xs_gate"]  = STATE(   GateVar, Xs_gate , 0.0087  );
      cellTypes[name].state["m_gate"]   = STATE(   GateVar, m_gate  , 0.00155 );
      cellTypes[name].state["h_gate"]   = STATE(   GateVar, h_gate  , 0.7573  );
      cellTypes[name].state["j_gate"]   = STATE(   GateVar, j_gate  , 0.7225  );
      cellTypes[name].state["r_gate"]   = STATE(   GateVar, r_gate  , 2.235e-8);
      cellTypes[name].state["d_gate"]   = STATE(   GateVar, d_gate  , 3.164e-5);
      cellTypes[name].state["f_gate"]   = STATE(   GateVar, f_gate  , 0.8009  );
      cellTypes[name].state["f2_gate"]  = STATE(   GateVar, f2_gate , 0.9778  );
      cellTypes[name].state["s_gate"]   = STATE(   GateVar, s_gate  , 0.3212  );
      cellTypes[name].state["jL_gate"]  = STATE(   GateVar, jL_gate , 0.066   );
   }
   {
      int indices[] = { 76,101} ;
      name="midCellML";
      cellTypes[name].name=name;
      cellTypes[name].s_switch=1; 
      cellTypes[name].P_NaK=2.724;
      cellTypes[name].g_Ks = 0.098; 
      cellTypes[name].g_Kr = 0.153; 
      cellTypes[name].g_to = 0.294; 
      cellTypes[name].g_NaL=0.0; 
      cellTypes[name].Vm = -85.423;
      cellTypes[name].state["K_i"]      = STATE(nonGateVar, K_i     , 138.52  );
      cellTypes[name].state["Na_i"]     = STATE(nonGateVar, Na_i    , 10.132  );
      cellTypes[name].state["Ca_i"]     = STATE(nonGateVar, Ca_i    , 0.000153);
      cellTypes[name].state["Ca_ss"]    = STATE(nonGateVar, Ca_ss   , 0.00042 );
      cellTypes[name].state["Ca_SR"]    = STATE(nonGateVar, Ca_SR   , 4.272   );
      cellTypes[name].state["R_prime"]  = STATE(nonGateVar, R_prime , 0.8978  );
      cellTypes[name].state["fCass"]    = STATE(nonGateVar, fCass   , 0.9942  );
      cellTypes[name].state["Xr1_gate"] = STATE(   GateVar, Xr1_gate, 0.0165  );
      cellTypes[name].state["Xr2_gate"] = STATE(   GateVar, Xr2_gate, 0.473   );
      cellTypes[name].state["Xs_gate"]  = STATE(   GateVar, Xs_gate , 0.0174  );
      cellTypes[name].state["m_gate"]   = STATE(   GateVar, m_gate  , 0.00165 );
      cellTypes[name].state["h_gate"]   = STATE(   GateVar, h_gate  , 0.749   );
      cellTypes[name].state["j_gate"]   = STATE(   GateVar, j_gate  , 0.6788  );
      cellTypes[name].state["r_gate"]   = STATE(   GateVar, r_gate  , 2.347e-8);
      cellTypes[name].state["d_gate"]   = STATE(   GateVar, d_gate  , 3.288e-5);
      cellTypes[name].state["f_gate"]   = STATE(   GateVar, f_gate  , 0.7026  );
      cellTypes[name].state["f2_gate"]  = STATE(   GateVar, f2_gate , 0.9526  );
      cellTypes[name].state["s_gate"]   = STATE(   GateVar, s_gate  , 0.999998);
      cellTypes[name].state["jL_gate"]  = STATE(   GateVar, jL_gate , 0.066   );
   }
   {
      int indices[] = {77,102} ;
      name="epiCellML";
      cellTypes[name].name=name;
      cellTypes[name].s_switch=1; 
      cellTypes[name].P_NaK=2.724;
      cellTypes[name].g_Ks = 0.392; 
      cellTypes[name].g_Kr = 0.153; 
      cellTypes[name].g_to = 0.294; 
      cellTypes[name].g_NaL=0.0; 
      cellTypes[name].Vm  = -85.23;
      cellTypes[name].state["K_i"]      = STATE(nonGateVar, K_i     , 136.89  );
      cellTypes[name].state["Na_i"]     = STATE(nonGateVar, Na_i    , 8.604   );
      cellTypes[name].state["Ca_i"]     = STATE(nonGateVar, Ca_i    , 0.000126);
      cellTypes[name].state["Ca_ss"]    = STATE(nonGateVar, Ca_ss   , 0.00036 );
      cellTypes[name].state["Ca_SR"]    = STATE(nonGateVar, Ca_SR   , 3.64    );
      cellTypes[name].state["R_prime"]  = STATE(nonGateVar, R_prime , 0.9073  );
      cellTypes[name].state["fCass"]    = STATE(nonGateVar, fCass   , 0.9953  );
      cellTypes[name].state["Xr1_gate"] = STATE(   GateVar, Xr1_gate, 0.00621 );
      cellTypes[name].state["Xr2_gate"] = STATE(   GateVar, Xr2_gate, 0.4712  );
      cellTypes[name].state["Xs_gate"]  = STATE(   GateVar, Xs_gate , 0.0095  );
      cellTypes[name].state["m_gate"]   = STATE(   GateVar, m_gate  , 0.00172 );
      cellTypes[name].state["h_gate"]   = STATE(   GateVar, h_gate  , 0.7444  );
      cellTypes[name].state["j_gate"]   = STATE(   GateVar, j_gate  , 0.7045  );
      cellTypes[name].state["r_gate"]   = STATE(   GateVar, r_gate  , 2.42e-8 );
      cellTypes[name].state["d_gate"]   = STATE(   GateVar, d_gate  , 3.373e-5);
      cellTypes[name].state["f_gate"]   = STATE(   GateVar, f_gate  , 0.7888  );
      cellTypes[name].state["f2_gate"]  = STATE(   GateVar, f2_gate , 0.9755  );
      cellTypes[name].state["s_gate"]   = STATE(   GateVar, s_gate  , 0.999998);
      cellTypes[name].state["jL_gate"]  = STATE(   GateVar, jL_gate , 0.066   );
    }
   {
      int indices[] = {30,31,100} ;
      name ="endoRRG"; 
      cellTypes[name].name=name;
      cellTypes[name].s_switch=0; 
      cellTypes[name].P_NaK=3.000;
      cellTypes[name].g_Ks = 0.392; 
      cellTypes[name].g_Kr = 0.153; 
      cellTypes[name].g_to = 0.073; 
      cellTypes[name].g_NaL=0.15; 
      cellTypes[name].Vm = -86.709;
      cellTypes[name].state["K_i"]      = STATE(nonGateVar, K_i     , 138.4   );
      cellTypes[name].state["Na_i"]     = STATE(nonGateVar, Na_i    , 10.355  );
      cellTypes[name].state["Ca_i"]     = STATE(nonGateVar, Ca_i    , 0.00013 );
      cellTypes[name].state["Ca_ss"]    = STATE(nonGateVar, Ca_ss   , 0.00036 );
      cellTypes[name].state["Ca_SR"]    = STATE(nonGateVar, Ca_SR   , 3.715   );
      cellTypes[name].state["R_prime"]  = STATE(nonGateVar, R_prime , 0.9068  );
      cellTypes[name].state["fCass"]    = STATE(nonGateVar, fCass   , 0.9953  );
      cellTypes[name].state["Xr1_gate"] = STATE(   GateVar, Xr1_gate, 0.00448 );
      cellTypes[name].state["Xr2_gate"] = STATE(   GateVar, Xr2_gate, 0.476   );
      cellTypes[name].state["Xs_gate"]  = STATE(   GateVar, Xs_gate , 0.0087  );
      cellTypes[name].state["m_gate"]   = STATE(   GateVar, m_gate  , 0.00155 );
      cellTypes[name].state["h_gate"]   = STATE(   GateVar, h_gate  , 0.7573  );
      cellTypes[name].state["j_gate"]   = STATE(   GateVar, j_gate  , 0.7225  );
      cellTypes[name].state["r_gate"]   = STATE(   GateVar, r_gate  , 2.235e-8);
      cellTypes[name].state["d_gate"]   = STATE(   GateVar, d_gate  , 3.164e-5);
      cellTypes[name].state["f_gate"]   = STATE(   GateVar, f_gate  , 0.8009  );
      cellTypes[name].state["f2_gate"]  = STATE(   GateVar, f2_gate , 0.9778  );
      cellTypes[name].state["s_gate"]   = STATE(   GateVar, s_gate  , 0.3212  );
      cellTypes[name].state["jL_gate"]  = STATE(   GateVar, jL_gate , 0.066   );
   }
   {
      int indices[] = { 76,101} ;
      name="midRRG";
      cellTypes[name].name=name;
      cellTypes[name].s_switch=1; 
      cellTypes[name].P_NaK=3.100;
      cellTypes[name].g_Ks = 0.098; 
      cellTypes[name].g_Kr = 0.153; 
      cellTypes[name].g_to = 0.294; 
      cellTypes[name].g_NaL=0.3; 
      cellTypes[name].Vm = -85.423;
      cellTypes[name].state["K_i"]      = STATE(nonGateVar, K_i     , 138.52  );
      cellTypes[name].state["Na_i"]     = STATE(nonGateVar, Na_i    , 10.132  );
      cellTypes[name].state["Ca_i"]     = STATE(nonGateVar, Ca_i    , 0.000153);
      cellTypes[name].state["Ca_ss"]    = STATE(nonGateVar, Ca_ss   , 0.00042 );
      cellTypes[name].state["Ca_SR"]    = STATE(nonGateVar, Ca_SR   , 4.272   );
      cellTypes[name].state["R_prime"]  = STATE(nonGateVar, R_prime , 0.8978  );
      cellTypes[name].state["fCass"]    = STATE(nonGateVar, fCass   , 0.9942  );
      cellTypes[name].state["Xr1_gate"] = STATE(   GateVar, Xr1_gate, 0.0165  );
      cellTypes[name].state["Xr2_gate"] = STATE(   GateVar, Xr2_gate, 0.473   );
      cellTypes[name].state["Xs_gate"]  = STATE(   GateVar, Xs_gate , 0.0174  );
      cellTypes[name].state["m_gate"]   = STATE(   GateVar, m_gate  , 0.00165 );
      cellTypes[name].state["h_gate"]   = STATE(   GateVar, h_gate  , 0.749   );
      cellTypes[name].state["j_gate"]   = STATE(   GateVar, j_gate  , 0.6788  );
      cellTypes[name].state["r_gate"]   = STATE(   GateVar, r_gate  , 2.347e-8);
      cellTypes[name].state["d_gate"]   = STATE(   GateVar, d_gate  , 3.288e-5);
      cellTypes[name].state["f_gate"]   = STATE(   GateVar, f_gate  , 0.7026  );
      cellTypes[name].state["f2_gate"]  = STATE(   GateVar, f2_gate , 0.9526  );
      cellTypes[name].state["s_gate"]   = STATE(   GateVar, s_gate  , 0.999998);
      cellTypes[name].state["jL_gate"]  = STATE(   GateVar, jL_gate , 0.066   );
   }
   {
      int indices[] = {77,102} ;
      name="epiRRG";
      cellTypes[name].name=name;
      cellTypes[name].s_switch=1; 
      cellTypes[name].P_NaK=3.000;
      cellTypes[name].g_Ks = 0.392; 
      cellTypes[name].g_Kr = 0.153; 
      cellTypes[name].g_to = 0.294; 
      cellTypes[name].g_NaL=0.15; 
      cellTypes[name].Vm  = -85.23;
      cellTypes[name].state["K_i"]      = STATE(nonGateVar, K_i     , 136.89  );
      cellTypes[name].state["Na_i"]     = STATE(nonGateVar, Na_i    , 8.604   );
      cellTypes[name].state["Ca_i"]     = STATE(nonGateVar, Ca_i    , 0.000126);
      cellTypes[name].state["Ca_ss"]    = STATE(nonGateVar, Ca_ss   , 0.00036 );
      cellTypes[name].state["Ca_SR"]    = STATE(nonGateVar, Ca_SR   , 3.64    );
      cellTypes[name].state["R_prime"]  = STATE(nonGateVar, R_prime , 0.9073  );
      cellTypes[name].state["fCass"]    = STATE(nonGateVar, fCass   , 0.9953  );
      cellTypes[name].state["Xr1_gate"] = STATE(   GateVar, Xr1_gate, 0.00621 );
      cellTypes[name].state["Xr2_gate"] = STATE(   GateVar, Xr2_gate, 0.4712  );
      cellTypes[name].state["Xs_gate"]  = STATE(   GateVar, Xs_gate , 0.0095  );
      cellTypes[name].state["m_gate"]   = STATE(   GateVar, m_gate  , 0.00172 );
      cellTypes[name].state["h_gate"]   = STATE(   GateVar, h_gate  , 0.7444  );
      cellTypes[name].state["j_gate"]   = STATE(   GateVar, j_gate  , 0.7045  );
      cellTypes[name].state["r_gate"]   = STATE(   GateVar, r_gate  , 2.42e-8 );
      cellTypes[name].state["d_gate"]   = STATE(   GateVar, d_gate  , 3.373e-5);
      cellTypes[name].state["f_gate"]   = STATE(   GateVar, f_gate  , 0.7888  );
      cellTypes[name].state["f2_gate"]  = STATE(   GateVar, f2_gate , 0.9755  );
      cellTypes[name].state["s_gate"]   = STATE(   GateVar, s_gate  , 0.999998);
      cellTypes[name].state["jL_gate"]  = STATE(   GateVar, jL_gate , 0.066   );
    }
    return cellTypes; 
}
}
void fv05General(void *fitIn, double Vm, double *fv)
{
   PADE *fit = (PADE *)fitIn;
   fv[0]=fit[0].afunc(Vm, fit[0].aparms); 
   fv[1]=fit[1].afunc(Vm, fit[1].aparms); 
   fv[2]=fit[2].afunc(Vm, fit[2].aparms); 
   fv[3]=fit[3].afunc(Vm, fit[3].aparms); 
   fv[4]=fit[4].afunc(Vm, fit[4].aparms); 
   fv[5]=fit[5].afunc(Vm, fit[5].aparms); 
}
double fv6General(void *fitIn, double dv)
{
   PADE *fit = (PADE *)fitIn;
   double fv6=fit[6].afunc(dv, fit[6].aparms); 
   return fv6; 
}
