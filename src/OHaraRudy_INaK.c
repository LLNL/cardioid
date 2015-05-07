// INaK  V,Nai,Ki::PNaK
#include <math.h>
#include "OHaraRudy.h"
#include "OHaraRudy_INaK.h"

static double k1Plus = 949.5;//Hz
static double k1Minus = 182.4;//mM^-1
static double k2Plus = 687.2;//Hz
static double k2Minus = 39.4 ;//Hz
static double k3Plus = 1899.;//Hz
static double k3Minus = 79300 ;//Hz*mM^-2
static double k4Plus = 639.0;//Hz
static double k4Minus = 40.0 ;//Hz
static double KoNai = 9.073 ; //mM
static double KoNao = 27.78 ; //mM
static double delta = -0.1550;  

static double FRT = -1; 

void OHaraRudy_INaKFunc(CELLPARMS *parmsPtr, STATE *state, int pOffset, DERIVED *derived, double dt)
{
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double V = state->Vm; 
   double Nai = state->Nai; 
   double Ki = state->Ki; 

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
   derived->I.NaK = cP->PNaK*(zNa*JNaKNa + zK*JNaKK); 
}
COMPONENTINFO OHaraRudy_INaKInit()
{
   if (FRT  < 0) FRT = F/(R*T);
   COMPONENTINFO info;
   info.nVar = nVar; 
   info.varInfo = varInfo;
   info.func = OHaraRudy_INaKFunc;
   info.access = OHaraRudy_INaKAccess;
   return info;
}

