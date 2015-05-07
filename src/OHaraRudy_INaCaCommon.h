#include <math.h>
static double kNa1 = 15; //mM
static double kNa2 =  5; //mM
static double kNa3 =  88.12; //mM
static double kasymm =  12.5; //mM
static double wNa = 6e4 ; // Hz
static double wCa = 6e4 ; // Hz
static double wNaCa = 5e3; //Hz
static double kCaOn = 1.5e6 ; //mM/ms
static double kCaOff=5e3; // Hz
static double qNa = 0.5224; 
static double qCa = 0.1670; 
static double FRT=-1; 

static inline double INaCaYFunc( double V,double CaY, double NaY  )
{
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
