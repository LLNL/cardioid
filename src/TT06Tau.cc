#include "TT06Tau.hh" 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SQ(x) ((x)*(x))

double mTauR(double Vm, void *) 
{ 
   double t1  = 1.0/(1.0+(exp(((- 60.0 - Vm)/5.0))));
   double t2  =  0.10000/(1.0+(exp(((Vm+35.0)/5.0))))+0.100000/(1.0+(exp(((Vm - 50.0)/200.0))));
   double tauR =  1.0/(t1*t2);
   return tauR;
}
double hTauRecip(double Vm, double *dtauRecip, double *ddtauRecip)
{
   double t1,t2; 
   double dt1,dt2; 
   //double t1_h  = (Vm<- 40.0 ?  0.0570000*(exp((- (Vm+80.0)/6.80000))) : 0.0);
   //double t2_h  = (Vm<- 40.0 ?  2.70000*(exp(( 0.0790000*Vm)))+ 310000.*(exp(( 0.348500*Vm))) : 0.770000/( 0.130000*(1.0+(exp(((Vm+10.6600)/- 11.1000))))));
   //double tau_h = 1.0/(t1_h+t2_h);
   if ( Vm < -40.0) 
   { 

   t1  = 0.0570000*exp(-(Vm+80.0)/6.80000) ;
   dt1 = -t1/6.8000; 
   t2  = 2.70000*exp(0.0790000*Vm)+ 310000.*exp(0.348500*Vm);
   dt2  =0.079* 2.70000*exp(0.0790000*Vm)+ 0.3485*310000.*exp(0.348500*Vm);
   } 
   else 
  {
   t1 = 0.0;    
   dt1 = 0.0;  
   double c =  (0.770000/0.130000);
   double e =  exp((Vm+10.6600)/-11.1000);
   double s = 1.0/(1.0+e);
   t2 =  c*s;
   dt2 = c*s*s*e/11.1;
   }
   *dtauRecip = dt1+dt2; 
   double  tauRecip = t1+t2; 
   return tauRecip;
}
double jTauRecip(double Vm, double *dtauRecip, double *ddtauRecip)
{
   double t1,t2; 
   double dt1,dt2; 
   double ddt1,ddt2; 
   if ( Vm < -40.0) 
   { 
   double a1  =  (-25428.0*exp(0.244400*Vm) -  6.94800e-06*exp(-0.0439100*Vm));
   double da1 = (-0.24400*25428.0*exp(0.244400*Vm) +  0.0439100*6.94800e-06*exp(-0.0439100*Vm));
   double dda1 = (-SQ(0.24400)*25428.0*exp(0.244400*Vm) -  SQ(0.0439100)*6.94800e-06*exp(-0.0439100*Vm));
   double a2  =(Vm+37.7800);
   double da2 = 1.0; 
   double dda2 = 0.0; 
   double c = 0.311; 
   double e =exp(c*(Vm+79.2300));
   double a3  = 1/(1.0+e) ;
   double da3 = -a3*a3*e*c;
   double dda3 = -2*da3*a3*e*c+da3*c;
   t1 = a1*a2*a3; 
   double q = da1/a1 + da2/a2 + da3/a3; 
   double dq = dda1/a1 -SQ(da1/a1) + dda2/a2 - SQ(da2/a2)+ dda3/a3 - SQ(da3/a3); 
   dt1 = q*t1; 
   ddt1 = dq *t1 + q * dt1; 
   //t1 = a3; dt1 = da3; ddt1 = dda3; 

   a1 = 0.0242400*exp(-0.0105200*Vm) ;
   da1 = -0.0105200*a1;
   dda1 = -0.0105200*da1;

   c = -0.137800; 
   e =exp(c*(Vm+40.1400)) ;
   a2 =1.0/(1.0+e) ;
   da2 = -c*a2*a2*e;
   dda2 = -2.0*c*e*a2*da2 + da2*c;

   t2 = a1*a2; 
   q = da1/a1 + da2/a2 ; 
   dq = dda1/a1 -SQ(da1/a1) + dda2/a2 - SQ(da2/a2); 
   dt2 = (da1/a1 + da2/a2)*t2; 
   ddt2 = dq *t2 + q * dt2; 

   //t1 = Vm*Vm*Vm; dt1 = 3*Vm*Vm; ddt1 = 6*Vm; 
   //t2 = 0.0; dt2 = 0.0; ddt2 =  0.0;
   //t1 = 0.0; dt1 = 0.0; ddt1 =  0.0;
   //t2 = a2; dt2 = da2; ddt2 = dda2; 

   } 
   else 
  {
   t1 = 0.0;    
   dt1 = 0.0;  
   ddt1 = 0.0;  
   double a1 = 0.6000000*exp(0.0570000*Vm) ;
   double da1 =  0.0570000*a1;
   double dda1 =  0.0570000*da1;
   double c = -0.1; 
   double e = exp(c*(Vm+32.0000));
   double a2 =1.0/(1.0+e);
   double da2 = -a2*a2*e*c;
   double dda2 = -2*da2*a2*e*c+da2*c;
   t2 = a1*a2; 
   double q = da1/a1 + da2/a2 ; 
   double dq = dda1/a1 -SQ(da1/a1) + dda2/a2 - SQ(da2/a2); 
   dt2 = (da1/a1 + da2/a2)*t2; 
   ddt2 = dq *t2 + q * dt2; 
  }
   *dtauRecip = dt1+dt2; 
   *ddtauRecip = ddt1+ddt2; 
   double  tauRecip = t1+t2; 
   return tauRecip;
}
void hermite(double x0, double p0, double dp0, double x1, double p1, double dp1, double *c)
{
    double dx = x1-x0; 
    c[0] = p0; 
    c[1] = dp0; 
    c[2] = 3*(p1-p0)/(dx*dx) -  (2*dp0+dp1)/dx; 
    c[3] =-2*(p1-p0)/(dx*dx*dx) + (dp0+dp1)/(dx*dx);
}
double hermitePoly(double x,double *c, double *d, double *dd)
{
        *d =   c[1] + x * (2*c[2] + x*3*c[3]); 
        *dd =  2*c[2] + x*6*c[3]; 
	return c[0] +  x*(c[1]+x*(c[2] + x*c[3])); 
}
TauRecipParms* makeTauRecipParms(double V0, double V1, double (*func)(double, double *, double *))
{
   TauRecipParms *parms = (TauRecipParms *) malloc(sizeof(TauRecipParms)); 
   double dR0,dR1; 
   double ddR0,ddR1; 
   double R0 = func(V0,&dR0,&ddR0); 
   double R1 = func(V1,&dR1,&ddR1); 
   hermite(V0,R0,dR0,V1,R1,dR1,parms->c);
   parms->x0=V0; 
   parms->x1=V1; 
   parms->func = func; 
   return parms; 
}
double TauRecipMod(double V, TauRecipParms *parms, double *dR, double *ddR)
{
   if (V < parms->x0 || V > parms->x1)  return  parms->func(V,dR,ddR); 
   return hermitePoly(V-parms->x0,parms->c,dR,ddR); 
}
#if 0 
int main()
{
   double tauR,dtauR,ddtauR;
   double a[] = { 1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0, 1.0/720.0, 1.0/5040.0}; 
   TauRecipParms *jParms =makeTauRecipParms(-48.85,-17.6,jTauRecip); 
   TauRecipParms *hParms =makeTauRecipParms(-64.20,-23.3,hTauRecip); 
   double Vm = -48.85; 
   //tauR = jTauRecip( Vm,  &dtauR,&ddtauR); printf("%e %e ",Vm,ddtauR); 
   //tauR = TauRecipMod(Vm,jParms,&dtauR,&ddtauR);  printf("%e\n",ddtauR); 
   Vm = -17.6; 
   //tauR = jTauRecip( Vm,  &dtauR,&ddtauR); printf("%e %e ",Vm,ddtauR); 
   //tauR = TauRecipMod(Vm,jParms,&dtauR,&ddtauR);  printf("%e\n",ddtauR); 
   for (double Vm=-91;Vm< 50;Vm+=1) 
   {
   //double tauRE = hTauRecip( Vm,  &dtauR,&ddtauR);
   //double tauRM = TauRecipMod(Vm,hParms,&dtauR,&ddtauR); 
   //printf("%16.12f %20.12e %20.12e %20.12e\n",Vm,1/tauRE,1/tauRM,ddtauR); 
   double tauR = mTauR(Vm,NULL); 
   double dt = 1.; 
   double dtTauR = dt*tauR; 
   double dtTauREff = 1-exp(-dtTauR); 
   double x = dtTauR/4; 
   double f =0.0; 
   double s=0; 
   for (int i=7;i>=0;i--) f = a[i]+s*x; 
   //double  f =   1+x*( 1 + x*(1.0/2 + x*(1.0/6 + x*(1.0/24+x*(1.0/120+x*0.0/720))))); 
   f *= f; 
   f *= f; 
   f =   1 - 1 / f ;
   //double f = (1-x/3) / (1+x*2/3.0+x*x/6.0); 
     //f = 1 -f ; 
   double err = dtTauR;
   
   printf("%16.12f %20.12e %20.12e %20.12e %e\n",Vm,tauR,dtTauREff,f,dtTauREff/f-1.0); 
   }
}
#endif
