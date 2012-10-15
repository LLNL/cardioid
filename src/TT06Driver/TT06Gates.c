#include <math.h>
void update_mGateNew(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex=0;
int mhu_l=10;
int mhu_m=5;
int tauR_l=1;
int tauR_m=18;
int exp_l=16;

   
 double  tauRdt_a[tauR_m];
 double exp_a[32]; 
 exp_a[0]=1.0; 
 for (int i=1;i<=exp_l;i++) exp_a[i] = exp_a[i-1]/i;
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   double mhu= sum1/sum2;
   for (int i=exp_l;i>=0;i--) sum4 = exp_a[i]+sum4*sum3; 
   double tauRdt =   1.0 - 1.0/sum4 ;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_mGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex=0;
int mhu_l=10;
int mhu_m=5;
int tauR_l=1;
int tauR_m=18;

   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   double mhu= sum1/sum2;
   double tauRdt= 1.0-exp(-sum3);                   //tauR*dt can be very large tauR ~=  114000.0 . UseRush_Larson to integrate. 
   g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
#define recipApprox(x) 1.0/(x)
void update_hGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex=1;
int mhu_l=7 ;
int mhu_m= 9;
int tauR_l=11;
int tauR_m=11;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double v = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + v*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + v*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j]+v*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + v*sum4; 
#if (1) 
    double tauRdt= sum3/sum4; 
    double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
#else
   double x = sum2*sum4; 
   double f = recipApprox(x); 
   double a = sum1-g[ii]*sum2; 
   double b = sum3*f*(2.0-f *x) ;
   g[ii] +=  a*b; 
#endif
  }
}
// hGate and jGate have the same mhu function.   There could be some advantage of combining these two updates
void update_jGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int  gateIndex=2;
int mhu_l= 7 ;
int mhu_m= 9;
int tauR_l= 1;
int tauR_m=13;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   for (int j=mhu_m-1;j>=0 ;j--)sum1   =      mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0 ;j--)sum3  =      tauRdt_a[j] + x*sum3; 
   //for (int j=tauR_k  ;j>=tauR_m;j--)sum4 = tauR_a[j] + x*sum4; 
   double tauRdt= sum3; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_Xr1Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex=3;
int mhu_l= 5 ;
int mhu_m= 8;
int tauR_l= 1;
int tauR_m=13;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0   ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0   ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_Xr2Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 4;
int mhu_l= 1 ;
int mhu_m= 10;
int tauR_l= 1;
int tauR_m=10;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum3=0;
   for (int j=mhu_m-1;j>=0   ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=tauR_m-1;j>=0   ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   double tauRdt= sum3; 
   double mhu= sum1;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_XsGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 5;
int mhu_l= 5 ;
int mhu_m= 5;
int tauR_l= 6;
int tauR_m= 9;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_rGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 6;
int mhu_l= 5 ;
int mhu_m= 8;
int tauR_l=5 ;
int tauR_m=7;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_dGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 7;
int mhu_l= 5 ;
int mhu_m= 7;
int tauR_l=7 ;
int tauR_m=10;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_fGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 8;
int mhu_l= 5 ;
int mhu_m= 8;
int tauR_l=13;
int tauR_m=11;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_f2Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 9;
int mhu_l= 8 ;
int mhu_m= 5;
int tauR_l=11;
int tauR_m=12;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_jLGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 10;
int mhu_l= 6 ;
int mhu_m= 12;
int tauR_l=1;
int tauR_m=1;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   double tauRdt= tauRdt_a[0]; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}

void update_s0Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 11;
int mhu_l= 7 ;
int mhu_m= 8;
int tauR_l=6 ;
int tauR_m=7 ;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
void update_s1Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a)
{
int gateIndex= 11;
int mhu_l= 7 ;
int mhu_m= 8;
int tauR_l=11;
int tauR_m=11;
   
 double  tauRdt_a[tauR_m];
 for (int j=tauR_m-1;j>=0;j--)    tauRdt_a[j] = tauR_a[j]*dt;
 int mhu_k  = mhu_m+mhu_l-1;  ;
 int tauR_k = tauR_m+tauR_l-1;  
 for (int ii=0;ii<nCells;ii++)
 { 
   double x = VM[ii];  
   double sum1=0;
   double sum2=0; 
   double sum3=0;
   double sum4=0;
   for (int j=mhu_m-1;j>=0    ;j--)sum1   =   mhu_a[j] + x*sum1;
   for (int j=mhu_k  ;j>=mhu_m;j--)sum2   =   mhu_a[j] + x*sum2;
   for (int j=tauR_m-1;j>=0    ;j--)sum3  =  tauRdt_a[j] + x*sum3;
   for (int j=tauR_k  ;j>=tauR_m;j--)sum4 =  tauR_a[j] + x*sum4; 
   double tauRdt= sum3/sum4; 
   double mhu= sum1/sum2;
    g[ii] +=  mhu*tauRdt - g[ii]*tauRdt;
  }
}
