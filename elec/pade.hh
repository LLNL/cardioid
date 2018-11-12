#ifndef FIT_H
#define FIT_H
#include <stdio.h>
#include <string>
struct PADE
{  
   std::string name; 
   double (*func)(double x, void *parms); 
   void *parms; 
   double (*afunc)(double x, void *aparms); 
   void *aparms; 
   int n; 
   double deltaX, x0,x1,ymin,ymax, tol, errMax, errRMS; 
   double *x, *y ; 
   int cost; 
   int l; 
   int m; 
   double *coef;
}; 
void padeApprox (PADE &pade, std::string name, double (*func)(double x, void *parms), void *parms, int size_parms,
                   double deltaX, double x0, double x1,
                   double tol, int lMax,int mMax,int maxCost,
                   int l, int m, double *coef);
void padeErrorInfo(PADE pade,int index);
void padeWrite(FILE *file,PADE pade);
#endif
