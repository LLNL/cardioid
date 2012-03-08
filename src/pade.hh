#ifndef FIT_H
#define FIT_H
#include <stdio.h>
typedef struct {  char *name; double (*func)(double x, void *parms); void *parms; double (*afunc)(double x, void *aparms); void *aparms; int n; double deltaX, x0,x1,ymin,ymax, tol, errMax, errRMS; double *x, *y ; int cost; int l; int m; double *coef;} PADE; 
PADE   *padeApprox (const char *name, double (*func)(double x, void *parms), void *parms, int size_parms,double tol, double deltaX, double x0, double x1);
void padeCalc(PADE *pade,int lMax, int  mMax, int maxCost);
void padeErrorInfo(PADE pade,int index);
void padeWrite(FILE *file,PADE pade);
double padeFunc(double x, PADE *parms) ;
double polyFunc(double x, PADE *parms) ;
#endif
