#include "pade.hh"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gsl.h"
#define C(i) (gsl_vector_get(c,(i)))

int costFunc(int l, int m)
{
    //if (l != m  ) -1; 
  //if (l > 1 ) -1 ; 
   int cost = (m-1)+(l-1)+4;
   if ( l==1 ) cost = (m-1) ; 
	return cost; 
}
double padeFunc(double x, PADE *pade) 
{
   int l=pade->l; 
   int m=pade->m; 
   double *a=pade->coef; 
   double sum1=0; 
   double sum2=0;    
   int k = m+l-1; 
   for (int j=m-1;j>=0;j--)sum1 =  a[j] + x*sum1; 
   if (l<2) return sum1; 
   for (int j=k;j>=m;j--)  sum2 =  a[j] + x*sum2; 
   return sum1/sum2; 
}

static void padeError(int l,int m,double *a,int n,double *x,double *y,double *errMax, double *errRMS)
{
   double eMax = 0.0; 
   double err2 = 0.0; 
   PADE pade; 
   pade.l=l; 
   pade.m=m; 
   pade.coef=a; 
   for (int i = 0; i<n;i++) 
   {
     double f = padeFunc(x[i],&pade); 
     double err = fabs(f-y[i]); 
     if (err > eMax) eMax = err; 
     err2 += err*err; 
   }
   *errMax = eMax; 
   *errRMS = sqrt(err2/n); 
}
void  findPadeApprox(int l, int m, int n, double *x, double *y, double *a )
{
   gsl_matrix *X=NULL, *cov=NULL;
   gsl_vector *yy=NULL, *w=NULL, *c=NULL;
   int k = m+l-1 ; 
      
    X = gsl_matrix_alloc(n, k);
    yy= gsl_vector_alloc(n);
    w = gsl_vector_alloc(n);
        
     c = gsl_vector_alloc(k);
     cov = gsl_matrix_alloc(k, k);
         
     for (int i = 0; i < n; i++)
     {
           double xj=1.0; 
           for (int j=0;j<m;j++) { gsl_matrix_set (X, i, j, xj); xj *= x[i]; }
           xj=y[i]*x[i]; 
           for (int j=m;j<k;j++) { gsl_matrix_set (X, i, j, xj); xj *= x[i]; }
           gsl_vector_set (yy, i, y[i]);
           gsl_vector_set (w, i, 1.0);
      }
      for (int i = 0; i < n; i++) gsl_vector_set(w, i, 1.0);
      double chisq; 
      {
          gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, k);
          gsl_multifit_wlinear (X, w, yy, c, cov, &chisq, work);
          gsl_multifit_linear_free (work);
       }
       for (int j=0;j<m;j++) a[j] = C(j);
       a[m]=1.0; 
       for (int j=m;j<k;j++) a[j+1] =-C(j);
      gsl_matrix_free (X);
      gsl_vector_free (yy);
      gsl_vector_free (w);
      gsl_vector_free (c);
      gsl_matrix_free (cov);
}
static void makeFunctionTable(PADE *pade) 
{
   double deltaX = pade->deltaX; 
   int nmax = (int) ((pade->x1-pade->x0)/deltaX  + 1.0001); 
   double x[nmax+1],y[nmax+1]; 
   double ymin=0.0,ymax=0.0; 
   int n =0; 
   for (int i = 0; i < nmax; i++)
   {
      double xn = pade->x0 + deltaX * i ; 
      if (  -50.0 < xn && xn < -30.0 ) continue; 
      x[n]=xn; 
      y[n]=pade->func(x[n], pade->parms); 

      if ( y[n] < ymin   || i==0) ymin    =  y[n]; 
      if ( y[n] > ymax   || i==0) ymax    =  y[n]; 
      n++; 
   }
   pade->x=(double *)malloc(n*sizeof(double));
   pade->y=(double *)malloc(n*sizeof(double));
   for (int i=0;i<n;i++) { pade->x[i]=x[i]; pade->y[i]=y[i];} 
   pade->n =n; 
   pade->ymin = ymin; 
   pade->ymax = ymax; 
}
void padeErrorInfo(PADE pade) 
{
   
   double dy = pade.ymax-pade.ymin; 
   char filename[256]; 
   sprintf(filename,"func_tt06_%s",pade.name); 
   FILE *file = fopen(filename,"w"); 
   fprintf(file,"#%10s: %6d %6d %6d ",pade.name,costFunc(pade.l,pade.m),pade.l,pade.m); fflush(stdout); 
   fprintf(file,"%10.2e %10.2e %10.2e ",pade.tol,pade.errMax,pade.errMax/dy); 
   fprintf(file,"%10.2e %10.2e\n",pade.errRMS,pade.errRMS/dy); 
   for (int i=0;i<pade.n;i++) 
   {
      double x = pade.x[i]; 
      double f = padeFunc(x, &pade); 
      double y=pade.func(x, pade.parms); 
      double err = f-y; 
      fprintf(file,"%e %e %e %e\n",x,y,f,err); fflush(file); 
   }
   fclose(file); 
}
static void  minimizeCost(PADE *pade,double tol, int maxCost, int lMax, int mMax)
{
   int n = pade->n; 
   double *x = pade->x; 
   double *y = pade->y; 
   pade->tol = tol; 
   double errTol = tol*(pade->ymax-pade->ymin); 
   int lmin=0,mmin=0; 
   int length; 
   double amin[lMax+mMax]; 
   double errMaxMin = -1.0; 
   double errRMSMin ; 
   for (int kk=4;kk<maxCost;kk++) 
   {
        for (int l=1;l<=lMax;l+=1)
        for (int m=1;m<=mMax;m+=1)
        {
           if (costFunc(l,m) != kk) continue; 
           double a[l+m], errMax, errRMS; 
	   findPadeApprox(l,m,n,x, y,a);
           padeError(l,m,a,n,x,y,&errMax,&errRMS);
           if (errMax < errMaxMin || errMaxMin==-1.0) 
           {
              lmin=l;
              mmin=m;
              for (int j=0;j<l+m;j++) amin[j] = a[j];
              errMaxMin = errMax; 
              errRMSMin = errRMS; 
           }
       }
       if (errMaxMin  < errTol && errMaxMin >= 0.0 ) break ; 
   }
   length = lmin+mmin; 
   pade->l = lmin; 
   pade->m = mmin; 
   pade->errMax=errMaxMin; 
   pade->errRMS=errRMSMin; 
   pade->coef = (double *)malloc(length*sizeof(double)); 
   for (int i =0;i<length;i++) pade->coef[i]=amin[i]; 
   pade->cost=costFunc(lmin,mmin); 
}
void padeWrite(FILE *file,PADE pade)
{
	fprintf(file,"%10s PADE { l=%d;m=%d;",pade.name,pade.l,pade.m);
	fprintf(file, " coef="); 
        for (int i=0;i<(pade.l+pade.m);i++) fprintf(file, "%21.14e ",pade.coef[i]); 
	fprintf(file, ";}\n"); 
}
PADE   *padeApprox (char *name, double (*func)(double x, void *parms), void *parms, int size_parms,double tol, double deltaX, double x0, double x1, int maxOrder, int maxCost)
{
   PADE *pade=(PADE *)malloc(sizeof(PADE)); 
   pade->name = strdup(name); 
   pade->func = func; 
   pade->parms = NULL ; 
   pade->afunc = (double (*)(double , void *))padeFunc; 
   pade->aparms = pade; 
   if (size_parms > 0 && parms != NULL) 
   {
     pade->parms = malloc(size_parms); 
     for (int i=0;i<size_parms;i++) ((char *)pade->parms)[i] = ((char*)parms)[i]; 
   }
   pade->deltaX = deltaX; 
   pade->n = 0; 
   pade->x0 = x0; 
   pade->x1 = x1; 
   pade->l = 0; 
   pade->m = 0; 
   pade->errMax=0.0; 
   pade->errRMS=0.0; 
   pade->coef = NULL; 
   pade->cost=0 ;
   makeFunctionTable(pade) ;
   if (tol > 0)  minimizeCost(pade, tol, maxCost, maxOrder+1, maxOrder+1);
   else { pade->afunc = func ; pade->aparms= pade->parms; }
   return pade; 
}
