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
   for (int j=k;j>=m;j--)  sum2 =  a[j] + x*sum2; 
   return sum1/sum2; 
}
double polyFunc(double x, PADE *pade) 
{
   int m=pade->m; 
   double *a=pade->coef; 
   double sum1=0; 
   for (int j=m-1;j>=0;j--)sum1 =  a[j] + x*sum1; 
   return sum1; 
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
   gsl_matrix *XX=NULL, *cov=NULL;
   gsl_vector *yy=NULL, *w=NULL, *cc=NULL;
   int k = m+l-1 ; 
      
   XX = gsl_matrix_alloc(n, k);
   yy= gsl_vector_alloc(n);
   w = gsl_vector_alloc(n);
        
   cc = gsl_vector_alloc(k);
   cov = gsl_matrix_alloc(k, k);
         
   for (int i = 0; i < n; i++)
   {
      double xj=1.0; 
      double xi = x[i]; 
      for (int j=0;j<m;j++) { gsl_matrix_set(XX, i, j, xj); xj *= xi; }
      xj=y[i]*xi; 
      for (int j=m;j<k;j++) { gsl_matrix_set(XX, i, j, xj); xj *= xi; }
      gsl_vector_set (yy, i, y[i]);
      gsl_vector_set (w, i, 1.0);
   }
   for (int i = 0; i < n; i++) gsl_vector_set(w, i, 1.0);
   double chisq; 
   {
      gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, k);
      gsl_multifit_wlinear(XX, w, yy, cc, cov, &chisq, work);
      gsl_multifit_linear_free (work);
   }
   //for (int j=0;j<m;j++) a[j] = C(j);
   for (int j=0;j<m;j++) a[j] = gsl_vector_get(cc,(j));
   a[m]=1.0; 
   //for (int j=m;j<k;j++) a[j+1] =-C(j);
   for (int j=m;j<k;j++) a[j+1] =-1.0*gsl_vector_get(cc,(j));
   gsl_matrix_free(XX);
   gsl_vector_free(yy);
   gsl_vector_free(w);
   gsl_vector_free(cc);
   gsl_matrix_free(cov);
}
void makeFunctionTable(PADE *pade) 
{
   double deltaX = pade->deltaX; 
   int n = (int) ((pade->x1-pade->x0)/deltaX  + 1.0001); 
   double *x = pade->x=(double *)malloc(n*sizeof(double));
   double *y = pade->y=(double *)malloc(n*sizeof(double));
   double ymin=0.0,ymax=0.0; 
   for (int i = 0; i < n; i++)
   {
      x[i] = pade->x0 + deltaX * i ; 
      y[i]=pade->func(x[i], pade->parms); 
      if ( y[i] < ymin   || i==0) ymin    =  y[i]; 
      if ( y[i] > ymax   || i==0) ymax    =  y[i]; 
   }
   pade->n =n; 
   pade->ymin = ymin; 
   pade->ymax = ymax; 
}
void padeErrorInfo(PADE pade,int index) 
{
   
   double dy = pade.ymax-pade.ymin; 
   char filename[256]; 
   sprintf(filename,"func_tt06_%s",pade.name.c_str()); 
   FILE *file = fopen(filename,"w"); 
   fprintf(file,"#%10s: %3d %6d %6d %6d ",pade.name.c_str(),index,costFunc(pade.l,pade.m),pade.l,pade.m); fflush(stdout); 
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
static void  minimizeCost(PADE *pade,int maxCost, int lMax, int mMax)
{
   int lMin = 1; 
   int mMin = 1; 
   int minCost=0; 
   if (lMax < 0 && mMax < 0)  minCost = costFunc(lMax,mMax); 
   if (lMax < 0) { lMax *= -1; lMin=lMax; }
   if (mMax < 0) { mMax *= -1; mMin=mMax; }
   int n = pade->n; 
   double *x = pade->x; 
   double *y = pade->y; 
   double tol = pade->tol; 
   double dy = (pade->ymax-pade->ymin); 
   double norm = (fabs(pade->ymax)+fabs(pade->ymin)); 
   double errTol = tol*dy; 
   if ( errTol < 1e-14*norm)  errTol = tol*norm; 
   int lmin=0,mmin=0; 
   int length; 
   double amin[lMax+mMax]; 
   double errMaxMin = -1.0; 
   double errRMSMin=0.0 ; 
   
   for (int kk=minCost;kk<=maxCost;kk++) 
   {
        for (int l=lMin;l<=lMax;l+=1)
        for (int m=mMin;m<=mMax;m+=1)
        {
           if (costFunc(l,m) != kk) continue; 
           double a[l+m], errMax, errRMS; 
           findPadeApprox(l,m,n,x,y,a);
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
	fprintf(file,"%10s FIT { tol=%e; deltaX=%21.14e; x0= %21.14e; x1= %21.14e; l=%d;m=%d;",pade.name.c_str(),pade.tol, pade.deltaX, pade.x0,pade.x1,pade.l,pade.m);
	fprintf(file, " coef="); 
        for (int i=0;i<(pade.l+pade.m);i++) fprintf(file, "%23.16e ",pade.coef[i]); 
	fprintf(file, ";}\n"); 
}
void padeCalc(PADE *pade, int lMax, int mMax, int maxCost)
{
   pade->aparms = pade; 
   pade->afunc = (double (*)(double , void *))padeFunc; 
   if (pade->tol > 0)  
   {
      minimizeCost(pade, maxCost, lMax, mMax);
      if (pade->l < 2  )  pade->afunc = (double (*)(double , void *))polyFunc; 
   }
   else { pade->afunc = pade->func ; pade->aparms= pade->parms; }
}
void padeSet(PADE *pade, int lMax, int mMax, int maxCost)
{
   if (pade->tol > 0)  
   {
      pade->afunc = (double (*)(double , void *))padeFunc; 
      pade->aparms = pade; 
      if (pade->coef == NULL) minimizeCost(pade, maxCost, lMax, mMax);
      if (pade->l < 2  )  pade->afunc = (double (*)(double , void *))polyFunc; 
   }
   else 
   { 
       pade->afunc = pade->func ; 
       pade->aparms= pade->parms; 
   }
}
void padeApprox(PADE &pade, std::string name, double (*func)(double x, void *parms), void *parms, int size_parms,
                double deltaX, double x0, double x1,
                double tol, int lMax,int mMax,int maxCost,
                int l, int m, double *coef)
{
   pade.name = name; 
   pade.func = func; 
   pade.parms = NULL ; 
   if (size_parms > 0 && parms != NULL) 
   {
     pade.parms = malloc(size_parms); 
     for (int i=0;i<size_parms;i++) ((char *)pade.parms)[i] = ((char*)parms)[i]; 
   }
   pade.tol=tol; 
   pade.deltaX = deltaX; 
   pade.x0 = x0; 
   pade.x1 = x1; 

   pade.n = 0; 
   makeFunctionTable(&pade) ;

   pade.afunc = (double (*)(double , void *))NULL; 
   pade.aparms = NULL; 

   pade.l = 0; 
   pade.m = 0; 
   pade.coef = NULL; 
   pade.cost=0 ;
   pade.errMax=0.0; 
   pade.errRMS=0.0; 
   if (pade.tol > 0)  
   {
      if (coef == NULL) 
      {
          minimizeCost(&pade, maxCost, lMax, mMax);
      }
      else
      {
         pade.l = l; 
         pade.m = m; 
         pade.coef = (double *)malloc((l+m)*sizeof(double)); 
         for (int j=0;j<(l+m);j++) pade.coef[j] = coef[j]; 
      }
      pade.afunc = (double (*)(double , void *))padeFunc; 
      pade.aparms = &pade; 
      if (pade.l < 2  )  pade.afunc = (double (*)(double , void *))polyFunc; 
   }
   else 
   { 
       pade.afunc = pade.func ; 
       pade.aparms= pade.parms; 
   }

}
