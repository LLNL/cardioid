#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "svd.h"
#define MAX(x,y) ((x)>(y) ? (x):(y))
#define SQ(x)   ((x)*(x))
static double length2d(double a, double b)
{
   double length = 0.0; 

   double at = fabs(a);
   double bt = fabs(b); 
   if (at > bt) { double ct = bt / at; length = at * sqrt(1.0 + ct * ct); }
   if (bt > at) { double ct = at / bt; length = bt * sqrt(1.0 + ct * ct); }
   return length;
}

/* accumulate the left-hand transformation */
void accumLeft(int m, int n, double **a, double *sigma, double **v)
{
   for (int i = n - 1; i >= 0; i--) 
   {
      int l = i + 1;
      double g = sigma[i];
      for (int j = l; j < n; j++) a[i][j] = 0.0;
      if (g) 
      {
         g = 1.0 / g;
         for (int j = l; j < n; j++) 
         {
            double s=0.0; 
            for (int k = l; k < m; k++) s += (a[k][i] * a[k][j]);
            double f = (s / a[i][i]) * g;
            for (int k = i; k < m; k++) a[k][j] += (f * a[k][i]);
         }
      }
      for (int j = i; j < m; j++) a[j][i] *= g;
      a[i][i]+=1.0;
   }
}
//multiply the right-hand transformations
void accumRight(int m, int n, double **a, double *sigma, double **v, double *vec)
{
   v[n-1][n-1] = 1.0;
   for (int i = n - 2; i >= 0; i--) 
   {
      int l = i+1; 
      double g = vec[l]; 
      if (g) 
      {
         for (int j = l; j < n; j++) v[j][i] = ((a[i][j] / a[i][l]) / g); /* double division to avoid underflow */
         for (int j = l; j < n; j++) 
         {
            double s=0.0; 
            for (int k = l; k < n; k++) s += (a[i][k] * v[k][j]);
            for (int k = l; k < n; k++) v[k][j] += (s * v[k][i]);
         }
      }
      for (int j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
      v[i][i] = 1.0;
   }
}
/* Householder reduction to bidiagonal form */
double  householder(int m, int n, double **a, double *sigma, double **v, double *vec)
{
   double anorm=0.0; 
   vec[0] = 0.0; 
   for (int i = 0; i < n; i++) 
   {
      /* left-hand reduction */
      int l = i + 1;
      double scale = 0.0;
      for (int k = i; k < m; k++) scale += fabs(a[k][i]);
      if (scale>0.0) 
      {
         double s =0.0; 
         for (int k = i; k < m; k++) 
         {
            a[k][i] = (a[k][i]/scale);
            s += (a[k][i] * a[k][i]);
         }
         double f = a[i][i];
         double g = -copysign(sqrt(s), f);
         double h = f * g - s;
         a[i][i] = (f - g);
         for (int j = l; j < n; j++) 
         {
            double s = 0.0; 
            for (int k = i; k < m; k++) s += (a[k][i] * a[k][j]);
            f = s / h;
            for (int k = i; k < m; k++) a[k][j] += (f * a[k][i]);
         }
         for (int k = i; k < m; k++) a[k][i] = (a[k][i]*scale);
         scale *= g; 
      }
      sigma[i] = scale;

      /* right-hand reduction */
      scale = 0.0;
      for (int k = l; k < n; k++) scale += fabs(a[i][k]);
      if (scale>0.0) 
      {
         double s =0.0; 
         for (int k = l; k < n; k++) 
         {
            a[i][k] = (a[i][k]/scale);
            s += (a[i][k] * a[i][k]);
         }
         double f = a[i][l];
         double g = -copysign(sqrt(s), f);
         double h = f * g - s;
         a[i][l] = (f - g);
         for (int k = l; k < n; k++) vec[k] = a[i][k] / h;
         for (int j = l; j < m; j++) 
         {
            double s = 0.0; 
            for (int k = l; k < n; k++) s += (a[j][k] * a[i][k]);
            for (int k = l; k < n; k++) a[j][k] += (s * vec[k]);
         }
         for (int k = l; k < n; k++) a[i][k] = (a[i][k]*scale);
         scale *=g; 
      }
      if (l<n) vec[l] = scale; 
      anorm = MAX(anorm, (fabs(sigma[i]) + fabs(vec[i])));
   }
   double eps = 8*DBL_EPSILON*anorm; 
   return eps; 
}
void svdQR(int k, int l, int m, int n, double **a, double *sigma, double **v, double *vec)
{
         /* shift from bottom 2 x 2 minor */
         double x = sigma[l];
         double y = sigma[k-1];
         double z = sigma[k];
         double g = vec[k-1];
         double h = vec[k];
         double f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
         g = length2d(f, 1.0);
         f = ((x - z) * (x + z) + h * ((y / (f + copysign(g, f))) - h)) / x;

         /* next QR transformation */
         double c = 1.0;
         double s = 1.0;
         for (int i = l+1; i <= k; i++) 
         {
            double g = vec[i];
            y = sigma[i];
            h = s * g;
            g = c * g;
            z = length2d(f, h);
            vec[i-1] = z;
            c = f / z;
            s = h / z;
            f = x * c + g * s;
            g = g * c - x * s;
            h = y * s;
            y = y * c;
            for (int j = 0; j < n; j++) 
            {
               x = v[j][i-1];
               z = v[j][i];
               v[j][i-1] = (x * c + z * s);
               v[j][i] = (z * c - x * s);
            }
            z = length2d(f, h);
            sigma[i-1] = z;
            if (z>0.0) 
            {
               z = 1.0 / z;
               c = f * z;
               s = h * z;
            }
            f = (c * g) + (s * y);
            x = (c * y) - (s * g);
            for (int j = 0; j < m; j++) 
            {
               y = a[j][i-1];
               z = a[j][i];
               a[j][i-1] = (y * c + z * s);
               a[j][i] = (z * c - y * s);
            }
         }
         vec[l] = 0.0;
         vec[k] = f;
         sigma[k] = x;
}
int  svdSplit(int k, int m,double **a, double *sigma, double *vec, double eps)
{
    int l; 
    for (l = k; l > 0; l--) // test for splitting 
    {                    
       if (fabs(vec[l]) < eps )  return l; 
       if (fabs(sigma[l-1]) < eps )break;  
    }
    if (fabs(vec[l]) < eps )  return l; 
    double c = 0.0;
    double s = 1.0;
    for (int i = l; i <= k; i++) 
    {
       double f = s * vec[i];
       if (fabs(f) < eps) 
       {
          double g = sigma[i];
          double h = length2d(f, g);
          sigma[i] = h; 
          h = 1.0 / h;
          c = g * h;
          s = (-f * h);
          for (int j = 0; j < m; j++) 
          {
             double y = a[j][l-1];
             double z = a[j][i];
             a[j][l-1] = (y * c + z * s);
             a[j][i] = (z * c - y * s);
          }
       }
    }
    return l; 
}
//Golub and Reinsch SVD algorithm. 
//Num. Math. 14, 403-420 (1970) by Golub 
//and Reinsch, Handbook for Auto. Comp., vol II-Linear 
//Algebra, 134-151 (1971)
int svd(int m, int n, double **a, double *sigma, double **u, double **v)
{
   assert (m>=n);  // if case m<n desired operate on a transpose. 
   for (int i=0;i<m;i++) for (int j=0;j<n;j++) u[i][j] = a[i][j]; 

   double vec[n]; 

// Householder reduction to bidiagonal form 
   double eps = householder(m, n, u, sigma, v,vec);
   accumRight(m, n, u, sigma, v,vec);
   accumLeft(m, n, u, sigma, v);

// diagonalize the bidiagonal form 
   for (int k = n - 1; k >= 0; k--) // loop over singular values 
   {              
      for (int loop = 0; loop <=30; loop++) // loop over allowed iterations
      {                         
         int l = svdSplit(k, m,u, sigma, vec, eps);
         if (l == k) break; // convergence 
         if (loop == 30) return(svdNotConverged);
         svdQR(k, l,  m, n, u, sigma, v, vec);
      }
      if (sigma[k]  < 0.0) /* make singular value nonnegative */
      {  
         sigma[k] *= -1.0;
         for (int j = 0; j < n; j++) v[j][k] *= -1.0;
      }
   }
   return(svdConverged);
}
void svdTest(int m, int n, double **a, double *sigma, double **u, double **v)
{
   double b[m][n]; 
   double diffU = 0; 
   for (int i=0;i<n;i++) 
   for (int j=0;j<n;j++) 
   {
       if (i==j) diffU-=1.0;   
       for (int k=0;k<m;k++) diffU += u[k][i]*u[k][j]; 
   }
   double diffV = 0; 
   for (int i=0;i<n;i++) 
   for (int j=0;j<n;j++) 
   {
       if (i==j) diffV-=1.0;   
       for (int k=0;k<n;k++) diffV += v[k][i]*v[k][j]; 
   }

   for (int i=0;i<m;i++) 
   for (int j=0;j<n;j++) 
   {
       b[i][j]=0.0; 
       for (int k=0;k<n;k++) b[i][j] += u[i][k]*sigma[k]*v[j][k]; 
   }
   //for (int i=0;i<m;i++) { for (int j=0;j<n;j++) printf("%12.8f ",b[i][j]); printf("\n"); }
   double diffA =0.0; 
   double sumA =0.0; 
   for (int i=0;i<m;i++)  for (int j=0;j<n;j++) 
   {
      diffA += fabs(b[i][j]-a[i][j]);
      sumA += fabs(a[i][j]); 
   }
   diffA /= sumA; 
   
   int rank=0; 
   for (rank=1;rank<n;rank++) if (sigma[rank] < DBL_EPSILON*sigma[0]) break; 
   double sigMin = sigma[0]; 
   double sigMax = sigma[0]; 
   for (int i=1;i<n;i++) 
   {
	if (sigma[i] < sigMin) sigMin = sigma[i]; 
	if (sigma[i] > sigMax) sigMax = sigma[i]; 
   }
   double conditionNumber = sigMax/sigMin; 
   printf("%dx%d  numerical Rank = %d conditionNunber = %e diffA = %e orthognality U=%e orthongality V = %e \n",m,n,rank,conditionNumber,diffA,diffU,diffV); 
}
void svdLinearLeastSquares(int m, int n, double **A, double *y, double *x) 
{
   double sigma[n]; 
   double Ubuffer[m*n]; 
   double Vbuffer[n*n]; 
   double  *U[m]; 
   double  *V[n]; 
   for(int i=0;i<m;i++)   U[i] = Ubuffer+i*n; 
   for(int i=0;i<n;i++)   V[i] = Vbuffer+i*n; 
   double xx[n]; 
   for (int i=0;i<n;i++) xx[i]=x[i]=1.0; 
   double diff=0.0; 
   double diffLast=0.0;
   int loop =0; 
   for (loop = 0; loop<16;loop++)
   {
      for (int i=0;i<m;i++) for (int j=0;j<n;j++) A[i][j] *=xx[j];
      int rc = svd(m,n, A,sigma,U,V);
      if (rc > 0) fprintf(stderr,"SVD column not converged\n");
      assert(rc == 0);
      double rr[n]; 
      for (int i=0;i<n;i++)
      {
         rr[i] =0.0; 
         for (int j=0;j<m;j++) rr[i]+= U[j][i]*y[j];     
      }
      diff=0.0; 
      for (int i=0;i<n;i++)
      {
         xx[i] = 0.0; 
         for (int j=0;j<n;j++) xx[i]+=V[i][j]*rr[j]/sigma[j];
         diff += SQ(xx[i]-1.0); 
      }
      if (loop > 6  && diffLast < 0.9*diff) break; 
      for (int i=0;i<n;i++)  x[i] *= xx[i]; 
      if (diff < n*1e-22) break; 
      diffLast = diff; 
   }
}


