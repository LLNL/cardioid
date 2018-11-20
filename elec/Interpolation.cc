

#include "Interpolation.hh"
#include "svd.h"
#include <set>
#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

#ifdef HAVE_LAPACK
extern "C" {
   int dgels_ ( char *trans , int *m , int *n , int *nrhs , double *a , int *lda , double *b , int *ldb , double *work , int *lwork , int *info );
}
#endif

static vector<double> hornFromCheby(const vector<double>& chebyCoeff,
                                    const double lb,
                                    const double ub)
{
   int inSize = chebyCoeff.size();

   /*
    * zCoeff = zCoeffFromCheby * chebyCoeff
    *
    * zCoeffFromCheby = [
    *  1  0 -1  0  1 ...
    *  0  1  0 -3  0 ...
    *  0  0  2  0 -8 ...
    *  0  0  0  4  0 ...
    *  0  0  0  0  8 ...
    *  .  .  .  .  .
    */
   vector<double> zCoeffFromCheby(inSize*inSize, 0);
   zCoeffFromCheby[0 + 0*inSize] = 1;
   if (1 < inSize)
   {
      zCoeffFromCheby[1 + 1*inSize] = 1;
   }
   
   for (int kk=2; kk<inSize; kk++)
   {
      for (int ll=0; ll<inSize; ll++)
      {
         double leadingTerm = 0;
         if (ll>0)
         {
            leadingTerm = 2*zCoeffFromCheby[ll-1 + (kk-1)*inSize];
         }
         zCoeffFromCheby[ll + kk*inSize] = leadingTerm - zCoeffFromCheby[ll + (kk-2)*inSize];
      }
   }

   vector<double> zCoeff(inSize);
   for (int ii=0; ii<inSize; ii++)
   {
      zCoeff[ii] = 0;
      for (int jj=0; jj<inSize; jj++)
      {
         zCoeff[ii] += zCoeffFromCheby[ii + jj*inSize]*chebyCoeff[jj];
      }
   }

   /*
    * x=lb -> z=-1 
    * x=ub -> z= 1
    * z = M*x + B
    * M = 2/(ub-lb)
    * B = -(ub+lb)/(ub-lb)
    *
    * z^2 = (M*x + B)*(M*x + B) = M^2*x^2 + 2*M*B*x + B^2
    * z^(i+1) = M*x*(z^i) + B*(z^i)
    *
    * xCoeff = xFromZ * zCoeff
    *
    * xFromZ = [
    * 1  B   B*B   B*B*B ...
    * 0  M 2*M*B 3*M*B*B ...
    * 0  0   M*M 3*M*M*B ...
    * 0  0     0   M*M*M ...
    * .  .     .       .
    */
   double MM=2/(ub-lb);
   double BB=-(ub+lb)/(ub-lb);

   vector<double> xCoeffFromZ(inSize*inSize, 0);
   xCoeffFromZ[0 + 0*inSize] = 1;
   for (int kk=1; kk<inSize; kk++)
   {
      for (int ll=0; ll<inSize; ll++)
      {
         double leadingTerm=0;
         if (ll>0)
         {
            leadingTerm = MM*xCoeffFromZ[ll-1 + (kk-1)*inSize];
         }
         xCoeffFromZ[ll + kk*inSize] = leadingTerm + BB*xCoeffFromZ[ll + (kk-1)*inSize];
      }
   }

   vector<double> xCoeff(inSize);
   for (int ii=0; ii<inSize; ii++)
   {
      xCoeff[ii] = 0;
      for (int jj=0; jj<inSize; jj++)
      {
         xCoeff[ii] += xCoeffFromZ[ii + jj*inSize]*zCoeff[jj];
      }
   }

   return xCoeff;
}

double Interpolation::create(const vector<double>& inputs,
                             const vector<double>& outputs,
                             const double tolerance,
                             const double rangeWindow)
{   
   double lb = inputs[0];
   double ub = inputs[inputs.size()-1];
   int inSize = inputs.size();
   vector<double> chebyPoly(MAX_TERM_COUNT*inSize);
   for (int ii=0; ii<inSize; ii++)
   {
      chebyPoly[ii + 0*inSize] = 1;
      chebyPoly[ii + 1*inSize] = (inputs[ii] - (ub+lb)/2) * 2/(ub-lb); //convert to z in [-1,1]
   }
   for (int jj=2; jj<MAX_TERM_COUNT; jj++)
   {
      for (int ii=0; ii<inSize; ii++)
      {
         chebyPoly[ii + jj*inSize] =
            2*chebyPoly[ii + 1*inSize]*chebyPoly[ii + (jj-1)*inSize]
            -
            chebyPoly[ii + (jj-2)*inSize]
            ;
      }
   }

   /* Fill up a vector of the dynamic range for the function.
    *
    * I'd like the tolerance to be a measure of the relative error.
    * To compute relative error, I need to have a measure of how "big"
    * the function is.
    *
    * Using relative error for each output point individually doesn't
    * work well, because some function evaluations might be zero or
    * near zero, which would lead to ridiculously stringent error
    * requirements.
    *
    * Using |max-min| for the function doesn't work for things like
    * exponential functions, where the extrema of the range are
    * completely ridiculous, leading to 100% loss of accuracy for the
    * majority of the range.
    *
    * Instead, I'm going to use a moving window, where I compute the
    * max-min for some subset of the range.  This solves the near zero
    * issue, but also allows for dynamic relative error fits for
    * different parts of the function.
    *
    */
   vector<double> range(inSize);
   if (0)
   {
      //range is the output for every point.
      for (int ii=0; ii<inSize; ii++)
      {
         range[ii] = fabs(outputs[ii]);
      }
   }
   else if (0)
   {
      //range is set within the whole domain
      double funcMax = -1e30;
      double funcMin =  1e30;
      for (int ii=0; ii<inSize; ii++)
      {
         funcMax = max(funcMax,outputs[ii]);
         funcMin = min(funcMin,outputs[ii]);
      }
      //handle trivial case, where function is a constant
      if (funcMax <= funcMin)
      {
         numNumer_ = 1;
         numDenom_ = 1;
         coeff_.resize(1);
         coeff_[0] = funcMax;
         return 0;
      }
      for (int ii=0; ii<inSize; ii++)
      {
         range[ii] = funcMax-funcMin;
      }
   }
   else
   {
      int filterSize = inSize*rangeWindow;
      //make sure filterSize is odd
      if (!(filterSize % 2)) { filterSize++; }
      if (filterSize > inSize)
      {
        filterSize = inSize;
        if (!(filterSize % 2)) { filterSize--; }
      }
      assert(filterSize <= inSize);
      assert(filterSize >= 2);
      multiset<double> window;
      for (int ii=0; ii<filterSize; ii++)
      {
         window.insert(outputs[ii]);
      }
      for (int ii=0; ii<inSize; ii++)
      {
         if (filterSize/2+1 <= ii && ii <= inSize-filterSize/2-1)
         {
            window.erase(window.find(outputs[ii-filterSize/2-1]));
            window.insert(outputs[ii+filterSize/2]);
         }
         range[ii] = *(window.rbegin())-*(window.begin());
      }     
   }   
   //check to make sure the range is reasonable.
   double rangeMax = -1e30;
   for (int ii=0; ii<inSize; ii++)
   {
      rangeMax = max(rangeMax,range[ii]);
   }
   double absMinRange = rangeMax*1e-7;
   for (int ii=0; ii<inSize; ii++)
   {
      range[ii] = max(absMinRange, range[ii]);
   }

   int totalMaxTerm = 2*MAX_TERM_COUNT;
   double bestError = 1e30;
   vector<double> bestCoeffs(totalMaxTerm);
   int bestNumer=-1;
   int bestDenom=-1;
   
   int cost = 2;
   while (bestError >= tolerance and cost <= 2*MAX_TERM_COUNT)
   {
      for (int denomTermCount=1; denomTermCount<cost; denomTermCount++)
      {
         //find the number of the numerators and denominators we're fitting here
         int numerTermCount = cost - denomTermCount;
         if (denomTermCount == 1)
         {
            //no division needed for a unitary numerator
            numerTermCount = cost;
         }
         if (numerTermCount > MAX_TERM_COUNT) { continue; }
         if (denomTermCount > MAX_TERM_COUNT) { continue; }
         
         int coeffCount = numerTermCount + denomTermCount - 1;

         //find the coefficient matrices
         vector<double> bbb = outputs;
         vector<double> currentCoeffs(coeffCount);
         #ifdef HAVE_LAPACK
         {
            vector<double> AAA(coeffCount*inSize);

            for (int nn=0; nn<numerTermCount; nn++)
            {
               for (int ii=0; ii<inSize; ii++)
               {
                  AAA[ii + nn*inSize] = chebyPoly[ii + nn*inSize];
               }
            }
            for (int dd=1; dd<denomTermCount; dd++)
            {
               for (int ii=0; ii<inSize; ii++)
               {
                  AAA[ii + (numerTermCount+dd-1)*inSize] = -outputs[ii]*chebyPoly[ii + dd*inSize];
               }
            }

            int info = -1;
            char trans = 'N';
            int lwork = -1;
            int one=1;
            double optimalWork;
            dgels_(&trans, &inSize, &coeffCount, &one, &AAA[0], &inSize, &bbb[0], &inSize, &optimalWork, &lwork, &info);
            assert(info == 0);
            info = -1;
            lwork = optimalWork;
            vector<double> work(lwork);
            dgels_(&trans, &inSize, &coeffCount, &one, &AAA[0], &inSize, &bbb[0], &inSize, &work[0], &lwork, &info);
            assert(info == 0);
            for (int cc=0; cc<coeffCount; cc++) currentCoeffs[cc] = bbb[cc];
         }
         #else
         assert(0 && "Disabling Jim's svd computation for now.");
         {
            vector<double> Abuffer(coeffCount*inSize);
            vector<double*> AAA(inSize);
            for (int ii=0; ii<inSize; ii++)
            {
               AAA[ii] = &Abuffer[ii*coeffCount];
            }
            
            for (int nn=0; nn<numerTermCount; nn++)
            {
               for (int ii=0; ii<inSize; ii++)
               {
                  AAA[ii][nn] = chebyPoly[ii + nn*inSize];
               }
            }
            for (int dd=1; dd<denomTermCount; dd++)
            {
               for (int ii=0; ii<inSize; ii++)
               {
                  AAA[ii][numerTermCount+dd-1] = -outputs[ii]*chebyPoly[ii + dd*inSize];
               }
            }
            svdLinearLeastSquares(inSize, coeffCount, &AAA[0], &bbb[0], &currentCoeffs[0]);
         }
         #endif

         //compute the approximate
         vector<double> approximate(inSize);
         for (int ii=0; ii<inSize; ii++)
         {
            approximate[ii] = 0;
            for (int nn=0; nn<numerTermCount; nn++)
            {
               approximate[ii] += chebyPoly[ii + nn*inSize]*currentCoeffs[nn];
            }
            double denom=1;
            for (int dd=1; dd<denomTermCount; dd++)
            {
               denom += chebyPoly[ii + dd*inSize]*currentCoeffs[numerTermCount+dd-1];
            }
            approximate[ii] /= denom;
         }

         //compute the error
         double currentError=0;
         if (0)
         {
            //maximum error
            for (int ii=0; ii<inSize; ii++)
            {
               currentError = max(currentError, fabs(outputs[ii]-approximate[ii])/range[ii]);
            }
         }
         else
         {
            /* l2 error, assuming uniform step
             *
             * Here, assume we're always integrating from 0-1 so the
             * domain doesn't throw things off.  This means that if
             * the domain of integration changes, we don't want to
             * have to change our tolerance.
             *
             * l2 = sqrt(int((output-approx)^2,lb,ub) / (ub-lb))
             * 
             * int(((output-approx)^2,lb,ub)) =
             *   sum(for i=1:n,(output_i-approx_i)^2 * (lb-ub)/n)
             *
             * therefore
             *
             * l2 = sqrt(
             *   sum(for i=1:n,(output_i-approx_i)^2) / n
             *   )
             */
                 
            for (int ii=0; ii<inSize; ii++)
            {
               double diff = outputs[ii]-approximate[ii];
               currentError += (diff/range[ii])*(diff/range[ii]);
            }
            
            currentError = sqrt(currentError/inSize);
         }

         //if this is the best answer we've found, save it.
         if (currentError < bestError)
         {

            bestNumer = numerTermCount;
            bestDenom = denomTermCount;
            bestCoeffs = currentCoeffs;
            bestError = currentError;
         }
      }
      cost++;
   }

   numNumer_ = bestNumer;
   numDenom_ = bestDenom;

   //convert chebychev coeff to horn coeff for better evaluation
   vector<double> chebyNumerCoeffs(numNumer_);
   vector<double> chebyDenomCoeffs(numDenom_);
   for (int nn=0; nn<numNumer_; nn++)
   {
      chebyNumerCoeffs[nn] = bestCoeffs[nn];
   }
   chebyDenomCoeffs[0] = 1;
   for (int dd=1; dd<numDenom_; dd++)
   {
      chebyDenomCoeffs[dd] = bestCoeffs[numNumer_+dd-1];
   }
   vector<double> hornNumer(hornFromCheby(chebyNumerCoeffs, lb, ub));
   vector<double> hornDenom(hornFromCheby(chebyDenomCoeffs, lb, ub));
      
   //save the final coefficients to be used with eval
   //remember to renormalize so that the denominator constant is always 1
   coeff_.resize(numNumer_+numDenom_-1);
   for (int nn=0; nn<numNumer_; nn++)
   {
      coeff_[nn] = hornNumer[nn]/hornDenom[0];
   }
   for (int dd=1; dd<numDenom_; dd++)
   {
      coeff_[numNumer_ + dd-1] = hornDenom[dd]/hornDenom[0];
   }

   return bestError;
}

