#ifndef INTERPOLATION_HH
#define INTERPOLATION_HH

#include <vector>

#define MAX_TERM_COUNT 32

class Interpolation {
 public:
   template <typename TTT>
   inline TTT eval(const TTT xx)
   {
      double * coeffCursor = &coeff_[0];
      TTT numer(coeffCursor[numNumer_-1]);
      switch (numNumer_)
      {
        case 32: numer = coeffCursor[30] + xx*numer;
        case 31: numer = coeffCursor[29] + xx*numer;
        case 30: numer = coeffCursor[28] + xx*numer;
        case 29: numer = coeffCursor[27] + xx*numer;
        case 28: numer = coeffCursor[26] + xx*numer;
        case 27: numer = coeffCursor[25] + xx*numer;
        case 26: numer = coeffCursor[24] + xx*numer;
        case 25: numer = coeffCursor[23] + xx*numer;
        case 24: numer = coeffCursor[22] + xx*numer;
        case 23: numer = coeffCursor[21] + xx*numer;
        case 22: numer = coeffCursor[20] + xx*numer;
        case 21: numer = coeffCursor[19] + xx*numer;
        case 20: numer = coeffCursor[18] + xx*numer;
        case 19: numer = coeffCursor[17] + xx*numer;
        case 18: numer = coeffCursor[16] + xx*numer;
        case 17: numer = coeffCursor[15] + xx*numer;
        case 16: numer = coeffCursor[14] + xx*numer;
        case 15: numer = coeffCursor[13] + xx*numer;
        case 14: numer = coeffCursor[12] + xx*numer;
        case 13: numer = coeffCursor[11] + xx*numer;
        case 12: numer = coeffCursor[10] + xx*numer;
        case 11: numer = coeffCursor[ 9] + xx*numer;
        case 10: numer = coeffCursor[ 8] + xx*numer;
        case  9: numer = coeffCursor[ 7] + xx*numer;
        case  8: numer = coeffCursor[ 6] + xx*numer;
        case  7: numer = coeffCursor[ 5] + xx*numer;
        case  6: numer = coeffCursor[ 4] + xx*numer;
        case  5: numer = coeffCursor[ 3] + xx*numer;
        case  4: numer = coeffCursor[ 2] + xx*numer;
        case  3: numer = coeffCursor[ 1] + xx*numer;
        case  2: numer = coeffCursor[ 0] + xx*numer;
        default:
         ;
      }
      TTT result;
      if (numDenom_ == 1)
      {
         result = numer;
      }
      else
      {
         coeffCursor += numNumer_;
         TTT denom(coeffCursor[numDenom_-2]);
         switch (numDenom_)
         {
           case 32: denom = coeffCursor[29] + xx*denom;
           case 31: denom = coeffCursor[28] + xx*denom;
           case 30: denom = coeffCursor[27] + xx*denom;
           case 29: denom = coeffCursor[26] + xx*denom;
           case 28: denom = coeffCursor[25] + xx*denom;
           case 27: denom = coeffCursor[24] + xx*denom;
           case 26: denom = coeffCursor[23] + xx*denom;
           case 25: denom = coeffCursor[22] + xx*denom;
           case 24: denom = coeffCursor[21] + xx*denom;
           case 23: denom = coeffCursor[20] + xx*denom;
           case 22: denom = coeffCursor[19] + xx*denom;
           case 21: denom = coeffCursor[18] + xx*denom;
           case 20: denom = coeffCursor[17] + xx*denom;
           case 19: denom = coeffCursor[16] + xx*denom;
           case 18: denom = coeffCursor[15] + xx*denom;
           case 17: denom = coeffCursor[14] + xx*denom;
           case 16: denom = coeffCursor[13] + xx*denom;
           case 15: denom = coeffCursor[12] + xx*denom;
           case 14: denom = coeffCursor[11] + xx*denom;
           case 13: denom = coeffCursor[10] + xx*denom;
           case 12: denom = coeffCursor[ 9] + xx*denom;
           case 11: denom = coeffCursor[ 8] + xx*denom;
           case 10: denom = coeffCursor[ 7] + xx*denom;
           case  9: denom = coeffCursor[ 6] + xx*denom;
           case  8: denom = coeffCursor[ 5] + xx*denom;
           case  7: denom = coeffCursor[ 4] + xx*denom;
           case  6: denom = coeffCursor[ 3] + xx*denom;
           case  5: denom = coeffCursor[ 2] + xx*denom;
           case  4: denom = coeffCursor[ 1] + xx*denom;
           case  3: denom = coeffCursor[ 0] + xx*denom;
           default:
           ;
         }
         denom = 1 + xx*denom;
         result = numer/denom;
      }
      return result;
   }
   double create(const std::vector<double>& inputs,
                 const std::vector<double>& outputs,
                 const double tolerance,
                 const double rangeWindow=0.1);

 public:
   int numNumer_;
   int numDenom_;
   
   std::vector<double> coeff_;
};

#endif
