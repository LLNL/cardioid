#ifndef FGR_UTILS_HH
#define FGR_UTILS_HH

#include <vector>


#ifdef Diff_Weight_Type_Single
#define WeightType float
#define WTSZ 4
const double weightSumTolerance = 1e-7;
#else
#define WeightType double
#define WTSZ 8
const double weightSumTolerance = 1e-14;
#endif

namespace FGRUtils
{
   enum NodeLocation
   { ZZZ, MZZ, ZMZ, PZZ, ZPZ, ZZM, ZZP, MMZ, PMZ, PPZ,
     MPZ, MZM, ZMM, PZM, ZPM, MZP, ZMP, PZP, ZPP
   };

   struct FGRDiffusionParms
   {
      bool   printBBox_;
      double diffusionScale_;
   };
   
   struct DiffWeight
   {
      WeightType A[19];
   };
   

   void setGradientWeights(double* grad, int* tissue,
                           NodeLocation n3, NodeLocation n2, NodeLocation n1,
                           NodeLocation n4, NodeLocation n5, NodeLocation n0);
   
   void f2(unsigned iFace, int tissue[19], double gradPhi[3][19]);
   
   inline int isTissue(int cellType){return cellType > 0;}
};

#endif
