#include "FGRUtils.hh"
#include <cassert>

using namespace std;



namespace FGRUtils
{
   void setGradientWeights(double* grad, int* tissue,
                           NodeLocation n3, NodeLocation n2, NodeLocation n1,
                           NodeLocation n4, NodeLocation n5, NodeLocation n0)
   {
      int caseNum = tissue[n4] + 2*tissue[n3] + 4*tissue[n0] + 8*tissue[n1];

      double& w0 = grad[n0];
      double& w1 = grad[n1];
      double& w2 = grad[n2];
      double& w3 = grad[n3];
      double& w4 = grad[n4];
      double& w5 = grad[n5];

      switch (caseNum)
      {
        case 0:
         // all weights are zero
         break;
        case 1:
        // w5 = 1.0; w4 = -1.0;
         break;
        case 2:
        // w2 = 1.0; w3 = -1.0;
         break;
        case 3:
    //     w2 = w5 = 0.5; w3 = w4 = -0.5;
         break;
        case 4:
     //    w0 = 1.0; w5 = -1.0;
         break;
        case 5:
      //   w0 = 0.5; w4 = -0.5;
         break;
        case 6:
       //  w0 = w2 = 0.5; w3 = w5 = -0.5;
         break;
        case 7:
        // w0 = 0.5; w2 = 0.25; w3 = w4 = w5 = -0.25;
         break;
        case 8:
        // w1 = 1.0; w2 = -1.0;
         break;
        case 9:
       //  w1 = w5 = 0.5; w2 = w4 = -0.5;
         break;
        case 10:
       //  w1 = 0.5; w3 = -0.5;
         break;
        case 11:
       //  w1 = 0.5; w5 = 0.25; w2 = w3 = w4 = -0.25;
         break;
        case 12:
        // w0 = w1 = 0.5; w2 = w5 = -0.5;
         break;
        case 13:
        // w4 = -0.5; w2 = -0.25; w0 = w1 = w5 = 0.25;
         break;
        case 14:
       //  w3 = -0.5; w5 = -0.25; w0 = w1 = w2 = 0.25;
         break;
        case 15:
         w0 = w1 = 0.25; w3 = w4 = -0.25;
         break;
        default:
         assert(false);
      }
   }
}

namespace FGRUtils
{
   void f2(unsigned iFace, int tissue[19], double gradPhi[3][19])
   {
      switch (iFace)
      {
        case 0:
         gradPhi[0][ZZZ] = 1.0;
         gradPhi[0][MZZ] = -1.0;
         setGradientWeights(gradPhi[1], tissue,
                            MMZ, MZZ, MPZ,
                            ZMZ, ZZZ, ZPZ);
         setGradientWeights(gradPhi[2], tissue,
                            ZZM, ZZZ, ZZP,
                            MZM, MZZ, MZP);
         break;
        case 1:
         setGradientWeights(gradPhi[0], tissue,
                            MZZ, ZZZ, PZZ,
                            MMZ, ZMZ, PMZ);
         gradPhi[1][ZZZ] = 1.0;
         gradPhi[1][ZMZ] = -1.0;
         setGradientWeights(gradPhi[2], tissue,
                            ZMM, ZMZ, ZMP,
                            ZZM, ZZZ, ZZP);
         break;
        case 2:
         gradPhi[0][PZZ] = 1.0;
         gradPhi[0][ZZZ] = -1.0;
         setGradientWeights(gradPhi[1], tissue,
                            ZMZ, ZZZ, ZPZ,
                            PMZ, PZZ, PPZ);
         setGradientWeights(gradPhi[2], tissue,
                            PZM, PZZ, PZP,
                            ZZM, ZZZ, ZZP);
         break;
        case 3:
         setGradientWeights(gradPhi[0], tissue,
                            MPZ, ZPZ, PPZ,
                            MZZ, ZZZ, PZZ);
         gradPhi[1][ZPZ] = 1.0;
         gradPhi[1][ZZZ] = -1.0;
         setGradientWeights(gradPhi[2], tissue,
                            ZZM, ZZZ, ZZP,
                            ZPM, ZPZ, ZPP);
         break;
        case 4:
         setGradientWeights(gradPhi[0], tissue,
                            MZM, ZZM, PZM,
                            MZZ, ZZZ, PZZ);
         setGradientWeights(gradPhi[1], tissue,
                            ZMZ, ZZZ, ZPZ,
                            ZMM, ZZM, ZPM);
         gradPhi[2][ZZZ] = 1.0;
         gradPhi[2][ZZM] = -1.0;
         break;
        case 5:
         setGradientWeights(gradPhi[0], tissue,
                            MZZ, ZZZ, PZZ,
                            MZP, ZZP, PZP);
         setGradientWeights(gradPhi[1], tissue,
                            ZMP, ZZP, ZPP,
                            ZMZ, ZZZ, ZPZ);
         gradPhi[2][ZZP] = 1.0;
         gradPhi[2][ZZZ] = -1.0;
         break;
        default:
         assert(false);
      }
   }
}
