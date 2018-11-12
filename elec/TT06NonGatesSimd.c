#include <math.h>
#include "portableSIMD.h"
#include "TT06Func.h" 
#include "TT06NonGates.h"
#include "fastLog.h"


void (*fv05Func)(double Vm, double *fv);
double (*fv6Func)(double dv);

double fvX_a[75]__attribute__((aligned(32))) = {
  6.02170007043617e-16, 2.60192856479526e-13, 6.04272220479134e-11, 1.46053164590980e-08, 2.41325169024597e-06, 1.57571330412697e-04,  
  -1.52370496566743e-12, -2.55280464729594e-10, 2.94210996176821e-08, 4.00772142804151e-07, -6.00926106818024e-04, 1.18653494080701e-01, -1.19531891610964e+01, 5.40504807620401e+02, 
  -1.69089871040937e-19, -1.57414794358693e-17, 2.10687633250732e-15, -3.07155144145422e-14, -6.15332407381053e-12, 1.86878131970215e-09, -2.93525992825876e-07, 2.87019818110585e-05, -2.24272406249374e-03, -1.45499489004044e+00, 
  5.64722239204646e-17, 2.40194725811093e-14, 4.32295749601499e-12, 4.31098192118690e-10, 2.61347101393462e-08, 9.83170991434406e-07, 2.16983094649991e-05, 2.19831706117241e-04, 
  2.24306976181154e-07, -2.65106715048650e-05, 2.00764738834124e-03, -6.61573242254081e-02, 1.00000000000000e+00,
  8.93564690058337e-12, -5.91201683806020e-09, 1.50747331544023e-06, -1.75996322715896e-04, 7.99339727741565e-03, -8.89837887994794e-03, 2.97942000558175e-01,
  1.95808469577737e-16, 2.07853354660018e-14, -4.46027161749137e-12, 8.84744201829959e-10, -5.36894894634127e-08, 3.38737191083109e-06, -7.48970504345503e-05, 1.59235097855721e-03, 2.17628774116084e-02, 2.17226317175310e-01, 1.00000000000000e+00, 
  -1.74506487282187e-20, -3.65847291145083e-18, 1.36854440783922e-16, 7.81666161976579e-14, 2.49915783972090e-12, -6.96639048514165e-10, -4.97022959453108e-08, 4.30241562974461e-06, 7.91075918239993e-04, 4.60580107291175e-02, 1.03965405517718e+00,
  1.89136344887477e-15, 4.55066993718199e-13, 1.18578052624143e-11, -5.22361846923760e-09, -3.79918615957248e-07, 3.41107297247441e-05, 6.31701538469039e-03, -7.81445921919197e-01, 2.55684105370893e+01};
 

double SP[40]__attribute__((aligned(32)));

void set_SP(struct nonGateCnst cnst)
{

 SP[0] =  cnst.c26 ;
 SP[1] =  cnst.c27 ;
 SP[2] =  cnst.c28 ;
 SP[3] =  cnst.c29 ;
 SP[4] =  cnst.c8 ;
 SP[5] =  cnst.c7 ;   
 SP[6] =  cnst.c24 ;
 SP[7] =  cnst.c43 ;
 SP[8] =  cnst.c44 ;
 SP[9] =  cnst.c23 ;
 SP[10] = cnst.c15 ;
 SP[11] = cnst.c16 ;
 SP[12] = cnst.c9 ;
 SP[13] = cnst.c25 ;
 SP[14] = cnst.c3 ;
 SP[15] = cnst.c5 ;
 SP[16] = cnst.c4 ;
 SP[17] = cnst.c2 ;
 SP[18] = cnst.c6 ;
 SP[19] = cnst.c11 ;
 SP[20] = cnst.c20 ;
 SP[21] = cnst.c21 ;
 SP[22] = cnst.c22 ;
 SP[23] = cnst.c30 ;
 SP[24] = cnst.c31 ;
 SP[25] = cnst.c32 ;
 SP[26] = cnst.c33 ;
 SP[27] = cnst.c34 ;
 SP[28] = cnst.c19 ;
 SP[29] = cnst.c18 ;
 SP[30] = cnst.c36 ;
 SP[31] = cnst.c17 ;
 SP[32] = 20.0 ;
 SP[33] = 0.6 ;
 SP[34] = 0.4 ;
 SP[35] = 80.0 ;
 SP[36] = cnst.c14 * cnst.c40; 
 SP[37] = cnst.c13 ;
 SP[38] = cnst.c40 ;
 SP[39] = cnst.c36 ;
}


#define sigm(x)   ((x)/(1+(x)))

#define vLog4 vLog4m
#define vLog8 vLog8m
#define vLogSeries4 vLogSeries4Oinf
#define updateFunction update_nonGateSimdM 
#include "update_nonGatesSimd.h"
 
#undef vLog4
#undef vLog8
#undef vLogSeries4
#undef updateFunction 

#define vLog4 vLog4f
#define vLog8 vLog8f
#define vLogSeries4 vLogSeries4O3
#define updateFunction update_nonGateSimdF 
#include "update_nonGatesSimd.h"

#undef vLog4
#undef vLog8
#undef vLogSeries4
#undef updateFunction 

#define vLog4 vLog4f
#define vLog8 vLog8f
#define vLogSeries4 vLogSeries4O3
#define AGRESSIVE_EXP
#define updateFunction update_nonGateSimdFA
#include "update_nonGatesSimd.h"
