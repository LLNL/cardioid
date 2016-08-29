#include <math.h>
#include "Grandi.h"
#include "Grandi_INCX.h"
void Grandi_INCXFunc(CELLPARMS *parmsPtr, double *cell, int pOffset, DERIVED *derived, double dt)
{

   VOLTAGE *voltage = (VOLTAGE *)cell; 
   CONCENTRATIONS   *concentrations = (CONCENTRATIONS*) (cell + CONCENTRATIONS_OFFSET); 
   PARAMETERS *cP  = (PARAMETERS *)parmsPtr; 
   double v = voltage->Vm; 
   double Caj = concentrations->Caj; 
   double Casl = concentrations->Casl; 
   double Naj = concentrations->Naj; 
   double Nasl = concentrations->Nasl; 

   double VFRT = v*FRT;
   double Qpow=(T-310.0)/10.0;

   double KmCai=3.59e-3;    // [mM]
   double KmCao=1.3;        // [mM]
   double KmNai=12.29;      // [mM]
   double KmNao=87.5;       // [mM]
   double ksat=0.27;        // [none]
   double nu=0.35;          // [none]
   double Kdact=0.384e-3;   // [mM] 0.256 rabbit
   double Q10NCX=1.57;      // [none]

   double Ka_junc=1.0/(1.0+pow(Kdact/Caj,2.0));
   double Ka_sl=1.0/(1.0+pow(Kdact/Casl,2.0));
   double s1_junc=exp(nu*VFRT)*Naj*Naj*Naj*Cao;
   double s1_sl=exp(nu*VFRT)*Nasl*Nasl*Nasl*Cao;
   double s2_junc=exp((nu-1.0)*VFRT)*Nao*Nao*Nao*Caj;
   double s3_junc=KmCai*Nao*Nao*Nao*(1+pow(Naj/KmNai,3.0))+KmNao*KmNao*KmNao*Caj*(1.0+Caj/KmCai)+KmCao*Naj*Naj*Naj+Naj*Naj*Naj*Cao+Nao*Nao*Nao*Caj;
   double s2_sl=exp((nu-1.0)*VFRT)*Nao*Nao*Nao*Casl;
   double s3_sl=KmCai*Nao*Nao*Nao*(1.0+pow(Nasl/KmNai,3.0))+KmNao*KmNao*KmNao*Casl*(1.0+Casl/KmCai)+KmCao*Nasl*Nasl*Nasl+Nasl*Nasl*Nasl*Cao+Nao*Nao*Nao*Casl;

   double phi = (1.0+0.4*cP->AF);
   double Fsl=1.0-Fjunc;

   derived->I.NCX_junc=Fjunc*cP->GNCX*pow(Q10NCX,Qpow)*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1.0+ksat*exp((nu-1.0)*VFRT));
   derived->I.NCX_sl=Fsl*cP->GNCX*pow(Q10NCX,Qpow)*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1.0+ksat*exp((nu-1.0)*VFRT));
}
