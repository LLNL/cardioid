#ifndef GRANDIINIT_H
#define GRANDIINIT_H
#ifdef __cplusplus
extern "C" 
{
#endif 
void Grandi_RevpotsInit();
COMPONENTINFO Grandi_VoltageInit();
COMPONENTINFO Grandi_ConcentInit();
COMPONENTINFO Grandi_FluxesInit();
COMPONENTINFO Grandi_INaInit();
COMPONENTINFO Grandi_INaLInit();
COMPONENTINFO Grandi_INabInit();
COMPONENTINFO Grandi_IKrInit();
COMPONENTINFO Grandi_IKsInit();
COMPONENTINFO Grandi_IKurInit();
COMPONENTINFO Grandi_IKpInit();
COMPONENTINFO Grandi_ItoInit();
COMPONENTINFO Grandi_IK1Init();
COMPONENTINFO Grandi_IClCaInit();
COMPONENTINFO Grandi_IClbInit();
COMPONENTINFO Grandi_ICaInit();
COMPONENTINFO Grandi_INCXInit();
COMPONENTINFO Grandi_INaKInit();
COMPONENTINFO Grandi_IpCaInit();
COMPONENTINFO Grandi_ICabInit();

COMPONENTINFO null_INullInit();
 
int GrandiGet_nComp(); 
COMPONENTINFO* GrandiGet_compInfo(); 

void Grandi_Revpots(double Naj, double Nasl, double Caj, double Casl, double Cli, double Ki, DERIVED *derived);

void GrandiCellular(); 
#ifdef __cplusplus
}
#endif 
#endif
