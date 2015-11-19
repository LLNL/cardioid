#ifndef OHARARUDYINIT_H
#define OHARARUDYINIT_H
#ifdef __cplusplus
extern "C" 
{
#endif 
void reversalPotentialsInit();
COMPONENTINFO OHaraRudy_ConcentInit();
COMPONENTINFO OHaraRudy_FluxesInit();
COMPONENTINFO OHaraRudy_CaMKtrapInit();
COMPONENTINFO OHaraRudy_ICaInit();
COMPONENTINFO OHaraRudy_ICabInit();
COMPONENTINFO OHaraRudy_IK1Init();
COMPONENTINFO OHaraRudy_IKbInit();
COMPONENTINFO OHaraRudy_IKrInit();
COMPONENTINFO OHaraRudy_IKsInit();
COMPONENTINFO OHaraRudy_INaCaiInit();
COMPONENTINFO OHaraRudy_INaCassInit();
COMPONENTINFO OHaraRudy_INaFastInit();
COMPONENTINFO OHaraRudy_INaKInit();
COMPONENTINFO OHaraRudy_INaLInit();
COMPONENTINFO OHaraRudy_INabInit();
COMPONENTINFO OHaraRudy_IpCaInit();
COMPONENTINFO OHaraRudy_ItoInit();

COMPONENTINFO OHaraRudyMod_INaFastInit();
COMPONENTINFO RTYSC14A_IKrInit();
COMPONENTINFO MYBGBKC_INaInit();
COMPONENTINFO null_INullInit();

int OHaraRudyGet_nComp(); 
COMPONENTINFO* OHaraRudyGet_compInfo(); 

void reversalPotentials(double Nai, double Ki, DERIVED *derived);
#ifdef __cplusplus
}
#endif 
#endif
