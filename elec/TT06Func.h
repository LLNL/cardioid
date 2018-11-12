#ifndef TT06FUNC_H
#define TT06FUNC_H
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#define gateFitOffset 7
#define nGateVar 12 
#define gateOffset 7
#define dVK_i   K_i 
enum TT06_DEV_STATE_INDEX
{ Ca_i, K_i, Na_i, Ca_ss, Ca_SR, R_prime, fCass, m_gate,                \
  h_gate, j_gate, Xr1_gate, Xr2_gate, Xs_gate, r_gate, d_gate,          \
  f_gate, f2_gate,  jL_gate, s_gate, nStateVar} ; 

enum STATETYPE { nonGateVar, GateVar }; 

struct CellTypeParms
{
   int cellType, s_switch;
   double P_NaK, g_Ks, g_Kr, g_to, g_NaL, minK_i, maxK_i, midK_i, minNa_i, maxNa_i, midNa_i;
};
// int mapCell2Dev[]                  {1,2,3,10,17,18,14,7,8,9,4,5,6,16,11,12,13,15};

#ifdef __cplusplus
extern "C" 
{
#endif 
void fv05General(void *fit, double Vm, double *fv);
double fv6General(void *fit, double dv);
#ifdef __cplusplus
}
#endif


#endif
