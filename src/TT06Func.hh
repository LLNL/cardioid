#ifndef TT06FUNC_HH
#define TT06FUNC_HH
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

enum TT06STATE { K_i, Na_i, Ca_i, Xr1_gate, Xr2_gate, Xs_gate, m_gate, h_gate, j_gate, r_gate, d_gate, f_gate, f2_gate, s_gate, fCass_gate, Ca_ss, Ca_SR, R_prime, nStateVar} ; 

void initState(double *state,int cellType);
void initCnst();
double computeUpdates(double dt, double Vm, double* state, int cellType);
double get_c9();
void makeFit(double tol, double V0, double V1, double deltaV, int maxOrder, int maxCost, double dt, int mod);
#endif
