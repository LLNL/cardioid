#ifndef TT06DEVFIT_HH
#define TT06DEVFIT_HH
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
enum { K_i,Na_i,Ca_i,Xr1_gate,Xr2_gate,Xs_gate,m_gate,h_gate,j_gate,r_gate,d_gate,f_gate,f2_gate,s_gate,fCass_gate,Ca_ss,Ca_SR,R_prime,nStateVar}; 
void initState(double *STATES,int cellType);
void initCnst();
double computeUpdates(double dt, double Vm, double* STATES, int cellType);
void  Approx (int n, int maxOrder, double tol, double dt);
double get_c9(); 
#endif
