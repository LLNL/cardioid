#ifndef TT06DEV_HH
#define TT06DEV_HH
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
enum { K_i,Na_i,Ca_i,Xr1_gate,Xr2_gate,Xs_gate,m_gate,h_gate,j_gate,r_gate,d_gate,f_gate,f2_gate,s_gate,fCass_gate,Ca_ss,Ca_SR,R_prime,nStateVar}; 
void initState(double *STATES,int cellType);
void initCnst();
double computeUpdates(double dt, double Vm, double* STATES, int cellType);
double get_c9(); 

void initState(double *STATES,int cellType);
void initCnst();
double computeUpdates(double dt, double Vm, double* STATES, int cellType);

class TT06Dev 
{
 public:
   TT06Dev(int);
   double calc(double dt, double Vm, double iStim, double states[nStateVar]);
   double defaultState(unsigned ii){return defaultState_[ii];}
 private:
   static double defaultState_[nStateVar];
   double c9; 
   double states_[nStateVar];
   int cellType_;
};

#endif
