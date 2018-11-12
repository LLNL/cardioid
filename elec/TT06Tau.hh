#ifndef TT06TAU_HH
#define TT06TAU_HH
struct TauRecipParms{ double x0, x1, c[4]; double (*func)(double, double *, double *) ; };
TauRecipParms* makeTauRecipParms(double V0, double V1, double (*func)(double, double *, double *));
double jTauRecip(double Vm, double *dtauRecip, double *ddtauRecip);
double hTauRecip(double Vm, double *dtauRecip, double *ddtauRecip);
double TauRecipMod(double V, TauRecipParms *parms, double *dR, double *ddR);
#endif 
