#ifndef OHARARUDYUPDATEGATE_H
#define OHARARUDYUPDATEGATE_H
#ifdef __cplusplus
extern "C" {
#endif
void initExp(void); 
void sGateInit(double *s1Mhu, double *s1TauR, int nCell_s0, int nCell_s1);

void update_mGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_hGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_jGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_Xr1Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_Xr2Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_XsGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_rGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_dGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_fGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_f2Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_jLGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_s0Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_s1Gate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_sGate(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);

#ifdef __cplusplus
}
#endif
#endif 

