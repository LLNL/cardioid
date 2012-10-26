#ifndef UPDATEGATE_H
#define UPDATEGATE_H
#ifdef __cplusplus
extern "C" {
#endif
void initExp(void); 
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

void update_mGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_hGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_jGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_Xr1Gate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_Xr2Gate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_XsGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_rGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_dGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_fGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_f2Gate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_jLGate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_s0Gate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
void update_s1Gate_v1(double dt, int nCells, double *VM, double *g, double *mhu_a, double *tauR_a);
#ifdef __cplusplus
}
#endif
#endif 

