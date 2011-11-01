/*  
This file has been modified by Matthias Reumann, IBM Research, to be stand alone
22. March 2010

    File: IBM_tenTusscher04.h
    automatically created by ExtractParameterClass.pl - done by dw (20.04.2007)
    Institute of Biomedical Engineering, Universität Karlsruhe (TH)
    send comments to dw@ibt.uka.de	*/

#define PARAMETER_INCLUDED

#include "IBM_tenTusscher04ParameterDefNLoad.hh"
#include "IBM_tenTusscher04_endoLUT.hh"
#include "IBM_tenTusscher04_midLUT.hh"
#include "IBM_tenTusscher04_epiLUT.hh"

class IBM_tenTusscher04
{
  public:

    // IBM parameter struct
    char EVfilename[255];// = "/bgpusr1/mr198/IBT/src/CellModel/ev_forceParameter/Flaim0910MutA/288_100.ev";
    TTParameterStruct ttps; //fps stands for TTParameterStruct
    double* y_TT;
                
#ifdef EULER
    double* dy_TT;
#endif // EULER
#ifdef RK4
    double* yk;
    double* k1;
    double* k2;
    double* k3;
    double* k4;
#endif // RK4
                                                        
    int my_cellType;

    // for values from lookup tables
    double rec_ipK;
    double d_inf;
    double f_inf;
    double tau_m;
    double tau_h;
    double tau_j;
    double tau_Xr1;
    double tau_Xr2;
    double tau_Xs;
    double tau_r;
    double tau_s;
    double m_inf;
    double h_inf;
    double j_inf;
    double Xr1_inf;
    double Xr2_inf;
    double Xs_inf;
    double r_inf;
    double s_inf;
    double tau_d;
    double tau_f;
    double NaCa_P1;
    double NaCa_P2;
    double NaK_P1;
    double CaL_P1;
    double CaL_P2;
    
    
    IBM_tenTusscher04(char *inputEVFilename, int cellType);
    ~IBM_tenTusscher04();

    inline double SurfaceToVolumeRatio(){return 1.0;}; 
    inline double Volume(){return 2.064403e-13*5.6;}; 
//    inline double GetVm(){return (v(VT_V_init));};
    inline double GetVm(){return (0.0);};
    inline double GetCai(){return (y_TT[tt_Ca_i]);};
//    inline double GetCao(){return (v(VT_Ca_o));};
    inline double GetCao(){return (0.0);};
    inline double GetNai(){return (y_TT[tt_Na_i]);};
//    inline double GetNao(){return (v(VT_Na_o));};
    inline double GetNao(){return (0.0);};
    inline double GetKi(){return (y_TT[tt_K_i]);};
//    inline double GetKo(){return (v(VT_K_o));};
    inline double GetKo(){return (0.0);};
    inline double GetIto(){return 0.0;};
    inline double GetIKr() {return 0.0;};
    inline double GetIKs() {return 0.0;};
//    inline int GetSize(void);
    inline double* GetBase(void){return &(y_TT[tt_Ca_i]);};
    inline double GetSpeedupMax(void){return .0;};
    double GetAmplitude(void){return 52.;};
    inline double GetStimTime() {return 0.001;};
    inline unsigned char getSpeed( double adVm);

    void parameterCalculate();

    void Init(int cellType);
   double Calc(double dt, double V, double iStim);
//    void Print(ostream &tempstr, double tArg,  double V);
    void Print();
//    void LongPrint(ostream &tempstr, double tArg,  double V);
    void LongPrint();
//    void GetParameterNames(vector<string> &getpara);
//    void GetLongParameterNames(vector<string> &getpara);
};
