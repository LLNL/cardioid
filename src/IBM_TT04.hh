/*  
Author: Matthias Reumann, IBM Research
*/
#include "IBM_TT04_LUT.hh"




class IBM_TT04
{
  public:

    char EVfilename[255];
    TT04ParameterStruct tt04ps; //tt04ps stands for TT04ParameterStruct
    double* y_TT04;
                
#ifdef EULER
    double* dy_TT04;
#endif // EULER
#ifdef RK4
    double* yk;
    double* k1;
    double* k2;
    double* k3;
    double* k4;
#endif // RK4
                                                        
    int my_cellType;

    double a_K1;
    double b_K1;
    double XK1_inf;

    double Xr1_inf;
    double tau_Xr1;
    double Xr2_inf;
    double tau_Xr2;

    double Xs_inf;
    double tau_Xs;    
    
    double m_inf;
    double tau_m;
    double h_inf;
    double tau_h;
    double j_inf;    
    double tau_j;

    double d_inf;
    double tau_d;
    double f_inf;
    double tau_f;
    double alpha_fCa;
    double beta_fCa;
    double gamma_fCa;
    double fCa_inf;
    double tau_fCa;
    double tau_fCa_reciproc;
    double CaL_P1;
    double CaL_P2;
    double g_inf;
    double tau_g;
    double tau_g_reciproc;

    double rec_ipK;

    double tau_r;
    double tau_s;
    double r_inf;
    double s_inf;
    double NaCa_P1;
    double NaCa_P2;
    double NaK_P1;


    double Kupsquare;
    double brelsquare;
    double Ca_ibufc;
    double Ca_srbufsr;


    double I_K1;
    double I_to;
    double I_Kr;
    double I_Ks;
    double I_CaL;
    double I_NaK;
    double I_Na;
    double I_bNa;
    double I_NaCa;
    double I_bCa;
    double I_pK;
    double I_pCa;
    double I_leak;
    double I_up;
    double I_rel;
    double I_ax;
    double Iion;
    
    double ENa;
    double EK;
    double ECa;
    double EKs;
    
    
    
    IBM_TT04(char *inputEVFilename, int cellType);
    ~IBM_TT04();

//    void parameterCalculate();
 
 
    inline double getCm(){return (tt04ps.Cm);};
    inline double getIK1(){return (I_K1);};
    inline double getIto(){return (I_to);};
    inline double getIKr(){return (I_Kr);};
    inline double getIKs(){return (I_Ks);};
    inline double getICaL(){return (I_CaL);};
    inline double getINaK(){return (I_NaK);};
    inline double getINa(){return (I_Na);};
    inline double getIbNa(){return (I_bNa);};
    inline double getINaCa(){return (I_NaCa);};
    inline double getIbCa(){return (I_bCa);};
    inline double getIpK(){return (I_pK);};
    inline double getIpCa(){return (I_pCa);};
    inline double getIleak(){return (I_leak);};
    inline double getIup(){return (I_up);};
    inline double getIrel(){return (I_rel);};
    inline double getIion(){return (Iion);};
    
    inline double getV(){return (y_TT04[tt04_V]);};
    inline double getNai(){return (y_TT04[tt04_Nai]);};
    inline double getKi(){return (y_TT04[tt04_Ki]);};
    inline double getCai(){return (y_TT04[tt04_Cai]);};
    inline double getXr1(){return (y_TT04[tt04_xr1]);};
    inline double getXr2(){return (y_TT04[tt04_xr2]);};
    inline double getXs(){return (y_TT04[tt04_xs]);};
    inline double getm(){return (y_TT04[tt04_m]);};
    inline double geth(){return (y_TT04[tt04_h]);};
    inline double getj(){return (y_TT04[tt04_j]);};
    inline double getd(){return (y_TT04[tt04_d]);};
    inline double getf(){return (y_TT04[tt04_f]);};
    inline double getfCa(){return (y_TT04[tt04_fCa]);};
    inline double gets(){return (y_TT04[tt04_s]);};
    inline double getr(){return (y_TT04[tt04_r]);};
    inline double getCaSR(){return (y_TT04[tt04_CaSR]);};
    inline double getg(){return (y_TT04[tt04_g]);};
      
    

    void Init(int cellType);
    double Calc(double tinc,  double V,  double i_external);
    void Print();
    void LongPrint();
};
