#include <math.h>
#include <sys/param.h>     //needed for endian stuff
#ifdef BYTE_ORDER
# if BYTE_ORDER ==LITTLE_ENDIAN
#define  LeadingBits 1
# else
#  if BYTE_ORDER == BIG_ENDIAN
#define LeadingBits 0
#  else
     Error: unknown byte order!
#  endif
# endif
#endif /* BYTE_ORDER */
static double c[] ={0.5*M_LN2,2.0,2/3.0,2/5.0,2.0/7.0,2.0/9.0,2.0/11.0,2/13.0}; // c[0] = 0.5*log(2); 
static double p[32];
static double d[32];
void fastLogInit()
{
    for (int n=0;n<32;n++) 
    {
       p[n] =  pow(0.5,n-16); 
       d[n] =  (n-16)*M_LN2; 
    }
}
static inline double fastLog(double x)
{
       
       union {unsigned int u[2]; double d ;} t;
       t.d=x; 
       int n = (t.u[LeadingBits]>>20)-1023; 
       t.u[LeadingBits] -= n<<20;
       double m =  t.d;  
       double z = (M_SQRT1_2*m-1.0)/(M_SQRT1_2*m+1.0);
       double z2 = z*z; 
       double f = n*M_LN2 +c[0] + z*(c[1] + z2*(c[2] + z2*(c[3]+z2*(c[4]+z2*(c[5]+z2*(c[6]+z2*c[7]))))));
       return f ; 
}
static inline double fastLogNew(double *x)
{
       register unsigned exp  = ((unsigned *)x)[LeadingBits];
       register int n = (exp >> 20)-1023+16; 
       double m =  *x*p[n];  
       double z = (M_SQRT1_2*m-1.0)/(M_SQRT1_2*m+1.0);
       double z2 = z*z; 
       double f = d[n] +c[0] + z*(c[1] + z2*(c[2] + z2*(c[3]+z2*(c[4]+z2*(c[5]+z2*(c[6]+z2*c[7]))))));
       return f ; 
}
/*
static inline double fastLog4(double *x)
{
       
       register unsigned expa  = ((unsigned *)x+0)[LeadingBits];
       register int na = (expa >> 20)-1023+16; 
       double ma =  x[0]*p[na];  
       double za = (M_SQRT1_2*ma-1.0)/(M_SQRT1_2*ma+1.0);
       double za2 = za*za; 

       register unsigned expb  = ((unsigned *)x+1)[LeadingBits];
       register int nb = (expb >> 20)-1023+16; 
       double mb =  x[1]*p[nb];  
       double zb = (M_SQRT1_2*mb-1.0)/(M_SQRT1_2*mb+1.0);
       double zb2 = zb*zb; 

       register unsigned expc  = ((unsigned *)x+2)[LeadingBits];
       register int nc = (expc >> 20)-1023+16; 
       double mc =  x[1]*p[nc];  
       double zc = (M_SQRT1_2*mc-1.0)/(M_SQRT1_2*mc+1.0);
       double zc2 = zc*zc; 

       register unsigned expd  = ((unsigned *)x+3)[LeadingBits];
       register int nd = (expd >> 20)-1023+16; 
       double md =  x[1]*p[nd];  
       double zd = (M_SQRT1_2*md-1.0)/(M_SQRT1_2*md+1.0);
       double zd2 = zd*zd; 

       double fa = d[na] +c[0] + za*(c[1] + za2*(c[2] + za2*(c[3]+za2*(c[4]+za2*(c[5]+za2*(c[6]+za2*c[7]))))));
       double fb = d[nb] +c[0] + zb*(c[1] + zb2*(c[2] + zb2*(c[3]+zb2*(c[4]+zb2*(c[5]+zb2*(c[6]+zb2*c[7]))))));
       double fc = d[nc] +c[0] + zc*(c[1] + zc2*(c[2] + zc2*(c[3]+zc2*(c[4]+zc2*(c[5]+zc2*(c[6]+zc2*c[7]))))));
       double fd = d[nd] +c[0] + zd*(c[1] + zd2*(c[2] + zd2*(c[3]+zd2*(c[4]+zd2*(c[5]+zd2*(c[6]+zd2*c[7]))))));
       return f ; 
}
*/

