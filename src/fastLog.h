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


/*
  Defining LEGACY_EXPFIX makes the code identical to the one used before.
  Not defining it, makes the fastLog functions use a macro to extract
  mantissa and exponent from a vector of double precision numbers.

  Defining VECTOR_PREFIX uses fast BGQ vectorized code.
  Not defining it uses integer code identical to the legacy code
  (except for being enclosed in a macro).

  Defining AGGRESSIVE uses an ultra-fast algorithm which loses
  11 bits of precision (the mantissa retains 42 out of 53 bits,
  equating 12 out of 15 digits). The aggressive code saves about
  2.2us on 16x16x14 blockPerformance test.
*/
//#define LEGACY_EXPFIX
#ifdef BGQ
#define VECTOR_EXPFIX
//#define AGGRESSIVE
#endif
#ifdef VECTOR_EXPFIX
/*
  Fast vectorized floating point pipeline code for extracting
  mantissa and exponent from vector of 4 double precision numbers.
*/
#ifdef AGGRESSIVE
/* Agressive: This version loses 11 bits of precision... */
#define expfix(x,m,nn)						\
  do {								\
    const double iscale = (1.0/(1<<26))/(1<<26);		\
    vector4double z,a,t;					\
    z = vec_cfidu(x);						\
    a = vec_msub(z,vec_splats(iscale),vec_splats(1023.0));	\
    t = vec_floor(a);						\
    nn = t;							\
    m = vec_add(vec_sub(a,t),vec_splats(1.0));			\
  } while(0)
#else
/* Conservative: This version has no accuracy loss... */
#define expfix(x,m,nn)							\
  do {									\
    const double scale = 1<<20, iscale = 1.0/scale;			\
    const double iscale2 = iscale/(1<<12);				\
    vector4double buf[1];						\
    double *p = (double *) buf;						\
    vector4double a,b,as,a2,b2,z1,z2;					\
									\
    vec_st(x,0,p);							\
    z1 = vec_cfidu(vec_ldiz(0,(unsigned int *) p));			\
    z2 = vec_cfidu(vec_ldiz(2*sizeof(double),(unsigned int *) p));	\
									\
    a = vec_perm(z1,z2, vec_gpci((0<<9) | (2<<6) | (4<<3) | 6) );	\
    b = vec_perm(z1,z2, vec_gpci((1<<9) | (3<<6) | (5<<3) | 7) );	\
									\
    as = vec_trunc(vec_mul(a,vec_splats(iscale)));			\
    nn = vec_sub(as,vec_splats(1023.0));				\
									\
    a2 = vec_nmsub(as,vec_splats(scale),a);				\
    b2 = vec_madd(b,vec_splats(iscale2),a2);				\
    m = vec_madd(b2,vec_splats(iscale),vec_splats(1.0));		\
  } while(0)
#endif

#else
/*
  Old serial integer pipeline code for extracting mantissa an
  exponent from vector of 4 double precision numbers.
*/
#define expfix(x,m,nn)					\
  do {							\
    int n0,n1,n2,n3;					\
    							\
    union {unsigned int u[4][2]; vector4double d ;} t;	\
    t.d=x;						\
    							\
    n0 = (t.u[0][LeadingBits]>>20)-1023;		\
    t.u[0][LeadingBits] -= n0<<20;			\
    							\
    n1 = (t.u[1][LeadingBits]>>20)-1023;		\
    t.u[1][LeadingBits] -= n1<<20;			\
    							\
    n2 = (t.u[2][LeadingBits]>>20)-1023;		\
    t.u[2][LeadingBits] -= n2<<20;			\
    							\
    n3 = (t.u[3][LeadingBits]>>20)-1023;		\
    t.u[3][LeadingBits] -= n3<<20;			\
    							\
    m = t.d;						\
    {							\
      vector4double nnx = {n0, n1, n2, n3};		\
      nn = nnx;						\
    }							\
  } while(0)
#endif

static inline vector4double fastLog4(vector4double x4)
{
  //vector4double y4 = {log(x4 get [0]), log(x4 get [1]), log(x4 get [2]), log(x4 get [3]) };
  //return y4; 

#ifdef LEGACY_EXPFIX
  int n0,n1,n2,n3;
#endif
  vector4double z, z2, m;

#ifdef LEGACY_EXPFIX
  union {unsigned int u[4][2]; vector4double d ;} t;
  t.d=x4; 
#endif

  vector4double c0 = vec_splats(c[0]);
  vector4double c1 = vec_splats(c[1]);
  vector4double c2 = vec_splats(c[2]);
  vector4double c3 = vec_splats(c[3]);
  vector4double c4 = vec_splats(c[4]);
  vector4double c5 = vec_splats(c[5]);
  vector4double c6 = vec_splats(c[6]);
  vector4double c7 = vec_splats(c[7]);

#ifdef LEGACY_EXPFIX
  n0 = (t.u[0][LeadingBits]>>20)-1023; 
  t.u[0][LeadingBits] -= n0<<20;

  n1 = (t.u[1][LeadingBits]>>20)-1023; 
  t.u[1][LeadingBits] -= n1<<20;

  n2 = (t.u[2][LeadingBits]>>20)-1023; 
  t.u[2][LeadingBits] -= n2<<20;

  n3 = (t.u[3][LeadingBits]>>20)-1023; 
  t.u[3][LeadingBits] -= n3<<20;

  m  = t.d;  
#else
  vector4double nn;
  expfix(x4,m,nn);
#endif

  vector4double m_sqrt1_2 = vec_splats(M_SQRT1_2);
  vector4double one = vec_splats(1.0);
  vector4double zn = vec_msub(m_sqrt1_2,m,one);
  vector4double zd = vec_madd(m_sqrt1_2,m,one);

  //vector4double zn = vec_msub(vec_splats(M_SQRT1_2), m, vec_splats(1.0));
  //vector4double zd = vec_madd(vec_splats(M_SQRT1_2), m, vec_splats(1.0));
  z = vec_swdiv(zn, zd);
  z2 = vec_mul(z,z);

#ifdef LEGACY_EXPFIX
  vector4double nn = {n0, n1, n2, n3};
#endif
  nn = vec_mul(nn, vec_splats(M_LN2));

  vector4double f;
  vector4double tmp;

#define vec_madd2(c,a,b) vec_madd(a,b,c)

  tmp = vec_madd2(c6, z2, c7);
  tmp = vec_madd2(c5, z2, tmp);
  tmp = vec_madd2(c4, z2, tmp);
  tmp = vec_madd2(c3, z2, tmp);
  tmp = vec_madd2(c2, z2, tmp);
  tmp = vec_madd2(c1, z2, tmp);
  tmp = vec_madd2(c0, z, tmp);
  f = vec_add(nn, tmp);

  return f;
}


//computes:
//*ret1 = fastLog4(y4);
//return fastLog4(x4);
static inline vector4double fastLog8(vector4double x4, vector4double y4, vector4double *ret1)
{
  //vector4double fx= {log(x4 get [0]), log(x4 get [1]), log(x4 get [2]), log(x4 get [3]) };
  //vector4double fy= {log(y4 get [0]), log(y4 get [1]), log(y4 get [2]), log(y4 get [3]) };
  //*ret1 = fy; 
  //return fx; 

#ifdef LEGACY_EXPFIX
  int n0,n1,n2,n3;
#endif

  vector4double zx, z2x, nnx;
  vector4double zy, z2y, nny;

#ifdef LEGACY_EXPFIX
  union {unsigned int u[4][2]; vector4double d ;} t;
#endif

  vector4double c0 = vec_splats(c[0]);
  vector4double c1 = vec_splats(c[1]);
  vector4double c2 = vec_splats(c[2]);
  vector4double c3 = vec_splats(c[3]);
  vector4double c4 = vec_splats(c[4]);
  vector4double c5 = vec_splats(c[5]);
  vector4double c6 = vec_splats(c[6]);
  vector4double c7 = vec_splats(c[7]);

  vector4double m_sqrt1_2 = vec_splats(M_SQRT1_2);
  vector4double one = vec_splats(1.0);

  {
#ifdef LEGACY_EXPFIX
    t.d=y4; 
  
    n0 = (t.u[0][LeadingBits]>>20)-1023; 
    t.u[0][LeadingBits] -= n0<<20;

    n1 = (t.u[1][LeadingBits]>>20)-1023; 
    t.u[1][LeadingBits] -= n1<<20;

    n2 = (t.u[2][LeadingBits]>>20)-1023; 
    t.u[2][LeadingBits] -= n2<<20;

    n3 = (t.u[3][LeadingBits]>>20)-1023; 
    t.u[3][LeadingBits] -= n3<<20;

    vector4double m  = t.d;  
#else
    vector4double m,nn;
    expfix(y4,m,nn);
#endif

    vector4double zn = vec_msub(m_sqrt1_2,m,one);
    vector4double zd = vec_madd(m_sqrt1_2,m,one);

    //vector4double zn = vec_msub(vec_splats(M_SQRT1_2), m, vec_splats(1.0));
    //vector4double zd = vec_madd(vec_splats(M_SQRT1_2), m, vec_splats(1.0));
    zx = vec_swdiv(zn, zd);
    z2x = vec_mul(zx,zx);


#ifdef LEGACY_EXPFIX
    vector4double nn = {n0, n1, n2, n3};
#endif
    nnx = vec_mul(nn, vec_splats(M_LN2));
  }

  {
#ifdef LEGACY_EXPFIX
    t.d=x4; 
  
    n0 = (t.u[0][LeadingBits]>>20)-1023; 
    t.u[0][LeadingBits] -= n0<<20;

    n1 = (t.u[1][LeadingBits]>>20)-1023; 
    t.u[1][LeadingBits] -= n1<<20;

    n2 = (t.u[2][LeadingBits]>>20)-1023; 
    t.u[2][LeadingBits] -= n2<<20;

    n3 = (t.u[3][LeadingBits]>>20)-1023; 
    t.u[3][LeadingBits] -= n3<<20;

    vector4double m  = t.d;  
#else
    vector4double m,nn;
    expfix(x4,m,nn);
#endif

    vector4double zn = vec_msub(m_sqrt1_2,m,one);
    vector4double zd = vec_madd(m_sqrt1_2,m,one);

    //vector4double zn = vec_msub(vec_splats(M_SQRT1_2), m, vec_splats(1.0));
    //vector4double zd = vec_madd(vec_splats(M_SQRT1_2), m, vec_splats(1.0));
    zy = vec_swdiv(zn, zd);
    z2y = vec_mul(zy,zy);

#ifdef LEGACY_EXPFIX
    vector4double nn = {n0, n1, n2, n3};
#endif
    nny = vec_mul(nn, vec_splats(M_LN2));

  }
  vector4double f;
  vector4double tmpx, tmpy;

#define vec_madd2(c,a,b) vec_madd(a,b,c)

  tmpx = vec_madd2(c6, z2x, c7);
  tmpy = vec_madd2(c6, z2y, c7);

  tmpx = vec_madd2(c5, z2x, tmpx);
  tmpy = vec_madd2(c5, z2y, tmpy);
  
  tmpx = vec_madd2(c4, z2x, tmpx);
  tmpy = vec_madd2(c4, z2y, tmpy);
  
  tmpx = vec_madd2(c3, z2x, tmpx);
  tmpy = vec_madd2(c3, z2y, tmpy);

  tmpx = vec_madd2(c2, z2x, tmpx);
  tmpy = vec_madd2(c2, z2y, tmpy);

  tmpx = vec_madd2(c1, z2x, tmpx);
  tmpy = vec_madd2(c1, z2y, tmpy);

  tmpx = vec_madd2(c0, zx, tmpx);
  tmpy = vec_madd2(c0, zy, tmpy);

  *ret1 = vec_add(nnx, tmpx);
  f = vec_add(nny, tmpy);

    return f;

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

