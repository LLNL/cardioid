/* $Id$ */
#ifndef THREE_ALGEBRA_H
#define THREE_ALGEBRA_H

#include <complex.h>
#ifdef __cplusplus
extern "C" {
#endif

#define TRACE(A) ( ((A).xx) + ((A).yy) + ((A).zz) )
#define ZEROVIR(A) ((A).xx=(A).xy=(A).xz=(A).yx=(A).yy=(A).yz=(A).zx=(A).zy=(A).zz=0)
#define MATSET(A,D) ((A).xx=(A).xy=(A).xz=(A).yx=(A).yy=(A).yz=(A).zx=(A).zy=(A).zz=(D))
#define MATNORM(A,D) do {(A).xx /= (D) ;(A).xy /= (D) ;(A).xz /= (D) ;(A).yx /= (D) ;(A).yy /= (D) ;(A).yz /= (D) ;(A).zx /= (D) ;(A).zy /= (D) ;(A).zz /= (D);} while(0)
#define SMATNORM(A,D) do {(A).xx /= (D) ;(A).xy /= (D) ;(A).xz /= (D) ;(A).yy /= (D) ;(A).yz /= (D) ;(A).zz /= (D);} while(0)
#define MATSCALE(A,D) do {(A).xx *= (D) ;(A).xy *= (D) ;(A).xz *= (D) ;(A).yx *= (D) ;(A).yy *= (D) ;(A).yz *= (D) ;(A).zx *= (D) ;(A).zy *= (D) ;(A).zz *= (D);} while(0)
#define SMATSCALE(A,D) do {(A).xx *= (D) ;(A).xy *= (D) ;(A).xz *= (D) ;(A).yy *= (D) ;(A).yz *= (D) ;(A).zz *= (D);} while(0)
#define SMATACUM(A,D)   do {(A).xx += (D).xx ;(A).xy += (D).xy ;(A).xz += (D).xz ; ;(A).yy += (D).yy ;(A).yz += (D).yz ; (A).zz += (D).zz;} while(0)
#define MATACUM(A,D)   do {(A).xx += (D).xx ;(A).xy += (D).xy ;(A).xz += (D).xz ; (A).yx += (D).yx ;(A).yy += (D).yy ;(A).yz += (D).yz ; (A).zx += (D).zx ;(A).zy += (D).zy ;(A).zz += (D).zz;} while(0)
#define ZER03D(A) ((A).x=(A).y=(A).z=0)
#define VECACUM(A,D)   do {(A).x += (D).x ;(A).y += (D).y ;(A).z += (D).z ;} while(0)
#define SQ(A)   ((A)*(A))
#define CUBE(A)   ((A)*(A)*(A))
#define DOT(A,B)   (((A).x)*((B).x)+((A).y)*((B).y)+((A).z)*((B).z))
#define DIFFSQ(A,B)   (SQ((A).x-(B).x)+SQ((A).y-(B).y)+SQ((A).z-(B).z))
#define VSQ(A)   ((SQ((A).x)+SQ((A).y)+SQ((A).z)))
#define QOP1(C,eq,A)  do {(C).v eq ((A).v) ;(C).x eq ((A).x) ; (C).y eq ((A).y) ; (C).z eq ((A).z);} while(0)
#define VOP1(C,eq,A)  do {(C.x) eq ((A).x) ;(C).y eq ((A).y) ; (C).z eq ((A).z);} while(0)
#define VOP2(C,eq,A,op,B)   do {(C.x) eq  (((A).x) op ((B).x)); C.y eq (((A).y) op ((B).y)); C.z eq (((A).z) op ((B).z));} while(0)
#define VSVOP(C,eq,A,op,B)  do {(C.x) eq  ( (A) op ((B).x)); C.y eq ((A) op ((B).y)); C.z eq ((A) op ((B).z));} while(0)
#define VSCALE(A,a)  do {((A).x) *= (a) ;  ((A).y) *= (a) ; ((A).z) *= (a);} while(0)
#define CROSS(c,a,b) do {(c).x = (a).y*(b).z-(a).z*(b).y; (c).y=(a).z*(b).x-(a).x*(b).z; (c).z=(a).x*(b).y-(a).y*(b).x;} while(0)
#define VSET(A,X,Y,Z)  do{(A).x=(X);(A).y=(Y);(A).z=(Z);} while(0)
#define VWRITE(fmt,A)  (printf(fmt,(A).x,(A).y,(A).z))
#define VREAD(fmt,A)  (scanf(fmt,&(A).x,&(A).y,&(A).z))
#define MFPRINT(handle,fmt,A)  (fprintf(handle,fmt,(A).xx,(A).xy,(A).xz,(A).yx,(A).yy,(A).yz,(A).zx,(A).zy,(A).zz))
#define SMFPRINT(handle,fmt,A)  (fprintf(handle,fmt,(A).xx,(A).xy,(A).xz,(A).xy,(A).yy,(A).yz,(A).xz,(A).yz,(A).zz))
enum EIGENVALUES { REAL_EIGENVALUES, COMPLEX_EIGENVALUES}; 

typedef struct { double x, y, z; } THREE_VECTOR;
typedef struct { double _Complex x, y, z; } CTHREE_VECTOR;
typedef struct { double v,x, y, z; } FOUR_VECTOR;
typedef struct { float  x, y, z; } THREE_VECTOR_FLOAT;
typedef struct
{
	double xx, xy, xz;
	double yx, yy, yz;
	double zx, zy, zz;
} THREE_MATRIX;

typedef struct
{
	double xx, yy, zz;
	double xy, xz, yz;
} THREE_SMATRIX;

typedef struct { int x, y, z; } THREE_INT;
typedef struct { short  x, y, z; } THREE_SHORT;
typedef struct { unsigned short  x, y, z; } THREE_USHORT;



extern const THREE_VECTOR vzero;
extern const FOUR_VECTOR qzero;
extern const THREE_MATRIX I_3x3;
extern const THREE_MATRIX mzero;
extern const THREE_SMATRIX szero;

double DET(THREE_MATRIX a );
double matinv(THREE_MATRIX*a0, THREE_MATRIX*ainv);
void matrix_set(THREE_MATRIX*a, double value, char *how);
int matrix_equal(THREE_MATRIX a1, THREE_MATRIX a2);
int matrix_equal_tol(THREE_MATRIX a1, THREE_MATRIX a2, double tol);
THREE_MATRIX matrix_sadd(double scale, THREE_MATRIX a1, THREE_MATRIX a2);
THREE_MATRIX matrix_vcadd(THREE_VECTOR m, THREE_MATRIX a1, THREE_MATRIX a2);
THREE_MATRIX matrix_mmadd(THREE_MATRIX m, THREE_MATRIX a1, THREE_MATRIX a2);
THREE_MATRIX matrix_matrix(THREE_MATRIX b, THREE_MATRIX c);
THREE_VECTOR matrix_vector(THREE_MATRIX m, THREE_VECTOR v);
void vector_matrix25X3(double *s, double *v, double *m);
void vector_matrix25X4(double *s, double *v, double *m);
void vector_matrix49X3(double *u, double *v, double *m);
void vector_matrix49X4(double *s, double *v, double *m);
double trace(double *a, double *b);
CTHREE_VECTOR eigenvalues(THREE_MATRIX m, int *status);
THREE_VECTOR cross(THREE_VECTOR*a, THREE_VECTOR*b);
THREE_VECTOR vector_sadd(double scale, THREE_VECTOR a1, THREE_VECTOR a2);
double dot1(THREE_VECTOR a, THREE_VECTOR b);
THREE_MATRIX transpose(THREE_MATRIX); 

#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
