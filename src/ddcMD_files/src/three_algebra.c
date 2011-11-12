/* $Id$ */
#include "three_algebra.h"
#include <string.h>
#include <complex.h>
#include <math.h>

const THREE_VECTOR vzero = { 0.0 , 0.0 , 0.0 }; 
const THREE_MATRIX I_3x3 = { 1.0 , 0.0 , 0.0 ,
			      0.0 , 1.0 , 0.0 ,
			      0.0 , 0.0 , 1.0 }; 
const THREE_MATRIX mzero = { 0.0 , 0.0 , 0.0 ,
			     0.0 , 0.0 , 0.0 ,
			     0.0 , 0.0 , 0.0 }; 

/*
* "matinv"; compute the determinent and inverse of a 3x3 matrix
 */
double DET(THREE_MATRIX a )
{
	double det;
	double d00, d01, d02;
	double d10, d11, d12;
	double d20, d21, d22;

	d00 = a.yy*a.zz - a.yz*a.zy;
	d11 = a.zz*a.xx - a.zx*a.xz;
	d22 = a.xx*a.yy - a.xy*a.yx;
	d01 = a.yz*a.zx - a.yx*a.zz;
	d12 = a.zx*a.xy - a.zy*a.xx;
	d20 = a.xy*a.yz - a.xz*a.yy;
	d02 = a.yx*a.zy - a.zx*a.yy;
	d10 = a.zy*a.xz - a.xy*a.zz;
	d21 = a.xz*a.yx - a.yz*a.xx;
	det = a.xx*d00 + a.xy*d01 + a.xz*d02;
	return det; 
}
double matinv(THREE_MATRIX*a0, THREE_MATRIX*ainv)
{
	double det;
	double d00, d01, d02;
	double d10, d11, d12;
	double d20, d21, d22;

	d00 = a0->yy*a0->zz - a0->yz*a0->zy;
	d11 = a0->zz*a0->xx - a0->zx*a0->xz;
	d22 = a0->xx*a0->yy - a0->xy*a0->yx;
	d01 = a0->yz*a0->zx - a0->yx*a0->zz;
	d12 = a0->zx*a0->xy - a0->zy*a0->xx;
	d20 = a0->xy*a0->yz - a0->xz*a0->yy;
	d02 = a0->yx*a0->zy - a0->zx*a0->yy;
	d10 = a0->zy*a0->xz - a0->xy*a0->zz;
	d21 = a0->xz*a0->yx - a0->yz*a0->xx;
	det = a0->xx*d00 + a0->xy*d01 + a0->xz*d02;
	ainv->xx = d00/det;
	ainv->yy = d11/det;
	ainv->zz = d22/det;
	ainv->xy = d10/det;
	ainv->yx = d01/det;
	ainv->xz = d20/det;
	ainv->zx = d02/det;
	ainv->yz = d21/det;
	ainv->zy = d12/det;
	return (det);		/* end of matinv */
}

void matrix_set(THREE_MATRIX*a, double value, char *how)
{
	if (strcasecmp(how, "ALL") == 0) a->xx = a->xy = a->xz = a->yx = a->yy = a->yz = a->zx = a->zy = a->zz = value;
	if (strcasecmp(how, "DIAGONAL") == 0) a->xx = a->yy = a->zz = value;
	if (strcasecmp(how, "LOWER") == 0) a->yx = a->zx = a->zy = value;
	if (strcasecmp(how, "UPPER") == 0) a->xy = a->xz = a->yz = value;
	if (strcasecmp(how, "COLX") == 0) a->xx = a->yx = a->zx = value;
	if (strcasecmp(how, "COLY") == 0) a->xy = a->yy = a->zy = value;
	if (strcasecmp(how, "COLZ") == 0) a->xz = a->yz = a->zz = value;
}

int matrix_equal(THREE_MATRIX a1, THREE_MATRIX a2)
{
	if (a1.xx != a2.xx) return 0;
	if (a1.xy != a2.xy) return 0;
	if (a1.xz != a2.xz) return 0;
	if (a1.yx != a2.yx) return 0;
	if (a1.yy != a2.yy) return 0;
	if (a1.yz != a2.yz) return 0;
	if (a1.zx != a2.zx) return 0;
	if (a1.zy != a2.zy) return 0;
	if (a1.zz != a2.zz) return 0;
	return 1;
}



int matrix_equal_tol(THREE_MATRIX a, THREE_MATRIX b, double tol)
{
   if (fabs(a.xx - b.xx) > tol) return 0;
   if (fabs(a.xy - b.xy) > tol) return 0;
   if (fabs(a.xz - b.xz) > tol) return 0;
   if (fabs(a.yx - b.yx) > tol) return 0;
   if (fabs(a.yy - b.yy) > tol) return 0;
   if (fabs(a.yz - b.yz) > tol) return 0;
   if (fabs(a.zx - b.zx) > tol) return 0;
   if (fabs(a.zy - b.zy) > tol) return 0;
   if (fabs(a.zz - b.zz) > tol) return 0;
   return 1;
}

THREE_MATRIX matrix_sadd(double scale, THREE_MATRIX a1, THREE_MATRIX a2)
{
	THREE_MATRIX a;
	a.xx = scale*a1.xx + a2.xx;
	a.xy = scale*a1.xy + a2.xy;
	a.xz = scale*a1.xz + a2.xz;
	a.yx = scale*a1.yx + a2.yx;
	a.yy = scale*a1.yy + a2.yy;
	a.yz = scale*a1.yz + a2.yz;
	a.zx = scale*a1.zx + a2.zx;
	a.zy = scale*a1.zy + a2.zy;
	a.zz = scale*a1.zz + a2.zz;
	return a;
}

THREE_MATRIX matrix_vcadd(THREE_VECTOR m, THREE_MATRIX a1, THREE_MATRIX a2)
{
	THREE_MATRIX a;
	a.xx = m.x*a1.xx + a2.xx;
	a.xy = m.y*a1.xy + a2.xy;
	a.xz = m.z*a1.xz + a2.xz;
	a.yx = m.x*a1.yx + a2.yx;
	a.yy = m.y*a1.yy + a2.yy;
	a.yz = m.z*a1.yz + a2.yz;
	a.zx = m.x*a1.zx + a2.zx;
	a.zy = m.y*a1.zy + a2.zy;
	a.zz = m.z*a1.zz + a2.zz;
	return a;
}

THREE_MATRIX matrix_mmadd(THREE_MATRIX m, THREE_MATRIX a1, THREE_MATRIX a2)
{
	THREE_MATRIX a;
	a.xx = m.xx*a1.xx + a2.xx;
	a.xy = m.xy*a1.xy + a2.xy;
	a.xz = m.xz*a1.xz + a2.xz;
	a.yx = m.yx*a1.yx + a2.yx;
	a.yy = m.yy*a1.yy + a2.yy;
	a.yz = m.yz*a1.yz + a2.yz;
	a.zx = m.zx*a1.zx + a2.zx;
	a.zy = m.zy*a1.zy + a2.zy;
	a.zz = m.zz*a1.zz + a2.zz;
	return a;
}
THREE_MATRIX matrix_matrix(THREE_MATRIX b, THREE_MATRIX c)
{
	THREE_MATRIX a;
	a.xx = b.xx*c.xx + b.xy*c.yx + b.xz*c.zx;
	a.xy = b.xx*c.xy + b.xy*c.yy + b.xz*c.zy;
	a.xz = b.xx*c.xz + b.xy*c.yz + b.xz*c.zz;

	a.yx = b.yx*c.xx + b.yy*c.yx + b.yz*c.zx;
	a.yy = b.yx*c.xy + b.yy*c.yy + b.yz*c.zy;
	a.yz = b.yx*c.xz + b.yy*c.yz + b.yz*c.zz;

	a.zx = b.zx*c.xx + b.zy*c.yx + b.zz*c.zx;
	a.zy = b.zx*c.xy + b.zy*c.yy + b.zz*c.zy;
	a.zz = b.zx*c.xz + b.zy*c.yz + b.zz*c.zz;
	return a;
}

THREE_VECTOR matrix_vector(THREE_MATRIX m, THREE_VECTOR v)
{
	THREE_VECTOR u;
	u.x = m.xx*v.x + m.xy*v.y + m.xz*v.z;
	u.y = m.yx*v.x + m.yy*v.y + m.yz*v.z;
	u.z = m.zx*v.x + m.zy*v.y + m.zz*v.z;
	return u;
}

void vector_matrix25X3(double *s, double *v, double *m)
{
	int i;
	register double s0, s1, s2;
	s0 = s1 = s2 = 0.0;
	for (i = 0; i < 25; i++)
	{
		s0 += m[i]*v[i];
		s1 += m[i + 25]*v[i];
		s2 += m[i + 50]*v[i];
	}
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
}

void vector_matrix25X4(double *s, double *v, double *m)
{
	int i;
	register double s0, s1, s2, s3;
	s0 = s1 = s2 = s3 = 0.0;
	for (i = 0; i < 25; i++)
	{
		s0 += m[i]*v[i];
		s1 += m[i + 25]*v[i];
		s2 += m[i + 50]*v[i];
		s3 += m[i + 75]*v[i];
	}
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

void vector_matrix49X3(double *u, double *v, double *m)
{
	int i;
	register double s0, s1, s2;
	s0 = s1 = s2 = 0.0;
	for (i = 0; i < 49; i++)
	{
		s0 += m[i]*v[i];
		s1 += m[i + 49]*v[i];
		s2 += m[i + 98]*v[i];
	}
	u[0] = s0;
	u[1] = s1;
	u[2] = s2;
}

void vector_matrix49X4(double *s, double *v, double *m)
{
	int i;
	register double s0, s1, s2, s3;
	s0 = s1 = s2 = s3 = 0.0;
	for (i = 0; i < 49; i++)
	{
		s0 += m[i]*v[i];
		s1 += m[i + 49]*v[i];
		s2 += m[i + 98]*v[i];
		s3 += m[i + 147]*v[i];
	}
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

double trace(double *a, double *b)
{
	double v;
	int i;
	v = 0.0;
	for (i = 0; i < 25; i++)
		v += a[i]*b[i];
	return v;
}

THREE_VECTOR cross(THREE_VECTOR*a, THREE_VECTOR*b)
{
	THREE_VECTOR c;
	c.x = a->y*b->z - a->z*b->y;
	c.y = a->z*b->x - a->x*b->z;
	c.z = a->x*b->y - a->y*b->x;
	return c;
}

double dot1(THREE_VECTOR a, THREE_VECTOR b)
{
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}
THREE_MATRIX transpose(THREE_MATRIX a) 
{
	THREE_MATRIX b; 
	b.xx = a.xx; 
	b.xy = a.yx; 
	b.xz = a.zx; 

	b.yx = a.xy; 
	b.yy = a.yy; 
	b.yz = a.zy; 

	b.zx = a.xz; 
	b.zy = a.yz; 
	b.zz = a.zz; 
	return b; 
	
}
THREE_VECTOR eigenvalues(THREE_MATRIX m)
{
	double a,b,c,d,q,w,t;
	double complex r,r2,s;
	THREE_VECTOR ev;

	if (m.xy == 0 && m.yx == 0 &&
	    m.xz == 0 && m.zx == 0 &&
	    m.zy == 0 && m.yz == 0 )
	{
	   ev.x = m.xx;
	   ev.y = m.yy;
	   ev.z = m.zz;
	   return ev;
	}

	a = 1.0;
    b = -TRACE(m);
	c =  m.yy*m.zz + m.xx*m.zz + m.yy*m.zz - m.xy*m.yx - m.yz*m.zy - m.xz*m.xz ; 
	d = DET(m); 
	q = (9*b*c-27*d-2*b*b*b)/54; 	
	w = (3*c-b*b)/9;
	r2 = (q*q + w*w*w);
	r = csqrt(r2);
	s = cpow(q+r,1.0/3.0);
	ev.x = 2.0*creal(s) - b/3.0;
	ev.y =-creal(s) - b/3.0 - sqrt(3.0)*cimag(s);
	ev.z =-creal(s) - b/3.0 + sqrt(3.0)*cimag(s);
	if (fabs(ev.x) < fabs(ev.y) ) { t=ev.x; ev.x = ev.y; ev.y = t; }
	if (fabs(ev.y) < fabs(ev.z) ) { t=ev.y; ev.y = ev.z; ev.z = t; }
	if (fabs(ev.x) < fabs(ev.y) ) { t=ev.x; ev.x = ev.y; ev.y = t; }
	return ev;
}
THREE_VECTOR vector_sadd(double scale, THREE_VECTOR a1, THREE_VECTOR a2)
{
	THREE_VECTOR a;
	a.x = scale*a1.x + a2.x;
	a.y = scale*a1.y + a2.y;
	a.z = scale*a1.z + a2.z;
        return a;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
