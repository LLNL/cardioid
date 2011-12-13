#include "cardiac_utilities.h"

//OOMPH headers
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include <generic.h>


typedef int integer;
typedef double doublereal;


// functions from LAPACK library
extern "C" void dgetri_(int &, double *, int &,
                        int *, double *,
                        int &, int &);
extern "C" void dgetrf_(int &, int &, double *, int &, int *, int &);

extern "C" int dsteqr_(char *compz, integer *n, doublereal *d, doublereal *e, doublereal *z,
                       integer *ldz, doublereal *work, integer *info);

extern "C" int dsytrd_(char *uplo, integer *n, doublereal *a, integer* lda, doublereal *d,
                       doublereal *e, doublereal *tau, doublereal *work, integer* lwork, integer *info);


namespace oomph
{


void get_eigen_vectors(double* lower_triangular_part, DenseMatrix<double> &output, double* eig)
{

// fix_later_tag, check for 3x3, may be expand to general case
  char c = 'V', uplo = 'L';
  int n = 3, lda = 3;
  doublereal d[3], e[2],
             z[9] = {lower_triangular_part[0], lower_triangular_part[1], lower_triangular_part[2],
                     0, lower_triangular_part[3], lower_triangular_part[4], 0, 0, lower_triangular_part[5]
                    }, work[100], tau[2];
  int ldz = 3, info, lwork = 100;

  dsytrd_(&uplo, &n, z, &lda, d, e, tau, work, &lwork, &info);

  if(info)
    throw OomphLibError(
      "dsytrd_ failed, probably zero determinant",
      "get_eigen_vectors()",
      OOMPH_EXCEPTION_LOCATION);

  double v = z[2];
  eig[0] = d[0];
  eig[1] = d[1];
  eig[2] = d[2];
  z[0] = 1;
  z[4] = 1 - tau[0];
  z[5] = -(1 - tau[1]) * v * tau[0];
  z[7] = -v * tau[0];
  z[8] = (1 - tau[1]) * (1 - v * v * tau[0]);
  dsteqr_(&c, &n, eig, e, z, &ldz, work, &info);

  if(info)
    throw OomphLibError(
      "dsteqr_ failed, probably zero determinant",
      "get_eigen_vectors()",
      OOMPH_EXCEPTION_LOCATION);

  // sorting by absolute value of eigenvalues
  std::pair<double, int> index[3] = {std::pair<double, int>(-fabs(eig[0]), 0), 
                                     std::pair<double, int>(-fabs(eig[1]), 1), 
                                     std::pair<double, int>(-fabs(eig[2]), 2)};
  std::sort(index, index + 3); //check_tag

  // converting from row-major to column-major
  for(int ii = 0; ii < 3; ii++)
    for(int jj = 0; jj < 3; jj++)
      output(ii, jj) = z[ii + index[jj].second * 3];

}

}
