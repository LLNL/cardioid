#ifndef OOMPH_CARDIAC_UTILITIES_HEADER
#define OOMPH_CARDIAC_UTILITIES_HEADER

#include <generic.h>

namespace oomph
{

/// \short Calculate eigenvectors and eigenvalues for 3x3 matrix
///
/// This funciton takes lower triangular part of the 3x3 matrix in the row-major
/// format and return a matrix
/// of eigen vectors (columns in the matrix), ordered in descending order of absolute values of the eigen
/// values. We use this function to obtain orthogonal matrix for rotation transformation
/// from tensors. For more details about vector-space of tensors see <A HREF="http://www.ncbi.nlm.nih.gov/pubmed/16788917"> Arsigny et al.
/// "Log‐Euclidean metrics for fast and simple calculus on diffusion tensors" </A>
void get_eigen_vectors(double* lower_triangular_part, DenseMatrix<double> &output, double* eig);
 
 /// \short C=A*B
 inline void mult_A_B(const oomph::DenseMatrix<double> &A, const oomph::DenseMatrix<double> &B,
                      oomph::DenseMatrix<double> &C){
   unsigned dim1 = A.nrow();
   unsigned dim2 = A.ncol();
   unsigned dim3 = B.ncol();

   C.initialise(0.0);
   for(unsigned ii = 0; ii < dim1; ii++)
      for(unsigned jj = 0; jj < dim3; jj++)
        for(unsigned kk = 0; kk < dim2; kk++)
          C(ii, jj)+= A(ii, kk)*B(kk, jj);
        
 }

 /// \short C=AT*B
 inline void mult_AT_B(const oomph::DenseMatrix<double> &A, const oomph::DenseMatrix<double> &B,
                      oomph::DenseMatrix<double> &C){
   unsigned dim1 = A.nrow();
   unsigned dim2 = A.ncol();
   unsigned dim3 = B.ncol();

   C.initialise(0.0);
   for(unsigned ii = 0; ii < dim2; ii++)
      for(unsigned jj = 0; jj < dim3; jj++)
        for(unsigned kk = 0; kk < dim1; kk++)
          C(ii, jj)+= A(kk, ii)*B(kk, jj);
        
 }

 /// \short C=A*BT
 inline void mult_A_BT(const oomph::DenseMatrix<double> &A, const oomph::DenseMatrix<double> &B,
                      oomph::DenseMatrix<double> &C){
   unsigned dim1 = A.nrow();
   unsigned dim2 = A.ncol();
   unsigned dim3 = B.nrow();

   C.initialise(0.0);
   for(unsigned ii = 0; ii < dim1; ii++)
      for(unsigned jj = 0; jj < dim3; jj++)
        for(unsigned kk = 0; kk < dim2; kk++)
          C(ii, jj)+= A(ii, kk)*B(jj, kk);
        
 }


}

#endif
