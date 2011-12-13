#include "cardiac_constitutive_laws.h"

namespace oomph
{


double UsykConstitutiveLaw::Default_coefficient_set[UsykConstitutiveLaw::size_of_coefficient_set] =
{0.88, 8.0, 6.0, 3.0, 24.0, 6.0, 6.0, 100.0};

const UsykConstitutiveLaw::Coefficient factor_index[] = {
  UsykConstitutiveLaw::fiber_fiber,
  UsykConstitutiveLaw::sheet_sheet,
  UsykConstitutiveLaw::cross_sheet_cross_sheet,
  UsykConstitutiveLaw::fiber_sheet,
  UsykConstitutiveLaw::fiber_cross_sheet,
  UsykConstitutiveLaw::sheet_cross_sheet
};

const UsykConstitutiveLaw::Axis_index first_strain_index[] = {
  UsykConstitutiveLaw::fiber,
  UsykConstitutiveLaw::sheet,
  UsykConstitutiveLaw::cross_sheet,
  UsykConstitutiveLaw::fiber,
  UsykConstitutiveLaw::fiber,
  UsykConstitutiveLaw::sheet
};

const UsykConstitutiveLaw::Axis_index second_strain_index[] = {
  UsykConstitutiveLaw::fiber,
  UsykConstitutiveLaw::sheet,
  UsykConstitutiveLaw::cross_sheet,
  UsykConstitutiveLaw::sheet,
  UsykConstitutiveLaw::cross_sheet,
  UsykConstitutiveLaw::cross_sheet
};

double UsykConstitutiveLaw::data_for_derivative_calculation[size_of_upper_tensor_triangle] = {0};
double UsykConstitutiveLaw::data_for_derivative_calculation_adj[size_of_upper_tensor_triangle] = {0};


void UsykConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG(
  const DenseMatrix<double> &g,
  const DenseMatrix<double> &G,
  const DenseMatrix<double> &sigma,
  RankFourTensor<double> &d_sigma_dG,
  const bool &symmetrize_tensor)
{

  AnisotropicConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG(
    g, G,
    sigma,
    d_sigma_dG,
    symmetrize_tensor);
  return;

//Initial error checking
#ifdef PARANOID

//Test that the matrices are of the same dimension
  if(!are_matrices_of_equal_dimensions(g, G)) {
    throw OomphLibError(
      "Matrices passed are not of equal dimension",
      "ConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG()",
      OOMPH_EXCEPTION_LOCATION);
  }

#endif

//Find the dimension of the matrix (assuming that it's square)
  const unsigned dim = G.ncol();

// Strain tensor
  DenseMatrix<double> strain(dim, dim), Gup(dim, dim);

// Upper triangle
  for (unsigned i = 0; i < dim; i++)
    for (unsigned j = i; j < dim; j++)
      strain(i, j) = 0.5 * (G(i, j) - g(i, j));


  double Q = 0;

  for(unsigned ii = 0; ii < size_of_upper_tensor_triangle; ii++) {
    double strain_value = strain( first_strain_index[ii], second_strain_index[ii]);
    Q += Coefficient_set[factor_index[ii]] * strain_value * strain_value;
  }

  double exp_Q = exp(Q);


  for(unsigned ii = 0; ii < size_of_upper_tensor_triangle; ii++) {
    data_for_derivative_calculation[ii] = Coefficient_set[scaling_factor] *
                                          Coefficient_set[factor_index[ii]] *
                                          exp_Q;

    for(unsigned jj = 0; jj < size_of_upper_tensor_triangle; jj++) {

      d_sigma_dG(first_strain_index[ii], second_strain_index[ii], first_strain_index[jj], second_strain_index[jj]) =
        data_for_derivative_calculation[ii] *
        strain(first_strain_index[ii], second_strain_index[ii]) *
        Coefficient_set[factor_index[jj]] *
        strain(first_strain_index[jj], second_strain_index[jj]);

      if(ii == jj)
        d_sigma_dG(first_strain_index[ii], second_strain_index[ii], first_strain_index[jj], second_strain_index[jj]) +=
          data_for_derivative_calculation[ii] * 0.5;
    }
  }

//incompressibility
  double detG = calculate_contravariant(G, Gup);
  double J = sqrt(detG);

  DenseMatrix<double> d_detG(dim, dim);
  d_detG(0, 0) = G(1, 1) * G(2, 2) - G(1, 2) * G(2, 1);
  d_detG(0, 1) = -G(1, 0) * G(2, 2) + G(2, 0) * G(1, 2);
  d_detG(0, 2) = G(1, 0) * G(2, 1) - G(2, 0) * G(1, 1);
  d_detG(1, 1) = G(0, 0) * G(2, 2) - G(2, 0) * G(0, 2);
  d_detG(1, 2) = -G(0, 0) * G(2, 1) + G(2, 0) * G(0, 1);
  d_detG(2, 2) = G(0, 0) * G(1, 1) - G(0, 1) * G(1, 0);

  //Cim = Bik*Ckj*Djm
  //Aij*Cjk = Akj*Cjk

  for(unsigned ii = 0; ii < size_of_upper_tensor_triangle; ii++)
    for(unsigned jj = 0; jj < size_of_upper_tensor_triangle; jj++) {
      double d_G = d_detG(first_strain_index[jj], second_strain_index[jj]);
      double d_J = 0.5 / J * d_G;
      d_sigma_dG(first_strain_index[ii], second_strain_index[ii], first_strain_index[jj], second_strain_index[jj]) +=
        Coefficient_set[penalty_factor] * (
          (2 * J - 1) * d_J * Gup(first_strain_index[ii], second_strain_index[ii]) -
          J * (J - 1) *
          Gup(first_strain_index[ii], second_strain_index[jj]) *
          Gup(first_strain_index[jj], second_strain_index[ii])
        );
    }


  for(unsigned i = 0; i < dim; i++)
    for(unsigned j = 0; j < i; j++)
      for(unsigned ii = 0; ii < dim; ii++)
        for(unsigned jj = 0; jj < ii; jj++)
          d_sigma_dG(ii, jj, i, j) = d_sigma_dG(jj, ii, j, i);

}


void UsykConstitutiveLaw::
calculate_second_piola_kirchhoff_stress_in_anisotropic_law(const DenseMatrix<double> &g,
    const DenseMatrix<double> &G,
    DenseMatrix<double> &sigma)
{
//Error checking
#ifdef PARANOID
  error_checking_in_input(g, G, sigma);
#endif

  //Find the dimension of the problem
  unsigned dim = G.nrow();

  if(dim != 3)
    throw OomphLibError(
      "Cardiac constitutive law is supported in 3D only",
      "UsykConstitutiveLaw::calculate_second_piola_kirchhoff_stress()",
      OOMPH_EXCEPTION_LOCATION);

  // Strain tensor
  DenseMatrix<double> strain(dim, dim);

  // Upper triangle
  for (unsigned i = 0; i < dim; i++)
    for (unsigned j = i; j < dim; j++)
      strain(i, j) = 0.5 * (G(i, j) - g(i, j));




  // Copy across
  for (unsigned i = 0; i < dim; i++)
    for (unsigned j = 0; j < i; j++)
      strain(i, j) = strain(j, i);


  double Q = 0;

  for(unsigned ii = 0; ii < size_of_upper_tensor_triangle; ii++) {
    double strain_value = strain( first_strain_index[ii], second_strain_index[ii]);
    Q += Coefficient_set[factor_index[ii]] * strain_value * strain_value;
  }

  double exp_Q = exp(Q);

  for(unsigned ii = 0; ii < size_of_upper_tensor_triangle; ii++) {
    data_for_derivative_calculation[ii] = Coefficient_set[scaling_factor] *
                                          Coefficient_set[factor_index[ii]] *
                                          exp_Q;

    sigma( first_strain_index[ii], second_strain_index[ii]) =
      data_for_derivative_calculation[ii] *
      strain( first_strain_index[ii], second_strain_index[ii]);
  }


  DenseMatrix<double> Gup(dim, dim);

  //incompressibility
  double detG = calculate_contravariant(G, Gup);

  double J = sqrt(detG);

  for(unsigned ii = 0; ii < size_of_upper_tensor_triangle; ii++) {
    sigma(first_strain_index[ii], second_strain_index[ii]) +=
      Coefficient_set[penalty_factor] *
      J * (J - 1) * Gup(first_strain_index[ii], second_strain_index[ii]);
  }

  for (unsigned ii = 0; ii < dim; ii++)
    for (unsigned jj = 0; jj < ii; jj++)
      sigma(ii, jj) = sigma(jj, ii);

}


} // oomph namespace
