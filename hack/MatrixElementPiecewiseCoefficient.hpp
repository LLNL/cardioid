#pragma once

#include "mfem.hpp"
#include <unordered_map>
#include <memory>

mfem::DenseMatrix quat2rot(const mfem::Vector& q);

class MatrixElementPiecewiseCoefficient : public mfem::MatrixCoefficient
{
public:
  MatrixElementPiecewiseCoefficient() : mfem::MatrixCoefficient(3) {}
  
  MatrixElementPiecewiseCoefficient(std::shared_ptr<mfem::ParGridFunction> x)
  : mfem::MatrixCoefficient(3), p_gf_(x) {}

  virtual void Eval(mfem::DenseMatrix &K,
		    mfem::ElementTransformation& T,
		    const mfem::IntegrationPoint &ip);

  std::shared_ptr<mfem::ParGridFunction> p_gf_;
  std::unordered_map<int,mfem::Vector> heartConductivities_;
  std::unordered_map<int,double> bathConductivities_;
};
