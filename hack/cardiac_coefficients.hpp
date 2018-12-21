#ifndef CARDIAC_COEF
#define CARDIAC_COEF

#include "mfem.hpp"
#include "Reaction.hh"
#include "ThreadServer.hh"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

/// Generic quadrature function coefficient class for using
/// coefficients which only live at integration points
class QuadratureFunctionCoefficient : public Coefficient
{
private:
   QuadratureFunction *QuadF;   
   
public:
   QuadratureFunctionCoefficient(QuadratureFunction *qf) { QuadF = qf; }

   void SetQuadratureFunction(QuadratureFunction *qf) { QuadF = qf; }

   QuadratureFunction * GetQuadFunction() const { return QuadF; }

   virtual double Eval(ElementTransformation &T,
                       const IntegrationPoint &ip);
   virtual double EvalQ(ElementTransformation &T,
                        const int num_ip);
};

/// Wrapper for interface to melodee-generated active tension
/// values that live only at quadrature points
class ReactionFunction : public QuadratureFunction
{
private:

   /// Underlying quadrature space
   QuadratureSpace *QuadS;

   /// Finite element space for the displacements
   FiniteElementSpace *fes;
   QuadratureFunction VmQuad;

   /// melodee compatible data members
   lazy_array<double> Vm;
   lazy_array<double> iStim;
   lazy_array<double> dVm;

   /// Number of integration points
   int nCells;
   double dt;
   std::string objectName;
   ThreadTeam threadGroup;

   /// Cell model object
   std::shared_ptr<Reaction> reaction;
   
public:
   ReactionFunction(QuadratureSpace *qs, FiniteElementSpace *f,const double new_dt, const std::string new_objectName, const ThreadTeam& group);

   /// Call to melodee generated initialization
   void Initialize();

   /// Call to melodee generated calc
   void Calc(const Vector& x);

 private:
   void CalcVm(const Vector& x);
};


#endif
