#ifndef CARDIAC_COEF
#define CARDIAC_COEF

#include "mfem.hpp"
#include "Lumens2009.hpp"
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
class ActiveTensionFunction : public QuadratureFunction
{
private:

   /// Underlying quadrature space
   QuadratureSpace *QuadS;

   /// Finite element space for the displacements
   FiniteElementSpace *fes;

   /// melodee compatible data members
   vector<double> stretch;
   vector<double> stretchVel;
   vector<double> tension;
   vector<double> dtension;
   vector<double> actTime;

   vector<string> outputNames;
   vector<int> outOrder;
   double *outArrays[2];

   vector<string> inputNames;
   vector<int> inOrder;
   double *inArrays[2];

   vector<double> nextStretch;

   /// Number of integration points
   int nCells;

   /// Cell model object
   Lumens2009::ThisModel tester;

   /// Fiber coefficient
   VectorCoefficient *Q;

   /// Stretch calculation 
   void CalcStretch(const Vector &x, const double dt);   
   
public:
   ActiveTensionFunction(QuadratureSpace *qs, FiniteElementSpace *f, VectorCoefficient &fib);

   /// Call to melodee generated initialization
   void Initialize();

   /// Try the active tension at each newton step
   void TryStep(const Vector &x, const double dt);

   /// Commit the active tension after a successful timestep
   void CommitStep(const double dt);

   /// Initialize the underlying vector with a constant (for testing only)
   ActiveTensionFunction &operator=(double value);
   
};


#endif
