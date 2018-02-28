#ifndef CARDIAC_COEF
#define CARDIAC_COEF

#include "mfem.hpp"
#include "Lumens2009.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

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


class ActiveTensionFunction : public QuadratureFunction
{
private:
   QuadratureSpace *QuadS;
   FiniteElementSpace *fes;
   
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

   int nCells;
   
   Lumens2009::ThisModel tester;
   
   void CalcStretch(const Vector &x, const double dt);
   
public:
   ActiveTensionFunction(QuadratureSpace *qs, FiniteElementSpace *f);
   
   void Initialize();
   
   void TryStep(const Vector &x, const double dt);
   
   void CommitStep(const double dt);

   /// Redefine '=' for GridFunction = constant.
   ActiveTensionFunction &operator=(double value);
   
};


#endif
