#include "mfem.hpp"
#include "cardiac_physics.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

void ReferenceConfiguration(const Vector &x, Vector &y)
{
   // set the reference, stress
   // free, configuration
   y = x;
}


void InitialDeformation(const Vector &x, Vector &y)
{
   // set the initial configuration. Having this different from the
   // reference configuration can help convergence
   y = x;
   y[1] = x[1] - 0.5*x[0];
}
                              
void BodyForceFunction(const Vector &x, Vector &y)
{
   y = 0.0;
   //y(1) = -1.0e0;
   //y(2) = -1.0e0;
}

void ActiveTensionFunction(const Vector &x, DenseMatrix &y)
{
   Vector dir(3);
   double tension_strength = 0.0;
   //double tension_strength = 1.0e0;
   y.SetSize(3);
   y = 0.0;

   FiberFunction(x, dir);
   MultVVt(dir, y);
   y *= tension_strength;
}

void FiberFunction(const Vector &x, Vector &y)
{
   y = 0.0;

   y(0) = 1.0;

   //y(1) = 1.0;
   //y(0) = 1.0;

   y /= y.Norml2();

}

void TractionFunction(const Vector &x, Vector &y)
{
   y = 0.0;

   
   //y(1) = -1.0e-3;

   //y(1) = (x(2) - 0.0005) * 0.8;
   //y(2) = (x(1) - 0.0005) * -0.8;


}

double PressureFunction(const Vector &x)
{
   return 4.0e-3;
   //return 0.0;
}
