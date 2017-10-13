#ifndef CARDIAC_PHYSICS
#define CARDIAC_PHYSICS

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

void ReferenceConfiguration(const Vector &x, Vector &y);
void InitialDeformation(const Vector &x, Vector &y);
void BodyForceFunction(const Vector &x, Vector &y);
void ActiveTensionFunction(const Vector &x, DenseMatrix &y);
void TractionFunction(const Vector &x, Vector &y);
double PressureFunction(const Vector &x);
void FiberFunction(const Vector &x, Vector &y);

#endif
