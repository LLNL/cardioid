#ifndef CARDIAC_PHYSICS
#define CARDIAC_PHYSICS

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

extern int run_mode;

/// Functions to set the loading and initial conditions. These are hard-coded with the run mode;
void ReferenceConfiguration(const Vector &x, Vector &y);
void InitialDeformation(const Vector &x, Vector &y);
double PressureFunction(const Vector &x, double t);
void FiberFunction(const Vector &x, Vector &y);
void VolumeFunction(const Vector &x, Vector &y);
void setSurfaces(Mesh *mesh);
void printSurfVTK(Mesh *mesh, std::ostream &out);
void findNeighbor(Element* ele, vector<Element*>& elements, int attr);
bool isPlanar(double *coor0, double *coor1, double *coor2, double cosTheta);



#endif
