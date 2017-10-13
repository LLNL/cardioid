#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace mfem;

void setSurfaces(Mesh *mesh);
void printSurfVTK(Mesh *mesh, std::ostream &out);
bool isPlanar(double *coor0, double *coor1, double *coor2, double cosTheta);
bool isTriInTet(vector<int>& tri, vector<int>& tet);
void findNeighbor(Element* ele, vector<Element*>& elements, int attr);
