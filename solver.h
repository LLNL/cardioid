/* 
 * File:   solver.h
 * Author: zhang30
 *
 * Created on February 22, 2017, 10:54 PM
 */

#ifndef SOLVER_H
#define	SOLVER_H

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include "io.h"
#include "option.h"

using namespace std;
using namespace mfem;

void setSurfaces(Mesh *mesh, vector<Vector>& boundingbox, double angle=20, int myid=0);
void setSurf4Surf(Mesh *surface, double angle=20);
void getVert2Elements(Mesh *mesh, vector<vector<int> >& vert2Elements);
GridFunction laplace(Mesh *mesh, Array<int> &all_ess_bdr, Array<int> &nonzero_ess_bdr, Array<int> &zero_ess_bdr, Option& options, int myid=0);
void getVetecesGradients(Mesh *mesh, GridFunction& x, vector<vector<int> >& vert2Elements, vector<double> &pot, vector<Vector> &gradients, string output, int myid=0);

#endif	/* SOLVER_H */

