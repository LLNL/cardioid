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

using namespace std;
using namespace mfem;

void setSurfaces(Mesh *mesh, double angle);
void getVert2Elements(Mesh *mesh, vector<vector<int> >& vert2Elements);
void laplace(Mesh *mesh, vector<vector<int> >& vert2Elements, vector<double> &pot, vector<Vector> &gradients, Array<int> &all_ess_bdr, Array<int> &nonzero_ess_bdr, Array<int> &zero_ess_bdr, string output, int order, bool static_cond);


#endif	/* SOLVER_H */

