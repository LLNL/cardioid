/* 
 * File:   io.h
 * Author: zhang30
 *
 * Created on February 22, 2017, 10:00 PM
 */

#ifndef IO_H
#define	IO_H

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace mfem;
using namespace std;

void printSurfVTK(Mesh *mesh, std::ostream &out);
void printFiberVTK(Mesh *mesh, vector<Vector>& fiber_vecs, std::ostream &out);

#endif	/* IO_H */

