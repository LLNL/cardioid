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

#include "cardfiber.h"

using namespace mfem;
using namespace std;

void printSurfVTK(Mesh *mesh, std::ostream &out);
void printSurf4SurfVTK(Mesh *mesh, std::ostream &out);
void printFiberVTK(Mesh *mesh, vector<Vector>& fiber_vecs, std::ostream &out);
void printAnatomy(vector<anatomy>& anatVectors, filerheader& header, std::ostream &out); 

#endif	/* IO_H */

