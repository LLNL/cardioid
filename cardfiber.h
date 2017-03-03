/* 
 * File:   cardfiber.h
 * Author: zhang30
 *
 * Created on March 2, 2017, 4:45 PM
 */

#ifndef CARDFIBER_H
#define	CARDFIBER_H

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include "triplet.h"

using namespace std;
using namespace mfem;

void buildKDTree(Mesh *mesh, tree_type& kdtree, vector<Vector>& boundingbox, double dd);

#endif	/* CARDFIBER_H */

