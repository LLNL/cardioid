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

struct filerheader{
    long nrecord;
    int nx;
    int ny;
    int nz;
    double dx;
    double dy;
    double dz;
};

struct anatomy{
    long  gid;
    int celltype;
    double sigma[6];
    
};
void buildKDTree(Mesh *mesh, tree_type& kdtree);
void getCardGradients(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, vector<Vector>& boundingbox, double dd, 
        Vector& conduct, Vector& fiberAngles);

#endif	/* CARDFIBER_H */

