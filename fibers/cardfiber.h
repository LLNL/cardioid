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
#include "option.h"

using namespace std;
using namespace mfem;

struct filerheader{
    int nrecord;
    int nx;
    int ny;
    int nz;
    double dx;
    double dy;
    double dz;
    double offset_x;
    double offset_y;
    double offset_z;
};

struct anatomy{
    long  gid;
    int celltype;
    double sigma[6];

};

struct Phi{
    double epi;
    double lv;
    double rv;
};

struct ThreeInts{
    int i;
    int j;
    int k;
};

void buildKDTree(Mesh *mesh, tree_type& kdtree);

double getMaxEdgeLen(Mesh *mesh);

double det4X4(DenseMatrix& matrix);
bool isInTetElement(const Vector& q, Mesh* mesh, int eleIndex);
void getCardEleGrads(GridFunction& x, const Vector& q, int eleIndex, Vector& grad_ele, double& xVal);
void calcSigma(DenseMatrix& Sigma, DenseMatrix& Q, Option& options);
int getCellType(Phi& phi);

void calcGradient(GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        Option& options, Vector& q, int eleIndex, DenseMatrix& QPfib, Phi& phi);

void getAnatomy(anatomy& anat, DenseMatrix& QPfib, Option& options, Phi& phi,
        ThreeInts& inds, ThreeInts& nns);

bool findPtEle(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        vector<vector<int> >& vert2Elements, Option& options, Vector& q, int vertex, std::string& elemnum, ostream& out);

bool findPtEleAnat(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        vector<vector<int> >& vert2Elements, Option& options,
        Vector& q, int vertex, ThreeInts& inds, ThreeInts& nns, vector<anatomy>& anatVectors);

#endif	/* CARDFIBER_H */
