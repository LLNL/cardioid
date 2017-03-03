/* 
 * File:   genfiber.h
 * Author: zhang30
 *
 * Created on February 22, 2017, 11:25 PM
 */

#ifndef GENFIBER_H
#define	GENFIBER_H

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

using namespace std;
using namespace mfem;

void bislerp(DenseMatrix& Q, DenseMatrix& Qa, DenseMatrix& Qb, double t);
void axis(DenseMatrix& Q,Vector &psi, Vector &phi);
void orient(DenseMatrix& Qp, DenseMatrix& Q, double a, double b);
void vectorEigen(Vector& psi_ab, DenseMatrix& QPfib);
void genfiber(Mesh *mesh, vector<DenseMatrix>& QPfibVectors,
        vector<double>& psi_ab, vector<Vector>& psi_ab_grads,
        vector<double>& phi_epi, vector<Vector>& phi_epi_grads,
        vector<double>& phi_lv, vector<Vector>& phi_lv_grads,
        vector<double>& phi_rv, vector<Vector>& phi_rv_grads,
        double a_endo, double a_epi, double b_endo, double b_epi
        );


#endif	/* GENFIBER_H */

