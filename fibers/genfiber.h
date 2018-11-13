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

#include "option.h"

using namespace std;
using namespace mfem;

void bislerp(DenseMatrix& Q, DenseMatrix& Qa, DenseMatrix& Qb, double t);
void axis(DenseMatrix& Q,Vector &psi, Vector &phi);
void orient(DenseMatrix& Qp, DenseMatrix& Q, double a, double b);
void vectorEigen(Vector& psi_ab, DenseMatrix& QPfib);
void biSlerpCombo(DenseMatrix& QPfib,
        double psi_ab, Vector& psi_ab_vec,
        double phi_epi, Vector& phi_epi_vec,
        double phi_lv, Vector& phi_lv_vec, 
        double phi_rv, Vector& phi_rv_vec,  
        Option& options);
void genfiber(vector<DenseMatrix>& QPfibVectors,
        vector<double>& psi_ab, vector<Vector>& psi_ab_grads,
        vector<double>& phi_epi, vector<Vector>& phi_epi_grads,
        vector<double>& phi_lv, vector<Vector>& phi_lv_grads,
        vector<double>& phi_rv, vector<Vector>& phi_rv_grads,
        Option& options);


#endif	/* GENFIBER_H */

