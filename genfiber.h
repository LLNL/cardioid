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

#endif	/* GENFIBER_H */

