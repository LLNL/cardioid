/* 
 * File:   utils.h
 * Author: zhang30
 *
 * Created on February 22, 2017, 11:08 PM
 */

#ifndef UTILS_H
#define	UTILS_H

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

double a_s_f(double a_endo, double a_epi, double d);
double a_w_f(double a_endo, double a_epi, double d);
double b_s_f(double b_endo, double b_epi, double d);
double b_w_f(double b_endo, double b_epi, double d);

bool vecisnonzero(Vector& vec);
bool vecdot(Vector &q1, Vector &q2);
bool vecisnorm(Vector &q);

#endif	/* UTILS_H */

