/* 
 * File:   cardgradientsp.h
 * Author: zhang30
 *
 * Created on March 6, 2017, 5:25 PM
 */

#ifndef CARDGRADIENTSP_H
#define	CARDGRADIENTSP_H

#include "option.h"

using namespace std;
using namespace mfem;

void getCardGradientsp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, vector<Vector>& boundingbox, Option& options, int num_procs, int myid);

void getRotMatrixp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, Option& options, int num_procs, int myid);

void getRotMatrixFastp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
                   vector<vector<int> >& vert2Elements, Option& options, int size, int rank);

void calcNodeFiberP(vector<DenseMatrix>& QPfibVectors, int num_procs, int myid);
#endif	/* CARDGRADIENTSP_H */

