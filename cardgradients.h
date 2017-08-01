/* 
 * File:   cardgradients.h
 * Author: zhang30
 *
 * Created on March 6, 2017, 5:15 PM
 */

#ifndef CARDGRADIENTS_H
#define	CARDGRADIENTS_H

using namespace std;
using namespace mfem;

void getCardGradients(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, vector<Vector>& boundingbox, double dd, 
        Vector& conduct, Vector& fiberAngles, double maxEdgeLen);

void getRotMatrix(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, Vector& fiberAngles, const char *fiblocs);

void getRotMatrixFast(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        vector<vector<int> >& vert2Elements, Vector& fiberAngles, const char *fiblocs);

void calcNodeFiber(vector<DenseMatrix>& QPfibVectors);

#endif	/* CARDGRADIENTS_H */

