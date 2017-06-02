/* 
 * File:   cardgradientsp.h
 * Author: zhang30
 *
 * Created on March 6, 2017, 5:25 PM
 */

#ifndef CARDGRADIENTSP_H
#define	CARDGRADIENTSP_H

using namespace std;
using namespace mfem;

void getCardGradientsp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, vector<Vector>& boundingbox, double dd, 
        Vector& conduct, Vector& fiberAngles, int num_procs, int myid);

void getRotMatrixp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, Vector& fiberAngles, const char *fiblocs, int num_procs, int myid);

void getRotMatrixFastp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
                   vector<vector<int> >& vert2Elements, Vector& fiberAngles, const char *fiblocs, int size, int rank);

#endif	/* CARDGRADIENTSP_H */

