//                                MFEM Example 1
//
// Compile with: make ex1
//
// Sample runs:  ex1 -m ../data/square-disc.mesh
//               ex1 -m ../data/star.mesh
//               ex1 -m ../data/escher.mesh
//               ex1 -m ../data/fichera.mesh
//               ex1 -m ../data/square-disc-p2.vtk -o 2
//               ex1 -m ../data/square-disc-p3.mesh -o 3
//               ex1 -m ../data/square-disc-nurbs.mesh -o -1
//               ex1 -m ../data/disc-nurbs.mesh -o -1
//               ex1 -m ../data/pipe-nurbs.mesh -o -1
//               ex1 -m ../data/star-surf.mesh
//               ex1 -m ../data/square-disc-surf.mesh
//               ex1 -m ../data/inline-segment.mesh
//               ex1 -m ../data/amr-quad.mesh
//               ex1 -m ../data/amr-hex.mesh
//               ex1 -m ../data/fichera-amr.mesh
//               ex1 -m ../data/mobius-strip.mesh
//               ex1 -m ../data/mobius-strip.mesh -o -1 -sc
//
// Description:  This example code demonstrates the use of MFEM to define a
//               simple finite element discretization of the Laplace problem
//               -Delta u = 1 with homogeneous Dirichlet boundary conditions.
//               Specifically, we discretize using a FE space of the specified
//               order, or if order < 1 using an isoparametric/isogeometric
//               space (i.e. quadratic for quadratic curvilinear mesh, NURBS for
//               NURBS mesh, etc.)
//
//               The example highlights the use of mesh refinement, finite
//               element grid functions, as well as linear and bilinear forms
//               corresponding to the left-hand side and right-hand side of the
//               discrete linear system. We also cover the explicit elimination
//               of essential boundary conditions, static condensation, and the
//               optional connection to the GLVis tool for visualization.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include "io.h"
#include "solver.h"
#include "constants.h"
#include "utils.h"
#include "genfiber.h"

using namespace std;
using namespace mfem;

int main(int argc, char *argv[]) {
    // 1. Parse command-line options.
    const char *mesh_file = "./human.vtk";
    int order = 1;
    bool static_cond = false;
    bool visualization = 1;

    double a_endo=40;
    double a_epi=-50;
    double b_endo=-65;
    double b_epi=25;
    
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
            "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
            "Finite element order (polynomial degree) or -1 for"
            " isoparametric space.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
            "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
            "--no-visualization",
            "Enable or disable GLVis visualization.");
    args.Parse();
    if (!args.Good()) {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
    //    the same code.
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    
    // Set the surfaces for the mesh: 0-Apex, 1-Base, 2-EPI, 3-LV, 4-RV.
    setSurfaces(mesh, 20); // use 30 degrees for determining the base surface. 
    ofstream surf_ofs("surfaces.vtk");
    printSurfVTK(mesh, surf_ofs);

    // 3. Solve the laplacian for four different boundary conditions.
    
    // get the vertex elements arrays.
    vector<vector<int> > vert2Elements;
    getVert2Elements(mesh, vert2Elements);
    ofstream v2e_ofs("vert2Elements.txt");
    for(int i=0; i<vert2Elements.size(); i++){
        vector<int> elements=vert2Elements[i];
        v2e_ofs << i  << " ";
        for(int j=0; j<elements.size(); j++){
            v2e_ofs << elements[j] << " ";
        }
        v2e_ofs << endl;
    }
    
    
    int bdr_attr_size=mesh->bdr_attributes.Max();
    Array<int> all_ess_bdr(bdr_attr_size);    
    Array<int> nonzero_ess_bdr(bdr_attr_size);
    Array<int> zero_ess_bdr(bdr_attr_size);
    int nv=mesh->GetNV();
    
    // 3a. Base → 1, Apex→ 0, Epi, LV, RV → no flux
     // Mark ALL boundaries as essential. This does not set what the actual Dirichlet
    // values are
    all_ess_bdr = 1;
    all_ess_bdr[2]=0;
    all_ess_bdr[3]=0;
    all_ess_bdr[4]=0;
    
    nonzero_ess_bdr = 0;    
    nonzero_ess_bdr[1] = 1;   

    zero_ess_bdr = 0;     
    zero_ess_bdr[0] = 1;
    
    string output="psi_ab";
    vector<double> psi_ab;
    vector<Vector> psi_ab_grads;
    laplace(mesh, vert2Elements, psi_ab, psi_ab_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
    MFEM_ASSERT(psi_ab.size()==nv, "size of psi_ab does not match number of vertices.");
    MFEM_ASSERT(psi_ab_grads.size()==nv, "size of psi_ab_grads does not match number of vertices.");
    
    
    // 3b. Apex, Epi → 1, LV, RV→ 0, Base→ no flux
    all_ess_bdr = 1;
    all_ess_bdr[1]=0;
    
    nonzero_ess_bdr = 0;    
    nonzero_ess_bdr[0] = 1;
    nonzero_ess_bdr[2] = 1;   

    zero_ess_bdr = 0;      
    zero_ess_bdr[3] = 1;
    zero_ess_bdr[4] = 1;
 
    output="phi_epi";
    vector<double> phi_epi;
    vector<Vector> phi_epi_grads;

    laplace(mesh, vert2Elements, phi_epi, phi_epi_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
    MFEM_ASSERT(phi_epi.size()==nv, "size of phi_epi does not match number of vertices.");
    MFEM_ASSERT(phi_epi_grads.size()==nv, "size of phi_epi_grads does not match number of vertices.");
    
    //3c. LV → 1, Apex, Epi, RV→ 0, Base→ no flux
    all_ess_bdr = 1;
    all_ess_bdr[1]=0;
    
    nonzero_ess_bdr = 0;    
    nonzero_ess_bdr[3] = 1;   

    zero_ess_bdr = 0;      
    zero_ess_bdr[0] = 1;
    zero_ess_bdr[2] = 1;
    zero_ess_bdr[4] = 1;
 
    output="phi_lv";
    vector<double> phi_lv;
    vector<Vector> phi_lv_grads;

    laplace(mesh, vert2Elements, phi_lv, phi_lv_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
    MFEM_ASSERT(phi_lv.size()==nv, "size of phi_lv does not match number of vertices.");
    MFEM_ASSERT(phi_lv_grads.size()==nv, "size of phi_lv_grads does not match number of vertices.");        
    
    //3d. RV → 1, Apex, Epi, LV→ 0, Base→ no flux
    all_ess_bdr = 1;
    all_ess_bdr[1]=0;
    
    nonzero_ess_bdr = 0;    
    nonzero_ess_bdr[4] = 1;   

    zero_ess_bdr = 0;      
    zero_ess_bdr[0] = 1;
    zero_ess_bdr[2] = 1;
    zero_ess_bdr[3] = 1;
 
    output="phi_rv";
    vector<double> phi_rv;
    vector<Vector> phi_rv_grads;

    laplace(mesh, vert2Elements, phi_rv, phi_rv_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
    MFEM_ASSERT(phi_rv.size()==nv, "size of phi_rv does not match number of vertices.");
    MFEM_ASSERT(phi_rv_grads.size()==nv, "size of phi_rv_grads does not match number of vertices.");
    
    ofstream psia_ofs("psi_ab_grads.vtk");
    ofstream phie_ofs("phi_epi_grads.vtk");
    ofstream phil_ofs("phi_lv_grads.vtk");
    ofstream phir_ofs("phi_rv_grads.vtk");

    printFiberVTK(mesh, psi_ab_grads, psia_ofs);
    printFiberVTK(mesh, phi_epi_grads, phie_ofs);
    printFiberVTK(mesh, phi_lv_grads, phil_ofs);    
    printFiberVTK(mesh, phi_rv_grads, phir_ofs); 
            
    vector<Vector> fvectors;
    vector<Vector> svectors;
    vector<Vector> tvectors;
    // Line 7 start for-loop
    for(int i=0; i <nv; i++){  
//        MFEM_ASSERT(phi_lv[i]>=0 && phi_lv[i] <=1, "phi_lv is not in range 0 to 1");
//        MFEM_ASSERT(phi_rv[i]>=0 && phi_rv[i] <=1, "phi_rv is not in range 0 to 1");
//        MFEM_ASSERT(phi_epi[i]>=0 && phi_epi[i] <=1, "phi_epi is not in range 0 to 1");
//        MFEM_ASSERT(psi_ab[i]>=0 && psi_ab[i] <=1, "psi_ab is not in range 0 to 1");
        //if(phi_lv[i] <0) phi_lv[i]=0;
        //if(phi_rv[i] <0) phi_rv[i]=0;
        //if(phi_epi[i] <0) phi_epi[i]=0;
        //if(psi_ab[i] <0) psi_ab[i]=0;

        double phi_v=phi_lv[i]+phi_rv[i];
        double frac=0.5;
        if(phi_v!=0){
            frac=phi_rv[i]/phi_v;
        }else{
            cout << "Warning: phi_v ==0" ;
            cout << " phi_lv[i]="<< phi_lv[i] << " phi_rv[i]=" << phi_rv[i]<< " phi_epi[i]=" << phi_epi[i] << " psi_ab[i]=" << psi_ab[i]<< endl;
        }
        double frac_epi=phi_epi[i];
        //stringstream ss;
        //ss << "i=" << i << " phi_rv[i]=" << phi_rv[i] << " phi_lv[i]=" << phi_lv[i] << " frac=" << frac;
        //MFEM_ASSERT(frac>=0 && frac<=1, ss.str());
        //MFEM_ASSERT(frac_epi>=0 && frac_epi<=1, "frac_epi is not in range 0 to 1");
        double as=a_s_f(a_endo, a_epi, frac);
        double bs=b_s_f(b_endo, b_epi, frac);
        double aw=a_w_f(a_endo, a_epi, frac_epi);
        double bw=b_w_f(b_endo, b_epi, frac_epi);
        

        Vector psi_ab_vec=psi_ab_grads[i];
        Vector phi_lv_vec=phi_lv_grads[i];
        Vector phi_rv_vec=phi_rv_grads[i];
        Vector phi_epi_vec=phi_epi_grads[i];
        
        bool phi_lv_isnonzero=vecisnonzero(phi_lv_vec);
        bool phi_rv_isnonzero=vecisnonzero(phi_rv_vec);
        bool phi_epi_isnonzero=vecisnonzero(phi_epi_vec);

        DenseMatrix QPendo(dim3,dim3);
        DenseMatrix QPfib(dim3,dim3);

        if (!vecisnonzero(psi_ab_vec)) {
            cout << "Warning psi_ab gradient " << i << "is zero" <<endl;
            Vector ten(3);
            ten = 10;
//            for(int i=0; i<dim3; i++){
//                QPfib.SetCol(i,ten);
//            }
            
            fvectors.push_back(ten);
            svectors.push_back(ten);
            tvectors.push_back(ten);
            continue;
        }
        
        
        DenseMatrix QPlv(dim3,dim3);        
        if(phi_lv_isnonzero){
            // Line 8
            Vector phi_lv_vec_neg=phi_lv_vec;
            phi_lv_vec_neg.Neg();
            DenseMatrix Qlv(dim3,dim3);
            if(vecdot(psi_ab_vec, phi_lv_vec_neg)){
                cout << "psi_ab_vec equal to phi_lv_vec_neg" << endl;
                phi_lv_isnonzero=false;
            }else{
                axis(Qlv, psi_ab_vec, phi_lv_vec_neg);
                orient(QPlv, Qlv, as, bs);
            }            
            // End of Line 8
        }
        
        DenseMatrix QPrv(dim3,dim3);         
        if(phi_rv_isnonzero){
            //Line 9
            DenseMatrix Qrv(dim3,dim3);
            if(vecdot(psi_ab_vec, phi_rv_vec)){
                cout << "psi_ab_vec equal to phi_rv_vec" << endl;
                phi_rv_isnonzero=false;
            }else{
                axis(Qrv, psi_ab_vec, phi_rv_vec);
                orient(QPrv, Qrv, as, bs); 
            }
        }
               
        DenseMatrix QPepi(dim3,dim3);
        if (phi_epi_isnonzero) {
            //Line 11
            DenseMatrix Qepi(dim3,dim3);
            if(vecdot(psi_ab_vec, phi_epi_vec)){
                cout << "psi_ab_vec equal to phi_epi_vec" << endl;
                phi_epi_isnonzero=false;
            }else{           
                axis(Qepi, psi_ab_vec, phi_epi_vec);
                orient(QPepi, Qepi, aw, bw);
            }
        }
                
        if(phi_lv_isnonzero){    
            
            if(phi_rv_isnonzero){
            
                if(phi_epi_isnonzero){
                    // if all three phi gradients are non-zero, use the original algorithm in paper. 
                    //Line 10
                    bislerp(QPendo, QPlv, QPrv, frac);
                    //QPendo=QPlv;
                    //Line 12 
                    bislerp(QPfib, QPendo, QPepi, frac_epi);
                    //QPfib=QPendo;
                }else{
                    // if phi_epi gradients are zero, phi_lv and phi are nonzero. use QPlv, QPrv, frac 
                    bislerp(QPfib, QPlv, QPrv, frac);
                    //QPfib=QPlv;
                }
                
            }else {
                if(phi_epi_isnonzero){
                    // if phi_rv is zero, phi_lv and phi_epi is nonzero
                    bislerp(QPfib, QPlv, QPepi, frac_epi);
                    //QPfib=QPlv;
                }else{
                    // if gradients of phi_lv, phi_rv are zero, then phi_epi is zero 
                    vectorEigen(psi_ab_vec, QPfib);
                }
            }
        }else{
            if(phi_rv_isnonzero && phi_epi_isnonzero){                
                // if phi_lv is zero, phi_rv and phi_epi is nonzero
                bislerp(QPfib, QPrv, QPepi, frac_epi);
                //QPfib=QPrv;
            }else{
                // if gradients of phi_lv, phi_rv are zero, then phi_epi is zero 
                vectorEigen(psi_ab_vec, QPfib);
            }
        }
                
        vector<Vector> qpVecs;
        for(int j=0; j<dim3; j++){
            Vector vec;
            QPfib.GetColumn(j, vec);
            qpVecs.push_back(vec);
        }
        fvectors.push_back(qpVecs[0]);
        svectors.push_back(qpVecs[1]);
        tvectors.push_back(qpVecs[2]);
    }
    
    ofstream f_ofs("fvectors.vtk");
    ofstream s_ofs("svectors.vtk");
    ofstream t_ofs("tvectors.vtk");

    printFiberVTK(mesh, fvectors, f_ofs);
    printFiberVTK(mesh, svectors, s_ofs);
    printFiberVTK(mesh, tvectors, t_ofs);
        
    delete mesh;

    return 0;
}

