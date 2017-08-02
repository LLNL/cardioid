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
#include <deque>
#include <limits>
#include <functional>
#include <set>

#include "kdtree++/kdtree.hpp"

#include "io.h"
#include "solver.h"
#include "genfiber.h"
#include "cardfiber.h"
#include "cardgradients.h"
#include "triplet.h"

using namespace std;
using namespace mfem;

int main(int argc, char *argv[]) {
    // 1. Parse command-line options.
    const char *mesh_file = "./human.vtk";
    int order = 1;
    bool static_cond = false;
    bool visualization = 1;
    // run omar's task
    bool omar_task=false;
    bool omar_fast=false;
    const char *fiblocs="";

    double a_endo=40;
    double a_epi=-50;
    double b_endo=-65;
    double b_epi=25;
    
    // grid spacing
    double dd=5;
    // conductivity
    double gL = 0.0001334177*1000; // mS/mm
    double gT = 0.0000176062*1000; // mS/mm
    double gN = 0.0000176062*1000; // mS/mm   

    // cutoff for kdtree point range search rangeCutoff=rcut*maxEdgeLen
    double rcut=1.0;
       
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
    args.AddOption(&omar_task, "-omar", "--omar_task", "-no-omar",
            "--no-omar_task",
            "Enable or disable Omar task.");  
    args.AddOption(&omar_fast, "-ofast", "--omar_fast", "-no-ofast",
            "--no-omar_fast",
            "Enable or disable Omar fast task.");    
    args.AddOption(&fiblocs, "-fl", "--fiblocs",
            "Fiber locagtion file to use.");    
    args.AddOption(&a_endo, "-ao", "--aendo", "Fiber angle alpha endo.");
    args.AddOption(&a_epi, "-ai", "--aepi", "Fiber angle alpha epi.");
    args.AddOption(&b_endo, "-bo", "--bendo", "Fiber angle beta endo.");
    args.AddOption(&b_epi, "-bi", "--bepi", "Fiber angle beta epi."); 
    args.AddOption(&dd, "-dd", "--dspacing", "Grid spacing for ddcMD gid.");
    args.AddOption(&gL, "-gl", "--gL", "Conductivity gL mS/mm.");
    args.AddOption(&gT, "-gt", "--gT", "Conductivity gT mS/mm.");
    args.AddOption(&gN, "-gn", "--gN", "Conductivity gN mS/mm.");
    args.AddOption(&rcut, "-rc", "--rcut", "rangeCutoff=rcut*maxEdgeLen (default rcut=1.0).");
    args.Parse();
    if (!args.Good()) {
        args.PrintUsage(cout);
        return 1;
    }
    if(omar_task && strlen(fiblocs)==0){
       std::cout << "The program runs Omar task but missing fiber location file!" << std::endl;
       args.PrintUsage(cout);
       return 1; 
    }
    
    if(omar_fast && strlen(fiblocs)==0){
       std::cout << "The program runs Omar fast task but missing fiber location file!" << std::endl;
       args.PrintUsage(cout);
       return 1;
    }    
    args.PrintOptions(cout);

    // Keep fiber angles in a Vector.
    Vector fiberAngles(4);
    fiberAngles(0)=a_endo;
    fiberAngles(1)=a_epi;
    fiberAngles(2)=b_endo;
    fiberAngles(3)=b_epi;    
    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
    //    the same code.
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    
    vector<Vector> boundingbox;
    // Set the surfaces for the mesh: 0-Apex, 1-Base, 2-EPI, 3-LV, 4-RV.
    setSurfaces(mesh, boundingbox, 20); // use 30 degrees for determining the base surface.
//    for(unsigned i = 0; i < boundingbox.size(); i++) {
//        Vector vec=boundingbox[i];     
//        cout << "Bounding Box " << i;
//        for(int j = 0; j < vec.Size(); j++) {
//            cout << " " << vec(j);
//        }
//        cout << endl;
//    }
    ofstream surf_ofs("surfaces.vtk");
    printSurfVTK(mesh, surf_ofs);

    // 3. Solve the laplacian for four different boundary conditions.
    
    // get the vertex elements arrays.    
    vector<vector<int> > vert2Elements;
    getVert2Elements(mesh, vert2Elements);
    ofstream v2e_ofs("vert2Elements.txt");
    for(unsigned i=0; i<vert2Elements.size(); i++){
        vector<int> elements=vert2Elements[i];
        v2e_ofs << i  << " ";
        for(unsigned j=0; j<elements.size(); j++){
            v2e_ofs << elements[j] << " ";
        }
        v2e_ofs << endl;
    }
    
    
    int bdr_attr_size=mesh->bdr_attributes.Max();
    Array<int> all_ess_bdr(bdr_attr_size);    
    Array<int> nonzero_ess_bdr(bdr_attr_size);
    Array<int> zero_ess_bdr(bdr_attr_size);
    unsigned nv=mesh->GetNV();
    
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
    GridFunction x_psi_ab=laplace(mesh, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, order, static_cond);
    getVetecesGradients(mesh, x_psi_ab, vert2Elements, psi_ab,psi_ab_grads, output);
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

    //laplace(mesh, vert2Elements, phi_epi, phi_epi_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
    GridFunction x_phi_epi=laplace(mesh, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, order, static_cond);
    getVetecesGradients(mesh, x_phi_epi, vert2Elements, phi_epi,phi_epi_grads, output);
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

    //laplace(mesh, vert2Elements, phi_lv, phi_lv_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
    GridFunction x_phi_lv=laplace(mesh, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, order, static_cond);
    getVetecesGradients(mesh, x_phi_lv, vert2Elements, phi_lv,phi_lv_grads, output);    
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

    //laplace(mesh, vert2Elements, phi_rv, phi_rv_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
    GridFunction x_phi_rv=laplace(mesh, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, order, static_cond);
    getVetecesGradients(mesh, x_phi_rv, vert2Elements, phi_rv,phi_rv_grads, output);
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
    
    vector<DenseMatrix> QPfibVectors;   
    genfiber(QPfibVectors, psi_ab, psi_ab_grads, phi_epi, phi_epi_grads, 
        phi_lv, phi_lv_grads, phi_rv, phi_rv_grads, fiberAngles);
    
    cout << "\nCalculate fiber on the nodes ...\n";
    calcNodeFiber(QPfibVectors); 
        
    vector<Vector> fvectors;
    vector<Vector> svectors;
    vector<Vector> tvectors;
    
    for(unsigned i=0; i< QPfibVectors.size(); i++){
        vector<Vector> qpVecs;
        for(int j=0; j<3; j++){
            Vector vec;
            QPfibVectors[i].GetColumn(j, vec);
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
    
    if(omar_fast){
      
      cout << "\nGet Omar's rotation matrix in fast way ...\n";
      getRotMatrixFast(mesh, x_psi_ab, x_phi_epi, x_phi_lv, x_phi_rv,
          vert2Elements, fiberAngles, fiblocs);  
      
      delete mesh;

      return 0;       
    }
    
    double maxEdgeLen=getMaxEdgeLen(mesh);
    cout << "\nThe maximum edge length in the mesh is "<< maxEdgeLen <<"\n";
    cout.flush();     
      
    cout << "\nStart to build k-D tree for the mesh...\n";
    cout.flush();
    tree_type kdtree(std::ptr_fun(tac));
    buildKDTree(mesh, kdtree);
    
    if(omar_task){
       cout << "\nGet Omar's rotation matrix ...\n";
       double rangeCutoff=rcut*maxEdgeLen;
       getRotMatrix(mesh, x_psi_ab, x_phi_epi, x_phi_lv, x_phi_rv,
          kdtree, vert2Elements, fiberAngles, rangeCutoff, fiblocs);
       
    }
    
    cout << "\nGet cardioid point gradients ...\n";
    cout.flush();
    Vector conduct(3);
    conduct(0)=gL;
    conduct(1)=gT;
    conduct(2)=gN;
    getCardGradients(mesh, x_psi_ab, x_phi_epi, x_phi_lv, x_phi_rv,
        kdtree, vert2Elements, boundingbox, dd,  
        conduct, fiberAngles, maxEdgeLen);
         
    delete mesh;

    return 0;
}

