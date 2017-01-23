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

using namespace std;
using namespace mfem;

void setSurfaces(Mesh *mesh);
void printSurfVTK(Mesh *mesh, std::ostream &out);

int main(int argc, char *argv[]) {
    // 1. Parse command-line options.
    const char *mesh_file = "./mechmesh.vtk";
    int order = 1;
    bool static_cond = false;
    bool visualization = 1;

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
    
    setSurfaces(mesh);
/*    
    ofstream surfh;
    surfh.open("surfaces.vtk");    
    printSurfVTK(mesh, surfh);
    //return 0;
      
    ofstream mfh;
    mfh.open("test.mesh");
    mesh->Print(mfh);
    return 0;
 
    ofstream vfh;
    vfh.open("test.vtk");
    mesh->PrintVTK(vfh);
    //return 0;
*/
    
    int dim = mesh->Dimension();
    cout << "Dimension =" << dim << endl;

    //delete mesh;

    //return 0;

    // 3. Refine the mesh to increase the resolution. In this example we do
    //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
    //    largest number that gives a final mesh with no more than 50,000
    //    elements.
   /*
    {
        int ref_levels =
                (int) floor(log(50000. / mesh->GetNE()) / log(2.) / dim);
        for (int l = 0; l < ref_levels; l++) {
            mesh->UniformRefinement();
        }
    }
    */
    
    // 4. Define a finite element space on the mesh. Here we use continuous
    //    Lagrange finite elements of the specified order. If order < 1, we
    //    instead use an isoparametric/isogeometric space.
    FiniteElementCollection *fec;
    if (order > 0) {
        fec = new H1_FECollection(order, dim);
    } else if (mesh->GetNodes()) {
        fec = mesh->GetNodes()->OwnFEC();
        cout << "Using isoparametric FEs: " << fec->Name() << endl;
    } else {
        fec = new H1_FECollection(order = 1, dim);
    }
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
    cout << "Number of finite element unknowns: "
            << fespace->GetTrueVSize() << endl;

    // 5. Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking all
    //    the boundary attributes from the mesh as essential (Dirichlet) and
    //    converting them to a list of true dofs.
    Array<int> ess_tdof_list;
    MFEM_ASSERT(mesh->bdr_attributes.Size()!=0, "Boundary size cannot be zero."); 


    // Mark ALL boundaries as essential. This does not set what the actual Dirichlet
    // values are
    Array<int> all_ess_bdr(mesh->bdr_attributes.Max());
    all_ess_bdr = 1;
    all_ess_bdr[1]=0;

    cout << "all_ess_bdr size=" << all_ess_bdr.Size() << endl;
    fespace->GetEssentialTrueDofs(all_ess_bdr, ess_tdof_list);
  

    // 6. Set up the linear form b(.) which corresponds to the right-hand side of
    //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
    //    the basis functions in the finite element fespace.
    LinearForm *b = new LinearForm(fespace);
    ConstantCoefficient zero(0.0);
    b->AddDomainIntegrator(new DomainLFIntegrator(zero));
    b->Assemble();

    // 7. Define the solution vector x as a finite element grid function
    //    corresponding to fespace. Initialize x with initial guess of zero,
    //    which satisfies the boundary conditions.
    GridFunction x(fespace);
    x = 0.0;
    
    cout << "x size " << x.Size() << endl;


    // 8. Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
    //    domain integrator.
    BilinearForm *a = new BilinearForm(fespace);

    // The diffusion integrator should have a coefficient of one, not zero
    ConstantCoefficient one(1.0);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));

    // 9. Assemble the bilinear form and the corresponding linear system,
    //    applying any necessary transformations such as: eliminating boundary
    //    conditions, applying conforming constraints for non-conforming AMR,
    //    static condensation, etc.
    if (static_cond) {
        a->EnableStaticCondensation();
    }
    a->Assemble();

    // Project the constant 14 value to all boundary attributes except 1
    Array<int> nonzero_ess_bdr(mesh->bdr_attributes.Max());
    nonzero_ess_bdr = 0;    
    nonzero_ess_bdr[4] = 1;
    //nonzero_ess_bdr[2] = 1;
    ConstantCoefficient nonzero_bdr(1.0);
    x.ProjectBdrCoefficient(nonzero_bdr, nonzero_ess_bdr);

    // Project the constant 0 value to boundary attribute 1
    Array<int> zero_ess_bdr(mesh->bdr_attributes.Max());
    zero_ess_bdr = 0;    
    //zero_ess_bdr[0] = 0;  
    zero_ess_bdr[0] = 1; 
    zero_ess_bdr[2] = 1;
    zero_ess_bdr[3] = 1;
    //zero_ess_bdr[2] = 0;
    ConstantCoefficient zero_bdr(0.0);
    x.ProjectBdrCoefficient(zero_bdr, zero_ess_bdr);

    
    int count=0;
    for(int i=0; i<x.Size(); i++){
        if(x[i] >0.01){
            //cout << "["<< i << "]=" <<x[i] << " ";
            count++;
        }
    }
    cout << "\nCount= " << count <<endl;
        
    //return 0;
    SparseMatrix A;
    Vector B, X;
    // Form the linear system using ALL of the essential boundary dofs (from all_ess_bdr)
    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

    cout << "Size of linear system: " << A.Height() << endl;

#ifndef MFEM_USE_SUITESPARSE
    // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to
    //     solve the system A X = B with PCG.
    GSSmoother M(A);
    PCG(A, M, B, X, 1, 1000, 1e-12, 0.0);
#else
    // 10. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
    UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(A);
    umf_solver.Mult(B, X);
#endif

    // 11. Recover the solution as a finite element grid function.
    a->RecoverFEMSolution(X, *b, x);

    // Create data collection for solution output: either VisItDataCollection for
    // ascii data files, or SidreDataCollection for binary data files.
    DataCollection *dc = NULL;
    dc = new VisItDataCollection("Fiber", mesh);
    dc->SetPrecision(8);
    dc->RegisterField("solution", &x);
    dc->SetCycle(0);
    dc->SetTime(0.0);
    dc->Save();

    // 12. Save the refined mesh and the solution. This output can be viewed later
    //     using GLVis: "glvis -m refined.mesh -g sol.gf".
    ofstream mesh_ofs("refined.mesh");
    mesh_ofs.precision(8);
    mesh->Print(mesh_ofs);
    
    ofstream vtk_ofs("refined.vtk");
    mesh->PrintVTK(vtk_ofs);
    
    ofstream sol_ofs("sol.gf");
    sol_ofs.precision(8);
    x.Save(sol_ofs);

    // 13. Send the solution by socket to a GLVis server.
    if (visualization) {
        char vishost[] = "localhost";
        int visport = 19916;
        socketstream sol_sock(vishost, visport);
        sol_sock.precision(8);
        sol_sock << "solution\n" << *mesh << x << flush;
    }

    // 14. Free the used memory.
    delete a;
    delete b;
    delete fespace;
    if (order > 0) {
        delete fec;
    }
    delete mesh;

    return 0;
}

bool isPlanar(double *coor0, double *coor1, double *coor2, double cosTheta){
    double u[3];
    double v[3];
    double w[3];
    u[0]=coor1[0]-coor0[0];
    u[1]=coor1[1]-coor0[1];
    u[2]=coor1[2]-coor0[2];
    v[0]=coor2[0]-coor0[0];
    v[1]=coor2[1]-coor0[1];
    v[2]=coor2[2]-coor0[2];
    
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];
    
    double r=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    double cosT=abs(w[2]/r);
    
    if(cosT>cosTheta){
        return true;
    }
    
    return false;
    
}

bool isTriInTet(vector<int>& tri, vector<int>& tet){
    MFEM_ASSERT(tri.size() == 3, "Wrong boundary size");
    MFEM_ASSERT(tet.size() == 4, "Wrong tetrahedral size");

    if(tri[0]==tet[0]){
        if(tri[1]==tet[1] && tri[2]==tet[2]){
            return true;
        }        
    }else if(tri[0]==tet[1]){
        if(tri[1]==tet[2] && tri[2]==tet[3]){
            return true;
        }          
    }    
    return false;
}

void findNeighbor(Element* ele, vector<Element*>& elements, int attr){
    const int *v = ele->GetVertices();
    const int nv = ele->GetNVertices(); 
    for(int i=0; i<elements.size(); i++){
        Element* queryEle=elements[i];
        // Only search for elements with unassigned attributes.
        if(queryEle->GetAttribute()==0){ 
            const int *qv = queryEle->GetVertices();
            const int nqv = queryEle->GetNVertices(); 
            bool isNeighbor=false;
            for (int j = 0; j < nv; j++) {
                for (int k = 0; k < nqv; k++) {
                    // If two elements share the same vertex they are neighbor. 
                    if(v[j]==qv[k]){                   
                        isNeighbor=true;
                        // Should break two loops can use lambda or function return.
                        break;  
                    }
                }
            }
            if(isNeighbor){
                queryEle->SetAttribute(attr);
                // recursively search for neighboring elements.
                findNeighbor(queryEle, elements, attr);
            }
        }
    }           
}

void setSurfaces(Mesh *mesh){
    // Attributes for different surface
    const int apexAttr=1;
    const int baseAttr=2;
    const int epiAttr=3;
    const int lvAttr=4;
    const int rvAttr=5;
       
    // Determine the max and min dimension of mesh and apex.
    double *coord;
    double coord_min[3];
    double coord_max[3];
    bool firstEle=true;
    int apexVet=0;
    int apexEleIndex=0;
    
    int nbe=mesh->GetNBE();
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        // The first loop has to initialize the min and max.
        if(firstEle){
            firstEle=false;
            coord=mesh->GetVertex(v[0]);
            for (int j = 0; j < 3; j++) {
                coord_min[j]=coord[j];
                coord_max[j]=coord[j];
            }            
        }
        
        for(int j=0; j<nv; j++){
            coord=mesh->GetVertex(v[j]);
            
            for (int k = 0; k < 3; k++) {
                if(coord[k]<coord_min[k]){
                    coord_min[k]=coord[k];
                    // Keep track vertex and element indeces for min in z-axis
                    if(k==2){  
                        apexVet=v[j];
                        apexEleIndex=i;
                    }
                }
                if(coord[k]>coord_max[k]){
                    coord_max[k]=coord[k];
                }            
            }                                    
        }
        
    }

    cout << "Min: " << coord_min[0] << " " << coord_min[1] << " " << coord_min[2] << endl;
    cout << "Max: " << coord_max[0] << " " << coord_max[1] << " " << coord_max[2] << endl;
    coord = mesh->GetVertex(apexVet);
    cout << "Apex: " << coord[0] << " " << coord[1] << " " << coord[2] << endl;
    
    // Top 5% of the z axis.
    double zTop5=coord_max[2]-(coord_max[2]-coord_min[2])*0.05;
    cout << "Top 5% z coordinate: " << zTop5 << endl;    

    // Initialization the attributes to 0 and set attribute of apex
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        // initialize the attribute for boundary.  
        ele->SetAttribute(0);
        
        //Found apex elements and set attribute.
        for (int j = 0; j < nv; j++) {
            if (v[j] ==apexVet){
                ele->SetAttribute(apexAttr);
                cout << "Element index = " << i << endl;
            }
        }        
    }
    
    // Base    
    // The base must be planar. Its norm must be within 20 degrees of z axis.
    double cosTheta = cos(20*3.14159265/180); 
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        MFEM_ASSERT(nv == 3, "Wrong boundary size");
        
        double *coord0 = mesh->GetVertex(v[0]);
        if(coord0[2]>zTop5){
            double *coord1 = mesh->GetVertex(v[1]);
            double *coord2 = mesh->GetVertex(v[2]);
            if(isPlanar(coord0, coord1, coord2, cosTheta)){
                ele->SetAttribute(baseAttr);
            }
        }
    }
    
    //EPI
    vector<Element *> elements;
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i); 
        if(ele->GetAttribute()==0){
            elements.push_back(ele);
        }
    }
    
    Element *apexEle=mesh->GetBdrElement(apexEleIndex);
    findNeighbor(apexEle, elements, epiAttr);
    
    // LV & RV
    vector<Element *> vElements;    
    for(int i=0; i<elements.size(); i++){
        Element *ele =elements[i];
        if(ele->GetAttribute()==0){
            vElements.push_back(ele);
        }
    }
    // pick one element in the container and assume it is in LV.
    // TODO: we need additional information to identify LV and RV.
    int last=vElements.size()-1;
    Element *lastEle=vElements[last];
    lastEle->SetAttribute(lvAttr);
    // get rid of last element in the container
    vElements.pop_back();
    findNeighbor(lastEle, vElements, lvAttr);
    
    for(int i=0; i<vElements.size(); i++){
        Element *ele =vElements[i];
        if(ele->GetAttribute()==0){
            ele->SetAttribute(rvAttr);
        }
    }

    // Check if there are unassigned elements left.
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i); 
        MFEM_ASSERT(ele->GetAttribute()!=0, "Unassigned element.");
    }  
    
    mesh->SetAttributes();
            
}

void setBaseOLD(Mesh *mesh, int attr){
    int nv = mesh->GetNV();
    double *coord;
    double coord_min[3];
    double coord_max[3];
    int apexVet=0;
    coord = mesh->GetVertex(0);
    if(coord!=NULL){
        for (int j = 0; j < 3; j++) {
            coord_min[j]=coord[j];
            coord_max[j]=coord[j];
        }
    }
    for (int i = 0; i < nv; i++) {
        coord = mesh->GetVertex(i);
        for (int j = 0; j < 3; j++) {
            if(coord[j]<coord_min[j]){
                coord_min[j]=coord[j];
                if(j==2){
                    apexVet=i;
                }
            }
            if(coord[j]>coord_max[j]){
                coord_max[j]=coord[j];
            }            
        }

    }

    cout << "Min: " << coord_min[0] << " " << coord_min[1] << " " << coord_min[2] << endl;
    cout << "Max: " << coord_max[0] << " " << coord_max[1] << " " << coord_max[2] << endl;
    coord = mesh->GetVertex(apexVet);
    cout << "Apex: " << coord[0] << " " << coord[1] << " " << coord[2] << endl;
    
    // Top 5% of the z axis.
    double zTop5=coord_max[2]-(coord_max[2]-coord_min[2])*0.05;
    cout << "Top 5% z coordinate: " << zTop5 << endl;
    
    
    const int apexAttr=0;
    int ne = mesh->GetNE();    
    for (int i = 0; i < ne; i++) {
        Element *ele = mesh->GetElement(i);
        ele->SetAttribute(5);
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        for (int j = 0; j < nv; j++) {
            if (v[j] ==apexVet){
                ele->SetAttribute(apexAttr);
                cout << "Element index = " << i << endl;
            }
        }
    }    
    
    // for debug
/*    ofstream befh;
    befh.open("boundary3.txt");
*/    //
    
    double cosTheta = cos(20*3.14159265/180); // within 20 degrees of z axis.
    vector<vector<int> > baseBoundary;
    int nbe=mesh->GetNBE();
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        if(nv!=3){
            cout << "Boundary element should be 3 but it is " << nv <<endl;
            return;
        }
        double *coord0 = mesh->GetVertex(v[0]);
        if(coord0[2]>zTop5){
            double *coord1 = mesh->GetVertex(v[1]);
            double *coord2 = mesh->GetVertex(v[2]);
            if(isPlanar(coord0, coord1, coord2, cosTheta)){
                //baseBoundary.Append(ele);
                vector<int> vertecies;                
                for (int j = 0; j < nv; j++) {
                    vertecies.push_back(v[j]);
                }
                sort(vertecies.begin(), vertecies.end());
                baseBoundary.push_back(vertecies);
                // for debug
/*                befh << nv;
                for (int j = 0; j < nv; j++) {                    
                    befh << ' ' << v[j];                    
                }
                befh << '\n'; 
                
                befh << nv;
                for (int j = 0; j < nv; j++) {                    
                    befh << ' ' << vertecies[j];                    
                }
                befh << '\n';                 
                //
 */
            }
        }
        
    }
      
    for (int i = 0; i < ne; i++) {
        Element *ele = mesh->GetElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        double *coord0 = mesh->GetVertex(v[0]);        
        if(coord0[2]>zTop5){
            vector<int> tet;
            for (int j = 0; j < nv; j++) {
                tet.push_back(v[j]);
            }
            sort(tet.begin(), tet.end());

            for(int j=0; j < baseBoundary.size(); j++){
                vector<int> tri=baseBoundary[j];
                if(isTriInTet(tri, tet)){
                    ele->SetAttribute(attr);
                    break;
                }
            }
        }
        
    }    
    
    
    
}

void printSurfVTK(Mesh *mesh, std::ostream &out){
   out <<
       "# vtk DataFile Version 3.0\n"
       "Generated by MFEM\n"
       "ASCII\n"
       "DATASET UNSTRUCTURED_GRID\n";
   
   int NumOfVertices=mesh->GetNV(); 
   int spaceDim=3;
   
    out << "POINTS " << NumOfVertices << " double\n";
    for (int i = 0; i < NumOfVertices; i++)
    {
       const double* coord=mesh->GetVertex(i);
       for(int j=0; j<spaceDim; j++){
           out << coord[j] << " ";
       }
       out << '\n';
    }
    
    int NumOfElements=mesh->GetNBE();
      int size = 0;
      for (int i = 0; i < NumOfElements; i++)
      {
         const Element *ele = mesh->GetBdrElement(i); 
         size += ele->GetNVertices() + 1;
      }
      
      out << "CELLS " << NumOfElements << ' ' << size << '\n';
      for (int i = 0; i < NumOfElements; i++)
      {
         const Element *ele = mesh->GetBdrElement(i);
         const int *v = ele->GetVertices();
         const int nv = ele->GetNVertices();
         out << nv;
         for (int j = 0; j < nv; j++)
         {
            out << ' ' << v[j];
         }
         out << '\n';
      } 
      
   out << "CELL_TYPES " << NumOfElements << '\n';
   for (int i = 0; i < NumOfElements; i++)
   {
      const Element *ele = mesh->GetBdrElement(i);
      int vtk_cell_type = 5;
      {
         switch (ele->GetGeometryType())
         {
            case Geometry::TRIANGLE:     vtk_cell_type = 5;   break;
            case Geometry::SQUARE:       vtk_cell_type = 9;   break;
            case Geometry::TETRAHEDRON:  vtk_cell_type = 10;  break;
            case Geometry::CUBE:         vtk_cell_type = 12;  break;
         }
      }

      out << vtk_cell_type << '\n';
   }
   
   // write attributes
   out << "CELL_DATA " << NumOfElements << '\n'
       << "SCALARS material int\n"
       << "LOOKUP_TABLE default\n";
   for (int i = 0; i < NumOfElements; i++)
   {
      const Element *ele = mesh->GetBdrElement(i);
      out << ele->GetAttribute() << '\n';
   }
   out.flush();   
      
}
