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

using namespace std;
using namespace mfem;

void setSurfaces(Mesh *mesh);
void printSurfVTK(Mesh *mesh, std::ostream &out);
void laplace(Mesh *mesh, vector<double> &pot, vector<Vector> &gradients, Array<int> &all_ess_bdr, Array<int> &nonzero_ess_bdr, Array<int> &zero_ess_bdr, string output, int order, bool static_cond);
void bislerp(vector<Vector>& Q, vector<Vector>& Qa, vector<Vector>& Qb, double t);
double a_s_f(double a_endo, double a_epi, double d);
double a_w_f(double a_endo, double a_epi, double d);
double b_s_f(double b_endo, double b_epi, double d);
double b_w_f(double b_endo, double b_epi, double d);
void axis(vector<Vector>& Q,Vector &psi, Vector &phi);
void orient(vector<Vector>& Qp, vector<Vector>& Q, double a, double b);

int main(int argc, char *argv[]) {
    // 1. Parse command-line options.
    const char *mesh_file = "./mechmesh.vtk";
    int order = 1;
    bool static_cond = false;
    bool visualization = 1;

    double a_endo=60;
    double a_epi=60;
    double b_endo=20;
    double b_epi=60;
    
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
    setSurfaces(mesh);

    // 3. Solve the laplacian for four different boundary conditions.
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
    laplace(mesh, psi_ab, psi_ab_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
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

    laplace(mesh, phi_epi, phi_epi_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
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

    laplace(mesh, phi_lv, phi_lv_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
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

    laplace(mesh, phi_rv, phi_rv_grads, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, output, order, static_cond);
    MFEM_ASSERT(phi_rv.size()==nv, "size of phi_rv does not match number of vertices.");
    MFEM_ASSERT(phi_rv_grads.size()==nv, "size of phi_rv_grads does not match number of vertices.");
            
    vector<Vector> fvectors;
    vector<Vector> svectors;
    vector<Vector> tvectors;
    // Line 7 start for-loop
    for(int i=0; i <nv; i++){        
        double frac=phi_rv[i]/(phi_lv[i]+phi_rv[i]);
        double frac_epi=phi_epi[i];
        double as=a_s_f(a_endo, a_epi, frac);
        double bs=b_s_f(b_endo, b_epi, frac);
        double aw=a_w_f(a_endo, a_epi, frac_epi);
        double bw=b_w_f(b_endo, b_epi, frac_epi);
        
        // Line 8
        Vector psi_ab_vec=psi_ab_grads[i];
        Vector phi_lv_vec=phi_lv_grads[i];
        Vector phi_lv_vec_neg(phi_lv_vec.Size());
        for(int i=0; i<phi_lv_vec.Size(); ++i){
            phi_lv_vec_neg(i)=phi_lv_vec(i);
        }
        vector<Vector> Qlv;
        axis(Qlv, psi_ab_vec, phi_lv_vec_neg);
        vector<Vector> QPlv;
        orient(QPlv, Qlv, as, bs);
        
        //Line 9
        Vector phi_rv_vec=phi_rv_grads[i];
        vector<Vector> Qrv;
        axis(Qrv, psi_ab_vec, phi_rv_vec);
        vector<Vector> QPrv;
        orient(QPrv, Qrv, as, bs); 
        
        //Line 10
        vector<Vector> QPendo;
        bislerp(QPendo, QPlv, QPrv, frac);
        
        //Line 11
        Vector phi_epi_vec=phi_epi_grads[i];
        vector<Vector> Qepi;
        axis(Qepi, psi_ab_vec, phi_epi_vec);
        vector<Vector> QPepi;
        orient(QPepi, Qepi, aw, bw);
        
        //Line 12
        vector<Vector> QPfib;
        bislerp(QPfib, QPendo, QPepi, frac_epi);
        MFEM_ASSERT(QPfib.size()==3, "Size of QPfilb should be 3");
        fvectors.push_back(QPfib[0]);
        svectors.push_back(QPfib[1]);
        tvectors.push_back(QPfib[2]);
    }
    
    ofstream f_ofs("fvectors.txt");
    ofstream s_ofs("svectors.txt");
    ofstream t_ofs("tvectors.txt");
    
    for(int i=0; i< fvectors.size(); ++i){
        fvectors[i].Print(f_ofs, 10);
        svectors[i].Print(s_ofs, 10);
        tvectors[i].Print(t_ofs, 10);
    }
    
        
    delete mesh;

    return 0;
}

double a_s_f(double a_endo, double a_epi, double d){
    return (a_endo*(1-d)-a_epi*d);
}

double a_w_f(double a_endo, double a_epi, double d){
    return (a_endo*(1-d)+a_epi*d);
}

double b_s_f(double b_endo, double b_epi, double d){
    return (b_endo*(1-d)-b_epi*d);
}

double b_w_f(double b_endo, double b_epi, double d){
    return (b_endo*(1-d)+b_epi*d);
}

void cross(Vector &e_0, Vector &e_1, Vector &e_2){
    MFEM_ASSERT(e_1.Size()==3, "size of e_1 should be 3");
    MFEM_ASSERT(e_2.Size()==3, "size of e_2 should be 3");
    e_0(0)=e_1(1)*e_2(2)-e_1(2)*e_2(1);
    e_0(1)=e_1(2)*e_2(0)-e_1(0)*e_2(2);
    e_0(2)=e_1(0)*e_2(1)-e_1(1)*e_2(0);    
}

void axis(vector<Vector>& Q,Vector &psi, Vector &phi){
    double norm=psi.Norml2();
    Vector e_1=psi;
    e_1/=norm;
    double psi_e1=e_1*psi;
    Vector e_2(3);
    add(phi, -psi_e1, e_1, e_2);
    norm=e_2.Norml2();
    e_2/=norm;
    Vector e_0(3);
    cross(e_0, e_1, e_2);
    Q.push_back(e_0);
    Q.push_back(e_1);
    Q.push_back(e_2);
}

void orient(vector<Vector>& Qp, vector<Vector>& Q, double a, double b){
    Vector e_0=Q[0];
    Vector e_1=Q[1];
    Vector e_2=Q[2];
    
    Vector ep_0(3);
    Vector ep_1(3);
    Vector ep_2(3);
    
    double sina=sin(a);
    double cosa=cos(a);
    double sinb=sin(b);
    double cosb=cos(b);
    
    add(cosa, e_0, sina, e_1, ep_0);
    
    Vector tmp_ep_1(3);
    add(-sina*cosb, e_0, cosa*cosb, e_1, tmp_ep_1);
    add(1, tmp_ep_1, -sinb, e_2, ep_1);
    
    Vector tmp_ep_2(3);
    add(-sina*sinb, e_0, cosa*sinb, e_1, tmp_ep_2);
    add(1, tmp_ep_2, cosb, e_2, ep_2);   
    
    Qp.push_back(ep_0);
    Qp.push_back(ep_1);
    Qp.push_back(ep_2);
    
}

void rot2quat(Vector &q,vector<Vector>& Q){
    MFEM_ASSERT(Q.size()==3, "Dimension of rotation matrix should be 3");
    q.SetSize(4); // quaternion q=w+x*i+y*j+z*k
    Vector e_0=Q[0];
    double M11=e_0(0);
    double M21=e_0(1);
    double M31=e_0(2);
    Vector e_1=Q[1];
    double M12=e_1(0);
    double M22=e_1(1);
    double M32=e_1(2);    
    Vector e_2=Q[2];
    double M13=e_2(0);
    double M23=e_2(1);
    double M33=e_2(2);  
    
    double w2=0.25*(1+M11+M22+M33);
    double error=0.001;
    if(w2>error){
        double w=sqrt(w2);
        q(0)=w;
        double wq=0.25*w;
        q(1)=(M23-M32)/wq;  //x
        q(2)=(M31-M13)/wq;  //y
        q(3)=(M12-M21)/wq;  //z
    }else{
        q(0)=0; //w
        double x2=-0.5*(M22+M33);
        if(x2>error){
            double x=sqrt(x2);
            q(1)=x;
            q(2)=M12/(2*x); //y
            q(3)=M13/(2*x); //z
        }else{
            q(1)=0; //x
            double y2=0.5*(1-M33);
            if(y2>error){
                double y=sqrt(y2);
                q(2)=y;
                q(3)=M23/(2*y); //z
            }else{
                q(2)=0; // y
                q(3)=1; // z
            }
            
        }
        
    }
}

void quat2rot(vector<Vector>& Q, Vector &q){
    MFEM_ASSERT(q.Size()==4, "quat2rot: Dimension of quaternion should be 4");
    double w=q(0);
    double x=q(1);
    double y=q(2);
    double z=q(3);
    
    double x2=x*x;
    double y2=y*y;
    double z2=z*z;
    double xy=x*y;
    double xz=x*z;
    double yz=y*z;
    double wx=w*x;
    double wy=w*y;
    double wz=w*z;
    
    Vector e_0(3);
    e_0(0)=1-2*y2-2*z2;
    e_0(1)=2*xy-2*wz;
    e_0(2)=2*xz+2*wy;
    
    Vector e_1(3);
    e_1(0)=2*xy+2*wz;
    e_1(1)=1-2*x2-2*z2;
    e_1(2)=2*yz-2*wx;
    
    Vector e_2(3);
    e_2(0)=2*xz-2*wy;
    e_2(1)=2*yz+2*wx;
    e_2(2)=1-2*x2-2*y2; 
    
    Q.push_back(e_0);
    Q.push_back(e_1);
    Q.push_back(e_2);
}

double quatdot(Vector &q1, Vector &q2){
    MFEM_ASSERT(q1.Size()==4, "quatdot Dimension of quaternion1 should be 4");
    MFEM_ASSERT(q2.Size()==4, "quatdot Dimension of quaternion2 should be 4");
    double sum=0.0;
    for(int i=0; i<q1.Size(); i++){
        sum=sum+q1(i)*q2(i);
    }    
    return sum;
}

void quatNormalized(Vector &q){
    double sum=0.0;
    for(int i=0; i<q.Size(); i++){
        sum=sum+q(i)*q(i);
    } 
    q/=(sqrt(sum));    
}


void slerp(Vector &q, Vector &q1, Vector &q2, double t) {
    double dot = quatdot(q1, q2);
    if (dot < 0) {
        dot = -dot;
        q2.Neg();
    } 
    
    if (dot < 0.95) {
        double angle = acos(dot);
        double a=sin(angle * (1 - t))/sin(angle);
        double b=sin(angle * t) / sin(angle);
        add(a, q1, b, q2, q); //q = (q1 * sin(angle * (1 - t)) + q2 * sin(angle * t)) / sin(angle);
    } else { // if the angle is small dot>0.95, use linear interpolation								
        add((1-t), q1, t, q2, q); //q=q1*(1-t)+q2*t;
        quatNormalized(q);
    }	
}

void bislerp(vector<Vector>& Q, vector<Vector>& Qa, vector<Vector>& Qb, double t){
    Vector qa;
    rot2quat(qa, Qa);
    Vector qb;
    rot2quat(qb, Qb);
    double a=qa(0);
    double b=qa(1);
    double c=qa(2);
    double d=qa(3);
   
    vector<Vector> qavec;
    Vector qa1(4);    
    qa1(0)=a;
    qa1(1)=b;
    qa1(2)=c;
    qa1(3)=d;
    qavec.push_back(qa1);
    
    Vector qa2(4);    
    qa2(0)=-b;
    qa2(1)=a;
    qa2(2)=-d;
    qa2(3)=c;
    qavec.push_back(qa2); 
 
    Vector qa3(4);    
    qa3(0)=-c;
    qa3(1)=d;
    qa3(2)=a;
    qa3(3)=-b;
    qavec.push_back(qa3); 
    
    Vector qa4(4);    
    qa4(0)=-d;
    qa4(1)=-c;
    qa4(2)=b;
    qa4(3)=a;
    qavec.push_back(qa4);   
    
    Vector qm(4);
    double maxdot=-1;
    for(int i=0; i<qavec.size(); i++){
        Vector qai=qavec[i];
        double dot=abs(quatdot(qai, qb));
        if(maxdot<dot){
            maxdot=dot;
            qm=qai;
        }     
    }
    
    Vector qab(4);
    slerp(qab, qm, qb, t);
    quat2rot(Q, qab);
    
}


void laplace(Mesh *mesh, vector<double> &pot, vector<Vector> &gradients, Array<int> &all_ess_bdr, Array<int> &nonzero_ess_bdr, Array<int> &zero_ess_bdr, string output, int order, bool static_cond){
    
    int dim = mesh->Dimension();
    cout << "Dimension =" << dim << endl;    
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
    ConstantCoefficient nonzero_bdr(1.0);
    x.ProjectBdrCoefficient(nonzero_bdr, nonzero_ess_bdr);

    // Project the constant 0 value to boundary attribute 1
    ConstantCoefficient zero_bdr(0.0);
    x.ProjectBdrCoefficient(zero_bdr, zero_ess_bdr);
        
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
      
    //double *x_data=x.GetData();
    for(int i=0; i<x.Size(); i++){         
        //pot.push_back(x_data[i]);
        pot.push_back(x(i));
        if(i<5) cout << "pot " << x(i) << endl;
    }
    
    const FiniteElementSpace *fes=x.FESpace();
    int nv=fes->GetNV();
    for(int i=0; i<nv; i++){   
        ElementTransformation *tr=fes->GetElementTransformation(i);
        Vector grad;
        x.GetGradient((*tr), grad);
        gradients.push_back(grad);
        if(i<5) {
            cout << "grad ";
            for (int j = 0; j < grad.Size(); j++) {
                cout << grad[j] << " ";
            }
            cout << endl;
        }
    }        
      

    // 12. Save the refined mesh and the solution. This output can be viewed later
    //     using GLVis: "glvis -m refined.mesh -g sol.gf".
    string fileName=output+".mesh";
    ofstream mesh_ofs(fileName.c_str());
    mesh_ofs.precision(8);
    mesh->Print(mesh_ofs);

    fileName=output+".gf";
    ofstream sol_ofs(fileName.c_str());
    sol_ofs.precision(8);
    x.Save(sol_ofs);

    // 14. Free the used memory.
    delete a;
    delete b;
    delete fespace;
    if (order > 0) {
        delete fec;
    }
    
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
    const int lvAttr=5; // TODO: Need something to indicate the LV and RV
    const int rvAttr=4;
       
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
