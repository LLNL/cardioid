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

const double PI=3.14159265;
const int dim3=3;

using namespace std;
using namespace mfem;

void setSurfaces(Mesh *mesh, double angle);
void printSurfVTK(Mesh *mesh, std::ostream &out);
void printFiberVTK(Mesh *mesh, vector<Vector>& fiber_vecs, std::ostream &out);
void getVert2Elements(Mesh *mesh, vector<vector<int> >& vert2Elements);
void laplace(Mesh *mesh, vector<vector<int> >& vert2Elements, vector<double> &pot, vector<Vector> &gradients, Array<int> &all_ess_bdr, Array<int> &nonzero_ess_bdr, Array<int> &zero_ess_bdr, string output, int order, bool static_cond);
void bislerp(DenseMatrix& Q, DenseMatrix& Qa, DenseMatrix& Qb, double t);
double a_s_f(double a_endo, double a_epi, double d);
double a_w_f(double a_endo, double a_epi, double d);
double b_s_f(double b_endo, double b_epi, double d);
double b_w_f(double b_endo, double b_epi, double d);
void axis(DenseMatrix& Q,Vector &psi, Vector &phi);
void orient(DenseMatrix& Qp, DenseMatrix& Q, double a, double b);
bool vecisnonzero(Vector& vec);
bool vecdot(Vector &q1, Vector &q2);
void vectorEigen(Vector& psi_ab, DenseMatrix& QPfib);

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


double a_s_f(double a_endo, double a_epi, double d){
    return (a_endo*(1-d)-a_endo*d);
}

double a_w_f(double a_endo, double a_epi, double d){
    return (a_endo*(1-d)+a_epi*d);
}

double b_s_f(double b_endo, double b_epi, double d){
    return (b_endo*(1-d)-b_endo*d);
}

double b_w_f(double b_endo, double b_epi, double d){
    return (b_endo*(1-d)+b_epi*d);
}

bool vecdot(Vector &q1, Vector &q2){
    Vector a=q1;
    Vector b=q2;
    
    double norm1=a.Norml2();
    double norm2=b.Norml2();
    
    if(norm1>0){
        a/=norm1;
    }
    
    if(norm2>0){
        b/=norm2;
    }
    
//    cout << "a=";
//    a.Print(cout);
//    cout << " b=";
//    b.Print(cout);
//    cout << " q1=";
//    q1.Print(cout);
//    cout << " q2=";
//    q2.Print(cout); 
//    cout << endl;
    
    double dot=a*b;
    if(dot<0){
        dot=-dot;
    }
    if(dot>0.9999){
        return true;
    }
    return false;
}

bool vecisnorm(Vector &q){
    double sum=0;
    for(int i=0; i<q.Size(); i++){
        sum=sum+q(i)*q(i);
    }
    if(sum>0.99 && sum <=1.01){
        return true;
    }
    return false;
}

bool vecisnonzero(Vector& vec){
    double sum=0.0;
    for(int i=0; i<vec.Size(); i++){
        sum=sum+vec(i)*vec(i);
    }
    if(sum>0){
        return true;
    }
    return false;
}


void cross(Vector &e_0, Vector &e_1, Vector &e_2){
    MFEM_ASSERT(e_1.Size()==3, "size of e_1 should be 3");
    MFEM_ASSERT(e_2.Size()==3, "size of e_2 should be 3");
    e_0(0)=e_1(1)*e_2(2)-e_1(2)*e_2(1);
    e_0(1)=e_1(2)*e_2(0)-e_1(0)*e_2(2);
    e_0(2)=e_1(0)*e_2(1)-e_1(1)*e_2(0);    
}

void axis(DenseMatrix& Q,Vector &psi, Vector &phi){
    double norm=psi.Norml2();
    Vector e_1=psi;
    e_1/=norm;

    norm=phi.Norml2();
    Vector phi_unit=phi;
    phi_unit/=norm;

    double phi_e1=e_1*phi_unit; 
    //cout << "dot " << phi_e1 << endl;
    Vector e_2(3);
    add(phi_unit, -phi_e1, e_1, e_2);
    norm=e_2.Norml2();    
    e_2/=norm;
    
    Vector e_0(3);
    cross(e_0, e_1, e_2);

//    if(!vecisnorm(e_0)){  
//        norm=e_0.Norml2();
//        e_0/=norm;
//    }
    
    Q.SetCol(0, e_0);
    Q.SetCol(1, e_1);
    Q.SetCol(2, e_2);    
          
    MFEM_ASSERT(vecisnorm(e_0), "axis: e_0 is not normalized");
//    if(!vecisnorm(e_0)){
//        cout << "Norm psi ";
//        for(int i=0; i <psi.Size(); i++){
//            cout << psi(i) << " ";
//        } 
//        cout << " phi ";
//        for(int i=0; i <phi.Size(); i++){
//            cout << phi(i) << " ";
//        } 
//        cout <<endl;        
//        for(int j=0; j<Q.size(); j++){
//            Vector e=Q[j];
//            cout << "Norm e" << j << " ";
//            for(int i=0; i <e.Size(); i++){
//                cout << e(i) << " ";
//            }
//            cout << endl;
//        }         
//    }
    MFEM_ASSERT(vecisnorm(e_1), "axis: e_1 is not normalized");
    MFEM_ASSERT(vecisnorm(e_2), "axis: e_2 is not normalized");
    
}

void orient(DenseMatrix& Qp, DenseMatrix& Q, double a, double b){
    
    double sina=sin(a*PI/180);
    double cosa=cos(a*PI/180);
    double sinb=sin(b*PI/180);
    double cosb=cos(b*PI/180);
    
    //| cosa   -sina   0  |
    //| sina    cosa   0  |
    //|  0        0    1  |
    DenseMatrix matrixA(dim3, dim3);
    matrixA=0.0;
    matrixA(0, 0)=cosa;
    matrixA(0, 1)=-sina;
    matrixA(1, 0)=sina;
    matrixA(1, 1)=cosa;
    matrixA(2, 2)=1;    
    //|  1     0      0   |
    //|  0    cosb   sinb |
    //|  0   -sinb   cosb |    
    DenseMatrix matrixB(dim3, dim3);
    matrixB=0.0;
    matrixB(0, 0)=1;
    matrixB(1, 1)=cosb;
    matrixB(1, 2)=sinb;
    matrixB(2, 1)=-sinb;
    matrixB(2, 2)=cosb;      
    
    DenseMatrix QA(dim3, dim3);
    
    Mult(Q, matrixA, QA);    
    Mult(QA, matrixB, Qp);
       
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

void rot2quat(Vector &q, DenseMatrix& Q){
    q.SetSize(4); // quaternion q=w+x*i+y*j+z*k
    double M11=Q(0,0);
    double M21=Q(1,0);
    double M31=Q(2,0);
    double M12=Q(0,1);
    double M22=Q(1,1);
    double M32=Q(2,1);    
    double M13=Q(0,2);
    double M23=Q(1,2);
    double M33=Q(2,2);  
        
    double w2=0.25*(1+M11+M22+M33);
    double error=0.000001;
    if(w2>error){
        double w=sqrt(w2);
        q(0)=w;
        double wq=4*w;
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
//    quatNormalized(q);
//    for(int i=0; i <q.Size(); i++){
//        cout << q(i) << endl;
//    }
//    MFEM_ASSERT(vecisnorm(q), "rot2quat: quaternion is not normalized");
}

void quat2rot(DenseMatrix& Q, Vector &q){
    MFEM_ASSERT(q.Size()==4, "quat2rot: Dimension of quaternion should be 4");
    //MFEM_ASSERT(vecisnorm(q), "quat2rot: quaternion is not normalized");
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
    
    Q(0,0)=1-2*y2-2*z2;
    Q(1,0)=2*xy-2*wz;
    Q(2,0)=2*xz+2*wy;
    
    Q(0,1)=2*xy+2*wz;
    Q(1,1)=1-2*x2-2*z2;
    Q(2,1)=2*yz-2*wx;
    
    Q(0,2)=2*xz-2*wy;
    Q(1,2)=2*yz+2*wx;
    Q(2,2)=1-2*x2-2*y2; 
    
}


void slerp(Vector &q, Vector &q1, Vector &q2, double t) {
    double dot = quatdot(q1, q2);
    Vector q3=q2;
    if (dot < 0) {
        dot = -dot;
        //MFEM_ASSERT(dot>=0 && dot<=1, "slerp: dot is not in the range from 0 to 1");
        q3.Neg();
    } 
    
//    if(dot>1){
//        cout << "slerp: Warning: dot is " << dot << ", use linear interpolation. "<< endl;
//    }
    
    if (dot < 0.9999) {
        double angle = acos(dot);
        double a=sin(angle * (1 - t))/sin(angle);
        double b=sin(angle * t) / sin(angle);
        add(a, q1, b, q3, q); //q = (q1 * sin(angle * (1 - t)) + q2 * sin(angle * t)) / sin(angle);
        //MFEM_ASSERT(vecisnorm(q), "slerp: quaternion is not normalized");
    } else { // if the angle is small dot>0.999, use linear interpolation								
        add((1-t), q1, t, q3, q); //q=q1*(1-t)+q2*t;
        
    }	
//    quatNormalized(q);
}

void bislerp(DenseMatrix& Q, DenseMatrix& Qa, DenseMatrix& Qb, double t){    
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
        //MFEM_ASSERT(vecisnorm(qai), "bislerp: quaternion qai is not normalized");
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

void vectorEigen(Vector& psi_ab, DenseMatrix& QPfib){
    
    Vector e=psi_ab;
    double norm=e.Norml2();
    e/=norm;
    
    DenseMatrix m(dim3,dim3);  
    MultVVt(psi_ab, m);
      
    double lamda[dim3];
    double vec[dim3*dim3];
    m.CalcEigenvalues(lamda, vec);
    double lamda_max=-1;
    int index_max=-1;
    for(int i=0; i<3; i++){
        double lamda_abs=abs(lamda[i]);
        if(lamda_abs>lamda_max) {
            lamda_max=lamda_abs;
            index_max=i;     
        } 
    }
    
    Vector fvec;
    vector<Vector> stvec;
    for(int i=0; i<3; i++){
        Vector vec_tmp(3);
        for(int j=0; j<3; j++){
            vec_tmp(j)=vec[i*3+j];
        }
        if(i==index_max){
            fvec=vec_tmp;
        }else{
            stvec.push_back(vec_tmp);
        }
    }
    QPfib.SetCol(0, fvec);
    for(int i=0; i<stvec.size(); i++){
        QPfib.SetCol(i+1, stvec[i]);
    }
    
}

void getVert2Elements(Mesh *mesh, vector<vector<int> >& vert2Elements) {

    int NumOfVertices = mesh->GetNV();
    for (int i = 0; i < NumOfVertices; i++) {
        vector<int> elements;
        vert2Elements.push_back(elements);
    }

    int NumOfElements = mesh->GetNE();
    for (int i = 0; i < NumOfElements; i++) {
        const Element *ele = mesh->GetElement(i);
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        for (int j = 0; j < nv; j++) {
            int vert = v[j];
            vert2Elements[vert].push_back(i);
        }
    }
    
    for (int i = 0; i < NumOfVertices; i++) {
        stringstream msg;
        msg << "getVert2Elements : vertex[" << i << "] size is zero"; 
        MFEM_ASSERT(vert2Elements[i].size()!=0, msg.str());
    }
    

}


void laplace(Mesh *mesh, vector<vector<int> >& vert2Elements, vector<double> &pot, vector<Vector> &gradients, Array<int> &all_ess_bdr, Array<int> &nonzero_ess_bdr, Array<int> &zero_ess_bdr, string output, int order, bool static_cond){
    
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
    // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to //ILU
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
    int nv = fes->GetNV();
    for (int i = 0; i < nv; i++) {
        vector<int> elements = vert2Elements[i];
        Vector grad(3);
        grad=0.0;
        stringstream msg;
        msg << "laplace : vertex[" << i << "] size is zero"; 
        MFEM_ASSERT(elements.size()!=0, msg.str());        
        for (int j = 0; j < elements.size(); j++) {
            ElementTransformation * tr = fes->GetElementTransformation(elements[j]);
            //Vector grad_ele;
            //x.GetGradient((*tr), grad_ele);
            //grad+=grad_ele; 
            x.GetGradient((*tr), grad);
            break;
        }
        //grad/=elements.size();
        gradients.push_back(grad);

        if (i < 5) {
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
    
    vector<Vector> xvectors;
    for(int i=0; i<x.Size(); i++){
        Vector vectmp(3);
        vectmp(0)=x(i);
        vectmp(1)=0.0;
        vectmp(2)=0.0;
        xvectors.push_back(vectmp);
    }
    
    fileName=output+"-x.vtk";
    ofstream x_ofs(fileName.c_str());
    printFiberVTK(mesh, xvectors, x_ofs);    

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

void setSurfaces(Mesh *mesh, double angle=20){
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
    double cosTheta = cos(angle*PI/180); 
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
    
    double cosTheta = cos(20*PI/180); // within 20 degrees of z axis.
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

void printFiberVTK(Mesh *mesh, vector<Vector>& fiber_vecs, std::ostream &out){
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
    
    int NumOfElements=mesh->GetNE();
      int size = 0;
      for (int i = 0; i < NumOfElements; i++)
      {
         const Element *ele = mesh->GetElement(i); 
         size += ele->GetNVertices() + 1;
      }
      
      out << "CELLS " << NumOfElements << ' ' << size << '\n';
      for (int i = 0; i < NumOfElements; i++)
      {
         const Element *ele = mesh->GetElement(i);
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
      const Element *ele = mesh->GetElement(i);
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
   
   // write point data
   MFEM_ASSERT(fiber_vecs.size()==NumOfVertices, "Size of fiber vector should equal to size of points");
   out << "POINT_DATA " << NumOfVertices<< "\n"
       << "VECTORS fiber double\n";
    for (int i = 0; i < NumOfVertices; i++) {
        fiber_vecs[i].Print(out, 10);
    }
   
   out << "\n";
               
   // write attributes
   out << "CELL_DATA " << NumOfElements << '\n'
       << "SCALARS material int\n"
       << "LOOKUP_TABLE default\n";
   for (int i = 0; i < NumOfElements; i++)
   {
      const Element *ele = mesh->GetElement(i);
      out << ele->GetAttribute() << '\n';
   }
   out.flush();   
      
}
