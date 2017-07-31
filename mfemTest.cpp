#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include "solver.h"
#include "genfiber.h"

using namespace std;
using namespace mfem;

bool vecisnorm(Vector &q) {
    double sum = 0;
    for (int i = 0; i < q.Size(); i++) {
        sum = sum + q(i) * q(i);
    }
    if (sum > 0.99 && sum <= 1.01) {
        return true;
    }
    return false;
}

void rot2quat(Vector &q, DenseMatrix& Q) {
    q.SetSize(4); // quaternion q=w+x*i+y*j+z*k
    double M11 = Q(0, 0);
    double M21 = Q(1, 0);
    double M31 = Q(2, 0);
    double M12 = Q(0, 1);
    double M22 = Q(1, 1);
    double M32 = Q(2, 1);
    double M13 = Q(0, 2);
    double M23 = Q(1, 2);
    double M33 = Q(2, 2);

    double w2 = 0.25 * (1 + M11 + M22 + M33);
    double error = 0.000001;
    if (w2 > error) {
        double w = sqrt(w2);
        q(0) = w;
        double wq = 4 * w;
        q(1) = (M23 - M32) / wq; //x
        q(2) = (M31 - M13) / wq; //y
        q(3) = (M12 - M21) / wq; //z
    } else {
        q(0) = 0; //w
        double x2 = -0.5 * (M22 + M33);
        if (x2 > error) {
            double x = sqrt(x2);
            q(1) = x;
            q(2) = M12 / (2 * x); //y
            q(3) = M13 / (2 * x); //z
        } else {
            q(1) = 0; //x
            double y2 = 0.5 * (1 - M33);
            if (y2 > error) {
                double y = sqrt(y2);
                q(2) = y;
                q(3) = M23 / (2 * y); //z
            } else {
                q(2) = 0; // y
                q(3) = 1; // z
            }

        }

    }
    //    quatNormalized(q);
    //    for(int i=0; i <q.Size(); i++){
    //        cout << q(i) << endl;
    //    }
    MFEM_ASSERT(vecisnorm(q), "rot2quat: quaternion is not normalized");
}

void quat2rot(DenseMatrix& Q, Vector &q) {
    MFEM_ASSERT(q.Size() == 4, "quat2rot: Dimension of quaternion should be 4");
    //MFEM_ASSERT(vecisnorm(q), "quat2rot: quaternion is not normalized");
    double w = q(0);
    double x = q(1);
    double y = q(2);
    double z = q(3);

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double xy = x*y;
    double xz = x*z;
    double yz = y*z;
    double wx = w*x;
    double wy = w*y;
    double wz = w*z;

    Q(0, 0) = 1 - 2 * y2 - 2 * z2;
    Q(1, 0) = 2 * xy - 2 * wz;
    Q(2, 0) = 2 * xz + 2 * wy;

    Q(0, 1) = 2 * xy + 2 * wz;
    Q(1, 1) = 1 - 2 * x2 - 2 * z2;
    Q(2, 1) = 2 * yz - 2 * wx;

    Q(0, 2) = 2 * xz - 2 * wy;
    Q(1, 2) = 2 * yz + 2 * wx;
    Q(2, 2) = 1 - 2 * x2 - 2 * y2;

}

//const double PI=3.14159265;

void testRotQat() {
    DenseMatrix m(3, 3);
    m = 0.0;
    //    double angle=30*PI/180.0;    
    //    m(0,0)=1.0;
    //    m(1,1)=cos(angle);
    //    m(1,2)=-sin(angle);
    //    m(2,2)=cos(angle);
    //    m(2,1)=sin(angle);
    ////    m(0,0)=0.36;
    ////    m(0,1)=0.48;
    ////    m(0,2)=-0.8;
    ////    m(1,0)=-0.8;
    ////    m(1,1)=0.60;
    ////    m(1,2)=0;    
    ////    m(2,0)=0.48;
    ////    m(2,1)=0.64;
    ////    m(2,2)=0.60;
    m(0, 2) = 1;
    m(1, 0) = 1;
    m(2, 1) = 1;

    cout << "Matrix 1" << endl;
    m.PrintMatlab(cout);

    Vector q;
    rot2quat(q, m);
    cout << "Quaternion " << endl;
    q.Print(cout);

    DenseMatrix m2(3, 3);
    quat2rot(m2, q);
    cout << "Matrix 2 should be the same as Matrix 1 within error." << endl;
    m2.PrintMatlab(cout);
}

void testMatrix() {
    DenseMatrix m1(3, 3);
    DenseMatrix m2(3, 3);
    DenseMatrix m3(3, 3);
    m1 = 0.0;
    Vector v0(3);
    v0 = 0.1;
    Vector v1(3);
    v1 = 0.2;
    Vector v2(3);
    v2 = 0.3;
    m1.SetCol(0, v0);
    m1.SetCol(1, v1);
    m1.SetCol(2, v2);
    m2 = 0.0;
    //    m1(0,0)=0.0;
    //    m1(1,1)=0.0;
    //    m1(2,2)=0.0;    
    m2(0, 0) = 1.0;
    m2(1, 1) = 1.0;
    m2(2, 2) = 1.0;

    Mult(m1, m2, m3);

    Vector v30;
    m3.GetColumn(0, v30);
    v30.Print(cout);
    Vector v31;
    m3.GetColumn(1, v31);
    v31.Print(cout);
    Vector v32;
    m3.GetColumn(2, v32);
    v32.Print(cout);

    Vector v32n = v32;
    v32n.Neg();
    v32n.Print(cout);


    //m3.Print(cout);
    cout << "Matrix 1" << endl;
    m1.PrintMatlab(cout);

    cout << "Matrix 2" << endl;
    m2.PrintMatlab(cout);

    cout << "Matrix 3" << endl;
    m3.PrintMatlab(cout);


}

void testMatrixEigen() {
    Vector psi_ab(3);
    //0.000192423 -0.000580408 -0.000580227
    // 
    //    psi_ab(0)=-0.00105952 ;
    //    psi_ab(1)=-0.00102312;
    //    psi_ab(2)=-0.000290406;

    psi_ab(0) = 0.34;
    psi_ab(1) = -0.22;
    psi_ab(2) = 0.33;

    Vector e = psi_ab;
    double norm = e.Norml2();
    e /= norm;
    e.Print(cout);

    //    psi_ab(0)=1;
    //    psi_ab(1)=0;
    //    psi_ab(2)=0;

    DenseMatrix m(3, 3);
    MultVVt(psi_ab, m);

    cout << "Matrix " << endl;
    m.PrintMatlab(cout);
    double lamda[3];
    double vec[9];
    m.CalcEigenvalues(lamda, vec);
    for (int i = 0; i < 3; i++) {
        cout << "lamda " << i << " = " << lamda[i] << endl;
    }
    for (int i = 0; i < 9; i++) {
        if (i % 3 == 0) {
            cout << endl;
        }
        cout << vec[i] << i << "    ";
    }
    cout << endl;
}

void cross(Vector &e_0, Vector &e_1, Vector &e_2) {
    MFEM_ASSERT(e_1.Size() == 3, "size of e_1 should be 3");
    MFEM_ASSERT(e_2.Size() == 3, "size of e_2 should be 3");
    e_0(0) = e_1(1) * e_2(2) - e_1(2) * e_2(1);
    e_0(1) = e_1(2) * e_2(0) - e_1(0) * e_2(2);
    e_0(2) = e_1(0) * e_2(1) - e_1(1) * e_2(0);
}

void testCross() {
    Vector e_0(3);
    Vector e_1(3);
    Vector e_2(3);
    e_1(0) = 2;
    e_1(1) = 3;
    e_1(2) = 4;
    e_2(0) = 5;
    e_2(1) = 6;
    e_2(2) = 7;

    cross(e_0, e_1, e_2);
    // should print out -3 6 -3
    e_0.Print(cout);
}

void testMesh(){
    const char *mesh_file = "examples/heart/heart.vtk";
    int order = 1;
    bool static_cond = false;
    bool visualization = 1;

    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    
    int vmax0=0;
    int vmax1=0;
    double* cMax0=0;
    double* cMax1=0;
    
    double maxEdgeLen=0;
    int NumOfElements=mesh->GetNE();
    for (int i = 0; i < NumOfElements; i++){
         const Element *ele = mesh->GetElement(i);
         const int *v = ele->GetVertices();
         const int ne = ele->GetNEdges();
         for (int j = 0; j < ne; j++){
            const int *e = ele->GetEdgeVertices(j);
            const double* coord0=mesh->GetVertex(v[e[0]]);
            const double* coord1=mesh->GetVertex(v[e[1]]);
            double dist2=0;
            for(int j=0; j<3; j++){
                dist2+=(coord0[j]-coord1[j])*(coord0[j]-coord1[j]);
            }            
//            std::cout << "e[0]=" << e[0] << " e[1]=" << e[1] 
//                    << " v[e[0]]=" << v[e[0]] << " v[e[1]]=" << v[e[1]] 
//                    << " dist2=" << dist2 <<std::endl; 
            if(maxEdgeLen<dist2){
                maxEdgeLen=dist2;
                vmax0=v[e[0]];
                vmax1=v[e[1]];
                cMax0=const_cast<double*>(coord0);
                cMax1=const_cast<double*>(coord1);
            }
         }
        //break;
    }
    std::cout << "maxEdgeLen=" << maxEdgeLen << std::endl;
    std::cout << "vmax0=" << vmax0  << " vmax1=" << vmax1 << std::endl;
    for(int j=0; j<3; j++){
        std::cout << cMax0[j]  << " " << cMax1[j] << std::endl;
    }
}

void testMFEM() {
    const char *mesh_file = "../fiber-test/human.vtk";
    int order = 1;
    bool static_cond = false;
    bool visualization = 1;

    Mesh *mesh = new Mesh(mesh_file, 1, 1);

    vector<Vector> boundingbox;
    // Set the surfaces for the mesh: 0-Apex, 1-Base, 2-EPI, 3-LV, 4-RV.
    setSurfaces(mesh, boundingbox, 20); // use 30 degrees for determining the base surface.

    int bdr_attr_size = mesh->bdr_attributes.Max();
    Array<int> all_ess_bdr(bdr_attr_size);
    Array<int> nonzero_ess_bdr(bdr_attr_size);
    Array<int> zero_ess_bdr(bdr_attr_size);

    // 3a. Base → 1, Apex→ 0, Epi, LV, RV → no flux
    // Mark ALL boundaries as essential. This does not set what the actual Dirichlet
    // values are
    all_ess_bdr = 1;
    all_ess_bdr[2] = 0;
    all_ess_bdr[3] = 0;
    all_ess_bdr[4] = 0;

    nonzero_ess_bdr = 0;
    nonzero_ess_bdr[1] = 1;

    zero_ess_bdr = 0;
    zero_ess_bdr[0] = 1;

    string output = "psi_ab";
    vector<double> psi_ab;
    vector<Vector> psi_ab_grads;
    GridFunction x_psi_ab = laplace(mesh, all_ess_bdr, nonzero_ess_bdr, zero_ess_bdr, order, static_cond);
    const FiniteElementSpace *fes = x_psi_ab.FESpace();

    ElementTransformation * tr = fes->GetElementTransformation(0);
    const FiniteElement * el= fes->GetFE(0);
    
    const IntegrationRule &ir = el->GetNodes(); // Get the parametric integration rule
    
    cout << "ir->GetNPoints() "<<ir.GetNPoints() <<endl;
    for (int k=0; k < ir.GetNPoints(); k++) {
        //ir.Print(cout);
        const IntegrationPoint &ip = ir.IntPoint(k);
        cout << "ip= " << k << " ip.x=" <<ip.x<< " ip.y=" <<ip.y<< " ip.z=" <<ip.z << " weight:" << ip.weight << endl; 
        tr->SetIntPoint(&ip); 
        cout << "tr  attr=" << tr->Attribute << " elementno" << tr->ElementNo << "weight" << tr->Weight() << endl;
        const IntegrationPoint &ipt=tr->GetIntPoint();
        cout << "ipt= " << k << " ipt.x=" <<ipt.x<< " ipt.y=" <<ipt.y<< " ipt.z=" <<ipt.z << " weight:" << ipt.weight << endl;
    }

}

void testMatrix2() {
    DenseMatrix m1(3, 3);
    m1 = 0.0;
    Vector v0(3);
    v0 = 1.0;
    Vector v1(3);
    v1 = 2.0;
    Vector v2(3);
    v2 = 3.0;
    m1.SetCol(0, v0);
    m1.SetCol(1, v1);
    m1.SetCol(2, v2);
    cout << "Matrix 1" << endl;
    m1.PrintMatlab(cout);

    double *d=m1.Data();
    for(int i=0; i<9; i++){
        cout << d[i] << " ";
    }
    cout <<endl;
    
    m1.SetRow(0, v0);
    m1.SetRow(1, v1);
    m1.SetRow(2, v2);
 
    cout << "Matrix 1 new" << endl;
    m1.PrintMatlab(cout);

    d=m1.Data();
    for(int i=0; i<9; i++){
        cout << d[i] << " ";
    }
    cout <<endl;    
}

double det4X4(DenseMatrix& matrix){
    double* m=matrix.Data();
    for(int i=0; i<16; i++){
        cout << m[i] << " ";
    }
    cout <<endl;     
     return
         m[12] * m[ 9] * m[ 6] * m[ 3] - m[ 8] * m[13] * m[ 6] * m[ 3] -
         m[12] * m[ 5] * m[10] * m[ 3] + m[ 4] * m[13] * m[10] * m[ 3] +
         m[ 8] * m[ 5] * m[14] * m[ 3] - m[ 4] * m[ 9] * m[14] * m[ 3] -
         m[12] * m[ 9] * m[ 2] * m[ 7] + m[ 8] * m[13] * m[ 2] * m[ 7] +
         m[12] * m[ 1] * m[10] * m[ 7] - m[ 0] * m[13] * m[10] * m[ 7] -
         m[ 8] * m[ 1] * m[14] * m[ 7] + m[ 0] * m[ 9] * m[14] * m[ 7] +
         m[12] * m[ 5] * m[ 2] * m[11] - m[ 4] * m[13] * m[ 2] * m[11] -
         m[12] * m[ 1] * m[ 6] * m[11] + m[ 0] * m[13] * m[ 6] * m[11] +
         m[ 4] * m[ 1] * m[14] * m[11] - m[ 0] * m[ 5] * m[14] * m[11] -
         m[ 8] * m[ 5] * m[ 2] * m[15] + m[ 4] * m[ 9] * m[ 2] * m[15] +
         m[ 8] * m[ 1] * m[ 6] * m[15] - m[ 0] * m[ 9] * m[ 6] * m[15] -
         m[ 4] * m[ 1] * m[10] * m[15] + m[ 0] * m[ 5] * m[10] * m[15]; 
}

void testTet(){
    Vector v1(4);
    v1(0)=1.0;
    v1(1)=2.0;
    v1(2)=3.0;
    v1(3)=1.0;
    
    Vector v2(4);
    v2(0)=2.0;
    v2(1)=2.0;
    v2(2)=3.0;
    v2(3)=1.0;    

    Vector v3(4);
    v3(0)=1.0;
    v3(1)=3.0;
    v3(2)=3.0;
    v3(3)=1.0;

    Vector v4(4);
    v4(0)=1.0;
    v4(1)=2.0;
    v4(2)=9.0;
    v4(3)=1.0;
    
    Vector q(4);
    q(0)=1.25;
    q(1)=2.25;
    q(2)=4.50;
    q(3)=1.0;  
    
    DenseMatrix D0(4,4);
    D0.SetRow(0, v1);
    D0.SetRow(1, v2);
    D0.SetRow(2, v3);
    D0.SetRow(3, v4);   
    D0.PrintMatlab(cout);
    double d0=det4X4(D0);
    
    DenseMatrix D1(4,4);
    D1.SetRow(0, q);
    D1.SetRow(1, v2);
    D1.SetRow(2, v3);
    D1.SetRow(3, v4); 
    D1.PrintMatlab(cout);    
    double d1=det4X4(D1);
    
    
    DenseMatrix D2(4,4);
    D2.SetRow(0, v1);
    D2.SetRow(1, q);
    D2.SetRow(2, v3);
    D2.SetRow(3, v4);  
    D2.PrintMatlab(cout);
    double d2=det4X4(D2);   

    DenseMatrix D3(4,4);
    D3.SetRow(0, v1);
    D3.SetRow(1, v2);
    D3.SetRow(2, q);
    D3.SetRow(3, v4); 
    D3.PrintMatlab(cout);
    double d3=det4X4(D3);    
  
    DenseMatrix D4(4,4);
    D4.SetRow(0, v1);
    D4.SetRow(1, v2);
    D4.SetRow(2, v3);
    D4.SetRow(3, q); 
    D4.PrintMatlab(cout);
    double d4=det4X4(D4);  
    cout << "d1= " << d1 << " d2= " << d2 << " d3= " << d3
            << " d4= " << d4 << " d0= " << d0 <<endl;    
    cout << "d1= " << d1/d0 << " d2= " << d2/d0 << " d3= " << d3/d0 
            << " d4= " << d4/d0 << " tot= " << (d1+d2+d3+d4)/d0 <<endl;
    
}
int main(int argc, char *argv[]) {
    
    //testTet();
    //testMatrix2();
    //testMFEM();
    testMesh();
    //testMatrix();
    //testMatrixEigen();
    //testRotQat();

    //testCross();
    return 0;
}