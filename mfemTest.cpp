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
    MFEM_ASSERT(vecisnorm(q), "rot2quat: quaternion is not normalized");
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

//const double PI=3.14159265;

void testRotQat() {
    DenseMatrix m(3,3);
    m=0.0;
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
    m(0,2)=1;
    m(1,0)=1;
    m(2,1)=1;
    
    cout << "Matrix 1" <<endl;
    m.PrintMatlab(cout);     
    
    Vector q;
    rot2quat(q, m);
    cout << "Quaternion " <<endl;
    q.Print(cout);
    
    DenseMatrix m2(3,3);
    quat2rot(m2, q);
    cout << "Matrix 2 should be the same as Matrix 1 within error." <<endl;
    m2.PrintMatlab(cout);     
}


void testMatrix() {
    DenseMatrix m1(3,3);    
    DenseMatrix m2(3,3);
    DenseMatrix m3(3,3);
    m1=0.0;
    Vector v0(3);
    v0=0.1;
    Vector v1(3);
    v1=0.2;
    Vector v2(3);
    v2=0.3;    
    m1.SetCol(0, v0);
    m1.SetCol(1, v1);
    m1.SetCol(2, v2);
    m2=0.0;
//    m1(0,0)=0.0;
//    m1(1,1)=0.0;
//    m1(2,2)=0.0;    
    m2(0,0)=1.0;
    m2(1,1)=1.0;
    m2(2,2)=1.0;
    
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
    
    Vector v32n=v32;
    v32n.Neg();
    v32n.Print(cout);
    
    
    //m3.Print(cout);
    cout << "Matrix 1" <<endl;
    m1.PrintMatlab(cout);     
    
    cout << "Matrix 2" <<endl;
    m2.PrintMatlab(cout);    
    
    cout << "Matrix 3" <<endl;
    m3.PrintMatlab(cout);
    
    
}

void testMatrixEigen(){
    Vector psi_ab(3);   
    //0.000192423 -0.000580408 -0.000580227
    // 
//    psi_ab(0)=-0.00105952 ;
//    psi_ab(1)=-0.00102312;
//    psi_ab(2)=-0.000290406;
    
    psi_ab(0)=0.34;
    psi_ab(1)=-0.22;
    psi_ab(2)=0.33;   
    
    Vector e=psi_ab;
    double norm=e.Norml2();
    e/=norm;
    e.Print(cout);
    
//    psi_ab(0)=1;
//    psi_ab(1)=0;
//    psi_ab(2)=0;
    
    DenseMatrix m(3,3);  
    MultVVt(psi_ab, m);
    
    cout << "Matrix " <<endl;
    m.PrintMatlab(cout);  
    double lamda[3];
    double vec[9];
    m.CalcEigenvalues(lamda, vec);
    for(int i=0; i<3; i++){
        cout << "lamda " << i << " = " << lamda[i] <<endl;       
    }
    for(int i=0; i<9; i++){
        if(i%3==0){
            cout <<endl;
        }
        cout << vec[i] << i << "    " ;           
    }    
    cout <<endl;
}

void cross(Vector &e_0, Vector &e_1, Vector &e_2){
    MFEM_ASSERT(e_1.Size()==3, "size of e_1 should be 3");
    MFEM_ASSERT(e_2.Size()==3, "size of e_2 should be 3");
    e_0(0)=e_1(1)*e_2(2)-e_1(2)*e_2(1);
    e_0(1)=e_1(2)*e_2(0)-e_1(0)*e_2(2);
    e_0(2)=e_1(0)*e_2(1)-e_1(1)*e_2(0);    
}

void testCross(){
    Vector e_0(3);
    Vector e_1(3);
    Vector e_2(3);
    e_1(0)=2;
    e_1(1)=3;
    e_1(2)=4;
    e_2(0)=5;
    e_2(1)=6;
    e_2(2)=7;

    cross(e_0, e_1, e_2);
    // should print out -3 6 -3
    e_0.Print(cout);
}

int main(int argc, char *argv[]) {
    //testMatrix();
    testMatrixEigen();
    //testRotQat();
      
    //testCross();
    return 0;
}