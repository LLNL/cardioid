#include "genfiber.h"
#include "constants.h"
#include "utils.h"

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


