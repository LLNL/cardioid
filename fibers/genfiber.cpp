#include "genfiber.h"
#include "constants.h"
#include "utils.h"
#include "io.h"
#include "option.h"

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
          
   // MFEM_ASSERT(vecisnorm(e_0), "axis: e_0 is not normalized");
    if(!vecisnorm(e_0)){
        cout << "Norm psi ";
        for(int i=0; i <psi.Size(); i++){
            cout << psi(i) << " ";
        } 
        cout << " phi ";
        for(int i=0; i <phi.Size(); i++){
            cout << phi(i) << " ";
        } 
        cout <<endl; 
        double *data=Q.Data();
        for(int j=0; j<9; j++){
            cout << " d " << j << " = " << data[j] << " ";
            
        }  
        cout <<endl; 
    }
    MFEM_ASSERT(vecisnorm(e_0), "axis: e_0 is not normalized");
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
    DenseMatrix matrixA(dim, dim);
    matrixA=0.0;
    matrixA(0, 0)=cosa;
    matrixA(0, 1)=-sina;
    matrixA(1, 0)=sina;
    matrixA(1, 1)=cosa;
    matrixA(2, 2)=1;    
    //|  1     0      0   |
    //|  0    cosb   sinb |
    //|  0   -sinb   cosb |    
    DenseMatrix matrixB(dim, dim);
    matrixB=0.0;
    matrixB(0, 0)=1;
    matrixB(1, 1)=cosb;
    matrixB(1, 2)=sinb;
    matrixB(2, 1)=-sinb;
    matrixB(2, 2)=cosb;      
    
    DenseMatrix QA(dim, dim);
    
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
    if(!vecisnorm(q)){
        //cout << "Warning rot2quat: quaternion is not normalized" << endl;
        // This matrix to quaternion method is numerical instable when w is small.
        quatNormalized(q);
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
    for(unsigned i=0; i<qavec.size(); i++){
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
    
    DenseMatrix m(dim,dim);  
    MultVVt(psi_ab, m);
      
    double lamda[dim];
    double vec[dim*dim];
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
    for(unsigned i=0; i<stvec.size(); i++){
        QPfib.SetCol(i+1, stvec[i]);
    }
    
}

void biSlerpCombo(DenseMatrix& QPfib,
        double psi_ab, Vector& psi_ab_vec,
        double phi_epi, Vector& phi_epi_vec,
        double phi_lv, Vector& phi_lv_vec,
        double phi_rv, Vector& phi_rv_vec,
        Option& options) {

   // Initialize big values for QPfib so we will know it is a wrong one.
   Vector nonVal(3);
   nonVal = 999;
   for (int i = 0; i < dim; i++) {
       QPfib.SetCol(i, nonVal);
   }   
   
    double a_endo=options.a_endo;
    double a_epi=options.a_epi;
    double b_endo=options.b_endo;
    double b_epi=options.b_epi;    

    double phi_v = phi_lv + phi_rv;
    double frac = 0.5;
    if (phi_v != 0) {
        frac = phi_rv / phi_v;
    } else {
        if(options.verbose){
            cout << "\tWarning: phi_v ==0";
            cout << " phi_lv[i]=" << phi_lv << " phi_rv[i]=" << phi_rv << " phi_epi[i]=" << phi_epi << " psi_ab[i]=" << psi_ab << endl;       
        }
    }
    double frac_epi = phi_epi;
    //stringstream ss;
    //ss << "i=" << i << " phi_rv[i]=" << phi_rv[i] << " phi_lv[i]=" << phi_lv[i] << " frac=" << frac;
    //MFEM_ASSERT(frac>=0 && frac<=1, ss.str());
    //MFEM_ASSERT(frac_epi>=0 && frac_epi<=1, "frac_epi is not in range 0 to 1");
    double as = a_s_f(a_endo, a_epi, frac);
    double bs = b_s_f(b_endo, b_epi, frac);
    double aw = a_w_f(a_endo, a_epi, frac_epi);
    double bw = b_w_f(b_endo, b_epi, frac_epi);

    bool phi_lv_isnonzero = vecisnonzero(phi_lv_vec);
    bool phi_rv_isnonzero = vecisnonzero(phi_rv_vec);
    bool phi_epi_isnonzero = vecisnonzero(phi_epi_vec);

    DenseMatrix QPendo(dim, dim);

    if (!vecisnonzero(psi_ab_vec)) {
        if(options.verbose){
            cout << "\tWarning psi_ab gradient is zero" << endl;
        }
        return;
    }


    DenseMatrix QPlv(dim, dim);
    if (phi_lv_isnonzero) {
        // Line 8
        Vector phi_lv_vec_neg = phi_lv_vec;
        phi_lv_vec_neg.Neg();
        DenseMatrix Qlv(dim, dim);
        if (vecdot(psi_ab_vec, phi_lv_vec_neg)) {
            if(options.verbose){
                cout << "\tpsi_ab_vec equal to phi_lv_vec_neg" << endl;
            }
            phi_lv_isnonzero = false;
        } else {
            axis(Qlv, psi_ab_vec, phi_lv_vec_neg);
            orient(QPlv, Qlv, as, bs);
        }
        // End of Line 8
    }

    DenseMatrix QPrv(dim, dim);
    if (phi_rv_isnonzero) {
        //Line 9
        DenseMatrix Qrv(dim, dim);
        if (vecdot(psi_ab_vec, phi_rv_vec)) {
            if(options.verbose){
                cout << "\tpsi_ab_vec equal to phi_rv_vec" << endl;
            }
            phi_rv_isnonzero = false;
        } else {
            axis(Qrv, psi_ab_vec, phi_rv_vec);
            orient(QPrv, Qrv, as, bs);
        }
    }

    DenseMatrix QPepi(dim, dim);
    if (phi_epi_isnonzero) {
        //Line 11
        DenseMatrix Qepi(dim, dim);
        if (vecdot(psi_ab_vec, phi_epi_vec)) {
            if(options.verbose){
                cout << "\tpsi_ab_vec equal to phi_epi_vec" << endl;
            }
            phi_epi_isnonzero = false;
        } else {
            axis(Qepi, psi_ab_vec, phi_epi_vec);
            orient(QPepi, Qepi, aw, bw);
        }
    }

    if (phi_lv_isnonzero && phi_rv_isnonzero && phi_epi_isnonzero) {
       // if all three phi gradients are non-zero, use the original algorithm in paper. 
       //Line 10
       bislerp(QPendo, QPlv, QPrv, frac);
       //QPendo=QPlv;
       //Line 12 
       bislerp(QPfib, QPendo, QPepi, frac_epi);
       //QPfib=QPendo;
    } else if (!phi_lv_isnonzero && phi_rv_isnonzero && phi_epi_isnonzero) {
       // if phi_lv is zero, phi_rv and phi_epi is nonzero
       bislerp(QPfib, QPrv, QPepi, frac_epi);
       //QPfib=QPrv;
    } else if (phi_lv_isnonzero && !phi_rv_isnonzero && phi_epi_isnonzero) {
       // if phi_rv is zero, phi_lv and phi_epi is nonzero
       bislerp(QPfib, QPlv, QPepi, frac_epi);
       //QPfib=QPlv;
    } else if (phi_lv_isnonzero && phi_rv_isnonzero && !phi_epi_isnonzero) {
       // if phi_epi gradients are zero, phi_lv and phi are nonzero. use QPlv, QPrv, frac 
       bislerp(QPfib, QPlv, QPrv, frac);
       //QPfib=QPlv;
    } else {
       // if gradients of phi_lv, phi_rv are zero, then phi_epi is zero 
       vectorEigen(psi_ab_vec, QPfib);
    }               
}

void genfiber(vector<DenseMatrix>& QPfibVectors,
        vector<double>& psi_ab, vector<Vector>& psi_ab_grads,
        vector<double>& phi_epi, vector<Vector>& phi_epi_grads,
        vector<double>& phi_lv, vector<Vector>& phi_lv_grads,
        vector<double>& phi_rv, vector<Vector>& phi_rv_grads,
        Option& options
        ){
    
    unsigned nv=psi_ab.size();
    // Line 7 start for-loop
    for(unsigned i=0; i <nv; i++){  
//        MFEM_ASSERT(phi_lv[i]>=0 && phi_lv[i] <=1, "phi_lv is not in range 0 to 1");
//        MFEM_ASSERT(phi_rv[i]>=0 && phi_rv[i] <=1, "phi_rv is not in range 0 to 1");
//        MFEM_ASSERT(phi_epi[i]>=0 && phi_epi[i] <=1, "phi_epi is not in range 0 to 1");
//        MFEM_ASSERT(psi_ab[i]>=0 && psi_ab[i] <=1, "psi_ab is not in range 0 to 1");
        //if(phi_lv[i] <0) phi_lv[i]=0;
        //if(phi_rv[i] <0) phi_rv[i]=0;
        //if(phi_epi[i] <0) phi_epi[i]=0;
        //if(psi_ab[i] <0) psi_ab[i]=0;

        Vector psi_ab_vec=psi_ab_grads[i];
        Vector phi_lv_vec=phi_lv_grads[i];
        Vector phi_rv_vec=phi_rv_grads[i];
        Vector phi_epi_vec=phi_epi_grads[i];
        DenseMatrix QPfib(dim,dim);
        
        biSlerpCombo(QPfib, psi_ab[i], psi_ab_vec, phi_epi[i], phi_epi_vec, phi_lv[i], phi_lv_vec, phi_rv[i], phi_rv_vec, options);

        QPfibVectors.push_back(QPfib);
//        vector<Vector> qpVecs;
//        for(int j=0; j<dim; j++){
//            Vector vec;
//            QPfib.GetColumn(j, vec);
//            qpVecs.push_back(vec);
//        }
//        fvectors.push_back(qpVecs[0]);
//        svectors.push_back(qpVecs[1]);
//        tvectors.push_back(qpVecs[2]);
    }
    
//    ofstream f_ofs("fvectors.vtk");
//    ofstream s_ofs("svectors.vtk");
//    ofstream t_ofs("tvectors.vtk");
//
//    printFiberVTK(mesh, fvectors, f_ofs);
//    printFiberVTK(mesh, svectors, s_ofs);
//    printFiberVTK(mesh, tvectors, t_ofs);
//    
    
}


