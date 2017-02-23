#include "utils.h"

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

