#define KDTREE_DEFINE_OSTREAM_OPERATORS   // Position sensitive.

#include "mfem.hpp"
#include "cardfiber.h"
#include "genfiber.h"
#include "constants.h"
#include "triplet.h"
#include "io.h"
#include "kdtree++/kdtree.hpp"

void buildKDTree(Mesh *mesh, tree_type& kdtree) {
    int NumOfVertices = mesh->GetNV();

    for (int i = 0; i < NumOfVertices; i++){
    //for (int i = 0; i < 10; i++) {
        const double* coord = mesh->GetVertex(i);
        triplet node(coord[0], coord[1], coord[2], i);
        kdtree.insert(node);
    }
    
//    triplet t(97.3604, 40.8014, 63.152, 0);
//      std::pair<tree_type::const_iterator,double> found = kdtree.find_nearest(t);
//      assert(found.first != kdtree.end());
//      std::cout << "Test with search type, found: " << *found.first << ", wanted " << t << std::endl;
    
//    cout << kdtree << endl;
//    for (tree_type::const_iterator target = kdtree.begin(); target != kdtree.end(); ++target) {
//        cout << *target << endl;
//    }
}


// MFEM doesn't have determinant function for 4X4 matrix. So write my own version.
double det4X4(DenseMatrix& matrix){
    double* m=matrix.Data();
//    for(int i=0; i<16; i++){
//        cout << m[i] << " ";
//    }
//    cout <<endl;     
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

bool isInTetElement(const Vector& q, Mesh* mesh, int eleIndex){
    vector<double> barycentric;
    vector<Vector> VV;    
    const Element* ele=mesh->GetElement(eleIndex);
    const int *v = ele->GetVertices();
    const int nv = ele->GetNVertices();
    MFEM_ASSERT(nv==4, "Tetrahedron Element should contain 4 vertex.");
    for (int i = 0; i < nv; i++) {
        Vector vec(4);
        int vert=v[i];
        const double* coord=mesh->GetVertex(vert);
       for(int j=0; j<dim3; j++){
           vec(j)=coord[j];
       }        
        vec(dim3)=1.0;
        VV.push_back(vec);
    }    

    DenseMatrix D0(4, 4);
    for (int j = 0; j < 4; j++) {
        D0.SetRow(j, VV[j]);
    }
    double d0=det4X4(D0);
       
    for (int i = 0; i < 4; i++) {
        DenseMatrix D(4, 4);
        for (int j = 0; j < 4; j++) {
            if (j == i) {
                D.SetRow(j, q);
            } else {
                D.SetRow(j, VV[j]);
            }
        }
        double d = det4X4(D);
        double bc=d/d0;
        barycentric.push_back(bc);
        if(bc<0 || bc>1){
            return false; // if any barycentric coord is negative, point is not in tet.
        }
    }
    
    double sumbc=0.0;
    for(unsigned i=0;i<barycentric.size(); i++){
        sumbc+=barycentric[i];
    }
    if(sumbc<0.999){
        cout << "Warning: isInTetElement summation of barycentric coords is not equal to 1.";
        return false;
    }
    
    return true;
    
}

void getCardEleGrads(GridFunction& x, const Vector& q, int eleIndex, Vector& grad_ele, double& xVal) {
    //initialize xVal grad_ele to 0.0;
    xVal=0.0;
    grad_ele=0.0;
    const FiniteElementSpace *fes = x.FESpace();
    ElementTransformation * tr = fes->GetElementTransformation(eleIndex);
    //Coordinates
    Vector coor(3);
    coor(0)=q(0);
    coor(1)=q(1);
    coor(2)=q(2);
    
    const IntegrationRule &ir = fes->GetFE(eleIndex)->GetNodes(); // Get the parametric integration rule

    for (int k = 0; k < ir.GetNPoints(); k++) {
        Vector grad_point(3);
        grad_point = 0.0;
        IntegrationPoint ip; // The current integration point is calculated from tr->TransformBack. 
        tr->TransformBack(coor, ip);
        //ip.weight=barycentric[k];
        
        //tr->SetIntPoint(&ip); // Set the integration point for the transformation
        x.GetGradient((*tr), grad_point);
        xVal+=x.GetValue(eleIndex, ip, 1);
        grad_ele += grad_point;
    }
    grad_ele /= ir.GetNPoints();
    xVal /= ir.GetNPoints();
    
}

void calcSigma(DenseMatrix& Sigma, DenseMatrix& Q, Vector& conduct){

    DenseMatrix diag(dim3, dim3);
    //Get Diag (conductivity)).
    for(int i=0; i<dim3; i++){
        Vector vec(3);
        vec=0.0;
        vec(i)=conduct(i);
        diag.SetCol(i, vec);
    }
    
    DenseMatrix tmp(dim3, dim3);
    Mult(Q, diag, tmp);
    Q.Transpose();
    Mult(tmp, Q, Sigma);
    
}

int getCellType(double phi_epi, double phi_lv, double phi_rv){
    if(phi_epi>0.66) return 102;                // EPI_CELL=102
    if(phi_lv>0.66 || phi_rv >0.66) return 100; // ENDO_CELL=100
    return 101;                                 // M_CELL=101
}



