//#define KDTREE_DEFINE_OSTREAM_OPERATORS   // Position sensitive.

#include "mfem.hpp"
#include "cardfiber.h"
#include "genfiber.h"
#include "constants.h"
#include "triplet.h"
#include "io.h"
#include "kdtree++/kdtree.hpp"

#include <cmath>

void buildKDTree(Mesh *mesh, tree_type& kdtree) {
    int NumOfVertices = mesh->GetNV();

    for (int i = 0; i < NumOfVertices; i++){
    //for (int i = 0; i < 10; i++) {
        const double* coord = mesh->GetVertex(i);
        triplet node(coord[0], coord[1], coord[2], i);
        kdtree.insert(node);
    }
    kdtree.optimise();
//    triplet t(97.3604, 40.8014, 63.152, 0);
//      std::pair<tree_type::const_iterator,double> found = kdtree.find_nearest(t);
//      assert(found.first != kdtree.end());
//      std::cout << "Test with search type, found: " << *found.first << ", wanted " << t << std::endl;
    
//    cout << kdtree << endl;
//    for (tree_type::const_iterator target = kdtree.begin(); target != kdtree.end(); ++target) {
//        cout << *target << endl;
//    }
}

double getMaxEdgeLen(Mesh *mesh){
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
            }
         }
        //break;
    }
    return std::sqrt(maxEdgeLen);
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
       for(int j=0; j<dim; j++){
           vec(j)=coord[j];
       }        
        vec(dim)=1.0;
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

void calcGradient(GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        Option& options, Vector& q, int eleIndex, DenseMatrix& QPfib, Phi& phi){
    
    Vector psi_ab_vec(3);                        
    double psi_ab=0.0;
    getCardEleGrads(x_psi_ab, q, eleIndex, psi_ab_vec, psi_ab);

    Vector phi_epi_vec(3);
    double phi_epi=0.0;
    getCardEleGrads(x_phi_epi, q, eleIndex, phi_epi_vec, phi_epi);                        

    Vector phi_lv_vec(3);
    double phi_lv=0.0;
    getCardEleGrads(x_phi_lv, q, eleIndex, phi_lv_vec, phi_lv);  

    Vector phi_rv_vec(3);
    double phi_rv=0.0;
    getCardEleGrads(x_phi_rv, q, eleIndex, phi_rv_vec, phi_rv);  

    phi.epi=phi_epi;
    phi.lv=phi_lv;
    phi.rv=phi_rv;
    
    biSlerpCombo(QPfib, psi_ab, psi_ab_vec, phi_epi, phi_epi_vec,
        phi_lv, phi_lv_vec, phi_rv, phi_rv_vec, options); 

 
}

void getAnatomy(anatomy& anat, DenseMatrix& QPfib, Option& options, Phi& phi,
        ThreeInts& inds, ThreeInts& nns){
    DenseMatrix Sigma(dim, dim);
    calcSigma(Sigma, QPfib, options);

    double *s=Sigma.Data();
    
    anat.gid=inds.i + inds.j*nns.i + inds.k*nns.i*nns.j;
    anat.celltype=getCellType(phi);
    anat.sigma[0]=s[0];
    anat.sigma[1]=s[3];
    anat.sigma[2]=s[6];
    anat.sigma[3]=s[4];
    anat.sigma[4]=s[7];
    anat.sigma[5]=s[8];      
}


bool findPtEle(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        vector<vector<int> >& vert2Elements, Option& options, Vector& q, int vertex, std::string& elemnum, ostream& out){
    
         vector<int> elements = vert2Elements[vertex];          
         bool findPt=false;        
         for (unsigned e = 0; e < elements.size(); e++)
         {
            int eleIndex = elements[e];

            if (isInTetElement(q, mesh, eleIndex))
            {
                DenseMatrix QPfib(dim, dim);
                Phi phi;
                calcGradient(x_psi_ab, x_phi_epi, x_phi_lv, x_phi_rv, options, q, eleIndex, QPfib, phi);  
                
               out << elemnum << " ";
               for(int ii=0; ii<dim; ii++){
                  for(int jj=0; jj<dim; jj++){
                     out << QPfib(ii, jj) << " ";
                  }
               }
               out << endl;
               
               findPt=true;
               break; // If the point is found in an element, don't need to check next one in the list. 

            }
                    
         }
         return findPt;
    
}

bool findPtEleAnat(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        vector<vector<int> >& vert2Elements, Option& options, 
        Vector& q, int vertex, ThreeInts& inds, ThreeInts& nns, vector<anatomy>& anatVectors){
    vector<int> elements = vert2Elements[vertex];
    bool foundPt=false;
    for (unsigned e = 0; e < elements.size(); e++) {
        int eleIndex=elements[e];                   
        if(isInTetElement(q, mesh, eleIndex)){
            DenseMatrix QPfib(dim, dim);
            Phi phi;
            calcGradient(x_psi_ab, x_phi_epi, x_phi_lv, x_phi_rv, options, q, eleIndex, QPfib, phi);  
            anatomy anat;
            
            getAnatomy(anat, QPfib, options, phi, inds, nns);
            anatVectors.push_back(anat);    
            foundPt=true;
            break; // If the point is found in an element, don't need to check next one in the list. 
        } 
    }
    
    return foundPt;
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

void calcSigma(DenseMatrix& Sigma, DenseMatrix& Q, Option& options){
    Vector conduct(3);
    conduct(0)=options.gL;
    conduct(1)=options.gT;
    conduct(2)=options.gN;
    DenseMatrix diag(dim, dim);
    //Get Diag (conductivity)).
    for(int i=0; i<dim; i++){
        Vector vec(3);
        vec=0.0;
        vec(i)=conduct(i);
        diag.SetCol(i, vec);
    }
    
    DenseMatrix tmp(dim, dim);
    Mult(Q, diag, tmp);
    DenseMatrix QT=Q;
    QT.Transpose();
    Mult(tmp, QT, Sigma);
    
}

int getCellType(Phi& phi){
    if(phi.epi>0.66) return 102;                // EPI_CELL=102
    if(phi.lv>0.66 || phi.rv >0.66) return 100; // ENDO_CELL=100
    return 101;                                 // M_CELL=101
}



