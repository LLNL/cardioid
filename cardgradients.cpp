#include "mfem.hpp"
#include "cardfiber.h"
#include "genfiber.h"
#include "constants.h"
#include "triplet.h"
#include "io.h"
#include "kdtree++/kdtree.hpp"
#include "cardgradients.h"

void getCardGradients(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, vector<Vector>& boundingbox, double dd,  
        Vector& conduct, Vector& fiberAngles, double maxEdgeLen){
    Vector min=boundingbox[0];
    Vector max=boundingbox[1];
    
    double xmin=min(0);
    double ymin=min(1);
    double zmin=min(2);
    
    double x_dim=max(0)-xmin;
    double y_dim=max(1)-ymin;
    double z_dim=max(2)-zmin;
    
    double dd10=dd*10; 
    
    double dx=x_dim/(int(x_dim/(dd10))*10-1);
    double dy=y_dim/(int(y_dim/(dd10))*10-1);
    double dz=z_dim/(int(z_dim/(dd10))*10-1);
    
    int nx=int(x_dim/dx)+1;
    int ny=int(y_dim/dy)+1;
    int nz=int(z_dim/dz)+1;
    
    long long totalCardPoints=0;
    vector<anatomy> anatVectors;
    
    double cutoff=maxEdgeLen*0.6124;  //Radius of circumsphere sqrt(6)/4 
    cout << "\nCutoff for nearest point is " << cutoff << std::endl;
       
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                double x=xmin+i*dx;
                double y=ymin+j*dy;
                double z=zmin+k*dz;
                //For k-D tree
                triplet pt(x, y, z, 0);                
                std::pair<tree_type::const_iterator,double> found = kdtree.find_nearest(pt);
                assert(found.first != kdtree.end());
                // Skip if the distance between pt and nearest is larger than cutoff
                if (found.second>cutoff) continue; 

                //For barycentric
                Vector q(4);
                q(0)=x;
                q(1)=y;
                q(2)=z;
                q(3)=1.0;                 
                
                triplet vetexNearPt=*found.first;
                int vertex=vetexNearPt.getIndex();
                vector<int> elements = vert2Elements[vertex];
                for (unsigned e = 0; e < elements.size(); e++) {
                    int eleIndex=elements[e];                   
                    if(isInTetElement(q, mesh, eleIndex)){
                        anatomy anat;
                        anat.gid=i+j*nx+k*nx*ny;
                        calcGradient(x_psi_ab, x_phi_epi, x_phi_lv, x_phi_rv, conduct, fiberAngles, q, eleIndex, anat);  
                        anatVectors.push_back(anat); 
                        totalCardPoints++;
                        if(totalCardPoints%10000==0){ 
                            cout << "Finish " << totalCardPoints << " points." << endl;
                            cout.flush();
                        }                          
                        break; // If the point is found in an element, don't need to check next one in the list. 
                    } 
                }
                
            }
        }
    }
    filerheader header;
    header.nx=nx;
    header.ny=ny;
    header.nz=nz;
    header.dx=dx;
    header.dy=dy;
    header.dz=dz;
    header.nrecord=totalCardPoints;
    
    ofstream anat_ofs("anatomy#000000");
    printAnatomy(anatVectors, header, anat_ofs);
    
    
}

template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters = " ", const bool trimEmpty = true) {
    typedef ContainerT Base;
    typedef typename Base::value_type ValueType;
    typedef typename ValueType::size_type SizeType;
    SizeType pos, lastPos = 0;
    while (true) {
        pos = str.find_first_of(delimiters, lastPos);
        if (pos == std::string::npos) {
            pos = str.length();

            if (pos != lastPos || !trimEmpty)
                tokens.push_back(ValueType(str.data() + lastPos, (SizeType)pos - lastPos));

            break;
        } else {
            if (pos != lastPos || !trimEmpty)
                tokens.push_back(ValueType(str.data() + lastPos, (SizeType)pos - lastPos));
        }

        lastPos = pos + 1;
    }
};



void getRotMatrix(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, Vector& fiberAngles, const char *fiblocs){

    long long totalCardPoints=0;

    ofstream f_ofs("rotmatrix.txt");
    f_ofs << "# elementnum mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33" << endl;
    
    ifstream f_ifs;
    try{
       f_ifs.open(fiblocs);
    }
    catch(...)
    {
       std::cout << "Cannot open file: " << fiblocs << std::endl;
       return;
    }

    std::string fileLine;
    
    const std::string comment="#";

   while (f_ifs)
   {
      std::getline(f_ifs, fileLine);
      if (fileLine.compare(0, 1, comment) == 0) continue;
      std::vector<std::string> tokens;
      tokenize(fileLine, tokens);
      if (tokens.size() > 3)
      {
         double x = atof(tokens[1].c_str());
         double y = atof(tokens[2].c_str());
         double z = atof(tokens[3].c_str());
         triplet pt(x, y, z, 0);
         std::pair<tree_type::const_iterator, double> found = kdtree.find_nearest(pt);
         assert(found.first != kdtree.end());
         triplet vetexNearPt = *found.first;
         int vertex = vetexNearPt.getIndex();
         vector<int> elements = vert2Elements[vertex];
         bool findPt=false;
         for (unsigned e = 0; e < elements.size(); e++)
         {
            int eleIndex = elements[e];
            //For barycentric
            Vector q(4);
            q(0) = x;
            q(1) = y;
            q(2) = z;
            q(3) = 1.0;
            vector<double> barycentric;
            if (isInTetElement(q, mesh, eleIndex))
            {
               //cout << "fiblocs element index=" << tokens[0] << "; k-D tree index=" << eleIndex << endl;
               Vector psi_ab_vec(3);
               double psi_ab = 0.0;
               getCardEleGrads(x_psi_ab, q, eleIndex, psi_ab_vec, psi_ab);

               Vector phi_epi_vec(3);
               double phi_epi = 0.0;
               getCardEleGrads(x_phi_epi, q, eleIndex, phi_epi_vec, phi_epi);

               Vector phi_lv_vec(3);
               double phi_lv = 0.0;
               getCardEleGrads(x_phi_lv, q, eleIndex, phi_lv_vec, phi_lv);

               Vector phi_rv_vec(3);
               double phi_rv = 0.0;
               getCardEleGrads(x_phi_rv, q, eleIndex, phi_rv_vec, phi_rv);

               DenseMatrix QPfib(dim3, dim3);
               biSlerpCombo(QPfib, psi_ab, psi_ab_vec, phi_epi, phi_epi_vec,
                            phi_lv, phi_lv_vec, phi_rv, phi_rv_vec, fiberAngles);
               
               f_ofs << tokens[0] << " ";
               for(int ii=0; ii<dim3; ii++){
                  for(int jj=0; jj<dim3; jj++){
                     f_ofs << QPfib(ii, jj) << " ";
                  }
               }
               f_ofs << endl;
               
               totalCardPoints++;
               if (totalCardPoints % 10000 == 0)
               {
                  cout << "Finish " << totalCardPoints << " points." << endl;
                  cout.flush();
               }
               findPt=true;
               break; // If the point is found in an element, don't need to check next one in the list. 

            }
                    
         }
         if(!findPt){

                f_ofs << tokens[0] << " ";
               for (int ii = 0; ii < dim3; ii++)
               {
                  for (int jj = 0; jj < dim3; jj++)
                  {
                     f_ofs << "999 ";
                  }
               }
               f_ofs << endl;            
                  
         }
      }
   }

          
}


void getRotMatrixFast(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        vector<vector<int> >& vert2Elements, Vector& fiberAngles, const char *fiblocs){

    long long totalCardPoints=0;

    ofstream f_ofs("rotmatrix.txt");
    f_ofs << "# elementnum mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33" << endl;
    
    ifstream f_ifs;
    try{
       f_ifs.open(fiblocs);
    }
    catch(...)
    {
       std::cout << "Cannot open file: " << fiblocs << std::endl;
       return;
    }

    std::string fileLine;
    
    const std::string comment="#";

   while (f_ifs)
   {
      std::getline(f_ifs, fileLine);
      if (fileLine.compare(0, 1, comment) == 0) continue;
      std::vector<std::string> tokens;
      tokenize(fileLine, tokens);
      if (tokens.size() > 3)
      {
         int eleIndex=atoi(tokens[0].c_str());
         double x = atof(tokens[1].c_str());
         double y = atof(tokens[2].c_str());
         double z = atof(tokens[3].c_str());

            //For barycentric
            Vector q(4);
            q(0) = x;
            q(1) = y;
            q(2) = z;
            q(3) = 1.0;
            vector<double> barycentric;
            if (isInTetElement(q, mesh, eleIndex))
            {
               //cout << "fiblocs element index=" << tokens[0] << "; k-D tree index=" << eleIndex << endl;
               Vector psi_ab_vec(3);
               double psi_ab = 0.0;
               getCardEleGrads(x_psi_ab, q, eleIndex, psi_ab_vec, psi_ab);

               Vector phi_epi_vec(3);
               double phi_epi = 0.0;
               getCardEleGrads(x_phi_epi, q, eleIndex, phi_epi_vec, phi_epi);

               Vector phi_lv_vec(3);
               double phi_lv = 0.0;
               getCardEleGrads(x_phi_lv, q, eleIndex, phi_lv_vec, phi_lv);

               Vector phi_rv_vec(3);
               double phi_rv = 0.0;
               getCardEleGrads(x_phi_rv, q, eleIndex, phi_rv_vec, phi_rv);

               DenseMatrix QPfib(dim3, dim3);
               biSlerpCombo(QPfib, psi_ab, psi_ab_vec, phi_epi, phi_epi_vec,
                            phi_lv, phi_lv_vec, phi_rv, phi_rv_vec, fiberAngles);
               
               f_ofs << tokens[0] << " ";
               for(int ii=0; ii<dim3; ii++){
                  for(int jj=0; jj<dim3; jj++){
                     f_ofs << QPfib(ii, jj) << " ";
                  }
               }
               f_ofs << endl;
               
               totalCardPoints++;
               if (totalCardPoints % 10000 == 0)
               {
                  cout << "Finish " << totalCardPoints << " points." << endl;
                  cout.flush();
               }

            }
                    
         

      }
   }

          
}


void calcNodeFiber(vector<DenseMatrix>& QPfibVectors){
    
    // Set up diag matrix
    // [ 3 0 0 
    //   0 2 0 
    //   0 0 1 ]    
    DenseMatrix diag(dim3, dim3);
    for(int i=0; i<dim3; i++){
        Vector vec(3);
        vec=0.0;
        vec(i)=3-i;
        diag.SetCol(i, vec);
    }
 
    ofstream f_ofs("heart.fiber");
    
     for(unsigned i=0; i< QPfibVectors.size(); i++){
        DenseMatrix Q=QPfibVectors[i]; 
        DenseMatrix tmp(dim3, dim3);
        Mult(Q, diag, tmp);    
        DenseMatrix QT=Q;
        QT.Transpose();  
        DenseMatrix Sigma(dim3, dim3);
        Mult(tmp, QT, Sigma);
        f_ofs << i << " " 
             << Sigma(0,0) << " " << Sigma(1,0) << " " << Sigma(2,0) << " " 
             << Sigma(1,1) << " " << Sigma(2,1) << " " << Sigma(2,2) 
             << std::endl;
        
    }   
    
}