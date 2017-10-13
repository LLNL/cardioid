#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace mfem;

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
    const int apexAttr=4;
    const int baseAttr=1;
    const int epiAttr=3;
    const int lvAttr=2;
    const int rvAttr=5;
       
    // Determine the max and min dimension of mesh and apex.
    double *coord;
    double coord_min[3];
    double coord_max[3];
    bool firstEle=true;
    int apexVet=0;
    int apexEleIndex=0;
    
    // First set the domain elements
    int ne=mesh->GetNE();
    for (int i=0; i<ne; i++) {
       Element * el = mesh->GetElement(i);
       el->SetAttribute(1);       
    }

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
