#include "mfem.hpp"
#include "cardiac_physics.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;

int run_mode = 1;

void ReferenceConfiguration(const Vector &x, Vector &y)
{
   // set the reference, stress
   // free, configuration
   y = x;
}


void InitialDeformation(const Vector &x, Vector &y)
{
   // set the initial configuration. Having this different from the
   // reference configuration can help convergence
   y = x;

   if (run_mode == 1) {
      y(1) = x(1) - 0.05 * x(0);
   }
   
}

// Define the fiber directions. This will probably become either
// a grid function or quadrature based coefficient depending on the
// output of fiber gen.
void FiberFunction(const Vector &x, Vector &y)
{
   y = 0.0;

   y(0) = 1.0;

   y /= y.Norml2();

}

double PressureFunction(const Vector &x, double t)
{

   // Assume final time of 1

   double pres = 0.0;

   if (run_mode == 1) {
      pres =  4.0e-3 * (t / 1.0);
   }
   else if (run_mode == 4) {
      pres =  4.0e-4 * (t / 1.0);
   }
   else {
      pres =  10.0 * (t / 1.0);
   }
   return pres;
}

// Surface integral function used to calculate volume
void VolumeFunction(const Vector &x, Vector &y)
{
   y(2) = 0.0;

   y(0) = 0.5 * x(0);
   y(1) = 0.5 * x(1);
}

// Mesh helper functions from fiber gen to correctly label a heart-like mesh
void findNeighbor(Element* ele, vector<Element*>& elements, int attr){
    const int *v = ele->GetVertices();
    const int nv = ele->GetNVertices(); 
    for(unsigned i=0; i<elements.size(); i++){
        Element* queryEle=elements[i];
        // Only search for elements with unassigned attributes.
        if(queryEle->GetAttribute()==999){ 
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
    const int apexAttr=2;
    const int baseAttr=1;
    const int epiAttr=2;    
    const int vAttr=4;
       
    // Determine the max and min dimension of mesh and apex.
    double *coord;
    Vector coord_min(3);
    Vector coord_max(3);
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
                coord_min(j)=coord[j];
                coord_max(j)=coord[j];
            }            
        }
        
        for(int j=0; j<nv; j++){
            coord=mesh->GetVertex(v[j]);
            
            for (int k = 0; k < 3; k++) {
                if(coord[k]<coord_min[k]){
                    coord_min(k)=coord[k];
                    // Keep track vertex and element indeces for min in z-axis
                    if(k==2){  
                        apexVet=v[j];
                        apexEleIndex=i;
                    }
                }
                if(coord[k]>coord_max[k]){
                    coord_max(k)=coord[k];
                }            
            }                                    
        }
        
    }
    

    coord = mesh->GetVertex(apexVet);
    
    // Top 5% of the z axis.
    double zTop5=coord_max[2]-(coord_max[2]-coord_min[2])*0.05;

    // Initialization the attributes to 0 and set attribute of apex
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        // initialize the attribute for boundary.  
        ele->SetAttribute(999);
        
        //Found apex elements and set attribute.
        for (int j = 0; j < nv; j++) {
            if (v[j] ==apexVet){
                ele->SetAttribute(apexAttr);
                //cout << "Element index = " << i << endl;
            }
        }        
    }
    
    // Base    
    // The base must be planar. Its norm must be within 20 degrees of z axis.
    double angle = 20.0;
    double cosTheta = cos(angle*M_PI/180); 
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
        if(ele->GetAttribute()==999){
            elements.push_back(ele);
        }
    }
    
    Element *apexEle=mesh->GetBdrElement(apexEleIndex);
    findNeighbor(apexEle, elements, epiAttr);
    
    // V
    vector<Element *> vElements;    
    for(unsigned i=0; i<elements.size(); i++){
        Element *ele =elements[i];
        if(ele->GetAttribute()==999){
            vElements.push_back(ele);
        }
    }
    // pick one element in the container and assigned it to a temporary attr value.
    int last=vElements.size()-1;
    Element *lastEle=vElements[last];
    lastEle->SetAttribute(vAttr);
    // get rid of last element in the container
    vElements.pop_back();
    findNeighbor(lastEle, vElements, vAttr);
    
    // Check if there are unassigned elements left.
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i); 
        MFEM_ASSERT(ele->GetAttribute()!=999, "Unassigned element.");
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
