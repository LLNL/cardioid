//file modified from example driver

//Generic routines
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"

using namespace std;
using namespace oomph;



//=======================start_mesh========================================
/// Triangle-based mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class VentricularTetMesh : public virtual TetgenMesh<ELEMENT>, 
                     public virtual SolidMesh 
{
 
public:
 
 /// Constructor: 
 VentricularTetMesh(const std::string& node_file_name,
              const std::string& element_file_name,
              const std::string& face_file_name,
              TimeStepper* time_stepper_pt=
              &Mesh::Default_TimeStepper) : 
  TetgenMesh<ELEMENT>(node_file_name, element_file_name,
                      face_file_name, time_stepper_pt)
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();
   
   // Find elements next to boundaries
   setup_boundary_element_info();
  }

 /// Empty Destructor
 virtual ~VentricularTetMesh() { }


};


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Parameters
{

 /// Poisson's ratio
 double Nu=0.3;
 
 /// Create constitutive law
 ConstitutiveLaw* Constitutive_law_pt=new GeneralisedHookean(&Nu);

 /// Non-dim gravity
 double Gravity=0.0;

 /// Non-dimensional gravity as body force
 void gravity(const double& time,
              const Vector<double> &xi,
              Vector<double> &b)
 {
  b[0]=0.0; 
  b[1]=0.0;
  b[2]=-Gravity;
 } // end gravity

 /// Uniform pressure at LV
  double LV_P = 0.0;
  double RV_P = 0.0;


 /// \short Constant pressure load. The arguments to this function are imposed
 /// on us by the SolidTractionElements which allow the traction to 
 /// depend on the Lagrangian and Eulerian coordinates x and xi, and on the 
 /// outer unit normal to the surface. Here we only need the outer unit
 /// normal.
  void constant_pressure_LV(const Vector<double> &xi, const Vector<double> &x,
                            const Vector<double> &n, Vector<double> &traction)
  {
    unsigned dim = traction.size();
    for(unsigned i=0;i<dim;i++)
      {
        traction[i] = -LV_P*n[i];
      }
  } // end traction
  
  void constant_pressure_RV(const Vector<double> &xi, const Vector<double> &x,
                            const Vector<double> &n, Vector<double> &traction)
  {
    unsigned dim = traction.size();
    for(unsigned i=0;i<dim;i++)
      {
        traction[i] = -RV_P*n[i];
      }
  } // end traction

 
 
} //end namespace






//=============start_problem===========================================
/// Unstructured solid problem
//=====================================================================
template<class ELEMENT>
class UnstructuredSolidProblem : public Problem
{

public:

 /// Constructor: 
 UnstructuredSolidProblem();

 /// Destructor (empty)
 ~UnstructuredSolidProblem(){}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Create traction elements
 void create_traction_elements();

 /// Bulk solid mesh
 VentricularTetMesh<ELEMENT>* Ventricular_mesh_pt;

 /// Meshes of traction elements
 Vector<SolidMesh*> Solid_traction_mesh_pt;

 /// IDs of solid mesh boundaries where displacements are pinned
 Vector<unsigned> Pinned_solid_boundary_id;

 /// \short IDs of solid mesh boundaries which make up the traction interface
 Vector<unsigned> Solid_traction_boundary_id;

};



//=============start_constructor==========================================
/// Constructor for unstructured solid problem
//========================================================================
template<class ELEMENT>
UnstructuredSolidProblem<ELEMENT>::UnstructuredSolidProblem()
{ 

  //Create solid bulk mesh
 string node_file_name="canine.node";
 string element_file_name="canine.ele";
 string face_file_name="canine.face";

 
 Ventricular_mesh_pt =  new VentricularTetMesh<ELEMENT>(node_file_name,
                                            element_file_name,
                                            face_file_name);
  // The following IDs corresponds to the boundary IDs specified in
 // the *.poly file from which tetgen generated the unstructured mesh.
 
 /// IDs of solid mesh boundaries where displacements are pinned
 Pinned_solid_boundary_id.resize(3);
 Pinned_solid_boundary_id[0]=0; // pinned z - coordinate  
 Pinned_solid_boundary_id[1]=1; // pinned x, z - coordinate  
 Pinned_solid_boundary_id[2]=2; // pinned y, z - coordinate  
 
  // The solid mesh boundaries where an internal pressure is applied
 Solid_traction_boundary_id.resize(2);
 Solid_traction_boundary_id[0]=4;
 Solid_traction_boundary_id[1]=5;
 
 
 // Apply BCs for solid
 //--------------------
 
 // Doc pinned solid nodes
 std::ofstream bc_file("RESLT/pinned_solid_nodes.dat");
 
 // Pin positions at inflow boundary (boundaries 0 and 1)
 unsigned n=Pinned_solid_boundary_id.size();
 for (unsigned i=0;i<n;i++)
  {
   // Get boundary ID
   unsigned b=Pinned_solid_boundary_id[i];
   unsigned num_nod= Ventricular_mesh_pt->nboundary_node(b);  
   for (unsigned inod=0;inod<num_nod;inod++)
    {    
     // Get node
     SolidNode* nod_pt=Ventricular_mesh_pt->boundary_node_pt(b,inod);
     
     switch(b){

     case 0:
       nod_pt->pin_position(2);
       bc_file << nod_pt->x(2) << " ";
       
     case 1:
       nod_pt->pin_position(0);       
       nod_pt->pin_position(2);       
       bc_file << nod_pt->x(0) << " " << nod_pt->x(2) << " " ;       
       break;
       
     case 2:
       nod_pt->pin_position(1);       
       nod_pt->pin_position(2);       
       bc_file << nod_pt->x(1) << " " << nod_pt->x(2) << " " ;       
       break;

     } 
     
     bc_file << std::endl;
    }
  }
 bc_file.close();
 
 
 
 // Complete the build of all elements so they are fully functional
 //----------------------------------------------------------------
 unsigned n_element = Ventricular_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(
    Ventricular_mesh_pt->element_pt(i));
   
   // Set the constitutive law   
   el_pt->constitutive_law_pt() =
    Global_Parameters::Constitutive_law_pt;
   
   //Set the body force
   el_pt->body_force_fct_pt() = Global_Parameters::gravity;
  }
 
 
 // Create traction elements
 //-------------------------
 
 // Create meshes of traction elements
 n=Solid_traction_boundary_id.size();
 Solid_traction_mesh_pt.resize(n);
 for (unsigned i=0;i<n;i++)
  {
   Solid_traction_mesh_pt[i]=new SolidMesh;
  }
 
 // Build the traction elements
 create_traction_elements();


 // Combine the lot
 //----------------
 
 // The solid bulk mesh
 add_sub_mesh(Ventricular_mesh_pt);

 // The solid traction meshes
 n=Solid_traction_boundary_id.size();
 for (unsigned i=0;i<n;i++)
  {
   add_sub_mesh(Solid_traction_mesh_pt[i]);
  }

 // Build global mesh
 build_global_mesh();

 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end constructor



//============start_of_create_traction_elements==========================
/// Create traction elements 
//=======================================================================
template<class ELEMENT>
void UnstructuredSolidProblem<ELEMENT>::create_traction_elements()
{

 // Loop over traction boundaries
 unsigned n=Solid_traction_boundary_id.size();
 for (unsigned i=0;i<n;i++)
  {
   // Get boundary ID
   unsigned b=Solid_traction_boundary_id[i];
   
   // How many bulk elements are adjacent to boundary b?
   unsigned n_element = Ventricular_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Ventricular_mesh_pt->boundary_element_pt(b,e));
     
     //What is the index of the face of the element e along boundary b
     int face_index = Ventricular_mesh_pt->face_index_at_boundary(b,e);
     
     // Create new element 
     SolidTractionElement<ELEMENT>* el_pt=
      new SolidTractionElement<ELEMENT>(bulk_elem_pt,face_index);
     
     // Add it to the mesh
     Solid_traction_mesh_pt[i]->add_element_pt(el_pt);
     
     //Set the traction function
     if(b==4)
       el_pt->traction_fct_pt() = Global_Parameters::constant_pressure_LV;
     else
       el_pt->traction_fct_pt() = Global_Parameters::constant_pressure_RV;
    }
  }
 
} // end of create_traction_elements



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredSolidProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 // Output solid solution
 //-----------------------
 sprintf(filename,"%s/solid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Ventricular_mesh_pt->output(some_file,npts);
 some_file.close();

    
 // Output traction
 //----------------
 sprintf(filename,"%s/traction%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned n=Solid_traction_boundary_id.size();
 for (unsigned i=0;i<n;i++)
  {
   Solid_traction_mesh_pt[i]->output(some_file,npts);
  }
 some_file.close();

} // end doc





//============================start_main==================================
/// Demonstrate how to solve an unstructured 3D solid problem
//========================================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 //Set up the problem
 UnstructuredSolidProblem<TPVDElement<3,3> > problem;
 
 //Output initial configuration
 problem.doc_solution(doc_info);
 doc_info.number()++;   

  // Parameter study
 Global_Parameters::LV_P=0.0; 
 Global_Parameters::RV_P=0.0; 
 double g_increment=1.0e-3/10.0; 
 double p_increment=1.0e-2; 

 unsigned nstep=6;
 if (CommandLineArgs::Argc!=1)
  {
   std::cout << "Validation -- only doing two steps" << std::endl;
   nstep=2;
  }
 
 // Do the parameter study
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.newton_solve();
   
   //Output solution
   problem.doc_solution(doc_info);
   doc_info.number()++;

   // Bump up load
   Global_Parameters::Gravity+=g_increment;
   Global_Parameters::LV_P+=p_increment; 
   Global_Parameters::RV_P+=p_increment; 
   
  }
 
} // end main




