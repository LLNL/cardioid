#include "cardiac_problem.h"
#include <string>


namespace oomph
{

// define template functions for AnisotropicTPVDElement<3, 3> in this object file
template class VentricularProblem<ActiveAnisotropicTPVDElement<3, 3, BasicActiveModel> >;
template VentricularProblem<ActiveAnisotropicTPVDElement<3, 3, BasicActiveModel> >::VentricularProblem();
template void VentricularProblem<ActiveAnisotropicTPVDElement<3, 3, BasicActiveModel> >::doc_solution(DocInfo& doc_info);
template void VentricularProblem<ActiveAnisotropicTPVDElement<3, 3, BasicActiveModel> >::set_pressure(double lv, double rv);


inline void gravity(const double& time,
                    const Vector<double> &xi,
                    Vector<double> &b)
{
  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
}

namespace ventricular_pressure
{

double LV_P = 0, RV_P = 0;

void constant_pressure_LV(const Vector<double> &xi, const Vector<double> &x,
                          const Vector<double> &n, Vector<double> &traction)
{
  unsigned dim = traction.size();

  for(unsigned i = 0; i < dim; i++) {
    traction[i] = -LV_P * n[i];
  }
} // end traction

void constant_pressure_RV(const Vector<double> &xi, const Vector<double> &x,
                          const Vector<double> &n, Vector<double> &traction)
{
  unsigned dim = traction.size();

  for(unsigned i = 0; i < dim; i++) {
    traction[i] = -RV_P * n[i];
  }
} // end traction

}


template <class ELEMENT>
void VentricularProblem<ELEMENT>::set_pressure(double lv, double rv)
{
  ventricular_pressure::LV_P = lv;
  ventricular_pressure::RV_P = rv;
}

template <class ELEMENT>
VentricularProblem<ELEMENT>::VentricularProblem()
{

  Simulation_time = 0;
  // set Usyk constitutive law with active tension
  constitutive_law_pt = new ActiveConstitutiveLaw<UsykConstitutiveLaw>();

  // fix_later_tag
  // Jacobian reuse does not implemented properly in OOMPH.
  // Newton_solve() does not check if residual decreases or not.
  // I will change this behavior soon by overriding default Newton_solve() function.
  // Jacobian_reuse_is_enabled=true;

  Max_newton_iterations = 1000;
  Max_residuals = 1000;
  Newton_solver_tolerance = 1.0e-4;

  // using SuperLUSolver
  linear_solver_pt() = new SuperLUSolver;

  // Create solid bulk mesh from file.
  // File names are hard coded for now. fix_it_later_tag
  std::string node_file_name = "canine.node";
  std::string element_file_name = "canine.ele";
  std::string face_file_name = "canine.face";


  Ventricular_mesh_pt =  new AnisotropicTetMesh<ELEMENT >(node_file_name,
      element_file_name,
      face_file_name);

  this->assign_time_to_elements();
  // The following IDs corresponds to the boundary IDs specified in
  // the *.poly file from which tetgen generated the unstructured mesh.

  // IDs of solid mesh boundaries where displacements are pinned
  enum Pinned_coordinates {pinned_z = 0, pinned_x_and_z = 1, pinned_y_and_z = 2};


  Pinned_solid_boundary_id.resize(3);
  Pinned_solid_boundary_id[0] = pinned_z;
  Pinned_solid_boundary_id[1] = pinned_x_and_z;
  Pinned_solid_boundary_id[2] = pinned_y_and_z;

  enum Ventricular_surfaces {lv_surface = 4, rv_surface = 5};

  // The solid mesh boundaries where an internal pressure is applied
  Solid_traction_boundary_id.resize(2);
  Solid_traction_boundary_id[0] = lv_surface;
  Solid_traction_boundary_id[1] = rv_surface;



// Doc pinned solid nodes
  std::ofstream bc_file("RESLT/pinned_solid_nodes.dat");

// Pin positions at inflow boundary (boundaries 0 and 1)
  unsigned n = Pinned_solid_boundary_id.size();

  for (unsigned i = 0; i < n; i++) {
    // Get boundary ID
    unsigned b = Pinned_solid_boundary_id[i];
    unsigned num_nod = Ventricular_mesh_pt->nboundary_node(b);

    for (unsigned inod = 0; inod < num_nod; inod++) {
      // Get node
      SolidNode* nod_pt = Ventricular_mesh_pt->boundary_node_pt(b, inod);

      switch(b) {

      case pinned_z:
        nod_pt->pin_position(2);
        bc_file << nod_pt->x(2) << " ";

      case pinned_x_and_z:
        nod_pt->pin_position(0);
        nod_pt->pin_position(2);
        bc_file << nod_pt->x(0) << " " << nod_pt->x(2) << " " ;
        break;

      case pinned_y_and_z:
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

  for(unsigned i = 0; i < n_element; i++) {
    //Cast to a solid element
    ELEMENT *el_pt = dynamic_cast<ELEMENT*>(
                       Ventricular_mesh_pt->element_pt(i));

    // Set the constitutive law
    el_pt->constitutive_law_pt() =
      this->constitutive_law_pt;

    //Set the body force
    el_pt->body_force_fct_pt() = gravity;
  }


// Create traction elements
//-------------------------

// Create meshes of traction elements
  n = Solid_traction_boundary_id.size();
  Solid_traction_mesh_pt.resize(n);

  for (unsigned i = 0; i < n; i++) {
    Solid_traction_mesh_pt[i] = new SolidMesh;
  }

// Build the traction elements
  create_traction_elements();


// Combine the lot
//----------------

// The solid bulk mesh
  add_sub_mesh(Ventricular_mesh_pt);

// The solid traction meshes
  n = Solid_traction_boundary_id.size();

  for (unsigned i = 0; i < n; i++) {
    add_sub_mesh(Solid_traction_mesh_pt[i]);
  }



// Build global mesh
  build_global_mesh();


// Setup equation numbering scheme
  std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

  // distribute() does not work in OOMPH. They have removing nodes from halo elements, and then used pointer to removed nodes in classify_...
  // I am working to fix it
  //distribute();

} // end constructor


template <class ELEMENT>
void VentricularProblem<ELEMENT>::create_traction_elements()
{

// Loop over traction boundaries
  unsigned n = Solid_traction_boundary_id.size();

  for (unsigned i = 0; i < n; i++) {
    // Get boundary ID
    unsigned b = Solid_traction_boundary_id[i];

    // How many bulk elements are adjacent to boundary b?
    unsigned n_element = Ventricular_mesh_pt->nboundary_element(b);

    // Loop over the bulk elements adjacent to boundary b
    for(unsigned e = 0; e < n_element; e++) {
      // Get pointer to the bulk element that is adjacent to boundary b
      AnisotropicTPVDElement<3, 3>* bulk_elem_pt = dynamic_cast<AnisotropicTPVDElement<3, 3>*>(
            Ventricular_mesh_pt->boundary_element_pt(b, e));

      //What is the index of the face of the element e along boundary b
      int face_index = Ventricular_mesh_pt->face_index_at_boundary(b, e);

      // Create new element
      SolidTractionElement<AnisotropicTPVDElement<3, 3> >* el_pt =
        new SolidTractionElement<AnisotropicTPVDElement<3, 3> >(bulk_elem_pt, face_index);

      // Add it to the mesh
      Solid_traction_mesh_pt[i]->add_element_pt(el_pt);

      //Set the traction function
      if(b == 4)
        el_pt->traction_fct_pt() = ventricular_pressure::constant_pressure_LV;
      else
        el_pt->traction_fct_pt() = ventricular_pressure::constant_pressure_RV;
    }
  }

}


template <class ELEMENT>
void VentricularProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{

  std::ofstream some_file;
  char filename[100];

// Number of points, which are taken to split edges, 3 is a minimal number for tetelements with nodes at edges
  unsigned npts;
  npts = 3;

  // fix_later_tag
  // we are using the same name of the file, only processor 0 writes
  // mesh is not distributed yet, OOPMPH distribute() have to be fixed
  sprintf(filename, "%s/solid_soln%i.dat", doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Ventricular_mesh_pt->output(some_file, npts);
  some_file.close();


  // fix_later_tag
  // the same problem as above
  sprintf(filename, "%s/traction%i.dat", doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  unsigned n = Solid_traction_boundary_id.size();

  for (unsigned i = 0; i < n; i++) {
    Solid_traction_mesh_pt[i]->output(some_file, npts);
  }

  some_file.close();

}


}
