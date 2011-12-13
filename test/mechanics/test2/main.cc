//This file is modified from example driver for solid mechanics

#include "cardiac_problem.h"
#include "cardiac_elements.h"
#include "cardiac_myofilaments.h"

#include <generic.h>
#include <solid.h>


using namespace std;
using namespace oomph;


int main(int argc, char **argv)
{

  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Problem can run in parallel
  MPI_Helpers::init(argc, argv);

  // Label for output
  DocInfo doc_info;

  // Set and create output directory
  doc_info.set_directory("RESLT");


  // Set up the problem, we are using active anisotropic elelements
  VentricularProblem<ActiveAnisotropicTPVDElement<3, 3, BasicActiveModel> > problem;


  // Output initial configuration
  problem.doc_solution(doc_info);
  doc_info.number()++;




  double pressure_increment = 1.0e-1; // kPa (units are the same as the
  //scaling factor in Usyk constitutive law)

  unsigned nstep = 100; //set pressure in LV 10kPa

  // pressure in left and right ventricles
  double LV_P = 0, RV_P = 0;

  // First, dilate ventricles by applying pressure to the left and right ventricles
  for (unsigned istep = 0; istep < nstep; istep++) {
    // Solve the problem
    problem.newton_solve();

    //Output solution
    problem.doc_solution(doc_info);
    doc_info.number()++;
    oomph_info << "Time step: " << istep << std::endl;

    // increase pressure
    LV_P += pressure_increment;
    RV_P += pressure_increment / 2.0;
    problem.set_pressure(LV_P, RV_P);
  }

  nstep = 200;
  // Now, under constant pressure run in time, so active model would generate tension
  // and ventricles contracts
  for (unsigned istep = 0; istep < nstep; istep++) {
    // Solve the problem
    problem.newton_solve();

    //Output solution
    problem.doc_solution(doc_info);
    doc_info.number()++;
    std::cout << "Step: " << istep << std::endl;
    //advance in time by 1ms
    problem.advance_in_time(1.0);
  }

  // finalize MPI
  MPI_Helpers::finalize();

} // end main


