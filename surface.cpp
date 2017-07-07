/* 
 * File:   surface.cpp
 * Author: zhang30
 *
 * Created on July 6, 2017, 4:45 PM
 */

#include "mfem.hpp"
#include <cstdlib>

#include "solver.h"
#include "io.h"

using namespace std;
using namespace mfem;

/*
 * 
 */
int main(int argc, char** argv)
{
    // 1. Parse command-line options.
    const char *mesh_file = "./human.vtk";
    int order = 1;
    bool static_cond = false;
    
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
            "Surface Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
            "Finite element order (polynomial degree) or -1 for"
            " isoparametric space.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
            "--no-static-condensation", "Enable static condensation.");    
    
    args.Parse();
    if (!args.Good()) {
        args.PrintUsage(cout);
        return 1;
    }
    
    args.PrintOptions(cout);
    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
    //    the same code.
    Mesh *surface = new Mesh(mesh_file, 0, 0, false);  
    setSurf4Surf(surface);
    
    ofstream surf_ofs("surf4surf.vtk");
    printSurf4SurfVTK(surface, surf_ofs);  
    
   return 0;
}

