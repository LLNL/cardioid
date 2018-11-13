# Fiber generation for heart mesh.

## 1. Algorithm

The fiber generation algorithm mainly comes from the following paper with some modification to take care of imperfect mesh.

Bayer, J. D.; Blake, R. C.; Plank, G.; Trayanova, N. A., A Novel Rule-Based Algorithm for Assigning Myocardial Fiber Orientation to Computational Heart Models. Ann. Biomed. Eng. 2012, 40, 2243-2254.

The whole work can be divided up into the following steps:
```
    1. Find all the requisite surfaces for the algorithm
    2. Solve the laplacian problem
    3. Interpolate the laplacian solutions to some number of points
    4. Use the interpolated laplacian solutions to compute fiber orientations.
    5. Repeat 2-4 for four different Dirichlet boundary conditions 
         a. Base → 1, Apex→ 0, Epi, LV, RV → no flux
         b. Apex, Epi → 1, LV, RV→ 0, Base→ no flux
         c. LV → 1, Apex, Epi, RV→ 0, Base→ no flux
         d. RV → 1, Apex, Epi, LV→ 0, Base→ no flux
    6. Working on fiber angles with the bislerp method
    7. Calculate fiber on the nodes
    8. Start to build k-D tree for the mesh
    9. Calcualte Omar's rotation matrix
    10.Calculate anatomy files
```

## 2. Dependent Libraries

### 2.1 Serial version

The code includes the libkdtree (https://github.com/nvmd/libkdtree) for the nearest and range points search.

The serial version of the code depends only on MFEM library (http://mfem.org/). It uses algorithm from following paper:

### 2.1 Parallel version

For the parallel version, the code needs hypre (https://github.com/LLNL/hypre) and metis (http://glaros.dtc.umn.edu/gkhome/metis/metis/download) besides libkdtree and MFEM library.

## 3. Installation

To install, follow instructions in the INSTALL directory for your particular system

Note: The installation of hypre and metis for parallel version follows the the MFEM MPI installation guide (http://mfem.org/building/). Particalular versions of hypre and metis libraries are required for the parallel MFEM library.

## 4. Running the program

### 4.1 Program options:

```
[zhang30@soifon fiberp]$ ./fiber -h

Usage: ./fiber [options] ...
Options:
   -h, --help
    Print this help message and exit.
   -m <string>, --mesh <string>, current value: ./heart.vtk
    Mesh file to use.
   -o <int>, --order <int>, current value: 1
    Finite element order (polynomial degree) or -1 for isoparametric space.
   -sc, --static-condensation, -no-sc, --no-static-condensation, current option: --no-static-condensation
    Enable static condensation.
   -vis, --visualization, -no-vis, --no-visualization, current option: --visualization
    Enable or disable GLVis visualization.
   -omar, --omar_task, -no-omar, --no-omar_task, current option: --no-omar_task
    Enable or disable Omar task.
   -ofast, --omar_fast, -no-ofast, --no-omar_fast, current option: --no-omar_fast
    Enable or disable Omar fast task.
   -fl <string>, --fiblocs <string>, current value: 
    Fiber locagtion file to use.
   -vv, --verbose, -novv, --no-verbose, current option: --no-verbose
    Enable verbose output.
   -ao <double>, --aendo <double>, current value: 40
    Fiber angle alpha endo.
   -ai <double>, --aepi <double>, current value: -50
    Fiber angle alpha epi.
   -bo <double>, --bendo <double>, current value: -65
    Fiber angle beta endo.
   -bi <double>, --bepi <double>, current value: 25
    Fiber angle beta epi.
   -dd <double>, --dspacing <double>, current value: 5
    Grid spacing for ddcMD gid.
   -gl <double>, --gL <double>, current value: 0.133418
    Conductivity gL mS/mm.
   -gt <double>, --gT <double>, current value: 0.0176062
    Conductivity gT mS/mm.
   -gn <double>, --gN <double>, current value: 0.0176062
    Conductivity gN mS/mm.
   -rc <double>, --rcut <double>, current value: 1
    rangeCutoff=rcut*maxEdgeLen.

```

Most options have default values. Remember to change them with your preferable values. if you want to run omar task to generate rotation matrix please remember to specify "-fl" option with your IPlocations.txt file. The "-ofast" option is deprecated since the element number in IPlocations.txt doesn't correpsond to that in the mesh.

### 4.1 serial version

To generate fiber and anatomy files without generating omar's rotation matrix (no -omar)
```
fiber -m heart.vtk -dd 1
```

With Omar's task:
```
fiber -omar -m heart.vtk -fl IPlocations.txt -dd 1
```

### 4.2 Parallel version

To generate fiber and anatomy files without generating omar's rotation matrix (no -omar)
```
srun -N16 -n256 fiberp -m heart.vtk -dd 0.1
```

With Omar's task:
```
run -N16 -n256 fiberp -omar -m heart.vtk -fl IPlocations.txt -dd 0.1
```

### 4.3 Tag surfaces

Surface is tagged separated since by default mfem re-order the verteces in the faces. It is turned off in the surface program to keep original order of faces.

```
ply2vtk heart.ply heart.surf.vtk
surface -m heart.surf.vtk
```

## 5. Results

The fiber or fiberp (for parallel program) will generate a series of the following files:

### 5.1 surface.vtk 

surface vtk file with element tagged: 0-Apex, 1-Base, 2-EPI, 3-LV, 4-RV. The face order is re-arranged in this file.

### 5.2 vert2Elements.txt

The vertex to elements table.

### 5.3 fiber files for four different Dirichlet boundary conditions

The refined mesh and the solution for glvis:
```
phi_epi.mesh  phi_lv.mesh  phi_rv.mesh  psi_ab.mesh
phi_epi.gf  phi_lv.gf  phi_rv.gf  psi_ab.gf
``` 

Gradients of fiber in VTK format:
```
phi_epi_grads.vtk  phi_lv_grads.vtk  phi_rv_grads.vtk  psi_ab_grads.vtk
```

X axis of fiber:
```
phi_epi-x.vtk  phi_lv-x.vtk  phi_rv-x.vtk  psi_ab-x.vtk
```

The F S T vtk files:
```
fvectors.vtk  svectors.vtk  tvectors.vtk
```

### 5.4 heart.fiber 

fiber on the nodes

### 5.5 rotmatrix.txt

Omar's rotation matrix

### 5.6 Anatomy file(s)

For the serial verison only one anatomy file will generated: anatomy#000000
For the parallel version a subdirectory (named "snapshot.initial") will be created.
within the subdirectory, the number of files equal to number of MPI tasks. 
They are named:
```
anatomy#000000  anatomy#000037  anatomy#000074  anatomy#000111  anatomy#000148  anatomy#000185  anatomy#000222
anatomy#000001  anatomy#000038  anatomy#000075  anatomy#000112  anatomy#000149  anatomy#000186  anatomy#000223
anatomy#000002  anatomy#000039  anatomy#000076  anatomy#000113  anatomy#000150  anatomy#000187  anatomy#000224
anatomy#000003  anatomy#000040  anatomy#000077  anatomy#000114  anatomy#000151  anatomy#000188  anatomy#000225
anatomy#000004  anatomy#000041  anatomy#000078  anatomy#000115  anatomy#000152  anatomy#000189  anatomy#000226
anatomy#000005  anatomy#000042  anatomy#000079  anatomy#000116  anatomy#000153  anatomy#000190  anatomy#000227
...
```

### 5.7 surf4surf.vtk 

Generate by the standalone "surface". Tagged surface vtk file with face order unchanged.



## 6 Use ParaView to visualize fiber and generate the figure.
In ParaBiew:
1. Open fvectors.vtk 
2. Apply → stream tracer (Filters -> Alphabetical -> Stream Tracer)
3. StreamTracer: Seed Type=Point Source; Radius=70; Number of Points=100000; 
4. Apply → tube  (Filters -> Alphabetical -> Tube)
5. Tube: Radius=0.3; Coloring=fiber, Z (or X, Y depend on axis heart algined against)
6. Save the image by File -> Save Screenshot...

