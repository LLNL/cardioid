# Fiber generation for heart mesh.

The code depends on MFEM library (mfem.org). It uses algorithm from following paper:

Bayer, J. D.; Blake, R. C.; Plank, G.; Trayanova, N. A., A Novel Rule-Based Algorithm for Assigning Myocardial Fiber Orientation to Computational Heart Models. Ann. Biomed. Eng. 2012, 40, 2243-2254.

Installation Guide

1. Download and install mfem library.
    a. for serial version, you only need to install mfem.
    b. for the MPI parallel version, you need to install hypre and metis library. follow the mfem MPI installation guide.

2. Unzip Fiber code to the mfem home directory. Compile the code by use make. The code depends on mfem configuration. If mfem is serial the code will be compiled serial. If paralell, the code will be compiled parallel. 




To get Omar's rotation matrices:
    fiber -omar -m DTI060904.1.vtk -fl fiberlocs 

The output is rotmatrix.txt 

# elementnum mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33
0 0.262398 0.64768 -0.715302 -0.620382 -0.454545 -0.639152 -0.739103 0.611473 0.282537
0 0.262398 0.64768 -0.715302 -0.620382 -0.454545 -0.639152 -0.739103 0.611473 0.282537
0 0.288829 0.615381 -0.733406 -0.606783 -0.474879 -0.63742 -0.740535 0.629124 0.236244
0 0.283507 0.621307 -0.73048 -0.609931 -0.470978 -0.637309 -0.740005 0.626225 0.245429
1 -0.69841 -0.531627 -0.479162 -0.315769 -0.371927 0.872903 -0.642272 0.760949 0.0918865
...
...
 

To run the parallel version:

mpirun -np 4 fiberp -omar -m DTI060904.1.vtk -fl fiberlocs


To get Omar's rotation matrices using the elementnum provided in the fiberlocs:

mpirun -np 4 fiberp -ofast -m DTI060904.1.vtk -fl fiberlocs


