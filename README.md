# Fiber generation for heart mesh.

The code depends on MFEM library (mfem.org). It uses algorithm from following paper:

Bayer, J. D.; Blake, R. C.; Plank, G.; Trayanova, N. A., A Novel Rule-Based Algorithm for Assigning Myocardial Fiber Orientation to Computational Heart Models. Ann. Biomed. Eng. 2012, 40, 2243-2254.

TO INSTALL: 
- Follow instructions in the INSTALL directory for your particular system

TO RUN:
Step 1: Generate fibers
- on LC account: srun -N1 -n4 -ppdebug $CARDIOID/fib-gen/fiberp -ofast -m heart.vtk -fl IPlocations.txt
- on personal Mac: mpirun -np 4 $CARDIOID/fib-gen/fiberp -ofast -m heart.vtk -fl IPlocations.txt
- on personal Ubuntu: mpirun -np 4 $CARDIOID/fib-gen/fiberp -ofast -m heart.vtk -fl IPlocations.txt

Step 2: Tag surfaces
$CARDIOID/fib-gen/ply2vtk heart.ply heart.surf.vtk
$CARDIOID/fib-gen/surface -m heart.surf.vtk
