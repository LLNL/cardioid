//#include "array3D.h"

//functions using pointers. Generally, the value at x, y, z is given. the neighboring values are then computed by indices

//as implemented in Saleheen et al IEEE TBME 1997;44(2):200 - 204 eq. 7
//in Saleheen et al, a distinction is made between intra-/interstitial-/extracellular potentials (tridomain)
//interstitial potential often refered to in literature as extracellular potential in bidomain formulation
//
// sig_ox, sig_oy, sig_oz denote the interstitial  conductivity tensor sigma in x, y, z direction
// sig_ix, sig_iy, sig_iz denote the intracellular conductivity tensor sigma in x, y, z direction
// sig_ex, sig_ey, sig_ez denote the extracellular conductivity tensor sigma in x, y, z direction
//
// Iso, Isi denote the interstitial and intracellular stimulus currents
// Iion denotes the sum of ionic currents in the electrophysiological cell model
//
//reciprocDx2, reciprocDy2, reciprocDz2 is given by 1/(dx*dx), 1/(dy*dy), 1/(dz*dz)
//which is the reciprocal of distance between two elements in x, y, z direction squared
//used here since multiplication is faster than division if not power 2

// as implemented in Saleheen et al IEEE TBME 1997;44(2):200 - 204 eq. 9
// NOTE: in our implementation the current term due to Vm, ie. FD3 is computed beforehand and added to the stimulus current Istim!!!
//double FDLaplacianHomoTissueSaleheen97(int*** &tissue, double*** &Vm, double*** &interstitialPot, sigmatensorMatrix*** &sigmaMintra, sigmatensorMatrix*** &sigmaMextra, double Istim, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal)
double FDLaplacianHomoTissueSaleheen97(int*** &tissue, double*** &interstitialPot, sigmatensorMatrix*** &sigmaMintra, sigmatensorMatrix*** &sigmaMextra, double Istim, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal)
{
  double result;                                // fct argument OR return value?
  double dxdxReciprocal, dydyReciprocal, dzdzReciprocal;
  double potXminus1, potXplus1;
  double potYminus1, potYplus1;
  double potZminus1, potZplus1;
          
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dzdzReciprocal = dzReciprocal * dzReciprocal;

  // eq. 11
  double FD2 = ((sigmaMintra[x][y][z].a11 + sigmaMextra[x][y][z].a11) * (interstitialPot[x+1][y][z] + interstitialPot[x-1][y][z]) * dxdxReciprocal)
             + ((sigmaMintra[x][y][z].a22 + sigmaMextra[x][y][z].a22) * (interstitialPot[x][y+1][z] + interstitialPot[x][y-1][z]) * dydyReciprocal)
             + ((sigmaMintra[x][y][z].a33 + sigmaMextra[x][y][z].a33) * (interstitialPot[x][y][z+1] + interstitialPot[x][y][z-1]) * dzdzReciprocal);
                                                                        
                                                                        
/********************************************************************************************************
 * the computation of FD3 is not needed since Istim already is 
 * Istim = FD3 + (Isi + Iso)                                                                        
 *
  double twoVm = (Vm[x][y][z])*2;
  // Neumann boundary conditions
  // if neighboring tissue is non-excitable, then use the value of center element
  // otherwise use value of element
  // interstitial Potential is equivalent extracellular Potential in non-excitable tissue
  
  potXminus1 = ((tissue[x-1][y][z]!=0)? Vm[x-1][y][z]:Vm[x][y][z]);
  potXplus1  = ((tissue[x+1][y][z]!=0)? Vm[x+1][y][z]:Vm[x][y][z]);
  potYminus1 = ((tissue[x][y-1][z]!=0)? Vm[x][y-1][z]:Vm[x][y][z]);
  potYplus1  = ((tissue[x][y+1][z]!=0)? Vm[x][y+1][z]:Vm[x][y][z]);
  potZminus1 = ((tissue[x][y][z-1]!=0)? Vm[x][y][z-1]:Vm[x][y][z]);
  potZplus1  = ((tissue[x][y][z+1]!=0)? Vm[x][y][z+1]:Vm[x][y][z]);
                                                                                                                
  // eq. 12
  double FD3 = (sigmaMintra[x][y][z].a11 * (potXplus1 + potXminus1 - twoVm) * dxdxReciprocal)
             + (sigmaMintra[x][y][z].a22 * (potYplus1 + potYminus1 - twoVm) * dydyReciprocal)
             + (sigmaMintra[x][y][z].a33 * (potZplus1 + potZminus1 - twoVm) * dzdzReciprocal);
 ********************************************************************************************************/
                                       
  // Istim = (Isi + Iso) which is the stimulus current
  //result = (FD2 + FD3 + (Istim))
  result = (FD2 + (Istim)) // Istim = FD3 + (Isi + Iso) 
         / (2 * (((sigmaMintra[x][y][z].a11 + sigmaMextra[x][y][z].a11) * dxdxReciprocal)
         +       ((sigmaMintra[x][y][z].a22 + sigmaMextra[x][y][z].a22) * dydyReciprocal)
         +       ((sigmaMintra[x][y][z].a33 + sigmaMextra[x][y][z].a33) * dzdzReciprocal)));
                                              
  return (result);
}

// as implemented in Saleheen et al IEEE TBME 1997;44(2):200 - 204 eq. 10
// arguments for sig_ex, sig_ey and sig_ez are missing
double FDLaplacianHomoBathSaleheen97(int*** &tissue, double*** &phi, double sigmaX, double sigmaY, double sigmaZ, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal)
{
  double result;                                // fct argument OR return value?
  double dxdxReciprocal, dydyReciprocal, dzdzReciprocal;
  double potXminus1, potXplus1;
  double potYminus1, potYplus1;
  double potZminus1, potZplus1;
  
  
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dzdzReciprocal = dzReciprocal * dzReciprocal;
      
  potXminus1 = ((tissue[x-1][y][z]!=0)? phi[x-1][y][z]:phi[x][y][z]);
  potXplus1  = ((tissue[x+1][y][z]!=0)? phi[x+1][y][z]:phi[x][y][z]);
  potYminus1 = ((tissue[x][y-1][z]!=0)? phi[x][y-1][z]:phi[x][y][z]);
  potYplus1  = ((tissue[x][y+1][z]!=0)? phi[x][y+1][z]:phi[x][y][z]);
  potZminus1 = ((tissue[x][y][z-1]!=0)? phi[x][y][z-1]:phi[x][y][z]);
  potZplus1  = ((tissue[x][y][z+1]!=0)? phi[x][y][z+1]:phi[x][y][z]);
                  
              
  // eq. 13
  double FD4 = (sigmaX * (potXplus1 + potXminus1) * dxdxReciprocal)
             + (sigmaY * (potYplus1 + potYminus1) * dydyReciprocal)
             + (sigmaZ * (potZplus1 + potZminus1) * dzdzReciprocal);
                                                    
  result = FD4/(2*((sigmaX * dxdxReciprocal)+(sigmaY * dydyReciprocal)+(sigmaZ * dzdzReciprocal)));
  return (result);
}
                                                              

// as implemented in Saleheen et al IEEE TBME 1997;44(2):200 - 204 eq. 8
double FDLaplacianHomoAnisotrop(int*** &tissue, double*** &phi, sigmatensorMatrix*** &sigmaM, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal)
{
  double result;                                // fct argument OR return value?
  double dxdxReciprocal, dydyReciprocal, dzdzReciprocal;
    
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dzdzReciprocal = dzReciprocal * dzReciprocal;

  double sigmaX = sigmaM[x][y][z].a11;
  double sigmaY = sigmaM[x][y][z].a22;
  double sigmaZ = sigmaM[x][y][z].a33;

  double phitimes2 = (2 * (phi[x][y][z]));
  double potXminus1, potXplus1;
  double potYminus1, potYplus1;
  double potZminus1, potZplus1;
          
  // Neumann boundary conditions
  // if neighboring tissue is non-excitable, then use the value of center element
  // otherwise use value of element
                

  potXminus1 = ((tissue[x-1][y][z]!=0)? phi[x-1][y][z]:phi[x][y][z]);
  potXplus1  = ((tissue[x+1][y][z]!=0)? phi[x+1][y][z]:phi[x][y][z]);
  potYminus1 = ((tissue[x][y-1][z]!=0)? phi[x][y-1][z]:phi[x][y][z]);
  potYplus1  = ((tissue[x][y+1][z]!=0)? phi[x][y+1][z]:phi[x][y][z]);
  potZminus1 = ((tissue[x][y][z-1]!=0)? phi[x][y][z-1]:phi[x][y][z]);
  potZplus1  = ((tissue[x][y][z+1]!=0)? phi[x][y][z+1]:phi[x][y][z]);

  result = (sigmaX * ((potXplus1 + potXminus1 -  phitimes2) * dxdxReciprocal))
         + (sigmaY * ((potYplus1 + potYminus1 -  phitimes2) * dydyReciprocal))
         + (sigmaZ * ((potZplus1 + potZminus1 -  phitimes2) * dzdzReciprocal));
                                                               

  return (result);
}

int boundaryFDLaplacian2DConstants(int*** &tissue, diffusion2D*** &diffConst, sigmatensorMatrix*** &sigmaMatrix, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal)
{
  int i, j;
  sigmatensorMatrix conductivityMatrix[3][3];
  double phi[3][3];
      
  double A1, A2, A3, A4, A5, A6, A7;
        
  double result;                                // fct argument OR return value?
  double dxdxReciprocal, dydyReciprocal, dxdyReciprocal;
            
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dxdyReciprocal = dxReciprocal * dyReciprocal;
                  
  for(i = 0; i < 3; i++)
       for(j = 0; j < 3; j++)
          {
            conductivityMatrix[i][j].a11 = ((tissue[x-1+i][y-1+j][z]!=0)? (sigmaMatrix[x-1+i][y-1+j][z].a11):0.0);
            conductivityMatrix[i][j].a12 = ((tissue[x-1+i][y-1+j][z]!=0)? (sigmaMatrix[x-1+i][y-1+j][z].a12):0.0);
            conductivityMatrix[i][j].a13 = ((tissue[x-1+i][y-1+j][z]!=0)? (sigmaMatrix[x-1+i][y-1+j][z].a13):0.0);
            conductivityMatrix[i][j].a22 = ((tissue[x-1+i][y-1+j][z]!=0)? (sigmaMatrix[x-1+i][y-1+j][z].a22):0.0);
          }
  i = j = 1;
                                                                                                                  
  // after Pullan et al. page 190
  // note that sigma.a12 = sigma.a21. Thus, we have a factor 2 for A3
  diffConst[x][y][z].A1 = conductivityMatrix[i][j].a11 * dxdxReciprocal;
  diffConst[x][y][z].A2 = conductivityMatrix[i][j].a22 * dydyReciprocal;
  diffConst[x][y][z].A3 = 2 * conductivityMatrix[i][j].a12 * 0.25 * dxdyReciprocal;
  diffConst[x][y][z].A4 = (conductivityMatrix[i+1][j].a11 - conductivityMatrix[i-1][j].a11) * 0.25 * dxdxReciprocal;
  diffConst[x][y][z].A5 = (conductivityMatrix[i+1][j].a12 - conductivityMatrix[i-1][j].a12) * 0.25 * dxdyReciprocal;
  diffConst[x][y][z].A6 = (conductivityMatrix[i][j+1].a12 - conductivityMatrix[i][j-1].a12) * 0.25 * dxdyReciprocal;
  diffConst[x][y][z].A7 = (conductivityMatrix[i][j+1].a22 - conductivityMatrix[i][j-1].a22) * 0.25 * dxdyReciprocal;

  return (0);                                                                                                                                                                                                                                                                      
}

double boundaryFDLaplacian2DSumPhi(double*** &phi, diffusion2D*** &diffConst, int &x, int &y, int &z)
{
  double result = diffConst[x][y][z].A1 * ( phi[x-1][y][z] - (2 * phi[x][y][z]) + phi[x+1][y][z] )
                + diffConst[x][y][z].A2 * ( phi[x][y-1][z] - (2 * phi[x][y][z]) + phi[x][y+1][z] )
                + diffConst[x][y][z].A3 * ( phi[x-1][y-1][z] - phi[x-1][y+1][z] - phi[x+1][y-1][z] + phi[x+1][y+1][z] )
                + diffConst[x][y][z].A4 * ( phi[x+1][y][z] - phi[x-1][y][z] )
                + diffConst[x][y][z].A5 * ( phi[x][y+1][z] - phi[x][y-1][z] )
                + diffConst[x][y][z].A6 * ( phi[x+1][y][z] - phi[x-1][y][z] )
                + diffConst[x][y][z].A7 * ( phi[x][y+1][z] - phi[x][y-1][z] );
  
  return(result);
}


double boundaryFDLaplacian2D(int*** &tissue, double*** &V, sigmatensorMatrix*** &sigmaMatrix, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal)
{
  int i, j;
  sigmatensorMatrix conductivityMatrix[3][3];
  double phi[3][3];
      
  double A1, A2, A3, A4, A5, A6, A7;
          
  double result;                                // fct argument OR return value?
  double dxdxReciprocal, dydyReciprocal, dxdyReciprocal;
                
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dxdyReciprocal = dxReciprocal * dyReciprocal;
    
  for(i = 0; i < 3; i++)
     for(j = 0; j < 3; j++)
           {
             conductivityMatrix[i][j].a11 = ((tissue[x-1+i][y-1+j][z]!=0)? (sigmaMatrix[x-1+i][y-1+j][z].a11):0.0);
             conductivityMatrix[i][j].a12 = ((tissue[x-1+i][y-1+j][z]!=0)? (sigmaMatrix[x-1+i][y-1+j][z].a12):0.0);
             conductivityMatrix[i][j].a13 = ((tissue[x-1+i][y-1+j][z]!=0)? (sigmaMatrix[x-1+i][y-1+j][z].a13):0.0);
             conductivityMatrix[i][j].a22 = ((tissue[x-1+i][y-1+j][z]!=0)? (sigmaMatrix[x-1+i][y-1+j][z].a22):0.0);
                                                                              
             phi[i][j] = (V[x-1+i][y-1+j][z]);
           }
                                                                                                      
  i = j = 1;


  // after Pullan et al. page 190
  // note that sigma.a12 = sigma.a21. Thus, we have a factor 2 for A3
  A1 = conductivityMatrix[i][j].a11 * dxdxReciprocal;
  A2 = conductivityMatrix[i][j].a22 * dydyReciprocal;
  A3 = 2 * conductivityMatrix[i][j].a12 * 0.25 * dxdyReciprocal;
  A4 = (conductivityMatrix[i+1][j].a11 - conductivityMatrix[i-1][j].a11) * 0.25 * dxdxReciprocal;
  A5 = (conductivityMatrix[i+1][j].a12 - conductivityMatrix[i-1][j].a12) * 0.25 * dxdyReciprocal;
  A6 = (conductivityMatrix[i][j+1].a12 - conductivityMatrix[i][j-1].a12) * 0.25 * dxdyReciprocal;
  A7 = (conductivityMatrix[i][j+1].a22 - conductivityMatrix[i][j-1].a22) * 0.25 * dxdyReciprocal;
        
  result  = A1 * ( phi[i-1][j] - (2 * phi[i][j]) + phi[i+1][j] )
          + A2 * ( phi[i][j-1] - (2 * phi[i][j]) + phi[i][j+1] )
          + A3 * ( phi[i-1][j-1] - phi[i-1][j+1] - phi[i+1][j-1] + phi[i+1][j+1] )
          + A4 * ( phi[i+1][j] - phi[i-1][j] )
          + A5 * ( phi[i][j+1] - phi[i][j-1] )
          + A6 * ( phi[i+1][j] - phi[i-1][j] )
          + A7 * ( phi[i][j+1] - phi[i][j-1] );
                    
  return(result);

}

int boundaryFDLaplacianSaleheen98Constants(int*** &tissue, diffusion*** &diffConst, sigmatensorMatrix*** &sigmaMatrix, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal)
{
  int i, j, k;
  sigmatensorMatrix conductivityMatrix[3][3][3];
  double sigmaX, sigmaY, sigmaZ;
  double dxdxReciprocal, dydyReciprocal, dzdzReciprocal;
  
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dzdzReciprocal = dzReciprocal * dzReciprocal;
                      
  for(i = 0; i < 3; i++)
     for(j = 0; j < 3; j++)
        for(k = 0; k < 3; k++)
           {
             conductivityMatrix[i][j][k].a11 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a11):0.0);
             conductivityMatrix[i][j][k].a12 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a12):0.0);
             conductivityMatrix[i][j][k].a13 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a13):0.0);
             conductivityMatrix[i][j][k].a22 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a22):0.0);
             conductivityMatrix[i][j][k].a23 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a23):0.0);
             conductivityMatrix[i][j][k].a33 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a33):0.0);
           }
                                                                                                                  
  i = j = k = 1;
  // original formulation - note that I replaced division by multiplication of reciprocal and division by 2 with multiplication by 0.5
  // I also used the central difference method for the differentiation. This results in a factor of 0.25 below
  // NOTE: sigmaX, sigmaY, sigmaZ could be precalculated if no deformation and no chance in conductivity
                                                                                                                    
  sigmaX  = 0.25 * dxReciprocal * (((conductivityMatrix[i+1][j][k].a11 - conductivityMatrix[i-1][j][k].a11) * dxReciprocal)
          + ((conductivityMatrix[i][j+1][k].a12 - conductivityMatrix[i][j-1][k].a12) * dyReciprocal)
          + ((conductivityMatrix[i][j][k+1].a13 - conductivityMatrix[i][j][k-1].a13) * dzReciprocal));
  sigmaY  = 0.25 * dyReciprocal * (((conductivityMatrix[i+1][j][k].a12 - conductivityMatrix[i-1][j][k].a12) * dxReciprocal)
          + ((conductivityMatrix[i][j+1][k].a22 - conductivityMatrix[i][j-1][k].a12) * dyReciprocal)
          + ((conductivityMatrix[i][j][k+1].a23 - conductivityMatrix[i][j][k-1].a23) * dzReciprocal));
  sigmaZ  = 0.25 * dzReciprocal * (((conductivityMatrix[i+1][j][k].a13 - conductivityMatrix[i-1][j][k].a13) * dxReciprocal)
          + ((conductivityMatrix[i][j+1][k].a23 - conductivityMatrix[i][j-1][k].a23) * dyReciprocal)
          + ((conductivityMatrix[i][j][k+1].a33 - conductivityMatrix[i][j][k-1].a33) * dzReciprocal));
  
  diffConst[x][y][z].A1      = conductivityMatrix[i+1][j][k].a11 * (dxdxReciprocal - (sigmaX/conductivityMatrix[i][j][k].a11));
  diffConst[x][y][z].A3      = conductivityMatrix[i-1][j][k].a11 * (dxdxReciprocal + (sigmaX/conductivityMatrix[i][j][k].a11));
  diffConst[x][y][z].A2      = conductivityMatrix[i][j+1][k].a22 * (dydyReciprocal - (sigmaY/conductivityMatrix[i][j][k].a22));
  diffConst[x][y][z].A4      = conductivityMatrix[i][j-1][k].a22 * (dydyReciprocal + (sigmaY/conductivityMatrix[i][j][k].a22));
  diffConst[x][y][z].A5 = diffConst[x][y][z].A6 = diffConst[x][y][z].A7 = diffConst[x][y][z].A8 = (0.5 * dxReciprocal * dyReciprocal);
  
  diffConst[x][y][z].A5 = conductivityMatrix[i+1][j+1][k].a12 * diffConst[x][y][z].A5;
  diffConst[x][y][z].A6 = -(conductivityMatrix[i-1][j+1][k].a12 * diffConst[x][y][z].A6);
  diffConst[x][y][z].A7 = conductivityMatrix[i-1][j-1][k].a12 * diffConst[x][y][z].A7;
  diffConst[x][y][z].A8 = -(conductivityMatrix[i+1][j-1][k].a12 * diffConst[x][y][z].A8);
  diffConst[x][y][z].A9      = conductivityMatrix[i][j][k+1].a33 * (dzdzReciprocal - (sigmaZ/conductivityMatrix[i][j][k].a33));
  diffConst[x][y][z].A10     = conductivityMatrix[i][j][k-1].a33 * (dzdzReciprocal + (sigmaZ/conductivityMatrix[i][j][k].a33));
  diffConst[x][y][z].A11 = diffConst[x][y][z].A12 = diffConst[x][y][z].A13 = diffConst[x][y][z].A14 = (0.5 * dyReciprocal * dzReciprocal);
  
  diffConst[x][y][z].A11 = conductivityMatrix[i][j+1][k+1].a23 * diffConst[x][y][z].A11;
  diffConst[x][y][z].A12 = -(conductivityMatrix[i][j+1][k-1].a23 * diffConst[x][y][z].A12);
  diffConst[x][y][z].A13 = conductivityMatrix[i][j-1][k-1].a23 * diffConst[x][y][z].A13;
  diffConst[x][y][z].A14 = -(conductivityMatrix[i][j-1][k+1].a23 * diffConst[x][y][z].A14);
  diffConst[x][y][z].A15 = diffConst[x][y][z].A16 = diffConst[x][y][z].A17 = diffConst[x][y][z].A18 = (0.5 * dxReciprocal * dzReciprocal);
  
  diffConst[x][y][z].A15 = conductivityMatrix[i+1][j][k+1].a13 * diffConst[x][y][z].A15;
  diffConst[x][y][z].A16 = -(conductivityMatrix[i-1][j][k+1].a13 * diffConst[x][y][z].A16);
  diffConst[x][y][z].A17 = conductivityMatrix[i-1][j][k-1].a13 * diffConst[x][y][z].A17;
  diffConst[x][y][z].A18 = -(conductivityMatrix[i+1][j][k-1].a13 * diffConst[x][y][z].A18);

  diffConst[x][y][z].sumA    = diffConst[x][y][z].A1 
                             + diffConst[x][y][z].A2 
                             + diffConst[x][y][z].A3 
                             + diffConst[x][y][z].A4 
                             + diffConst[x][y][z].A5 
                             + diffConst[x][y][z].A6 
                             + diffConst[x][y][z].A7 
                             + diffConst[x][y][z].A8 
                             + diffConst[x][y][z].A9 
                             + diffConst[x][y][z].A10 
                             + diffConst[x][y][z].A11 
                             + diffConst[x][y][z].A12 
                             + diffConst[x][y][z].A13 
                             + diffConst[x][y][z].A14 
                             + diffConst[x][y][z].A15 
                             + diffConst[x][y][z].A16 
                             + diffConst[x][y][z].A17 
                             + diffConst[x][y][z].A18;
  
  
  return(0);
}

double boundaryFDLaplacianSaleheen98SumPhi(double*** &phi, diffusion*** &diffConst, int &x, int &y, int &z)
{

  double SumAphi = (diffConst[x][y][z].A1  * (phi[x+1][y][z]))
                 + (diffConst[x][y][z].A2  * (phi[x][y+1][z]))
                 + (diffConst[x][y][z].A3  * (phi[x-1][y][z]))
                 + (diffConst[x][y][z].A4  * (phi[x][y-1][z]))
                 + (diffConst[x][y][z].A5  * (phi[x+1][y+1][z]))
                 + (diffConst[x][y][z].A6  * (phi[x-1][y+1][z]))
                 + (diffConst[x][y][z].A7  * (phi[x-1][y-1][z]))
                 + (diffConst[x][y][z].A8  * (phi[x+1][y-1][z]))
                 + (diffConst[x][y][z].A9  * (phi[x][y][z+1]))
                 + (diffConst[x][y][z].A10 * (phi[x][y][z-1]))
                 + (diffConst[x][y][z].A11 * (phi[x][y+1][z+1]))
                 + (diffConst[x][y][z].A12 * (phi[x][y+1][z-1]))
                 + (diffConst[x][y][z].A13 * (phi[x][y-1][z-1]))
                 + (diffConst[x][y][z].A14 * (phi[x][y-1][z+1]))
                 + (diffConst[x][y][z].A15 * (phi[x+1][y][z+1]))
                 + (diffConst[x][y][z].A16 * (phi[x-1][y][z+1]))
                 + (diffConst[x][y][z].A17 * (phi[x-1][y][z-1]))
                 + (diffConst[x][y][z].A18 * (phi[x+1][y][z-1]));
  
  double result = SumAphi - (diffConst[x][y][z].sumA * (phi[x][y][z]));
  return(result);
}


  
double boundaryFDLaplacianSaleheen98(int*** &tissue, double*** &V, sigmatensorMatrix*** &sigmaMatrix, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal)
{
  int i, j, k;
  sigmatensorMatrix conductivityMatrix[3][3][3];
  double phi[3][3][3];
      
  double sigmaX, sigmaY, sigmaZ;
  double A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18;
            
  double SumA, SumAphi;
  double result;                                // fct argument OR return value?
  double dxdxReciprocal, dydyReciprocal, dzdzReciprocal;
                
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dzdzReciprocal = dzReciprocal * dzReciprocal;
                      
  for(i = 0; i < 3; i++)
     for(j = 0; j < 3; j++)
        for(k = 0; k < 3; k++)
           {
/*
             conductivityMatrix[i][j][k].a11 = (sigmaMatrix[x-1+i][y-1+j][z-1+k].a11);
             conductivityMatrix[i][j][k].a12 = (sigmaMatrix[x-1+i][y-1+j][z-1+k].a12);
             conductivityMatrix[i][j][k].a13 = (sigmaMatrix[x-1+i][y-1+j][z-1+k].a13);
             conductivityMatrix[i][j][k].a22 = (sigmaMatrix[x-1+i][y-1+j][z-1+k].a22);
             conductivityMatrix[i][j][k].a23 = (sigmaMatrix[x-1+i][y-1+j][z-1+k].a23);
             conductivityMatrix[i][j][k].a33 = (sigmaMatrix[x-1+i][y-1+j][z-1+k].a33);
*/                                                                                         
           
             //((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a11):0.0);             
             conductivityMatrix[i][j][k].a11 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a11):0.0);
             conductivityMatrix[i][j][k].a12 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a12):0.0);
             conductivityMatrix[i][j][k].a13 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a13):0.0);
             conductivityMatrix[i][j][k].a22 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a22):0.0);
             conductivityMatrix[i][j][k].a23 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a23):0.0);
             conductivityMatrix[i][j][k].a33 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? (sigmaMatrix[x-1+i][y-1+j][z-1+k].a33):0.0);
             
             phi[i][j][k] = (V[x-1+i][y-1+j][z-1+k]);
           }
                                                                                                                                                      
  i = j = k = 1;

// original formulation - note that I replaced division by multiplication of reciprocal and division by 2 with multiplication by 0.5
// I also used the central difference method for the differentiation. This results in a factor of 0.25 below
// NOTE: sigmaX, sigmaY, sigmaZ could be precalculated if no deformation and no chance in conductivity


  sigmaX  = 0.25 * dxReciprocal * (((conductivityMatrix[i+1][j][k].a11 - conductivityMatrix[i-1][j][k].a11) * dxReciprocal)
                                 + ((conductivityMatrix[i][j+1][k].a12 - conductivityMatrix[i][j-1][k].a12) * dyReciprocal)
                                 + ((conductivityMatrix[i][j][k+1].a13 - conductivityMatrix[i][j][k-1].a13) * dzReciprocal));
  sigmaY  = 0.25 * dyReciprocal * (((conductivityMatrix[i+1][j][k].a12 - conductivityMatrix[i-1][j][k].a12) * dxReciprocal)
                                 + ((conductivityMatrix[i][j+1][k].a22 - conductivityMatrix[i][j-1][k].a12) * dyReciprocal)
                                 + ((conductivityMatrix[i][j][k+1].a23 - conductivityMatrix[i][j][k-1].a23) * dzReciprocal));
  sigmaZ  = 0.25 * dzReciprocal * (((conductivityMatrix[i+1][j][k].a13 - conductivityMatrix[i-1][j][k].a13) * dxReciprocal)
                                 + ((conductivityMatrix[i][j+1][k].a23 - conductivityMatrix[i][j-1][k].a23) * dyReciprocal)
                                 + ((conductivityMatrix[i][j][k+1].a33 - conductivityMatrix[i][j][k-1].a33) * dzReciprocal));
  
  A1      = conductivityMatrix[i+1][j][k].a11 * (dxdxReciprocal - (sigmaX/conductivityMatrix[i][j][k].a11));
  A3      = conductivityMatrix[i-1][j][k].a11 * (dxdxReciprocal + (sigmaX/conductivityMatrix[i][j][k].a11));
  
  A2      = conductivityMatrix[i][j+1][k].a22 * (dydyReciprocal - (sigmaY/conductivityMatrix[i][j][k].a22));
  A4      = conductivityMatrix[i][j-1][k].a22 * (dydyReciprocal + (sigmaY/conductivityMatrix[i][j][k].a22));
  
  A5 = A6 = A7 = A8 = (0.5 * dxReciprocal * dyReciprocal);
  A5 = conductivityMatrix[i+1][j+1][k].a12 * A5;
  A6 = -(conductivityMatrix[i-1][j+1][k].a12 * A6);
  A7 = conductivityMatrix[i-1][j-1][k].a12 * A7;
  A8 = -(conductivityMatrix[i+1][j-1][k].a12 * A8);
  
  A9      = conductivityMatrix[i][j][k+1].a33 * (dzdzReciprocal - (sigmaZ/conductivityMatrix[i][j][k].a33));
  A10     = conductivityMatrix[i][j][k-1].a33 * (dzdzReciprocal + (sigmaZ/conductivityMatrix[i][j][k].a33));
  
  A11 = A12 = A13 = A14 = (0.5 * dyReciprocal * dzReciprocal);
  A11 = conductivityMatrix[i][j+1][k+1].a23 * A11;
  A12 = -(conductivityMatrix[i][j+1][k-1].a23 * A12);
  A13 = conductivityMatrix[i][j-1][k-1].a23 * A13;
  A14 = -(conductivityMatrix[i][j-1][k+1].a23 * A14);
  
  A15 = A16 = A17 = A18 = (0.5 * dxReciprocal * dzReciprocal);
  A15 = conductivityMatrix[i+1][j][k+1].a13 * A15;
  A16 = -(conductivityMatrix[i-1][j][k+1].a13 * A16);
  A17 = conductivityMatrix[i-1][j][k-1].a13 * A17;
  A18 = -(conductivityMatrix[i+1][j][k-1].a13 * A18);
  
  SumA    = A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10 + A11 + A12 + A13 + A14 + A15 + A16 + A17 + A18;
  SumAphi = (A1  * (phi[i+1][j][k]))
          + (A2  * (phi[i][j+1][k]))
          + (A3  * (phi[i-1][j][k]))
          + (A4  * (phi[i][j-1][k]))
          + (A5  * (phi[i+1][j+1][k]))
          + (A6  * (phi[i-1][j+1][k]))
          + (A7  * (phi[i-1][j-1][k]))
          + (A8  * (phi[i+1][j-1][k]))
          + (A9  * (phi[i][j][k+1]))
          + (A10 * (phi[i][j][k-1]))
          + (A11 * (phi[i][j+1][k+1]))
          + (A12 * (phi[i][j+1][k-1]))
          + (A13 * (phi[i][j-1][k-1]))
          + (A14 * (phi[i][j-1][k+1]))
          + (A15 * (phi[i+1][j][k+1]))
          + (A16 * (phi[i-1][j][k+1]))
          + (A17 * (phi[i-1][j][k-1]))
          + (A18 * (phi[i+1][j][k-1]));
                                                                                                                                                                            
  result = SumAphi - (SumA * (phi[i][j][k]));
  return(result);
}


double FDLaplacianHeteroTissuePhiExtraSaleheen98Constants(int*** &tissue, diffusion*** &diffConst, sigmatensorMatrix*** &sigmaMintra, sigmatensorMatrix*** &sigmaMextra,  int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal)
{
  int i, j, k;
  sigmatensorMatrix conductivityMatrix[3][3][3];
  double sigmaX, sigmaY, sigmaZ;
  double dxdxReciprocal, dydyReciprocal, dzdzReciprocal;
                
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dzdzReciprocal = dzReciprocal * dzReciprocal;
                      
  for(i = 0; i < 3; i++)
     for(j = 0; j < 3; j++)
        for(k = 0; k < 3; k++)
           {
             conductivityMatrix[i][j][k].a11 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a11) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a11)):0.0);
             conductivityMatrix[i][j][k].a12 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a12) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a12)):0.0);
             conductivityMatrix[i][j][k].a13 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a13) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a13)):0.0);
             conductivityMatrix[i][j][k].a22 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a22) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a22)):0.0);
             conductivityMatrix[i][j][k].a23 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a23) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a23)):0.0);
             conductivityMatrix[i][j][k].a33 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a33) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a33)):0.0);
           }
  
  i = j = k = 1;
  
  // original formulation - note that I replaced division by multiplication of reciprocal and division by 2 with multiplication by 0.5
  // I also used the central difference method for the differentiation. This results in a factor of 0.25 below
  // NOTE: sigmaX, sigmaY, sigmaZ could be precalculated if no deformation and no chance in conductivity
  
  sigmaX  = 0.25 * dxReciprocal * (((conductivityMatrix[i+1][j][k].a11 - conductivityMatrix[i-1][j][k].a11) * dxReciprocal)
          + ((conductivityMatrix[i][j+1][k].a12 - conductivityMatrix[i][j-1][k].a12) * dyReciprocal)
          + ((conductivityMatrix[i][j][k+1].a13 - conductivityMatrix[i][j][k-1].a13) * dzReciprocal));
  sigmaY  = 0.25 * dyReciprocal * (((conductivityMatrix[i+1][j][k].a12 - conductivityMatrix[i-1][j][k].a12) * dxReciprocal)
          + ((conductivityMatrix[i][j+1][k].a22 - conductivityMatrix[i][j-1][k].a12) * dyReciprocal)
          + ((conductivityMatrix[i][j][k+1].a23 - conductivityMatrix[i][j][k-1].a23) * dzReciprocal));
  sigmaZ  = 0.25 * dzReciprocal * (((conductivityMatrix[i+1][j][k].a13 - conductivityMatrix[i-1][j][k].a13) * dxReciprocal)
          + ((conductivityMatrix[i][j+1][k].a23 - conductivityMatrix[i][j-1][k].a23) * dyReciprocal)
          + ((conductivityMatrix[i][j][k+1].a33 - conductivityMatrix[i][j][k-1].a33) * dzReciprocal));
          

  diffConst[x][y][z].A1      = conductivityMatrix[i+1][j][k].a11 * (dxdxReciprocal - (sigmaX/conductivityMatrix[i][j][k].a11));
  diffConst[x][y][z].A3      = conductivityMatrix[i-1][j][k].a11 * (dxdxReciprocal + (sigmaX/conductivityMatrix[i][j][k].a11));
  diffConst[x][y][z].A2      = conductivityMatrix[i][j+1][k].a22 * (dydyReciprocal - (sigmaY/conductivityMatrix[i][j][k].a22));
  diffConst[x][y][z].A4      = conductivityMatrix[i][j-1][k].a22 * (dydyReciprocal + (sigmaY/conductivityMatrix[i][j][k].a22));
  
  diffConst[x][y][z].A5 = diffConst[x][y][z].A6 = diffConst[x][y][z].A7 = diffConst[x][y][z].A8 = (0.5 * dxReciprocal * dyReciprocal);
  diffConst[x][y][z].A5 = conductivityMatrix[i+1][j+1][k].a12 * diffConst[x][y][z].A5;
  diffConst[x][y][z].A6 = -(conductivityMatrix[i-1][j+1][k].a12 * diffConst[x][y][z].A6);
  diffConst[x][y][z].A7 = conductivityMatrix[i-1][j-1][k].a12 * diffConst[x][y][z].A7;
  diffConst[x][y][z].A8 = -(conductivityMatrix[i+1][j-1][k].a12 * diffConst[x][y][z].A8);
  
  diffConst[x][y][z].A9      = conductivityMatrix[i][j][k+1].a33 * (dzdzReciprocal - (sigmaZ/conductivityMatrix[i][j][k].a33));
  diffConst[x][y][z].A10     = conductivityMatrix[i][j][k-1].a33 * (dzdzReciprocal + (sigmaZ/conductivityMatrix[i][j][k].a33));
                      
  diffConst[x][y][z].A11 = diffConst[x][y][z].A12 = diffConst[x][y][z].A13 = diffConst[x][y][z].A14 = (0.5 * dyReciprocal * dzReciprocal);
  diffConst[x][y][z].A11 = conductivityMatrix[i][j+1][k+1].a23 * diffConst[x][y][z].A11;
  diffConst[x][y][z].A12 = -(conductivityMatrix[i][j+1][k-1].a23 * diffConst[x][y][z].A12);
  diffConst[x][y][z].A13 = conductivityMatrix[i][j-1][k-1].a23 * diffConst[x][y][z].A13;
  diffConst[x][y][z].A14 = -(conductivityMatrix[i][j-1][k+1].a23 * diffConst[x][y][z].A14);
                                
  diffConst[x][y][z].A15 = diffConst[x][y][z].A16 = diffConst[x][y][z].A17 = diffConst[x][y][z].A18 = (0.5 * dxReciprocal * dzReciprocal);
  diffConst[x][y][z].A15 = conductivityMatrix[i+1][j][k+1].a13 * diffConst[x][y][z].A15;
  diffConst[x][y][z].A16 = -(conductivityMatrix[i-1][j][k+1].a13 * diffConst[x][y][z].A16);
  diffConst[x][y][z].A17 = conductivityMatrix[i-1][j][k-1].a13 * diffConst[x][y][z].A17;
  diffConst[x][y][z].A18 = -(conductivityMatrix[i+1][j][k-1].a13 * diffConst[x][y][z].A18);
                                          
  diffConst[x][y][z].sumA    = diffConst[x][y][z].A1 
                             + diffConst[x][y][z].A2 
                             + diffConst[x][y][z].A3 
                             + diffConst[x][y][z].A4 
                             + diffConst[x][y][z].A5 
                             + diffConst[x][y][z].A6 
                             + diffConst[x][y][z].A7 
                             + diffConst[x][y][z].A8 
                             + diffConst[x][y][z].A9 
                             + diffConst[x][y][z].A10 
                             + diffConst[x][y][z].A11 
                             + diffConst[x][y][z].A12 
                             + diffConst[x][y][z].A13 
                             + diffConst[x][y][z].A14 
                             + diffConst[x][y][z].A15 
                             + diffConst[x][y][z].A16 
                             + diffConst[x][y][z].A17 
                             + diffConst[x][y][z].A18;
                             
  return(0);                             
}

double FDLaplacianHeteroTissuePhiExtraSaleheen98SumPhiExtra(double*** &phi, diffusion*** diffConst,  double*** &Istim, int &x, int &y, int &z)
{

  double SumAphi = (diffConst[x][y][z].A1  * (phi[x+1][y][z]))
                 + (diffConst[x][y][z].A2  * (phi[x][y+1][z]))
                 + (diffConst[x][y][z].A3  * (phi[x-1][y][z]))
                 + (diffConst[x][y][z].A4  * (phi[x][y-1][z]))
                 + (diffConst[x][y][z].A5  * (phi[x+1][y+1][z]))
                 + (diffConst[x][y][z].A6  * (phi[x-1][y+1][z]))
                 + (diffConst[x][y][z].A7  * (phi[x-1][y-1][z]))
                 + (diffConst[x][y][z].A8  * (phi[x+1][y-1][z]))
                 + (diffConst[x][y][z].A9  * (phi[x][y][z+1]))
                 + (diffConst[x][y][z].A10 * (phi[x][y][z-1]))
                 + (diffConst[x][y][z].A11 * (phi[x][y+1][z+1]))
                 + (diffConst[x][y][z].A12 * (phi[x][y+1][z-1]))
                 + (diffConst[x][y][z].A13 * (phi[x][y-1][z-1]))
                 + (diffConst[x][y][z].A14 * (phi[x][y-1][z+1]))
                 + (diffConst[x][y][z].A15 * (phi[x+1][y][z+1]))
                 + (diffConst[x][y][z].A16 * (phi[x-1][y][z+1]))
                 + (diffConst[x][y][z].A17 * (phi[x-1][y][z-1]))
                 + (diffConst[x][y][z].A18 * (phi[x+1][y][z-1]));
  
  double result = (SumAphi + (Istim[x][y][z]))/diffConst[x][y][z].sumA;
  return(result);
}



//double FDLaplacianHeteroTissuePhiExtraSaleheen98(int*** &tissue, double*** &phiExtra, sigmatensorMatrix*** &sigmaMintra, sigmatensorMatrix*** &sigmaMextra, double Istim, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal, int rank)
double FDLaplacianHeteroTissuePhiExtraSaleheen98(int*** &tissue, double*** &phiExtra, sigmatensorMatrix*** &sigmaMintra, sigmatensorMatrix*** &sigmaMextra,  double*** &Istim, int &x, int &y, int &z, double &dxReciprocal, double &dyReciprocal, double &dzReciprocal)
{
  int i, j, k;
  sigmatensorMatrix conductivityMatrix[3][3][3];
  double phi[3][3][3];
    
  double sigmaX, sigmaY, sigmaZ;
  double A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18;
    
  double SumA, SumAphi;
  double result = 0.0;         
  double dxdxReciprocal, dydyReciprocal, dzdzReciprocal;
                  
  dxdxReciprocal = dxReciprocal * dxReciprocal;
  dydyReciprocal = dyReciprocal * dyReciprocal;
  dzdzReciprocal = dzReciprocal * dzReciprocal;

  for(i = 0; i < 3; i++)
     for(j = 0; j < 3; j++) 
        for(k = 0; k < 3; k++)
           {
             conductivityMatrix[i][j][k].a11 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a11) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a11)):0.0);
             conductivityMatrix[i][j][k].a12 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a12) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a12)):0.0);
             conductivityMatrix[i][j][k].a13 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a13) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a13)):0.0);
             conductivityMatrix[i][j][k].a22 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a22) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a22)):0.0);
             conductivityMatrix[i][j][k].a23 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a23) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a23)):0.0);
             conductivityMatrix[i][j][k].a33 = ((tissue[x-1+i][y-1+j][z-1+k]!=0)? ((sigmaMintra[x-1+i][y-1+j][z-1+k].a33) + (sigmaMextra[x-1+i][y-1+j][z-1+k].a33)):0.0);

             phi[i][j][k] = (phiExtra[x-1+i][y-1+j][z-1+k]);
           }

  i = j = k = 1;           

  // original formulation - note that I replaced division by multiplication of reciprocal and division by 2 with multiplication by 0.5
  // I also used the central difference method for the differentiation. This results in a factor of 0.25 below
  // NOTE: sigmaX, sigmaY, sigmaZ could be precalculated if no deformation and no chance in conductivity
                                            
  sigmaX  = 0.25 * dxReciprocal * (((conductivityMatrix[i+1][j][k].a11 - conductivityMatrix[i-1][j][k].a11) * dxReciprocal)
                                 + ((conductivityMatrix[i][j+1][k].a12 - conductivityMatrix[i][j-1][k].a12) * dyReciprocal)
                                 + ((conductivityMatrix[i][j][k+1].a13 - conductivityMatrix[i][j][k-1].a13) * dzReciprocal));
  sigmaY  = 0.25 * dyReciprocal * (((conductivityMatrix[i+1][j][k].a12 - conductivityMatrix[i-1][j][k].a12) * dxReciprocal)
                                 + ((conductivityMatrix[i][j+1][k].a22 - conductivityMatrix[i][j-1][k].a12) * dyReciprocal)
                                 + ((conductivityMatrix[i][j][k+1].a23 - conductivityMatrix[i][j][k-1].a23) * dzReciprocal));
  sigmaZ  = 0.25 * dzReciprocal * (((conductivityMatrix[i+1][j][k].a13 - conductivityMatrix[i-1][j][k].a13) * dxReciprocal)
                                 + ((conductivityMatrix[i][j+1][k].a23 - conductivityMatrix[i][j-1][k].a23) * dyReciprocal)
                                 + ((conductivityMatrix[i][j][k+1].a33 - conductivityMatrix[i][j][k-1].a33) * dzReciprocal));
                                                                                                                                                                                                                                                        
  A1      = conductivityMatrix[i+1][j][k].a11 * (dxdxReciprocal - (sigmaX/conductivityMatrix[i][j][k].a11));
  A3      = conductivityMatrix[i-1][j][k].a11 * (dxdxReciprocal + (sigmaX/conductivityMatrix[i][j][k].a11));
    
  A2      = conductivityMatrix[i][j+1][k].a22 * (dydyReciprocal - (sigmaY/conductivityMatrix[i][j][k].a22));
  A4      = conductivityMatrix[i][j-1][k].a22 * (dydyReciprocal + (sigmaY/conductivityMatrix[i][j][k].a22));
        
  A5 = A6 = A7 = A8 = (0.5 * dxReciprocal * dyReciprocal);
  A5 = conductivityMatrix[i+1][j+1][k].a12 * A5;
  A6 = -(conductivityMatrix[i-1][j+1][k].a12 * A6);
  A7 = conductivityMatrix[i-1][j-1][k].a12 * A7;
  A8 = -(conductivityMatrix[i+1][j-1][k].a12 * A8);
                  
  A9      = conductivityMatrix[i][j][k+1].a33 * (dzdzReciprocal - (sigmaZ/conductivityMatrix[i][j][k].a33));
  A10     = conductivityMatrix[i][j][k-1].a33 * (dzdzReciprocal + (sigmaZ/conductivityMatrix[i][j][k].a33));
                      
  A11 = A12 = A13 = A14 = (0.5 * dyReciprocal * dzReciprocal);
  A11 = conductivityMatrix[i][j+1][k+1].a23 * A11;
  A12 = -(conductivityMatrix[i][j+1][k-1].a23 * A12);
  A13 = conductivityMatrix[i][j-1][k-1].a23 * A13;
  A14 = -(conductivityMatrix[i][j-1][k+1].a23 * A14);
  
  A15 = A16 = A17 = A18 = (0.5 * dxReciprocal * dzReciprocal);
  A15 = conductivityMatrix[i+1][j][k+1].a13 * A15;
  A16 = -(conductivityMatrix[i-1][j][k+1].a13 * A16);
  A17 = conductivityMatrix[i-1][j][k-1].a13 * A17;
  A18 = -(conductivityMatrix[i+1][j][k-1].a13 * A18);
                                              
  SumA    = A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10 + A11 + A12 + A13 + A14 + A15 + A16 + A17 + A18;

  SumAphi = (A1  * (phi[i+1][j][k]))
          + (A2  * (phi[i][j+1][k]))
          + (A3  * (phi[i-1][j][k]))
          + (A4  * (phi[i][j-1][k]))
          + (A5  * (phi[i+1][j+1][k]))
          + (A6  * (phi[i-1][j+1][k]))
          + (A7  * (phi[i-1][j-1][k]))
          + (A8  * (phi[i+1][j-1][k]))
          + (A9  * (phi[i][j][k+1]))
          + (A10 * (phi[i][j][k-1]))
          + (A11 * (phi[i][j+1][k+1]))
          + (A12 * (phi[i][j+1][k-1]))
          + (A13 * (phi[i][j-1][k-1]))
          + (A14 * (phi[i][j-1][k+1]))
          + (A15 * (phi[i+1][j][k+1]))
          + (A16 * (phi[i-1][j][k+1]))
          + (A17 * (phi[i-1][j][k-1]))
          + (A18 * (phi[i+1][j][k-1]));
          
  //Istim[x][y][z] = SumAphi - (SumA * (phi[i][j][k]));
  
  result = (SumAphi + (Istim[x][y][z]))/SumA;
  return(result);

  }
                                                                                                                                                                                                                                  
  
  
  
  
                                                                                                                                                                                                                                                                               


// from Saleheen et al. IEEE TBME 1998;45(1):15-25
// Nodes - 18 neighborhood
// 0  - [x][y][z]
// 1  - [x+1][y][z]
// 2  - [x][y+1][z]
// 3  - [x-1][y][z]
// 4  - [x][y-1][z]
// 5  - [x+1][y+1][z]
// 6  - [x-1][y+1][z]
// 7  - [x-1][y-1][z]
// 8  - [x+1][y-1][z]
// 9  - [x][y][z+1]
// 10 - [x][y][z-1]
// 11 - [x][y+1][z+1]
// 12 - [x][y+1][z-1]
// 13 - [x][y-1][z-1]
// 14 - [x][y-1][z+1]
// 15 - [x+1][y][z+1]
// 16 - [x-1][y][z+1]
// 17 - [x-1][y][z-1]
// 18 - [x+1][y][z-1]


// original formulation - note that I replaced division by multiplication of reciprocal and division by 2 with multiplication by 0.5
// I also used the central difference method for the differentiation. This results in a factor of 0.25 below
// NOTE: sigmaX, sigmaY, sigmaZ could be precalculated if no deformation and no chance in conductivity

  // sigmaX  = 0.25 * dxReciprocal * (((conductivityMatrix[x+1][y][z].a11 - conductivityMatrix[x-1][y][z].a11) * dxReciprocal)
  //                               + ((conductivityMatrix[x][y+1][z].a12 - conductivityMatrix[x][y-1][z].a12) * dyReciprocal)
  //                               + ((conductivityMatrix[x][y][z+1].a13 - conductivityMatrix[x][y][z-1].a13) * dzReciprocal));
  // sigmaY  = 0.25 * dyReciprocal * (((conductivityMatrix[x+1][y][z].a12 - conductivityMatrix[x-1][y][z].a12) * dxReciprocal)
  //                               + ((conductivityMatrix[x][y+1][z].a22 - conductivityMatrix[x][y-1][z].a12) * dyReciprocal)
  //                               + ((conductivityMatrix[x][y][z+1].a23 - conductivityMatrix[x][y][z-1].a23) * dzReciprocal));
  // sigmaZ  = 0.25 * dzReciprocal * (((conductivityMatrix[x+1][y][z].a13 - conductivityMatrix[x-1][y][z].a13) * dxReciprocal)
  //                               + ((conductivityMatrix[x][y+1][z].a23 - conductivityMatrix[x][y-1][z].a23) * dyReciprocal)
  //                               + ((conductivityMatrix[x][y][z+1].a33 - conductivityMatrix[x][y][z-1].a33) * dzReciprocal));
                                                                                                                                                                                                            
  // original formulation - note that I replaced division by multiplication of reciprocal and division by 2 with multiplication by 0.5
  //  A1      = conductivityMatrix[x+1][y][z].a11 * (dxdxReciprocal - (sigmaX/conductivityMatrix[x][y][z].a11));
  //  A2      = conductivityMatrix[x][y+1][z].a22 * (dydyReciprocal - (sigmaY/conductivityMatrix[x][y][z].a22));
  //  A3      = conductivityMatrix[x-1][y][z].a11 * (dxdxReciprocal - (sigmaX/conductivityMatrix[x][y][z].a11));
  //  A4      = conductivityMatrix[x][y-1][z].a22 * (dydyReciprocal - (sigmaY/conductivityMatrix[x][y][z].a22));
  //  A5      = conductivityMatrix[x+1][y+1][z].a12 * (0.5 * dxReciprocal * dyReciprocal);
  //  A6      = -(conductivityMatrix[x-1][y+1][z].a12 * (0.5 * dxReciprocal * dyReciprocal));
  //  A7      = conductivityMatrix[x-1][y-1][z].a12 * (0.5 * dxReciprocal * dyReciprocal);
  //  A8      = -(conductivityMatrix[x+1][y-1][z].a12 * (0.5 * dxReciprocal * dyReciprocal));
  //  A9      = conductivityMatrix[x][y][z+1].a33 * (dzdzReciprocal - (sigmaZ/conductivityMatrix[x][y][z].a33));
  //  A10     = conductivityMatrix[x][y][z-1].a33 * (dzdzReciprocal - (sigmaZ/conductivityMatrix[x][y][z].a33));
  //  A11     = conductivityMatrix[x][y+1][z+1].a23 * (0.5 * dyReciprocal * dzReciprocal);
  //  A12     = -(conductivityMatrix[x][y+1][z-1].a23 * (0.5 * dyReciprocal * dzReciprocal));
  //  A13     = conductivityMatrix[x][y-1][z-1].a23 * (0.5 * dyReciprocal * dzReciprocal);
  //  A14     = -(conductivityMatrix[x][y-1][z+1].a23 * (0.5 * dyReciprocal * dzReciprocal));
  //  A15     = conductivityMatrix[x+1][y][z+1].a13 * (0.5 * dxReciprocal * dzReciprocal);
  //  A16     = -(conductivityMatrix[x-1][y][z+1].a13 * (0.5 * dxReciprocal * dzReciprocal));
  //  A17     = conductivityMatrix[x-1][y][z-1].a13 * (0.5 * dxReciprocal * dzReciprocal);
  //  A18     = -(conductivityMatrix[x+1][y][z-1].a13 * (0.5 * dxReciprocal * dzReciprocal));
  