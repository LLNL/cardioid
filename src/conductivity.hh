//#include "setupParameters.h"

#ifndef PI
 #define PI 3.14159265358979323846
#endif
 

typedef struct
   {
     double x;
     double y;
     double z;
   }sigmatensor;

typedef struct
   {
     double a11, a12, a13;
     double      a22, a23; // symmetric matrix. Thus a21 = a12
     double           a33; // symmetric matrix. Thus a31 = a13 and a32 = a23
   }sigmatensorMatrix;

typedef struct
   {
     double A1;
     double A2;
     double A3;
     double A4;
     double A5;
     double A6;
     double A7;
   }diffusion2D;
                                                                                                     



typedef struct
   {
     double sumA;
     double A1;
     double A2;
     double A3;
     double A4;
     double A5;
     double A6;
     double A7;
     double A8;
     double A9;
     double A10;
     double A11;
     double A12;
     double A13;
     double A14;
     double A15;
     double A16;
     double A17;
     double A18;
   }diffusion;


template<class T, int count> class kaTabularizedFunction 
{
   public:
      T val[count];
      kaTabularizedFunction(T (*func)(T))
         {
           for(int i=0; i<count; i++)
              val[i]=func((T)i*(T)PI/(T)count);
         }
};
                           
//extern kaTabularizedFunction<double, 255> SinTable;
//extern kaTabularizedFunction<double, 255> CosTable;
                           

int TablularizeSin(double* &SinTable, int count)
{
  for(int i=0; i<count; i++)
     {
       SinTable[i] = sin(i*PI/count);
     }

  return(0);
}  

int TablularizeCos(double* &CosTable, int count)
{
  for(int i=0; i<count; i++)
     {
       CosTable[i] = cos(i*PI/count);
     }
                   
  return(0);
}
                     



/*               
// This is according to my understanding
// See Co-ordinate transform theory and my documentation notes
// Currently co-ordinate transform 3 is used similar to Seemann G PhD Thesis
// Maybe should just implement rotation matrix (http://genex.hgu.mrc.ac.uk/Software/paint/paint/node5.html) and multiply with local co-ordinates
int CalcConductivityMatrix(sigmatensorMatrix *conductivity, double *sigmal, double *sigmat, unsigned char *cphi, unsigned char *ctheta)
{
  double sinPhi, cosPhi, sinPhiSquared, cosPhiSquared;//according to Colli-Franzone et al.
  double sinTheta, cosTheta, sinThetaSquared, cosThetaSquared;
  double tmp;

  sinPhi = SinTable.val[(*cphi)];
  cosPhi = CosTable.val[(*cphi)];
  sinTheta = SinTable.val[(*ctheta)];
  cosTheta = CosTable.val[(*ctheta)];
              
  sinPhiSquared = sinPhi * sinPhi;
  cosPhiSquared = cosPhi * cosPhi;
  sinThetaSquared = sinTheta * sinTheta;
  cosThetaSquared = cosTheta * cosTheta;
                      
  // Equations can still be optimized computationally!!!
                        
  //M(0,0)
  conductivity->a11 = ((*sigmal) * cosPhiSquared * cosThetaSquared) + ((*sigmat) * ((cosPhiSquared * cosThetaSquared) + sinThetaSquared));
  //M(1,1)
  conductivity->a22 = ((*sigmal) * sinPhiSquared * cosThetaSquared)+((*sigmat) * (cosPhiSquared + (sinThetaSquared * sinPhiSquared)));
  //M(2,2)
  conductivity->a33 = ((*sigmal) * sinThetaSquared)+((*sigmat) *  cosThetaSquared);
                                    
                                    
  //M(0,1) = M(1,0)
  conductivity->a12 = ((*sigmat)  - (*sigmal)) * cosPhi * cosTheta * sinPhi;
  //M(0,2) = M(2,0)
  conductivity->a13 = ((*sigmal) - (*sigmat)) * cosTheta * cosPhi * sinTheta;
  //M(1,2) = M(2,1)
  conductivity->a23 = ((*sigmat) - (*sigmal)) * sinPhi * cosTheta * sinTheta;
                                                
  return 0;
}

                                                 
// This is according to my understanding
// See Co-ordinate transform theory and my documentation notes
// Currently co-ordinate transform 3 is used similar to Seemann G PhD Thesis
// Maybe should just implement rotation matrix (http://genex.hgu.mrc.ac.uk/Software/paint/paint/node5.html) and multiply with local co-ordinates
int CalcXYZConductivities(sigmatensor *conductivity, double *sigmal, double *sigmat, unsigned char *cphi, unsigned char *ctheta)
{
  double sinPhi, cosPhi, sinPhiSquared, cosPhiSquared;//according to Colli-Franzone et al.  
  double sinTheta, cosTheta, sinThetaSquared, cosThetaSquared;
  double tmp;

  sinPhi = SinTable.val[(*cphi)];
  cosPhi = CosTable.val[(*cphi)];
  sinTheta = SinTable.val[(*ctheta)];
  cosTheta = CosTable.val[(*ctheta)];
     
  sinPhiSquared = sinPhi * sinPhi;
  cosPhiSquared = cosPhi * cosPhi;
  sinThetaSquared = sinTheta * sinTheta;
  cosThetaSquared = cosTheta * cosTheta;
  
  // Equations can still be optimized computationally!!!
  
  //M(0,1) = M(1,0)
  conductivity->y = ((*sigmat)  - (*sigmal)) * cosPhi * cosTheta * sinPhi;
  //M(0,2) = M(2,0)
  conductivity->z = ((*sigmal) - (*sigmat)) * cosTheta * cosPhi * sinTheta;
  //M(1,2) = M(2,1)
  tmp = ((*sigmat) - (*sigmal)) * sinPhi * cosTheta * sinTheta;

  //M(0,0)
  conductivity->x = ((*sigmal) * cosPhiSquared * cosThetaSquared) + ((*sigmat) * ((cosPhiSquared * cosThetaSquared) + sinThetaSquared));      
  //sigma_x = M(0,0)+M(0,1)+M(0,2)
  conductivity->x += (conductivity->y + conductivity->z); // this is the new sigma_x;
  //sigma_y = M(1,0) + M(1,2) + M(1,2)
  conductivity->y += (((*sigmal) * sinPhiSquared * cosThetaSquared)+((*sigmat) * (cosPhiSquared + (sinThetaSquared * sinPhiSquared))) + tmp);
  //sigma_z = M(2,0) + M(2,1) + M(2,2)  
  conductivity->z += (((*sigmal) * sinThetaSquared)+((*sigmat) *  cosThetaSquared) + tmp);

//  cout << "Transform 3: sigma_x " << conductivity->x;
//  cout << " sigma_y " << conductivity->y;
//  cout << " sigma_z " << conductivity->z << endl;

  return 0;
}

*/

// This is according to the IBT transformation/rotation Matrix R as implemented in the code
int CalcConductivityMatrixIBT(sigmatensorMatrix *conductivity, double sigmal, double sigmat, int phi, int theta)
{
  double sinPhi, cosPhi, sinPhiSquared, cosPhiSquared;
  double sinTheta, cosTheta, sinThetaSquared, cosThetaSquared;
    
  double tmp;
  
/*    
  double* SinTable = (double*) malloc(sizeof(double)*256);
  double* CosTable = (double*) malloc(sizeof(double)*256);
  
  TablularizeSin(SinTable, 256);
  TablularizeCos(CosTable, 256);


  if (  ((phi >=0) && (phi < 256)) 
      &&((phi >=0) && (phi < 256)) )
     {
       phi = 0;
       theta = 0; 
     }
  else
     {
       cout << "phi/theta value not within scope" << endl;
     }     
          
          
  sinPhi = SinTable[phi];
  cosPhi = CosTable[phi];
  sinTheta = SinTable[theta];
  cosTheta = CosTable[theta];
                              
      
//  sinPhi = SinTable.val[(phi)];
//  cosPhi = CosTable.val[(phi)];
//  sinTheta = SinTable.val[(theta)];
//  cosTheta = CosTable.val[(theta)];

*/


//SinTable[i] = sin(i*PI/count);
//CosTable[i] = cos(i*PI/count);
                                         
  sinPhi   = sin(  phi*PI/256); // SinTable[phi];
  cosPhi   = cos(  phi*PI/256); // CosTable[phi];
  sinTheta = sin(theta*PI/256); // SinTable[theta];
  cosTheta = cos(theta*PI/256); // CosTable[theta];
                                                 
              
  sinPhiSquared = sinPhi * sinPhi;
  cosPhiSquared = cosPhi * cosPhi;
  sinThetaSquared = sinTheta * sinTheta;
  cosThetaSquared = cosTheta * cosTheta;
                      
  // Equations can still be optimized computationally!!!
                        
  //M(0,0)
  conductivity->a11 = ((sigmal) * cosPhiSquared * sinThetaSquared) + ((sigmat) * ((cosPhiSquared * cosThetaSquared) + sinPhiSquared));
  //M(1,1)
  conductivity->a22 = ((sigmal) * sinPhiSquared * sinThetaSquared)+((sigmat) * (cosPhiSquared + (sinPhiSquared*cosThetaSquared)));
  //M(2,2)
  conductivity->a33 = ((sigmal) * cosThetaSquared)+((sigmat) *  sinThetaSquared);
                                    
  //M(0,1) = M(1,0)
  conductivity->a12 = ((sigmal)  - (sigmat)) * cosPhi * sinThetaSquared * sinPhi;
  //M(0,2) = M(2,0)
  conductivity->a13 = ((sigmal) - (sigmat)) * cosTheta * cosPhi * sinTheta;
  //M(1,2) = M(2,1)
  conductivity->a23 = ((sigmal) - (sigmat)) * sinPhi * cosTheta * sinTheta;

  if (conductivity->a11 > 0.0)
     {
     }
  else
     {
  cout << sinPhi << " " << cosPhi << " " << sinTheta << " " << cosTheta << "\n"
       << conductivity->a11 << "\t"
       << conductivity->a12 << "\t"
       << conductivity->a13 << "\n"
       << conductivity->a22 << "\t"
       << conductivity->a23 << "\n"
       << conductivity->a33 << "\t"
       << endl;
     }                         

  return 0;
}
                                                  

/*
// This is according to the IBT transformation/rotation Matrix R as implemented in the code
int CalcXYZConductivitiesIBT(sigmatensor *conductivity, double *sigmal, double *sigmat, unsigned char *cphi, unsigned char *ctheta)
{
  double sinPhi, cosPhi, sinPhiSquared, cosPhiSquared;
  double sinTheta, cosTheta, sinThetaSquared, cosThetaSquared;
    
  double tmp;
  
  sinPhi = SinTable.val[(*cphi)];
  cosPhi = CosTable.val[(*cphi)];
  sinTheta = SinTable.val[(*ctheta)];
  cosTheta = CosTable.val[(*ctheta)];
              
  sinPhiSquared = sinPhi * sinPhi;
  cosPhiSquared = cosPhi * cosPhi;
  sinThetaSquared = sinTheta * sinTheta;
  cosThetaSquared = cosTheta * cosTheta;
                      
  // Equations can still be optimized computationally!!!
                        
  //M(0,1) = M(1,0)
  conductivity->y = ((*sigmal)  - (*sigmat)) * cosPhi * sinThetaSquared * sinPhi;
  //M(0,2) = M(2,0)
  conductivity->z = ((*sigmal) - (*sigmat)) * cosTheta * cosPhi * sinTheta;
  //M(1,2) = M(2,1)
  tmp = ((*sigmal) - (*sigmat)) * sinPhi * cosTheta * sinTheta;
                                    
  //M(0,0)
  conductivity->x = ((*sigmal) * cosPhiSquared * sinThetaSquared) + ((*sigmat) * ((cosPhiSquared * cosThetaSquared) + sinPhiSquared));
  //sigma_x = M(0,0)+M(0,1)+M(0,2)
  conductivity->x += (conductivity->y + conductivity->z); // this is the new sigma_x;
  //sigma_y = M(1,0) + M(1,2) + M(1,2)
  conductivity->y += (((*sigmal) * sinPhiSquared * sinThetaSquared)+((*sigmat) * (cosPhiSquared + (sinPhiSquared*cosThetaSquared))) + tmp);
  //sigma_z = M(2,0) + M(2,1) + M(2,2)
  conductivity->z += (((*sigmal) * cosThetaSquared)+((*sigmat) *  sinThetaSquared) + tmp);

//  cout << "IBT: sigma_x " << conductivity->x;
//  cout << " sigma_y " << conductivity->y;
//  cout << " sigma_z " << conductivity->z << endl;

  return 0;
}

*/

/*
 * Compute conductivity tensor in a block of tissue
 * Original code from Tabish Almas at JHU
 * Modified by Matthias Reumann at IBM T. J. Watson Research Center
 * Date Oct. 1st 2008
 * The code is mostly left unmodified. Only small changes were added to the function calls
 * i.e. the interface to the IBM code was established
 */
       
       
double getFiberAngleJHU(int x, int y, int z, int transmuralAxis)
{
  double FiberAngle;
  double twist_per_cell;
  double Endo_angle = 60.0; // +60 degrees
  double Epi_angle = -60.0; // -60 degrees
  double total_twist;
                 
                 
#ifdef HOMOGENEOUS_FIBER_ORIENTATION
  // Homogeneous Wedge
  // Constant fiber angle, i.e. fibers are parallel to each other
  FiberAngle = -60;             // Set this to the desired value
#endif
                       
#ifdef NONHOMOGENEOUS_FIBER_ORIENTATION
  // Tranmural wedge
  // Fibers twist linearly from +60 degrees at Endo to -60 degrees at Epi. twist angle = 120 degrees
  // fiber twist = twist angle/fiber length
                             
  total_twist = Endo_angle - Epi_angle;    // 120 degrees
  twist_per_cell =  total_twist/transmuralAxis;      // YDIM is TM axis (in JHU code) defined by max y coord
                                 
  FiberAngle = Endo_angle - (twist_per_cell * (y+1));
#endif
                                   
/*******************
#ifdef LAYERS_OF_FIBERS
  // 10 Layers of Transmural Tissue with constant intra-layer fiber orientation
  // but different inter-layer fiber orientation. Total fiber twist = 120 degrees
  // From Layer 1 (+60 degrees at Endo) to Layer 10 (-60 degrees at EPI).
  // For 1.5cm TM length at 0.01cm spatial resolution (i.e. 150 nodes in TM direction in the model wedge)
  // We have 15 layers with inter-layer shoft of 8 degrees to give total 120 degrees twist.
                                             
  // LOOK IN ORIGINAL CODE.
                                               
#endif
*******************/
                                               
  return FiberAngle;
}


double getImbriAngleJHU(int x, int y, int z)
{
  double ImbrAng;
  ImbrAng =0;
  return ImbrAng;
}
      
double getSheetAngleJHU(int x, int y, int z)
{
  double SheetAngle;
  SheetAngle = SHEET_ANGLE;
  return SheetAngle;
}
            
int GetDiffusionTensorJHU(sigmatensorMatrix *conductivity, int x, int y, int z, double G_L, double G_T, double G_N, int transmuralAxis)
{
  int i;
  double FiberAngle, ImbriAngle, SheetAngle;
  double FibAng, ImbrAng, ShetAng;
  double determinant_E, determinant_D;
  double** D;
                      
  D = new double*[3];
  for(i=0;i<3;i++)
     {
       D[i]= new double[3];
     }
  double E1[3];
  double E3[3];
  double E2[3];
                                                 
  //D = R* [G_L 0 0; 0 G_T 0; 0 0 G_N]*R' where R = [E1 E2 E3] = Euler Vecotr Rotation Matrix
                                                   
  FiberAngle = getFiberAngleJHU(x, y, z, transmuralAxis);
  ImbriAngle = getImbriAngleJHU(x, y, z);
  SheetAngle = getSheetAngleJHU(x, y, z);
                                                         
  // Convert from Degrees to radians
  FibAng = FiberAngle*PI/180;
  ImbrAng = ImbriAngle*PI/180;
  ShetAng = SheetAngle*PI/180;
                                                                 
  // !!!!!!!!!!!!!!!
  // NEEDS MORE WORK - comment by Tabish
  // !!!!!!!!!!!!!!!!
                                                                       

#if ROTATION_MATRIX==0
  // Use Katherine's Rotation Matrix
  
  E1[0] = cos(FibAng);               // E_1x
  E1[1] = 0;                         // E_1y
  E1[2] = sin(FibAng);               // E_1z
  E2[0] = -sin(FibAng)*sin(ShetAng); // E_2x
  E2[1] = cos(ShetAng);              // E_2y
  E2[2] = cos(FibAng)*sin(ShetAng);  // E_2z
  E3[0] = sin(FibAng)*cos(ShetAng);  // E_3x
  E3[1] = sin(ShetAng);              // E_3y
  E3[2] = -cos(FibAng)*cos(ShetAng); // E_3z
  
#elif ROTATION_MATRIX==1
  E1[0] = cos(FibAng);               // E_1x
  E1[1] = 0;                         // E_1y
  E1[2] = -sin(FibAng);              // E_1z
  E2[0] = sin(FibAng)*sin(ShetAng);  // E_2x
  E2[1] = cos(ShetAng);              // E_2y
  E2[2] = cos(FibAng)*sin(ShetAng);  // E_2z
  E3[0] = sin(FibAng)*cos(ShetAng);  // E_3x
  E3[1] = -sin(ShetAng);             // E_3y
  E3[2] = cos(FibAng)*cos(ShetAng);  // E_3z
                                      
#else
  // DO NOTHING
#endif  // ROTATION_MATRIX
                                        
  // Check: Determinant of Rotation matrix 'E' should be 1
  determinant_E = E1[0]*E2[1]*E3[2] + E2[0]*E3[1]*E1[2] + E3[0]*E1[1]*E2[2] -
                  E1[0]*E2[2]*E3[1] - E2[0]*E1[1]*E3[2] - E3[0]*E2[1]*E1[2];
                                                              
  D[0][0] = G_L*pow(E1[0],2) + G_T*pow(E2[0],2) + G_N*pow(E3[0],2);
  D[1][0] = G_L*E1[0]*E1[1] + G_T*E2[0]*E2[1] + G_N*E3[0]*E3[1];
  D[2][0] = G_L*E1[0]*E1[2] + G_T*E2[0]*E2[2] + G_N*E3[0]*E3[2];
  D[0][1] = D[1][0];
  D[1][1] = G_L*pow(E1[1],2) + G_T*pow(E2[1],2) + G_N*pow(E3[1],2);
  D[2][1] = G_L*E1[1]*E1[2] + G_T*E2[1]*E2[2] + G_N*E3[1]*E3[2];
  D[0][2] = D[2][0];
  D[1][2] = D[2][1];
  D[2][2] = G_L*pow(E1[2],2) + G_T*pow(E2[2],2) + G_N*pow(E3[2],2);
                                                                                
  // Check: Determinant of Diffusion Tensor 'D' should be 1
  determinant_D = D[0][0]*D[1][1]*D[2][2] + D[0][1]*D[1][2]*D[2][0] + D[0][2]*D[1][0]*D[2][1] -
                  D[0][0]*D[1][2]*D[2][1] - D[0][2]*D[1][1]*D[2][0] - D[0][1]*D[1][0]*D[2][2];
                                                                                                      

  // setting the conductivity matrix of the IBM code
  conductivity->a11 = D[0][0];
  conductivity->a12 = D[1][0];
  conductivity->a13 = D[2][0];
  conductivity->a22 = D[1][1];
  conductivity->a23 = D[1][2];
  conductivity->a33 = D[2][2];
              
  return(0);
}
