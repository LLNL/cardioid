// FlaimParameterType.h
// First written June 26th, 2009
// Author: Matthias K Reumann
// Copy Right: IBM
//
// This file contains the type definition of the parameter struct for parameters of the Flaim model

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;



#define STIMULATION
#ifdef STIMULATION
const int RangeTab=5000;
const int DivisionTab=10;
#else
const int RangeTab=400;
const int DivisionTab=10;
#endif

const int V_RESOLUTION=RangeTab*DivisionTab;
const int RangeTabhalf=RangeTab/2;
const double dDivisionTab=1.0/DivisionTab;


#define TT04_STATE_VARIABLES 17

enum TTStateVariables
{
  tt_V,
  tt_m,
  tt_h,
  tt_j,
  tt_xr1,
  tt_xr2,
  tt_xs,
  tt_r,
  tt_s,
  tt_d,
  tt_f,
  tt_fCa,
  tt_g,
  tt_Ca_i,
  tt_CaSR,
  tt_Na_i,
  tt_K_i
};


typedef struct
{
  double Rgas;
  double Tx;
  double Faraday;
  double K_o;
  double Ca_o;
  double Na_o;
  double Vc;
  double Vsr;
  double Bufc;
  double Kbufc;
  double Bufsr;
  double Kbufsr;
  double taufca;
  double taug;
  double Vmaxup;
  double Kup;
  double C;
  double g_Kr;
  double pKNa;
  double g_Ks;
  double g_K1;
  double g_to;
  double g_Na;
  double g_bNa;
  double KmK;
  double KmNa;
  double knak;
  double g_CaL;
  double g_bCa;
  double kNaCa;
  double KmNai;
  double KmCa;
  double ksat;
  double n;
  double g_pCa;
  double KpCa;
  double g_pK;
  double s_inf_vHalf;
  double tau_s_f1;
  double tau_s_slope1;
  double tau_s_vHalf1;
  double tau_s_f2;
  double tau_s_f3;

  double RToverF;
  double inverseRToverF;
  double inverseVcF2C;
  double inverseVcFC;
  double VcdVsr;
  double Kupsquare;
  double BufcPKbufc;
  double Kbufcsquare;
  double Kbufc2;
  double BufsrPKbufsr;
  double Kbufsrsquare;
  double Kbufsr2;
  double KopKNaNao;
  double KmNai3;
  double Nao3;
                              
  // not in parameter file
  double m_Xr1_1;
  double m_Xr1_2;
  double a_Xr1_1;
  double a_Xr1_2;
  double b_Xr1_1;
  double b_Xr1_2;
  double K_Q10Xr1;
  double m_Xr2_1;
  double m_Xr2_2;
  double a_Xr2_1;
  double a_Xr2_2;
  double b_Xr2_1;
  double b_Xr2_2;
  double K_Q10Xr2;


} TTParameterStruct;

static int loadTTParameter(TTParameterStruct &parameter, double* &y, char *parameterFile)
{

  FILE *inputFP;
  char readline[255];
  char delims[] = " ";
  char *ptmpstring = NULL;
  int numparameters = 0;
  int nparam = 0;
  bool readall = false;      

  // setting default values

    
  inputFP = fopen (parameterFile, "r");
  if (inputFP != NULL)
     {
       int lines = 0;
       //while ((!feof(inputFP)) && (readall == false))
       while (!feof(inputFP))
          {
            lines++;
            fgets(readline, 255, inputFP);
            
            if (readline != NULL)
               {
                 if ((atoi(readline)) != 0)
                    {
                      numparameters = atoi(readline);
                      //cout << "number of parameters " << numparameters << endl;
                    }
                 else   
                    {
                      ptmpstring=strtok (readline, delims);
                      //cout << "line " << lines << " parameter " << nparam << "/" << numparameters << " first ptmpstring " << ptmpstring << endl;               
                    }
               }     
            else
               {
                 //cout << "ERROR: fgets(readline) returned NULL at line " << (lines + 1) << endl;
               }

            // IBT Magic Numer stuff is disregarded because it will not match a string
            // the same for the cell name
            // the same for the parameter number            
            
            if (ptmpstring!=NULL)
               {
                 nparam++;
                 if (strcmp("R", ptmpstring) == 0)
                    {
                      ptmpstring=strtok (NULL, delims);
                      if (ptmpstring!=NULL)
                           parameter.Rgas = (double) atof(ptmpstring);
                      //cout << "parameter.Rgas " << parameter.Rgas << endl;
                    }
                    
                 else if (strcmp("Tx", ptmpstring) == 0)
                    {
                      ptmpstring=strtok (NULL, delims);
                      if (ptmpstring!=NULL)
                           parameter.Tx = (double) atof(ptmpstring);
                      //cout << "parameter.Tx " << parameter.Tx << endl;
                    }
                 else if (strcmp("F", ptmpstring) == 0)
                    {
                      ptmpstring=strtok (NULL, delims);
                      if (ptmpstring!=NULL)
                           parameter.Faraday = (double) atof(ptmpstring);
                    }

                 else
                    {
                      nparam--;
                    }
                    
               } // end if ()
          } // end while
          
//       cout << "Closing file" << endl;
       fclose(inputFP);        
     } // end if (inputFP != NULL)
  else
     {
       cout << "Error opening file " << parameterFile << endl;
     }
  
  return(0);
}