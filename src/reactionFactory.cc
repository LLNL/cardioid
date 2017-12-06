#include "reactionFactory.hh"

#include <cassert>
#include <cstdlib>
#include <mpi.h>
#include "object_cc.hh"
#include "ioUtils.h"
#include "mpiUtils.h"

#include "TT04_CellML_Reaction.hh" // TT04 implementation from CellML (Nov 2011)
#include "TT06_CellML_Reaction.hh" // TT06 implementation from CellML (Nov 2011)
#include "TT06Dev_Reaction.hh"
#include "pade.hh"
#include "TT06_RRG_Reaction.hh"    // TT06 with modifications from Rice et al.
#include "OHaraRudy_Reaction.hh"    // OHara Rudy .
//#include "Grandi_Reaction.hh"      // Grandi
#include "ReactionFHN.hh"
#include "NullReaction.hh"
#include "TestReaction.hh"
#include "ConstantReaction.hh"
#include "ThreadServer.hh"
#include "Anatomy.hh"
#include "string.h"

#include <iostream>
using namespace std;

namespace  scanReaction
{
   Reaction* scanTT04_CellML(OBJECT* obj, const int numPoints);
   Reaction* scanTT06_CellML(OBJECT* obj, const int numPoints);
   Reaction* scanTT06Dev(OBJECT* obj, double dt, const int numPoints, const ThreadTeam& group, const vector<string>& scaleCurrents);
   Reaction* scanTT06_RRG(OBJECT* obj, const int numPoints);
   Reaction* scanOHaraRudy(OBJECT* obj, const int numPoints);
   Reaction* scanSimpleOHaraRudy(OBJECT* obj, const int numPoints);
   Reaction* scanGrandi(OBJECT* obj, const int numPoints);
   Reaction* scanSimpleGrandi(OBJECT* obj, const int numPoints);
   Reaction* scanFHN(OBJECT* obj, const int numPoints);
   Reaction* scanPassive(OBJECT* obj, const int numPoints);
   Reaction* scanNull(OBJECT* obj);
   Reaction* scanTest(OBJECT* obj);
}


Reaction* reactionFactory(const string& name, double dt, const int numPoints,
                          const ThreadTeam& group, const vector<string>& scaleCurrents)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   
   OBJECT* obj = objectFind(name, "REACTION");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "TT04_CellML" || method == "tenTusscher04_CellML")
      return scanReaction::scanTT04_CellML(obj, numPoints);
   else if (method == "TT06_CellML" || method == "tenTusscher06_CellML")
      return scanReaction::scanTT06_CellML(obj, numPoints);
   else if (method == "TT06Dev" )
      return scanReaction::scanTT06Dev(obj, dt, numPoints, group, scaleCurrents);
   else if (method == "TT06Opt" )
      return scanReaction::scanTT06Dev(obj, dt, numPoints, group, scaleCurrents);
   else if (method == "TT06_RRG" )
      return scanReaction::scanTT06_RRG(obj, numPoints);
   else if (method == "OHaraRudy" )
      return scanReaction::scanOHaraRudy(obj, numPoints);
   else if (method == "SimpleOHaraRudy" )
      return scanReaction::scanSimpleOHaraRudy(obj, numPoints);
   else if (method == "Grandi" )
      return scanReaction::scanGrandi(obj, numPoints);
   else if (method == "SimpleGrandi" )
      return scanReaction::scanSimpleGrandi(obj, numPoints);
   else if (method == "FHN" || method == "FitzhughNagumo")
      return scanReaction::scanFHN(obj, numPoints);
   else if (method == "Passive")
      return scanReaction::scanPassive(obj, numPoints);
   else if (method == "null")
      return scanReaction::scanNull(obj);
   else if (method == "test")
      return scanReaction::scanTest(obj);
   
   if (myRank == 0)
      cerr<<"ERROR: Undefined reaction model in reactionFactory"<<endl;
   assert(false); // reachable only due to bad input
   return 0;
}

namespace  scanReaction
{
   Reaction* scanTT04_CellML(OBJECT* obj, const int numPoints)
   {
      TT04_CellML_Reaction::IntegratorType integrator;
      string tmp;
      objectGet(obj, "integrator", tmp, "rushLarsen");
      if      (tmp == "rushLarsen")   integrator = TT04_CellML_Reaction::rushLarsen;
      else if (tmp == "rushLarson")   integrator = TT04_CellML_Reaction::rushLarsen;
      else if (tmp == "forwardEuler") integrator = TT04_CellML_Reaction::forwardEuler;
      else    assert(false);
      int ttType;
      objectGet(obj, "ttType", ttType, "0");
      
      return new TT04_CellML_Reaction(numPoints, ttType, integrator);
   }
}

namespace  scanReaction
{
   Reaction* scanTT06_CellML(OBJECT* obj, const int numPoints)
   {
      TT06_CellML_Reaction::IntegratorType integrator;
      string tmp;
      objectGet(obj, "integrator", tmp, "rushLarsen");
      if      (tmp == "rushLarsen")   integrator = TT06_CellML_Reaction::rushLarsen;
      else if (tmp == "forwardEuler") integrator = TT06_CellML_Reaction::forwardEuler;
      else    assert(false);    
      int ttType;
      objectGet(obj, "ttType", ttType, "0");
      return new TT06_CellML_Reaction(numPoints, ttType, integrator);
   }
}

namespace  scanReaction
{
   PADE *scanFit()
   {
      if (!object_exists("functions", "FIT") ) return NULL;
      OBJECT* obj = object_find("functions", "FIT");
      vector<string>fitName; 
      objectGet(obj,"functions",fitName); 
      PADE *fit =  new PADE  [fitName.size()] ;
      for (int i=0;i<fitName.size();i++)
      {
         vector<double> coef; 
         double tol,deltaV,V0,V1;
         int l,m;
         string name = fitName[i]; 
         OBJECT* obj = object_find(name.c_str(), "FIT");
         objectGet(obj,"tol", tol,"1.0"); 
         objectGet(obj,"deltaX",deltaV,"1"); 
         objectGet(obj,"x0"    ,    V0,"0"); 
         objectGet(obj,"x1"    ,    V1,"0"); 
         objectGet(obj,"l"     ,     l,"0"); 
         objectGet(obj,"m"     ,     m,"0"); 
         objectGet(obj,"coef"  ,coef); 
         padeApprox(fit[i],name,fitFuncMap(name),NULL,0,deltaV,V0,V1,tol,0,0,0,l,m,&coef[0]); 
      //   fit[i].aparms=fit+i; 
      }
      return fit; 
   }
    Reaction* scanTT06Dev(OBJECT* obj, double dt, const int numPoints, const ThreadTeam& group, const vector<string>& scaleCurrents)
   {
      TT06Dev_ReactionParms parms;
      map<string,TT06Func::CellTypeParmsFull> cellTypeParms = TT06Func::getStandardCellTypes(); 
      int fastGate =-1; 
      int fastNonGate =-1; 
      TT06Func::initCnst(); 
      string fitFit; 
/*
      if (object_testforkeyword(obj,"mod") )
      {
          //if (getRank(0) == 0) printf("keyword <mod> is deprecated by <jhTauSmooth>.  Replace keyword <mod> in the reaction object <%s>  with <jhTauSmooth>\n",obj->name); 
          
          //exit(1); 
         objectGet(obj, "mod",          parms.jhTauSmooth, "0") ;
      }
*/
      objectGet(obj, "mod",          parms.jhTauSmooth, "0") ;
      objectGet(obj, "tolerance",    parms.tolerance, "0.0") ;
      objectGet(obj, "fastReaction", parms.fastReaction, "-1") ;
      objectGet(obj, "fastGate",     fastGate, "-1") ;
      objectGet(obj, "fastNonGate",  fastNonGate, "-1") ;
      if (object_testforkeyword(obj, "gateThreadMap") ) objectGet(obj,"gateThreadMap",parms.gateThreadMap);
      parms.fit = NULL; 
      parms.fitFile = "fit.data"; 
      if (object_testforkeyword(obj,"fitFile") )
      {
        objectGet(obj, "fitFile",     parms.fitFile,"fit.data"); 
        int fileFitExists=0; 
        if (getRank(0) ==0) 
        {
          if (filetest(parms.fitFile.c_str(), S_IFREG) == 0) fileFitExists=1; 
          if (fileFitExists) object_compilefile(parms.fitFile.c_str()); 
        }
        MPI_Bcast(&fileFitExists, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (fileFitExists) 
        {
           object_Bcast(0, MPI_COMM_WORLD);
           parms.fit = scanFit();
        }
      }
      if (parms.fastReaction == -1) 
      {
         if (fastGate    ==  -1 )  fastGate   =0; 
         if (fastNonGate ==  -1 )  fastNonGate=0; 
         parms.fastReaction = fastGate+256*fastNonGate; 
      }
      string cellTypeName;
      objectGet(obj, "cellTypeName", cellTypeName, "endoCellML");

      OBJECT* cellobj = object_find2(cellTypeName.c_str(), "CELLTYPE", IGNORE_IF_NOT_FOUND);
      if (cellobj) {
         string clone; 
         objectGet(cellobj, "clone", clone, "") ;
         if (clone != "")
         {
            assert(cellTypeName != clone); 
            assert(cellTypeParms.count(clone) != 0);
            cellTypeParms[cellTypeName]=cellTypeParms[clone]; 
            cellTypeParms[cellTypeName].name=cellTypeName;
         }
         if (object_testforkeyword(cellobj, "s_switch")       ) objectGet(cellobj,"s_switch"      ,cellTypeParms[cellTypeName].s_switch,"0"); 
         if (object_testforkeyword(cellobj, "P_NaK")          ) objectGet(cellobj,"P_NaK"         ,cellTypeParms[cellTypeName].P_NaK,"0.0"); 
         if (object_testforkeyword(cellobj, "g_NaL")          ) objectGet(cellobj,"g_NaL"         ,cellTypeParms[cellTypeName].g_NaL,"0.0"); 
         if (object_testforkeyword(cellobj, "g_Ks")           ) objectGet(cellobj,"g_Ks"        ,cellTypeParms[cellTypeName].g_Ks,"0.0"); 
         if (object_testforkeyword(cellobj, "g_Kr")           ) objectGet(cellobj,"g_Kr"        ,cellTypeParms[cellTypeName].g_Kr,"0.0"); 
         if (object_testforkeyword(cellobj, "g_to")           ) objectGet(cellobj,"g_to"        ,cellTypeParms[cellTypeName].g_to,"0.0"); 
      } else {
         assert(cellTypeName == cellTypeParms[cellTypeName].name);
      }
      parms.cellTypeParm = cellTypeParms[cellTypeName];

      parms.currentNames = scaleCurrents;
      
      Reaction *reaction = new TT06Dev_Reaction(dt, numPoints, parms, group);
      return  reaction; 
   }
}

namespace  scanReaction
{
   Reaction* scanTT06_RRG(OBJECT* obj, const int numPoints)
   {
      TT06_RRG_ReactionParms parms;
      objectGet(obj, "Ko",    parms.Ko, "-1") ;
      int ttType;
      objectGet(obj, "ttType", ttType, "0");
      Reaction *reaction = new TT06_RRG_Reaction(numPoints,ttType,parms);
      return  reaction; 
   }
}


namespace  scanReaction
{
   Reaction* scanFHN(OBJECT* obj, const int numPoints)
   {
      // None of the FHN model parameters are currently wired to the
      // input deck.
      return new ReactionFHN(numPoints);
   }
}

namespace  scanReaction
{
   Reaction* scanNull(OBJECT* obj)
   {
      NullReactionParms parms;
      objectGet(obj, "V0",    parms.initialVoltage, "-85",  "voltage");
      return new NullReaction(parms);
   }
}

namespace  scanReaction
{
   Reaction* scanTest(OBJECT* obj)
   {
      TestReactionParms parms;
      objectGet(obj, "V0",    parms.initialVoltage, "-85",  "voltage");
      objectGet(obj, "delta", parms.delta,          "1e-3", "voltage");
      return new TestReaction(parms);
   }
}

namespace  scanReaction
{
/*
   Reaction* scanConstant(OBJECT* obj, const int numPoints)
   {
      vector<double> eta;
      objectGet(obj, "eta", eta);
      
      double alpha = 0.0; 
      objectGet(obj, "alpha", alpha, "0.0") ;
      double beta = 0.0; 
      objectGet(obj, "beta", beta, "0.0") ;
      double gamma = 0.0; 
      objectGet(obj, "gamma", gamma, "0.0") ;
      int printRate;
      objectGet(obj, "printRate", printRate, "0") ;

      SymmetricTensor sigma1, sigma2, sigma3;
      double* buffer;
      int nn=object_getv(obj, "sigma1", (void*)&buffer, DOUBLE, ABORT_IF_NOT_FOUND);
      assert( nn==6 );
      sigma1.a11=buffer[0];
      sigma1.a12=buffer[1];
      sigma1.a13=buffer[2];
      sigma1.a22=buffer[3];
      sigma1.a23=buffer[4];
      sigma1.a33=buffer[5];
      free(buffer);
      nn=object_getv(obj, "sigma2", (void*)&buffer, DOUBLE, ABORT_IF_NOT_FOUND);
      assert( nn==6 );
      sigma2.a11=buffer[0];
      sigma2.a12=buffer[1];
      sigma2.a13=buffer[2];
      sigma2.a22=buffer[3];
      sigma2.a23=buffer[4];
      sigma2.a33=buffer[5];
      free(buffer);
      nn=object_getv(obj, "sigma3", (void*)&buffer, DOUBLE, ABORT_IF_NOT_FOUND);
      assert( nn==6 );
      sigma3.a11=buffer[0];
      sigma3.a12=buffer[1];
      sigma3.a13=buffer[2];
      sigma3.a22=buffer[3];
      sigma3.a23=buffer[4];
      sigma3.a33=buffer[5];
      free(buffer);

      return new ConstantReaction(numPoints,
                                  eta,
                                  sigma1, sigma2, sigma3,
                                  alpha, beta, gamma,
                                  printRate);
   }
*/
}
