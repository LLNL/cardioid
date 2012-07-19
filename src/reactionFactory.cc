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
#include "TT06Func.hh"
#include "ReactionFHN.hh"
#include "NullReaction.hh"
#include "TestReaction.hh"
#include "ConstantReaction.hh"
#include "ThreadServer.hh"

#include <iostream>
using namespace std;

namespace
{
   Reaction* scanTT04_CellML(OBJECT* obj, const Anatomy& anatomy);
   Reaction* scanTT06_CellML(OBJECT* obj, const Anatomy& anatomy);
   Reaction* scanTT06Dev(OBJECT* obj, Anatomy& anatomy, const ThreadTeam& group);
   Reaction* scanTT06_RRG(OBJECT* obj, const Anatomy& anatomy);
   Reaction* scanFHN(OBJECT* obj, const Anatomy& anatomy);
   Reaction* scanNull(OBJECT* obj);
   Reaction* scanTest(OBJECT* obj);
   Reaction* scanConstant(OBJECT* obj, const Anatomy& anatomy);
}


Reaction* reactionFactory(const string& name, Anatomy& anatomy,
                          const ThreadTeam& group)
{
   OBJECT* obj = objectFind(name, "REACTION");
   string method; objectGet(obj, "method", method, "undefined");

   if (method == "TT04_CellML" || method == "tenTusscher04_CellML")
      return scanTT04_CellML(obj, anatomy);
   else if (method == "TT06_CellML" || method == "tenTusscher06_CellML")
      return scanTT06_CellML(obj, anatomy);
   else if (method == "TT06Dev" )
      return scanTT06Dev(obj, anatomy, group);
   else if (method == "TT06Opt" )
      return scanTT06Dev(obj, anatomy, group);
   else if (method == "TT06_RRG" )
      return scanTT06_RRG(obj, anatomy);
   else if (method == "FHN" || method == "FitzhughNagumo")
      return scanFHN(obj, anatomy);
   else if (method == "null")
      return scanNull(obj);
   else if (method == "test")
      return scanTest(obj);
   else if (method == "constant")
      return scanConstant(obj, anatomy);
   
   cerr<<"ERROR: Undefined reaction model in reactionFactory"<<endl;
   assert(false); // reachable only due to bad input
   return 0;
}

namespace
{
   Reaction* scanTT04_CellML(OBJECT* obj, const Anatomy& anatomy)
   {
      TT04_CellML_Reaction::IntegratorType integrator;
      string tmp;
      objectGet(obj, "integrator", tmp, "rushLarsen");
      if      (tmp == "rushLarsen")   integrator = TT04_CellML_Reaction::rushLarsen;
      else if (tmp == "rushLarson")   integrator = TT04_CellML_Reaction::rushLarsen;
      else if (tmp == "forwardEuler") integrator = TT04_CellML_Reaction::forwardEuler;
      else    assert(false);    
      
      return new TT04_CellML_Reaction(anatomy, integrator);
   }
}

namespace
{
   Reaction* scanTT06_CellML(OBJECT* obj, const Anatomy& anatomy)
   {
      TT06_CellML_Reaction::IntegratorType integrator;
      string tmp;
      objectGet(obj, "integrator", tmp, "rushLarsen");
      if      (tmp == "rushLarsen")   integrator = TT06_CellML_Reaction::rushLarsen;
      else if (tmp == "forwardEuler") integrator = TT06_CellML_Reaction::forwardEuler;
      else    assert(false);    
      return new TT06_CellML_Reaction(anatomy, integrator);
   }
}

namespace
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
   Reaction* scanTT06Dev(OBJECT* obj, Anatomy& anatomy, const ThreadTeam& group)
   {
      TT06Dev_ReactionParms parms;
      parms.cellTypeParms=TT06Func::getStandardCellTypes(); 
      int fastGate =-1; 
      int fastNonGate =-1; 
      TT06Func::initCnst(); 
      string fitFit; 
      objectGet(obj, "tolerance",    parms.tolerance, "0.0") ;
      objectGet(obj, "mod",          parms.mod, "0") ;
      objectGet(obj, "fastReaction", parms.fastReaction, "-1") ;
      objectGet(obj, "fastGate",     fastGate, "-1") ;
      objectGet(obj, "fastNonGate",  fastNonGate, "-1") ;
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
         if (fastGate     > -1 || fastNonGate > -1 )  
         {
            if (fastGate    ==  -1 )  fastGate   =0; 
            if (fastNonGate ==  -1 )  fastNonGate=0; 
            parms.fastReaction = fastGate+256*fastNonGate; 
         } 
      }
      objectGet(obj, "cellTypes", parms.cellTypeNames) ;
      if (parms.cellTypeNames.size() == 0)
      {
         parms.cellTypeNames.push_back("endoCellML");
         parms.cellTypeNames.push_back("midCellML");
         parms.cellTypeNames.push_back("epiCellML");
      }
      for (int ii=0;ii<parms.cellTypeNames.size();ii++) 
      {
         string name = parms.cellTypeNames[ii]; 
         if (parms.cellTypeParms.count(name) == 0)
            parms.cellTypeParms[name] = parms.cellTypeParms["endoCellML"]; 

         OBJECT* cellobj = object_find2(name.c_str(), "CELLTYPE", IGNORE_IF_NOT_FOUND);
         if (! cellobj)
            continue;
        
         string clone; 
         objectGet(cellobj, "clone", clone, "") ;
         if (clone != "")
         {
            assert(name != clone); 
            assert(parms.cellTypeParms.count(clone) != 0); 
            parms.cellTypeParms[name]=parms.cellTypeParms[clone]; 
            parms.cellTypeParms[name].name=name;
         }
         if (object_testforkeyword(cellobj, "anatomyIndices") ) objectGet(cellobj,"anatomyIndices",parms.cellTypeParms[name].anatomyIndices); 
         if (object_testforkeyword(cellobj, "s_switch")       ) objectGet(cellobj,"s_switch"      ,parms.cellTypeParms[name].s_switch,"0"); 
         if (object_testforkeyword(cellobj, "P_NaK")          ) objectGet(cellobj,"P_NaK"         ,parms.cellTypeParms[name].P_NaK,"0.0"); 
         if (object_testforkeyword(cellobj, "g_NaL")          ) objectGet(cellobj,"g_NaL"         ,parms.cellTypeParms[name].g_NaL,"0.0"); 
         if (object_testforkeyword(cellobj, "g_Ks")           ) objectGet(cellobj,"P_g_Ks"        ,parms.cellTypeParms[name].g_Ks,"0.0"); 
         if (object_testforkeyword(cellobj, "g_to")           ) objectGet(cellobj,"P_g_to"        ,parms.cellTypeParms[name].g_to,"0.0"); 
         
      }
      
      Reaction *reaction = new TT06Dev_Reaction(anatomy, parms, group);
      return  reaction; 
   }
}

namespace
{
   Reaction* scanTT06_RRG(OBJECT* obj, const Anatomy& anatomy)
   {
      Reaction *reaction = new TT06_RRG_Reaction(anatomy);
      return  reaction; 
   }
}

namespace
{
   Reaction* scanFHN(OBJECT* obj, const Anatomy& anatomy)
   {
      // None of the FHN model parameters are currently wired to the
      // input deck.
      return new ReactionFHN(anatomy);
   }
}

namespace
{
   Reaction* scanNull(OBJECT* obj)
   {
      NullReactionParms parms;
      objectGet(obj, "V0",    parms.initialVoltage, "-85",  "voltage");
      return new NullReaction(parms);
   }
}

namespace
{
   Reaction* scanTest(OBJECT* obj)
   {
      TestReactionParms parms;
      objectGet(obj, "V0",    parms.initialVoltage, "-85",  "voltage");
      objectGet(obj, "delta", parms.delta,          "1e-3", "voltage");
      return new TestReaction(parms);
   }
}

namespace
{
   Reaction* scanConstant(OBJECT* obj, const Anatomy& anatomy)
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

      return new ConstantReaction(anatomy,
                                  eta,
                                  sigma1, sigma2, sigma3,
                                  alpha, beta, gamma,
                                  printRate);
   }
}
