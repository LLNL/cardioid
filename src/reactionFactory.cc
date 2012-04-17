#include "reactionFactory.hh"

#include <cassert>
#include <cstdlib>
#include "object_cc.hh"

#include "TT04_CellML_Reaction.hh" // TT04 implementation from CellML (Nov 2011)
#include "TT06_CellML_Reaction.hh" // TT06 implementation from CellML (Nov 2011)
#include "TT06Dev_Reaction.hh"
#include "TT06_RRG_Reaction.hh"    // TT06 with modifications from Rice et al.
#include "ReactionFHN.hh"
#include "NullReaction.hh"
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
   else if (method == "TT06_RRG" )
      return scanTT06_RRG(obj, anatomy);
   else if (method == "FHN" || method == "FitzhughNagumo")
      return scanFHN(obj, anatomy);
   else if (method == "null")
      return new NullReaction();
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
   Reaction* scanTT06Dev(OBJECT* obj, Anatomy& anatomy, const ThreadTeam& group)
   {
#include "TT06Func.hh"
      map<string,TT06Func::CellTypeParmsFull>cellTypeParms=TT06Func::getStandardCellTypes(); 
      double tolerance=0.0 ; 
      int fastReaction =-1; 
      vector<string> cellTypeNames; 
      int mod=0; 
      objectGet(obj, "tolerance", tolerance, "0.0") ;
      objectGet(obj, "mod", mod, "0") ;
      objectGet(obj, "fastReaction", fastReaction, "-1") ;
      objectGet(obj, "cellTypes", cellTypeNames) ;
      if (cellTypeNames.size() == 0)
      {
         cellTypeNames.push_back("endoCellML");
         cellTypeNames.push_back("midCellML");
         cellTypeNames.push_back("epiCellML");
      }
      for (int ii=0;ii<cellTypeNames.size();ii++) 
      {
         string name = cellTypeNames[ii]; 
         if (cellTypeParms.count(name) == 0)
            cellTypeParms[name] = cellTypeParms["endoCellML"]; 

         OBJECT* cellobj = object_find2(name.c_str(), "CELLTYPE", IGNORE_IF_NOT_FOUND);
         if (! cellobj)
            continue;
        
         string clone; 
         objectGet(cellobj, "clone", clone, "") ;
         if (clone != "")
         {
            assert(name != clone); 
            assert(cellTypeParms.count(clone) != 0); 
            cellTypeParms[name]=cellTypeParms[clone]; 
            cellTypeParms[name].name=name;
         }
         vector<int> iparm; 
         vector<double> dparm; 
         objectGet(cellobj,"anatomyIndices",iparm); 
         if (iparm.size() > 0) cellTypeParms[name].anatomyIndices = iparm; 
         objectGet(cellobj,"s_switch",iparm); 
         if (iparm.size() == 1) cellTypeParms[name].s_switch = iparm[0]; 
         objectGet(cellobj,"P_NaK",dparm); 
         if (dparm.size() == 1) cellTypeParms[name].P_NaK = dparm[0]; 
         printf("%s P_NaK=%f s_switch=%d\n",name.c_str(),cellTypeParms[name].P_NaK,cellTypeParms[name].s_switch); 
         
      }
      
      Reaction *reaction = new TT06Dev_Reaction(anatomy, cellTypeParms, cellTypeNames, tolerance, mod, fastReaction, group);
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
