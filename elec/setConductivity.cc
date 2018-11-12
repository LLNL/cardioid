#include "setConductivity.hh"

#include <cassert>
#include <cmath>
#include <iostream>

#include "object_cc.hh"
#include "units.h"
#include "mpiUtils.h"
#include "BucketOfBits.hh"
#include "Tuple.hh"
#include "AnatomyCell.hh"

#include "FibreConductivity.hh"
#include "JHUConductivity.hh"

using namespace std;

/**
 *   1. method = pio (default)
 *      Requires the BucketOfBits to contain tensor conductivity.  Each
 *      of the six required fields must be defined.
 *   2. method = fibre (or fiber)
 *      Theta and Phi are obtained from the BucketOfBits, or, if the
 *      Bucket has no data are assumed to be zero for all cells
 *   3. method = uniform.
 *      This one is easy.  Values come from object.data
 *   4. method = JHU
 *      Conductivity is calculated from the global coordinates of each
 *      cell.
 */
namespace
{
   void pioConductivity(const BucketOfBits& data,
                        vector<AnatomyCell>& cell);
   void fibreConductivity(OBJECT* obj,
                          const BucketOfBits& data,
                          vector<AnatomyCell>& cell);
   void jhuConductivity(OBJECT* obj,
                        const Tuple& globalGridSize,
                        vector<AnatomyCell> cell);
   void uniformConductivity(OBJECT* obj,
                            vector<AnatomyCell>& cell);

   void unitConsistencyError();
   void missingFieldError();
}


void setConductivity(const string& name,
                     const BucketOfBits& data,
                     const Tuple& globalGridSize,
                     vector<AnatomyCell>& cell)
{
   OBJECT* obj = 0;
   if (object_exists(name.c_str(), "CONDUCTIVITY"))
      obj = objectFind(name, "CONDUCTIVITY");
   string method = "pio";
   if (obj)
      objectGet(obj, "method", method, "pio");

   if (method == "pio")
      pioConductivity(data, cell);
   else if (method == "fiber" || method == "fibre")
      fibreConductivity(obj, data, cell);
   else if (method == "JHU")
      jhuConductivity(obj, globalGridSize, cell);
   else if (method == "uniform")
      uniformConductivity(obj, cell);
   else
      assert(false);
}

namespace
{
   void pioConductivity(const BucketOfBits& data,
                        vector<AnatomyCell>& cell)
   {
      unsigned nFields = data.nFields();
      
      unsigned index11 = data.getIndex("sigma11");
      unsigned index12 = data.getIndex("sigma12");
      unsigned index13 = data.getIndex("sigma13");
      unsigned index22 = data.getIndex("sigma22");
      unsigned index23 = data.getIndex("sigma23");
      unsigned index33 = data.getIndex("sigma33");
      unsigned indexGid = data.getIndex("gid");
      unsigned indexType = data.getIndex("cellType");
      
      if (index11 == nFields || index12 == nFields ||
          index13 == nFields || index22 == nFields ||
          index23 == nFields || index33 == nFields ||
          indexGid == nFields || indexType == nFields)
         missingFieldError();

      // get factor to convert input resisitivities to internal units.
      const char* to = "resistivity_internal/length_internal";
      if ( units_check(data.units(index11).c_str(), to) == 1 ||
           units_check(data.units(index12).c_str(), to) == 1 ||
           units_check(data.units(index13).c_str(), to) == 1 ||
           units_check(data.units(index22).c_str(), to) == 1 ||
           units_check(data.units(index23).c_str(), to) == 1 ||
           units_check(data.units(index33).c_str(), to) == 1 )
         unitConsistencyError();
      double uc11 = units_convert(1.0, data.units(index11).c_str(), to);
      double uc12 = units_convert(1.0, data.units(index12).c_str(), to);
      double uc13 = units_convert(1.0, data.units(index13).c_str(), to);
      double uc22 = units_convert(1.0, data.units(index22).c_str(), to);
      double uc23 = units_convert(1.0, data.units(index23).c_str(), to);
      double uc33 = units_convert(1.0, data.units(index33).c_str(), to);
      
      unsigned iCell = 0;
      for (unsigned ii=0; ii<data.nRecords(); ++ii)
      {
         BucketOfBits::Record rr = data.getRecord(ii);
         int cellType;
         rr.getValue(indexType, cellType);
         rr.getValue(index11, cell[iCell].sigma_.a11);
         rr.getValue(index12, cell[iCell].sigma_.a12);
         rr.getValue(index13, cell[iCell].sigma_.a13);
         rr.getValue(index22, cell[iCell].sigma_.a22);
         rr.getValue(index23, cell[iCell].sigma_.a23);
         rr.getValue(index33, cell[iCell].sigma_.a33);

         cell[iCell].sigma_.a11 *= uc11;
         cell[iCell].sigma_.a12 *= uc12;
         cell[iCell].sigma_.a13 *= uc13;
         cell[iCell].sigma_.a22 *= uc22;
         cell[iCell].sigma_.a23 *= uc23;
         cell[iCell].sigma_.a33 *= uc33;
         
         Long64 gid;
         rr.getValue(indexGid, gid);
         assert(gid == cell[iCell].gid_);
         ++iCell;
      }
   }
}

namespace
{
   void fibreConductivity(OBJECT* obj,
                          const BucketOfBits& data,
                          vector<AnatomyCell>& cell)
   {
      FibreConductivityParms p;
      objectGet(obj, "sigmaTi", p.sigmaTi, "0.0315e-3", "resistivity/l");
      objectGet(obj, "sigmaLi", p.sigmaLi, "0.3e-3",    "resistivity/l");
      FibreConductivity fibre(p);

      vector<double> theta(cell.size()); // angles default to zero
      vector<double> phi(cell.size());
      
      unsigned nFields = data.nFields();
      unsigned thetaIndex = data.getIndex("theta");
      unsigned phiIndex =   data.getIndex("phi");
      unsigned gidIndex =   data.getIndex("gid");
      unsigned typeIndex =  data.getIndex("cellType");
      if (thetaIndex != nFields || phiIndex != nFields)
      {
         // Extract data read from pio file
         assert(thetaIndex != nFields &&
                phiIndex != nFields &&
                gidIndex != nFields &&
                typeIndex != nFields);

         BucketOfBits::DataType dataType = data.dataType(thetaIndex);
         assert (data.dataType(phiIndex) == dataType);

         unsigned iCell=0;
         for (unsigned ii=0; ii<data.nRecords(); ++ii)
         {
            BucketOfBits::Record rr = data.getRecord(ii);
            int cellType;
            rr.getValue(typeIndex, cellType);
            
            Long64 gid;
            rr.getValue(gidIndex, gid);
            assert(gid == cell[iCell].gid_);
            if (dataType == BucketOfBits::intType)
            {
               int tmp;
               rr.getValue(thetaIndex, tmp);
               theta[iCell] = tmp*(M_PI/256);
               rr.getValue(phiIndex, tmp);
               phi[iCell] = tmp*(M_PI/256);
            }
            else
            {
               rr.getValue(thetaIndex, theta[iCell]);
               rr.getValue(phiIndex, phi[iCell]);
            }
            ++iCell;
         }
      }
      
      for (unsigned ii=0; ii<cell.size(); ++ii)
         fibre.compute(theta[ii], phi[ii], cell[ii].sigma_);
   }
}


namespace
{
   /** See JHUConductivity.hh for an explanation of the parameters and
    *  their default values. */
   void jhuConductivity(OBJECT* obj,
                        const Tuple& globalGridSize,
                        vector<AnatomyCell> cell)
   {
      JHUConductivityParms p;
      p.nx = globalGridSize.x();
      p.ny = globalGridSize.y();
      p.nz = globalGridSize.z();
      p.transmuralAxis = globalGridSize.y() - 2;

      objectGet(obj, "sigmaTi", p.sigmaTi, "0.0315e-3", "resistivity/l");
      objectGet(obj, "sigmaLi", p.sigmaLi, "0.3e-3",    "resistivity/l");

      objectGet(obj, "sheetAngle", p.sheetAngle, "45");
      objectGet(obj, "rotatationMatrix", p.rotatationMatrix, "1");
      objectGet(obj, "homogeneousFiber", p.homogeneousFiber, "0");

      JHUConductivity jhu(p);

      for (unsigned ii=0; ii<cell.size(); ++ii)
         jhu.compute(cell[ii]);

   }
}

namespace 
{
   void uniformConductivity(OBJECT* obj, vector<AnatomyCell>& cell)
   {
      SymmetricTensor sigma;
      objectGet(obj, "sigma11", sigma.a11, "0.1e-3", "resistivity/l");
      objectGet(obj, "sigma22", sigma.a22, "0.1e-3", "resistivity/l");
      objectGet(obj, "sigma33", sigma.a33, "0.1e-3", "resistivity/l");
      objectGet(obj, "sigma12", sigma.a12, "0.0e-3", "resistivity/l");
      objectGet(obj, "sigma13", sigma.a13, "0.0e-3", "resistivity/l");
      objectGet(obj, "sigma23", sigma.a23, "0.0e-3", "resistivity/l");

      for (unsigned ii=0; ii<cell.size(); ++ii)
         cell[ii].sigma_ = sigma;
   }
}

namespace
{
   void unitConsistencyError()
   {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == 0)
         cout << "Fatal Error in anatomy file.\n"
            "Anatomy files that contain conductivity tensor data must specify\n"
            "the field_units keyword in the header.  Units for conductivity\n"
            "must be resistivity per length." << endl;
      abortAll(1);
   }
}

namespace
{
   void missingFieldError()
   {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == 0)
         cout << "Fatal Error in anatomy file.\n"
            "You have chosen to read conductivty tensor data from the anatomy\n"
            "file.  However at least one component of the tensor is missing.\n"
            "Your anatomy file must contain sigma11, sigma12, sigma13, \n"
            "sigma22, sigma23, sigma33, gid, and cellType." << endl;
      abortAll(1);
   }
}
