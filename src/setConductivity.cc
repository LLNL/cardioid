#include "setConductivity.hh"

#include <cassert>
#include <cmath>

#include "object_cc.hh"
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
                        const set<int>& typeSet,
                        vector<AnatomyCell> cell);
   void fibreConductivity(OBJECT* obj,
                          const BucketOfBits& data,
                          const set<int>& typeSet,
                          vector<AnatomyCell>& cell);
   void jhuConductivity(OBJECT* obj,
                        const Tuple& globalGridSize,
                        vector<AnatomyCell> cell);
   void uniformConductivity(OBJECT* obj,
                            vector<AnatomyCell>& cell);
}


void setConductivity(const string& name,
                     const BucketOfBits& data,
                     const Tuple& globalGridSize,
                     const set<int> typeSet,
                     vector<AnatomyCell>& cell)
{
   OBJECT* obj = 0;
   if (object_exists(name.c_str(), "CONDUCTIVITY"))
      obj = objectFind(name, "CONDUCTIVITY");
   string method = "pio";
   if (obj)
      objectGet(obj, "method", method, "pio");

   if (method == "pio")
      pioConductivity(data, typeSet, cell);
   else if (method == "fiber" || method == "fibre")
      fibreConductivity(obj, data, typeSet, cell);
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
                        const set<int>& typeSet,
                        vector<AnatomyCell> cell)
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
      
      assert (index11 != nFields &&
              index12 != nFields &&
              index13 != nFields &&
              index22 != nFields &&
              index23 != nFields &&
              index33 != nFields &&
              indexGid != nFields &&
              indexType != nFields);

      unsigned iCell = 0;
      for (unsigned ii=0; ii<data.nRecords(); ++ii)
      {
         BucketOfBits::Record rr = data.getRecord(ii);
         int cellType;
         rr.getValue(indexType, cellType);
         if (typeSet.count(cellType) == 0)
            continue;
         rr.getValue(index11, cell[iCell].sigma_.a11);
         rr.getValue(index12, cell[iCell].sigma_.a12);
         rr.getValue(index13, cell[iCell].sigma_.a13);
         rr.getValue(index22, cell[iCell].sigma_.a22);
         rr.getValue(index23, cell[iCell].sigma_.a23);
         rr.getValue(index33, cell[iCell].sigma_.a33);

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
                          const set<int>& typeSet,
                          vector<AnatomyCell>& cell)
   {
      FibreConductivityParms p;
      objectGet(obj, "sigmaTi", p.sigmaTi, "0.0315");
      objectGet(obj, "sigmaLi", p.sigmaLi, "0.3");
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
            if (typeSet.count(cellType) == 0)
               continue;
            
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

      objectGet(obj, "sigmaTi", p.sigmaTi, "0.0315");
      objectGet(obj, "sigmaLi", p.sigmaLi, "0.3");

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
      objectGet(obj, "sigma11", sigma.a11, "0.1");
      objectGet(obj, "sigma22", sigma.a22, "0.1");
      objectGet(obj, "sigma33", sigma.a33, "0.1");
      objectGet(obj, "sigma12", sigma.a12, "0.0");
      objectGet(obj, "sigma13", sigma.a13, "0.0");
      objectGet(obj, "sigma23", sigma.a23, "0.0");

      for (unsigned ii=0; ii<cell.size(); ++ii)
         cell[ii].sigma_ = sigma;
   }
}
