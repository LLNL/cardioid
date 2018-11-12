#include "StateVariableSensor.hh"
#include "Anatomy.hh"
#include "ioUtils.h"
#include "Simulate.hh"
#include "ReactionManager.hh"
#include "pio.h"
#include "BoundingBox.hh"
#include "TupleToIndex.hh"
#include "IndexToTuple.hh"
#include "readCellList.hh"

#include "stringUtils.hh"

#include <mpi.h>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

// maps names to varHandle and units
typedef map<string, pair<int, string> > FieldMap;

namespace
{
   set<Long64> mkCellSet(const Anatomy& anatomy,
                         const vector<Long64>& cells,
                         string cellListFilename,
                         double radius);
   FieldMap mkFieldMap(const ReactionManager* reaction);
}


StateVariableSensor::StateVariableSensor(
   const SensorParms& sp, const StateVariableSensorParms& p, 
   const Simulate& sim)
    : Sensor(sp),
      binaryOutput_(p.binaryOutput),
      filename_(p.filename),
      sim_(sim)
{

   gidFormat_ = "%12llu ";
   varFormat_ = " %21.13e";
   int gidRecLen = 13;
   int varRecLen = 22;
   
   //    processFieldList(p.fieldList);

   FieldMap availableFields = mkFieldMap(sim_.reaction_);
   vector<string> fieldNames;
   vector<string> fieldUnits;

   if (p.allFields)
   {
      for (FieldMap::iterator iter = availableFields.begin(); iter != availableFields.end(); iter++)
      {
         fieldNames.push_back(iter->first);
         handles_.push_back(iter->second.first);
         fieldUnits.push_back(iter->second.second);
      }         
   }
   else
   {
      assert(p.fieldList.size() > 0);
      for (unsigned ii=0; ii<p.fieldList.size(); ++ii)
      {
         FieldMap::iterator iter = availableFields.find(p.fieldList[ii]);
         if (iter == availableFields.end())
            continue; // silently ignore typos.
         fieldNames.push_back(iter->first);
         handles_.push_back(iter->second.first);
         fieldUnits.push_back(iter->second.second);
      }
   }
      
   set<Long64> requestedCells = mkCellSet(
      sim.anatomy_, p.cells, p.cellListFilename, p.radius);

   for (unsigned ii=0; ii<sim.anatomy_.nLocal(); ++ii)
   {
      Long64 gid = sim.anatomy_.gid(ii) ;
      if (requestedCells.count(gid) > 0 || p.allCells)
         localCells_[gid] = ii;
   }

   int localRecords = localCells_.size();
   int nRecords;
   MPI_Allreduce(&localRecords, &nRecords, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   
   lRec_ = gidRecLen + handles_.size()*varRecLen +1;
   if (binaryOutput_)
      lRec_ = 8*(1+handles_.size());
   
   header_.objectName_ = "stateVariable";
   header_.className_ = "FILEHEADER";
   header_.dataType_ = PioHeaderData::ASCII;
   header_.nRecords_ = nRecords;
   header_.lRec_ = lRec_;
   header_.nFields_ = 1 + handles_.size();
   header_.fieldNames_ = "gid " + concat(fieldNames);
   header_.fieldTypes_ = "u " + concat(vector<string>(handles_.size(), "f"));
   header_.fieldUnits_ = "1 " + concat(fieldUnits);
   header_.addItem("nx", sim.anatomy_.nx());
   header_.addItem("ny", sim.anatomy_.ny());
   header_.addItem("nz", sim.anatomy_.nz());

   if (binaryOutput_)
   {
      header_.dataType_ = PioHeaderData::BINARY;
      header_.fieldTypes_ = "u8 " + concat(vector<string>(handles_.size(), "f8"));
   }
   
   
}

StateVariableSensor::~StateVariableSensor()
{
}

/** It can be argued that looping over the handles_ and getting each
 *  value individually incurs more branching and function call overhead
 *  than necessary.  However I doubt it will be an issue, and it allows
 *  us to write the data in the file in the same order that the user
 *  specified the fields.  We can always change it if the performance is
 *  an issue.
 */
void StateVariableSensor::print(double time, int loop)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   int myRank;
   MPI_Comm_rank(comm, &myRank);

   stringstream name;
   name << "snapshot."<<setfill('0')<<setw(12)<<loop;
   string dirname = name.str();
   if (myRank == 0)
      DirTestCreate(dirname.c_str());
   MPI_Barrier(comm); // none shall pass before task 0 creates directory
   string pFilename = dirname + "/" + filename_;
   PFILE* file = Popen(pFilename.c_str(), "w", comm);
   if (myRank == 0)
      header_.writeHeader(file, loop, time);

   vector<double> values(handles_.size());
   char buf[lRec_+1];
   for (MapType::const_iterator iter = localCells_.begin();
        iter != localCells_.end(); ++iter)
   {
      for (unsigned ii=0; ii<handles_.size(); ++ii)
         if (handles_[ii] >= 0)
            values[ii] = sim_.reaction_->getValue(iter->second, handles_[ii]);
         else
            values[ii] = getSimValue(iter->second, handles_[ii]);

      if (binaryOutput_)
      {
         copyBytes(buf, &iter->first, 8);
         int bufPos = 8;
         for (unsigned ii=0; ii<values.size(); ++ii)
         {
            copyBytes(buf+bufPos, &values[ii], 8);
            bufPos += 8;
         }
      }
      else
      {
         int bufPos = sprintf(buf, gidFormat_, iter->first);
         for (unsigned ii=0; ii<values.size(); ++ii)
            bufPos += sprintf(buf+bufPos, varFormat_, values[ii]);
         sprintf(buf+bufPos, "\n");
      }
      Pwrite(buf, lRec_, 1, file);
   }
   Pclose(file);
}

double StateVariableSensor::getSimValue(int iCell, int varHandle)
{
   double value;
   switch (varHandle)
   {
     case -1:
      {
         ro_array_ptr<double> VmArray = sim_.vdata_.VmTransport_.useOn(CPU);
         value = VmArray[iCell];
      }
      break;
     case -2:
      {
         ro_array_ptr<double> dVmDiffusion = sim_.vdata_.dVmDiffusionTransport_.useOn(CPU);
         value = dVmDiffusion[iCell];
      }
      break;
     case -3:
      {
         ro_array_ptr<double> dVmReaction = sim_.vdata_.dVmReactionTransport_.useOn(CPU);
         value = dVmReaction[iCell];
      }
      break;
     default:
      assert(false);
   }
   return value;
}


/** Construct a locally relevant set of cell gids for which the user has
 *  requested output.  Locally relevant means that the set contains only
 *  gid that might possibly exist on this task.
 *
 *  Algorithm:
 *  - Form a vector, requestedCells, that contains the union of the
 *    gids that the user specified in the input file, and those that
 *    were read from the file specified by cellListFilename.
 *  - Form the bounding box around the local cells, then push it out by
 *    radius in each direction.
 *  - Discard any requested cell that is not contained within this
 *    bounding box.  Such requests can't be relevant on this task.
 *    The cells that pass this screening are stored in screenedTuples.
 *  - Form a vector of all Tuple displacements, inRange, that are
 *    shorter than radius.  At a minimum this includes (0, 0, 0).
 *  - Iterate the screenedTuples.  Form a set that includes every cell
 *    within the radius of one of the screenedTuples by adding all of
 *    the inRange Tuple displacements.  This is the locally relevant set.
*/
namespace
{
   set<Long64> mkCellSet(const Anatomy& anatomy,
                         const vector<Long64>& cells,
                         string cellListFilename,
                         double radius)
   {
      // Union of user specification
      vector<Long64> requestedCells;
      if (!cellListFilename.empty())
         readCellList(cellListFilename, requestedCells);
      requestedCells.insert(requestedCells.end(), cells.begin(), cells.end());
   
      // Form bounding box
      Tuple minCorner = anatomy.globalTuple(0);
      Tuple maxCorner = anatomy.globalTuple(0);
      for (unsigned ii=1; ii<anatomy.nLocal(); ++ii)
      {
         Tuple tt = anatomy.globalTuple(ii);
         minCorner.x() = min(tt.x(), minCorner.x());
         minCorner.y() = min(tt.y(), minCorner.y());
         minCorner.z() = min(tt.z(), minCorner.z());
         maxCorner.x() = max(tt.x(), maxCorner.x());
         maxCorner.y() = max(tt.y(), maxCorner.y());
         maxCorner.z() = max(tt.z(), maxCorner.z());
      }
      int intRadius = ceil(radius);
      Tuple radiusTuple(intRadius, intRadius, intRadius);
      minCorner -= radiusTuple;
      maxCorner += radiusTuple;
      BoundingBox bb(minCorner, maxCorner);
   
      // screen
      IndexToTuple indexToTuple(anatomy.nx(), anatomy.ny(), anatomy.nz());
      vector<Tuple> screenedTuples;
      for (unsigned ii=0; ii<requestedCells.size(); ++ii)
      {
         Tuple tt = indexToTuple(requestedCells[ii]);
         if (  bb.contains(tt) )
            screenedTuples.push_back(tt);
      }
   
      // Form tuple displacements that are in range
      double radiusSq = radius*radius;
      vector<Tuple> inRange;
      for (int ix=-intRadius; ix<=intRadius; ++ix)
         for (int iy=-intRadius; iy<=intRadius; ++iy)
            for (int iz=-intRadius; iz<=intRadius; ++iz)
               if ( (ix*ix + iy*iy + iz*iz) <= radiusSq)
                  inRange.push_back(Tuple(ix, iy, iz));
      assert(inRange.size() > 0);
   
      // Add all inRange displacements to all screened Tuples.
      TupleToIndex tupleToIndex(anatomy.nx(), anatomy.ny(), anatomy.nz());
      set<Long64> locallyRelevantSet;
      for (unsigned ii=0; ii<screenedTuples.size(); ++ii)
         for (unsigned jj=0; jj<inRange.size(); ++jj)
         {
            Tuple tt = screenedTuples[ii] + inRange[jj];
            locallyRelevantSet.insert(tupleToIndex(tt));
         }
   
      return locallyRelevantSet;
   }
}


namespace
{
   FieldMap mkFieldMap(const ReactionManager* reaction)
   {
      vector<string> reactionFields;
      vector<string> reactionFieldUnits;
      reaction->getCheckpointInfo(reactionFields, reactionFieldUnits);
      FieldMap availableFields;
      for (unsigned ii=0; ii<reactionFields.size(); ++ii)
      {
         const string& iName = reactionFields[ii];
         availableFields[iName] =
            make_pair(reaction->getVarHandle(iName), reactionFieldUnits[ii]);
      }
      availableFields["Vm"]   = make_pair(-1, "mV");
      availableFields["dVmD"] = make_pair(-2, "mV/fs");
      availableFields["dVmR"] = make_pair(-3, "mV/fs");
      return availableFields;
   }
}
