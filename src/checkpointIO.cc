#include "checkpointIO.hh"

#include <cassert>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include "pio.h"
#include "ioUtils.h"
#include "Simulate.hh"
#include "Anatomy.hh"
#include "ReactionManager.hh"
#include "Long64.hh"
#include "BucketOfBits.hh"
#include "stateLoader.hh"
#include "units.h"
#include "Version.hh"
#include "utilities.h"
#include "stringUtils.hh"
#include <cstring>

using namespace std;


namespace
{
   /** One stop shopping for all of the data that needs to go into the
    *  header of a checkpoint file.  Items such as nfiles and endian_key
    *  that are entirely under the control of the pio module are
    *  intentionally omitted. */
   struct CheckpointHeaderData
   {
      enum DataType {ASCII, BINARY};

      string simulateName_;
      string stateFileName_;
      DataType dataType_;
      Long64 nRecord_;
      unsigned lRec_;
      double time_;
      int loop_;
      std::string reactionMethod_;
      unsigned nFields_;
      std::string fieldNames_;
      std::string fieldTypes_;
      std::string fieldUnits_;
      int nx_;
      int ny_;
      int nz_;
      std::string comment_;
   };
}

namespace
{
   void writeHeader(const CheckpointHeaderData& headerData, PFILE* file)
   {
      string dataType;
      switch (headerData.dataType_)
      {
        case CheckpointHeaderData::ASCII:
         dataType = "FIXRECORDASCII";
         break;
        case CheckpointHeaderData::BINARY:
         dataType = "FIXRECORDBINARY";
         break;
        default:
         assert(false);
      }
      int endianKey;
      memcpy(&endianKey, "1234", 4);
      int nFiles = file->ngroup;
      
      const Version& vv = Version::getInstance();
      
      Pprintf(file, "state FILEHEADER {\n");
      Pprintf(file, "   create_time = %s;\n",  timestamp_string());
      Pprintf(file, "   user = %s;",           vv.user().c_str());
      Pprintf(file, "   host = %s;\n",         vv.host().c_str());
      Pprintf(file, "   exe_version = %s;",    vv.version().c_str());
      Pprintf(file, "   srcpath = %s;\n",      vv.srcPath().c_str());
      Pprintf(file, "   compile_time = %s;\n", vv.compileTime().c_str());
      Pprintf(file, "   datatype = %s;\n", dataType.c_str());
      Pprintf(file, "   nfiles = %d;\n", nFiles);
      Pprintf(file, "   nrecord = %d;\n", headerData.nRecord_);
      Pprintf(file, "   lrec = %d;\n", headerData.lRec_);
      Pprintf(file, "   endian_key = %d;\n", endianKey);
      Pprintf(file, "   time = %f;\n", headerData.time_);
      Pprintf(file, "   loop = %d;\n", headerData.loop_);
      Pprintf(file, "   reactionMethod = %s;\n", headerData.reactionMethod_.c_str());
      Pprintf(file, "   nfields = %d;\n", headerData.nFields_);
      Pprintf(file, "   field_names = %s;\n", headerData.fieldNames_.c_str());
      Pprintf(file, "   field_types = %s;\n", headerData.fieldTypes_.c_str());
      Pprintf(file, "   field_units = %s;\n", headerData.fieldUnits_.c_str());
      Pprintf(file, "   nx = %u; ny = %u; nz = %u;\n",
              headerData.nx_, headerData.ny_, headerData.nz_);
      Pprintf(file, "}\n\n");
   }
}

namespace
{
   void writeRestart(const CheckpointHeaderData& headerData, const string& dirName)
   {
      string filename = dirName + "/restart";
      FILE* file = fopen(filename.c_str(), "w");
      fprintf(file, "%s SIMULATE {\n"
              "   loop=%d; time=%f; stateFile=%s;\n"
              "}\n",
              headerData.simulateName_.c_str(),
              headerData.loop_,
              units_convert(headerData.time_, NULL, "t"),
              headerData.stateFileName_.c_str());
      fflush(file); 
      fclose(file);
   }
}


void writeCheckpoint(const Simulate& sim, MPI_Comm comm)
{
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   const Anatomy& anatomy = sim.anatomy_;

   stringstream name;
   name << "snapshot."<<setfill('0')<<setw(12)<<sim.loop_;
   string dirName = name.str();
   if (myRank == 0)
      DirTestCreate(dirName.c_str());
   
   vector<string> fieldNames;
   vector<string> fieldUnits;
   sim.reaction_->getCheckpointInfo(fieldNames, fieldUnits);
   vector<int> handle = sim.reaction_->getVarHandle(fieldNames);

   const char* gidVmFormat = "%12llu %21.13e";
   const char* itemFormat = " %21.13e";
   int lRec = 34 + 22*fieldNames.size() + 1;
   
   CheckpointHeaderData headerData;
   headerData.stateFileName_ = dirName + "/state";
   headerData.simulateName_ = sim.name_;
   headerData.dataType_ = CheckpointHeaderData::ASCII;
   headerData.nRecord_ = anatomy.nGlobal();
   headerData.lRec_ = lRec;
   headerData.loop_ = sim.loop_;
   headerData.time_ = sim.time_;
   headerData.reactionMethod_ = sim.reaction_->stateDescription();
   headerData.nFields_ = 2+ fieldNames.size();
   headerData.fieldNames_ = "gid Vm " + concat(fieldNames);
   headerData.fieldTypes_ = "u f " + concat(vector<string>(fieldNames.size(), "f"));
   headerData.fieldUnits_ = "1 mV " + concat(fieldUnits);
   headerData.nx_ = anatomy.nx();
   headerData.ny_ = anatomy.ny();
   headerData.nz_ = anatomy.nz();

   // header was just setup for ASCII checkpoints.  If user asked for
   // BINARY we need a couple of tweaks
   if (!sim.asciiCheckpoints_)
   {
      headerData.dataType_ = CheckpointHeaderData::BINARY;
      lRec = 8 * (fieldNames.size() + 2);
      headerData.lRec_ = lRec;
      headerData.fieldTypes_ = "u8 f8 " + concat(vector<string>(fieldNames.size(), "f8"));
   }

   PFILE* file = Popen(headerData.stateFileName_.c_str(), "w", comm);
   if (myRank == 0)
   {
      writeHeader(headerData, file);
      writeRestart(headerData, dirName);
   }
   
   char buf[lRec+1];
   vector<double> value(handle.size(), 0.0);
   ro_array_ptr<double> vmarray = sim.vdata_.VmTransport_.useOn(CPU);
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      if (sim.asciiCheckpoints_)
      {
         int bufPos = sprintf(buf, gidVmFormat, anatomy.gid(ii), vmarray[ii]);
         sim.reaction_->getValue(ii, handle, value);
         for (unsigned jj=0; jj<value.size(); ++jj)
            bufPos += sprintf(buf+bufPos, itemFormat, value[jj]);
         sprintf(buf+bufPos, "\n");
      }
      else
      {
         Long64 gid = anatomy.gid(ii);
         copyBytes(buf, &gid, 8);
         copyBytes(buf+8, &vmarray[ii], 8);
         sim.reaction_->getValue(ii, handle, value);
         for (unsigned jj=0; jj<value.size(); ++jj)
            copyBytes(buf+16+jj*8, &(value[0])+jj, 8);
      }
      Pwrite(buf, lRec, 1, file);
   }
   int rc = Pclose(file);
   if (rc == 0) 
   {
      unlink("restart");
      string restartName = dirName + "/restart";
      symlink(restartName.c_str(), "restart");
   }
}

void readCheckpoint(const string& filename, Simulate& sim, MPI_Comm comm)
{
   BucketOfBits* data = 
      loadAndDistributeState(filename, sim.anatomy_);
   assert(data->nRecords() == sim.anatomy_.nLocal());

   vector<double> unitConvert(data->nFields(), 1.0);
   typedef map<int, int> FieldMap;
   FieldMap fieldMap;

   for (unsigned ii=0; ii<data->nFields(); ++ii)
   {
      int handle = sim.reaction_->getVarHandle(data->fieldName(ii));
      if (handle >= 0 ) // handle < 0 -> undefined name
      {
         fieldMap[ii] = handle;
         string from = data->units(ii);
         string to = sim.reaction_->getUnit(data->fieldName(ii));
         unitConvert[ii] = units_convert(1.0, from.c_str(), to.c_str());
      }
   }
   
   for (unsigned ii=0; ii<sim.anatomy_.nLocal(); ++ii)
   {
      BucketOfBits::Record iRec = data->getRecord(ii);
      for (FieldMap::const_iterator iter=fieldMap.begin();
           iter!=fieldMap.end(); ++iter)
      {
         int iField = iter->first;
         int handle = iter->second;
         double value;
         switch (data->dataType(iField))
         {
           case BucketOfBits::floatType:
           case BucketOfBits::f8Type:
           case BucketOfBits::f4Type:
            iRec.getValue(iField, value);
            break;
           case BucketOfBits::intType:
           case BucketOfBits::u8Type:
            int tmp;
            iRec.getValue(iField, tmp);
            value = double(tmp);
            break;
           default:
            assert(false);
         }
         value *= unitConvert[iField];
         sim.reaction_->setValue(ii, handle, value);
      }
   }

   // Load membrane voltage from checkpoint file into VmArray.
   wo_array_ptr<double> vmarray = sim.vdata_.VmTransport_.useOn(CPU); 
   unsigned vmIndex = data->getIndex("Vm");
   if (vmIndex != data->nFields())
      for (unsigned ii=0; ii<data->nRecords(); ++ii)
         data->getRecord(ii).getValue(vmIndex, vmarray.raw()[ii]);
   delete data;
}
