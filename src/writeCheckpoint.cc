#include "writeCheckpoint.hh"

#include <cassert>
#include <vector>

#include "pio.h"
#include "Simulate.hh"
#include "Anatomy.hh"
#include "Reaction.hh"
#include "Long64.hh"
#include "units.h"

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
      
      
      Pprintf(file, "state FILEHEADER {\n");
      Pprintf(file, "   datatype = %s;\n", dataType.c_str());
      Pprintf(file, "   nfiles = %d;\n", nFiles);
      Pprintf(file, "   nrecord = %d;\n", headerData.nRecord_);
      Pprintf(file, "   lrec = %d;\n", headerData.lRec_);
      Pprintf(file, "   endidan_key = %d;\n", endianKey);
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





namespace
{
   /** Concatenates the strings in vv into a single string with a space
    * separating each element.  No space is added to the beginning or
    * end. */
   string concat(const vector<string> vv)
   {
      if (vv.size() == 0)
         return "";
      string tmp = vv[0];
      for (unsigned ii=1; ii<vv.size(); ++ii)
         tmp += " " + vv[ii];
      return tmp;
   }
}


   
void writeCheckpoint(const Simulate& sim, const string& dirName, MPI_Comm comm)
{
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   const Anatomy& anatomy = sim.anatomy_;

   vector<string> fieldNames;
   vector<string> fieldUnits;
   sim.reaction_->getCheckpointInfo(fieldNames, fieldUnits);
   vector<int> handle = sim.reaction_->getHandle(fieldNames);

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
   headerData.reactionMethod_ = sim.reaction_->methodName();
   headerData.nFields_ = 2+ fieldNames.size();
   headerData.fieldNames_ = "gid Vm " + concat(fieldNames);
   headerData.fieldTypes_ = "u f " + concat(vector<string>(fieldNames.size(), "f"));
   headerData.fieldUnits_ = "1 mV/ms " + concat(fieldUnits);
   headerData.nx_ = anatomy.nx();
   headerData.ny_ = anatomy.ny();
   headerData.nz_ = anatomy.nz();
//   headerData.comment_ = ""; assert(false);
   
      
   
   
   

   PFILE* file = Popen(headerData.stateFileName_.c_str(), "w", comm);
   if (myRank == 0)
   {
      writeHeader(headerData, file);
      writeRestart(headerData, dirName);
   }
   
   
   char buf[lRec+1];
   vector<double> value(handle.size(), 0.0);
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      int bufPos = sprintf(buf, gidVmFormat, anatomy.gid(ii), sim.VmArray_[ii]);
      sim.reaction_->getValue(ii, handle, value);
      for (unsigned jj=0; jj<value.size(); ++jj)
         bufPos += sprintf(buf+bufPos, itemFormat, value[jj]);
      sprintf(buf+bufPos, "\n");
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

