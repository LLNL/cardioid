#include "CaAverageSensor.hh"

#include "PerformanceTimers.hh"
#include "Anatomy.hh"
#include "pio.h"
#include "ioUtils.h"
#include "ReactionManager.hh"
#include "CommTable.hh"
#include "stringUtils.hh"

#include <vector>
#include <iomanip>
#include <sstream>

using namespace std;
using namespace PerformanceTimers;

/////////////////////////////////////////////////////////////////////

CaAverageSensor::CaAverageSensor(const SensorParms& sp,
                                 string filename,
                                 unsigned nFiles,
                                 const Anatomy& anatomy,
                                 vector<Long64>& gid,
                                 const ReactionManager& reaction,
                                 const CommTable* commtable)
   :Sensor(sp),
    coarsening_(anatomy, gid, 10000.0, commtable),
    filename_(filename),
    nFiles_(nFiles),
    reaction_(reaction),
    comm_(commtable->_comm),
    loop_buffer_(-1)
{
   
   nx_=anatomy.nx();
   ny_=anatomy.ny();
   nz_=anatomy.nz();
   
   nlocal_=anatomy.nLocal();

   buffer_val_.resize(nlocal_);

   ca_handle_=reaction_.getVarHandle("Ca_i");
}

void CaAverageSensor::computeColorAverages(ro_array_ptr<double> val)
{
   // calculate local sums
   coarsening_.accumulateValues(val,avg_valcolors_);
   
   coarsening_.exchangeAndSum(avg_valcolors_);
}

void CaAverageSensor::writeAverages(const string& filename,
                             const double current_time,
                             const int current_loop)const
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);

   PFILE* file = Popen(filename.c_str(), "w", comm_);
   if (nFiles_ > 0)
     PioSet(file, "ngroup", nFiles_);
   
   const int nfields = 4+(int)times_.size(); 
   const int lrec    = 20+13*(int)times_.size() + 1;

   const std::set<int>& owned_colors=coarsening_.getOwnedColors();

   static Long64 nSnapSub = -1;
   static bool first_time = true;
   if( first_time )
   {
      Long64 nSnapSubLoc = owned_colors.size();
      MPI_Allreduce(&nSnapSubLoc, &nSnapSub, 1, MPI_LONG_LONG, MPI_SUM, comm_);
   }

   if (myRank == 0)
   {
      // write header
      int nfiles;
      Pget(file,"ngroup",&nfiles);
      Pprintf(file, "coarsenedCa FILEHEADER {\n");
      Pprintf(file, "  lrec = %d;\n", lrec);
      Pprintf(file, "  datatype = FIXRECORDASCII;\n");
      Pprintf(file, "  nrecords = %llu;\n", nSnapSub);
      Pprintf(file, "  nfields = %d;\n", nfields);
      string fieldNames="gid nvals " + concat(vector<string>(times_.size(), "avgCa"));
      Pprintf(file, "  field_names = %s;\n", fieldNames.c_str());
      string fieldTypes="u d " + concat(vector<string>(times_.size(), "f"));
      Pprintf(file, "  field_types = %s;\n", fieldTypes.c_str());
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  time = %f; loop = %u;\n",times_[0], current_loop);
      if( times_.size()>1 )
         Pprintf(file, "  nsteps = %d; dt = %f\n", times_.size(), times_[1]-times_[0]);
      Pprintf(file, "  h = %4u  0    0\n", nx_);
      Pprintf(file, "        0    %4u  0\n", ny_);
      Pprintf(file, "        0    0    %4u;\n", nz_);
      Pprintf(file, "}\n\n");
   }
   
   for(set<int>::const_iterator it = owned_colors.begin();
                                it!= owned_colors.end();
                              ++it)
   {
      const int color=(*it);
      
      const map< int, vector<float> >::const_iterator itn=averages_.find(color);
      const vector<float>& color_avg(itn->second);
      
      stringstream ss;
      ss << setw(12)<< right << coarsening_.getCenterGid(color) <<" ";
      ss << setw(7)<< right << avg_valcolors_.nValues(color);
      
      ss << setprecision(8);
      for(int it=0;it<times_.size();++it)
      {
         ss <<" "<< setw(12)<< right<< color_avg[it];
      }
      ss << endl;
      string line(ss.str());
      Pwrite(line.c_str(), line.size(), 1, file);
   }
   
   Pclose(file);
}

void CaAverageSensor::bufferReactionData(const int loop)
{
   bufferReactionData(0,nlocal_,loop);
}

void CaAverageSensor::bufferReactionData(const int begin, const int end, const int loop)
{
   if(begin==0)loop_buffer_=loop;
   
   for (unsigned ii=begin; ii<end; ++ii)
      buffer_val_[ii]=reaction_.getValue(ii, ca_handle_);
}

void CaAverageSensor::eval(double time, int loop)
{
   if( loop_buffer_<0 )return;
   
   assert( loop==loop_buffer_ );
   
   startTimer(sensorEvalTimer);
   
   times_.push_back(time);
   
   computeColorAverages(buffer_val_);
   
   const std::set<int>& owned_colors(coarsening_.getOwnedColors());
   for(set<int>::const_iterator it = owned_colors.begin();
                                it!= owned_colors.end();
                              ++it)
   {
      const int color=(*it);
      
      vector<float>& color_avg(averages_[color]);
      color_avg.push_back( avg_valcolors_.averageValue(color) );
   }
   
   stopTimer(sensorEvalTimer);
}

void CaAverageSensor::print(double time, int loop)
{
   startTimer(sensorPrintTimer);
   
   if( loop_buffer_<0 )return;

   int myRank;
   MPI_Comm_rank(comm_, &myRank);

   stringstream name;
   name << "snapshot."<<setfill('0')<<setw(12)<<loop;
   string fullname = name.str();
   if (myRank == 0)
      DirTestCreate(fullname.c_str());
   fullname += "/" + filename_;

   writeAverages(fullname,time, loop);   

   times_.clear();
   for(map<int,std::vector<float> >::iterator itg =averages_.begin();
                                              itg!=averages_.end();
                                            ++itg)
   {
      (itg->second).clear();
   }
   
   stopTimer(sensorPrintTimer);
}
