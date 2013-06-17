#include "DataVoronoiCoarsening.hh"
#include "PerformanceTimers.hh"
#include "pio.h"
#include "ioUtils.h"
#include "Simulate.hh"
#include "CommTable.hh"
#include "stringUtils.hh"

using namespace PerformanceTimers;

#include <iostream>
#include <sstream>
#include <iomanip>
using namespace std;

/////////////////////////////////////////////////////////////////////

DataVoronoiCoarsening::DataVoronoiCoarsening(const SensorParms& sp,
                                             string filename,
                                             unsigned nFiles,
                                             const Anatomy& anatomy,
                                             vector<Long64>& gid,
                                             const PotentialData& vdata,
                                             const CommTable* commtable,
                                             const double maxDistance)
   :Sensor(sp),
    coarsening_(anatomy, gid, maxDistance, commtable),
    filename_(filename),
    nFiles_(nFiles),
    anatomy_(anatomy),
    vdata_(vdata),
    comm_(commtable->_comm)
{
}

void DataVoronoiCoarsening::computeColorAverages(const VectorDouble32& val)
{
   // calculate local sums
   coarsening_.accumulateValues(val,avg_valcolors_);
   
   coarsening_.exchangeAndSum(avg_valcolors_);
}

void DataVoronoiCoarsening::writeAverages(const string& filename,
                                          const double current_time,
                                          const int current_loop)const
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);

   PFILE* file = Popen(filename.c_str(), "w", comm_);
   if (nFiles_ > 0)
     PioSet(file, "ngroup", nFiles_);

   const int nfields = 4+(int)times_.size(); 
   const int lrec    = 25+13*(int)times_.size();

   const std::set<int>& owned_colors(coarsening_.getOwnedColors());

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
      Pprintf(file, "voronoi FILEHEADER {\n");
      Pprintf(file, "  lrec = %d;\n", lrec);
      Pprintf(file, "  datatype = FIXRECORDASCII;\n");
      Pprintf(file, "  nrecords = %llu;\n", nSnapSub);
      Pprintf(file, "  nfields = %d;\n", nfields);
      string fieldNames="rx ry rz nvals " + concat(vector<string>(times_.size(), "Vm"));
      Pprintf(file, "  field_names = %s;\n", fieldNames.c_str());
      string fieldTypes="d d d d " + concat(vector<string>(times_.size(), "f"));
      Pprintf(file, "  field_types = %s;\n", fieldTypes.c_str());
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  time = %f; loop = %u;\n", times_[0], current_loop);
      if( times_.size()>1 )
         Pprintf(file, "  nsteps = %d; dt = %f\n", times_.size(), times_[1]-times_[0]);
      Pprintf(file, "  h = %4u  0    0\n", anatomy_.nx());
      Pprintf(file, "        0    %4u  0\n", anatomy_.ny());
      Pprintf(file, "        0    0    %4u;\n", anatomy_.nz());
      Pprintf(file, "}\n\n");
   }
   
   const int halfNx = anatomy_.nx()/2;
   const int halfNy = anatomy_.ny()/2;
   const int halfNz = anatomy_.nz()/2;
   
   for(set<int>::const_iterator it = owned_colors.begin();
                                it!= owned_colors.end();
                              ++it)
   {
      const int color=(*it);
      const Vector& v = coarsening_.center(color);
      int ix = int(v.x()) - halfNx;
      int iy = int(v.y()) - halfNy;
      int iz = int(v.z()) - halfNz;
      
      const map< int, vector<float> >::const_iterator itn=averages_.find(color);
      const vector<float>& color_avg(itn->second);
      
      stringstream ss;
      ss << setw(5)<< right << ix<<" ";
      ss << setw(5)<< right << iy<<" ";
      ss << setw(5)<< right << iz<<" ";
      ss << setw(7)<< right << avg_valcolors_.nValues(color);
      
      ss << setprecision(8);
      for(int it=0;it<times_.size();++it)
      {
         ss<< " " << setw(12)<< color_avg[it];
      }
      ss << endl;
      string line(ss.str());
      Pwrite(line.c_str(), line.size(), 1, file);
   }
   
   Pclose(file);
}

void DataVoronoiCoarsening::eval(double time, int loop)
{
   startTimer(sensorEvalTimer);
   
   times_.push_back(time);
   
   computeColorAverages(vdata_.VmArray_);
   
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

void DataVoronoiCoarsening::print(double time, int loop)
{
   startTimer(sensorPrintTimer);
   
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

