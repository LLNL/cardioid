#include "DataVoronoiCoarsening.hh"
#include "pio.h"
#include "ioUtils.h"

#include <iostream>
#include <sstream>
#include <iomanip>
using namespace std;

/////////////////////////////////////////////////////////////////////

DataVoronoiCoarsening::DataVoronoiCoarsening(const SensorParms& sp,
                                     string filename,
                                     const Anatomy& anatomy,
                                     const vector<Long64>& gid,
                                     MPI_Comm comm)
   :Sensor(sp),
    coarsening_(anatomy,gid,comm),
    filename_(filename),
    anatomy_(anatomy),
    comm_(comm)
{
   // color local cells
   int ret=coarsening_.bruteForceColoring();
   assert( ret>=0 );
   
   coarsening_.computeRemoteTasks();
}

void DataVoronoiCoarsening::computeColorAverages(const vector<double>& val)
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

   char fmt[] = "%5d %5d %5d %5d %18.12f";
   int lrec = 55;
   int nfields = 5; 

   Long64 nSnapSub = -1;
   const std::set<int>& owned_colors=coarsening_.getOwnedColors();
   Long64 nSnapSubLoc = owned_colors.size();
   MPI_Allreduce(&nSnapSubLoc, &nSnapSub, 1, MPI_LONG_LONG, MPI_SUM, comm_);

   if (myRank == 0)
   {
      // write header
      int nfiles;
      Pget(file,"ngroup",&nfiles);
      Pprintf(file, "cellViz FILEHEADER {\n");
      Pprintf(file, "  lrec = %d;\n", lrec);
      Pprintf(file, "  datatype = FIXRECORDASCII;\n");
      Pprintf(file, "  nrecords = %llu;\n", nSnapSub);
      Pprintf(file, "  nfields = %d;\n", nfields);
      Pprintf(file, "  field_names = rx ry rz nvals avgVm;\n");
      Pprintf(file, "  field_types = u u u u f;\n" );
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  time = %f; loop = %u;\n", current_time, current_loop);
      Pprintf(file, "  h = %4u  0    0\n", anatomy_.nx());
      Pprintf(file, "        0    %4u  0\n", anatomy_.ny());
      Pprintf(file, "        0    0    %4u;\n", anatomy_.nz());
      Pprintf(file, "}\n\n");
   }
   
   char line[lrec+1];
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
      
      int l = snprintf(line, lrec, fmt,
                       ix, iy, iz,
                       avg_valcolors_.nValues(color),
                       avg_valcolors_.averageValue(color));
      
      if (myRank == 0 && l>=lrec ){
         cerr<<"ERROR: printed record truncated in file "<<filename<<endl;
         cerr<<"This could be caused by out of range values"<<endl;
         cerr<<"record="<<line<<endl;
         break;
      }
      for (; l < lrec - 1; l++) line[l] = (char)' ';
      line[l++] = (char)'\n';
      assert (l==lrec);
      Pwrite(line, lrec, 1, file);
   }
   
   Pclose(file);
}

void DataVoronoiCoarsening::eval(double time, int loop,
                             const vector<double>& Vm,
                             const vector<double>&,
                             const vector<double>&)
{
   computeColorAverages(Vm);
}

void DataVoronoiCoarsening::print(double time, int loop,
                              const vector<double>& Vm,
                              const vector<double>&,
                              const vector<double>&)
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);

   stringstream name;
   name << "snapshot."<<setfill('0')<<setw(12)<<loop;
   string fullname = name.str();
   if (myRank == 0)
      DirTestCreate(fullname.c_str());
   fullname += "/" + filename_;

   writeAverages(fullname,time, loop);   
}

