#include "CaAverageSensor.hh"

#include "Anatomy.hh"
#include "pio.h"
#include "ioUtils.h"
#include "Reaction.hh"

#include <vector>
#include <iomanip>
#include <sstream>

CaAverageSensor::CaAverageSensor(const SensorParms& sp,
                   string filename,
                   const Anatomy& anatomy,
                   const vector<Long64>& gid,
                   const Reaction& reaction,
                   MPI_Comm comm)
   :Sensor(sp),
    coarsening_(anatomy,gid,comm),
    filename_(filename),
    reaction_(reaction),
    comm_(comm),
    loop_buffer_(-1)
{
   // color local cells
   int ret=coarsening_.bruteForceColoring();
   assert( ret>=0 );
   
   coarsening_.computeRemoteTasks();
   
   nx_=anatomy.nx();
   ny_=anatomy.ny();
   nz_=anatomy.nz();
   
   nlocal_=anatomy.nLocal();

   buffer_val_.resize(nlocal_);
}

void CaAverageSensor::computeColorAverages(const VectorDouble32& val)
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
      Pprintf(file, "  field_names = rx ry rz nvals avgCa;\n");
      Pprintf(file, "  field_types = u u u u f;\n" );
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  time = %f; loop = %u;\n", current_time, current_loop);
      Pprintf(file, "  h = %4u  0    0\n", nx_);
      Pprintf(file, "        0    %4u  0\n", ny_);
      Pprintf(file, "        0    0    %4u;\n", nz_);
      Pprintf(file, "}\n\n");
   }
   
   char line[lrec+1];
   const int halfNx = nx_/2;
   const int halfNy = ny_/2;
   const int halfNz = nz_/2;
   
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
      const short m=coarsening_.multiplicity(color);
      for(short ii=0;ii<m;ii++)
         Pwrite(line, lrec, 1, file);
   }
   
   Pclose(file);
}

void CaAverageSensor::bufferReactionData(const int loop)
{
   bufferReactionData(0,nlocal_,loop);
}

void CaAverageSensor::bufferReactionData(const int begin, const int end, const int loop)
{
   loop_buffer_=loop;
   const int handle=reaction_.getVarHandle("Ca_i");
   
   for (unsigned ii=begin; ii<end; ++ii)
      buffer_val_[ii]=reaction_.getValue(ii, handle);
}

void CaAverageSensor::eval(double time, int loop)
{
   assert( loop==loop_buffer_ );
   
   computeColorAverages(buffer_val_);
}

void CaAverageSensor::print(double time, int loop)
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
