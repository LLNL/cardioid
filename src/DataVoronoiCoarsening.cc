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

//Includes hacks from JPC to get Coarsened Activation Times out (changes marked with //AT-HACK), would eventually like to make it less hacky and interface with object.data (one bad thing I am doing is disabling the capability that was in previous version of DataVoronoiCoarsening code that allows you to print out data from multiple time points when evalrate<printrate, would be better if this was an object.data option.  I'm disabling this b/c to find Activation Times you'd want to evaluate often but print less frequently, which triggers printing out data from multiple time points in coarsened_anatomy#* files, which crashes the post-processing anatomy2ensight visualization tool as well as ECG code)

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
   // Setting up some variables to keep track of activationTime and activation status for each gid in sensor.txt, activationTime is normalized to initiation of simulation at time=0 ms	//AT-HACK
   const std::set<int>& owned_colors(coarsening_.getOwnedColors());	//AT-HACK
   activationTime_.resize(owned_colors.size());	//AT-HACK, activation time in ms (normalized to initation of simulation at t=0 ms) of select gids in sensor.txt
   active_.resize(owned_colors.size(), false);	//AT-HACK, active status of select gids in sensor.txt?, true or false
   clear();	//AT-HACK, function that sets active to false and AT to -1000 ms for all select gids in sensor.txt
}

void DataVoronoiCoarsening::computeColorAverages(ro_array_ptr<double> val)
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

//   const int nfields = 2+(int)times_.size(); 	//AT-HACK, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data
 //  const int lrec    = 20+13*(int)times_.size() + 1;	//AT-HACK, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data

   const int nfields = 3; 	//AT-HACK, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data
   const int lrec    = 20+13 + 1;	//AT-HACK, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data
 
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
//      string fieldNames="gid nvals " + concat(vector<string>(times_.size(), "Vm"));	//AT-HACK, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data
      string fieldNames="gid nvals Vm";	//AT-HACK, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data
      Pprintf(file, "  field_names = %s;\n", fieldNames.c_str());
      string fieldTypes="u d f";
//      string fieldTypes="u d " + concat(vector<string>(times_.size(), "f"));	//AT-HACK, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data
      Pprintf(file, "  field_types = %s;\n", fieldTypes.c_str());	//AT-HACK, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  time = %f; loop = %u;\n", times_[0], current_loop);
      if( times_.size()>1 )
         Pprintf(file, "  nsteps = %d; dt = %f;\n", times_.size(), times_[1]-times_[0]);	//AT-HACK, adding semicolon after dt
      Pprintf(file, "  h = %4u  0    0\n", anatomy_.nx());
      Pprintf(file, "        0    %4u  0\n", anatomy_.ny());
      Pprintf(file, "        0    0    %4u;\n", anatomy_.nz());
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
      
      ss << setprecision(5) << scientific;
      
//      for(int it=0;it<times_.size();++it)			//AT-HACK, lines 122 o 126, disabling ability to print out data from multiple time points when evalrate<printrate in object.data, and just prints most recently evaluated data
//      {
//           ss<< " " << setw(12)<< color_avg[it];
//      }
      ss<< " " << setw(12)<< color_avg[times_.size()-1];
      ss << endl;
      string line(ss.str());
      Pwrite(line.c_str(), line.size(), 1, file);
   }
   
   Pclose(file);
}

void DataVoronoiCoarsening::writeAveragesAT(const string& filename,		
                                          const double current_time,
                                          const int current_loop)const	//AT-HACK, this whole function is new, and is a sister function of writeAverages, but this one tells cardioid how to write out coarsened AT, instead of coarsened Vm.  For confusing parts I've added further AT-HACK comments within this function.
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);

   PFILE* fileAT = Popen(filename.c_str(),"w",comm_);
   if (nFiles_ > 0)
   {
     PioSet(fileAT, "ngroup", nFiles_);
   }


//   const int nfields = 2+(int)times_.size(); 
//   const int lrec    = 20+13*(int)times_.size() + 1;

   const int nfields = 3; 
   const int lrec    = 20+13 + 1;
 

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
//      string fieldNames="gid nvals " + concat(vector<string>(times_.size(), "Vm"));
//      string fieldTypes="u d " + concat(vector<string>(times_.size(), "f"));

      string fieldNames="gid nvals Vm";
      string fieldTypes="u d f";
      
      //write header for activation times files
      Pget(fileAT,"ngroup",&nfiles);
      Pprintf(fileAT, "voronoi FILEHEADER {\n");
      Pprintf(fileAT, "  lrec = %d;\n", lrec);
      Pprintf(fileAT, "  datatype = FIXRECORDASCII;\n");
      Pprintf(fileAT, "  nrecords = %llu;\n", nSnapSub);
      Pprintf(fileAT, "  nfields = %d;\n", nfields);
      Pprintf(fileAT, "  field_names = %s;\n", fieldNames.c_str());
      Pprintf(fileAT, "  field_types = %s;\n", fieldTypes.c_str());
      Pprintf(fileAT, "  nfiles = %u;\n", nfiles);
      Pprintf(fileAT, "  time = %f; loop = %u;\n", times_[0], current_loop);
      if( times_.size()>1 )
         Pprintf(fileAT, "  nsteps = %d; dt = %f;\n", times_.size(), times_[1]-times_[0]);
      Pprintf(fileAT, "  h = %4u  0    0\n", anatomy_.nx());
      Pprintf(fileAT, "        0    %4u  0\n", anatomy_.ny());
      Pprintf(fileAT, "        0    0    %4u;\n", anatomy_.nz());
      Pprintf(fileAT, "}\n\n");
   }
   
   int ii=0;		//AT-HACK, a counter to index into activationTime array, follows the iterator "it" as it iterates through owned_colors, i.e., the colors that are averaged together to get final averaged color in Voronoi Coarsening
   for(set<int>::const_iterator it = owned_colors.begin();
                                it!= owned_colors.end();
                              ++it)
   {
      const int color=(*it);
      
      const map< int, vector<float> >::const_iterator itn=averages_.find(color);
      const vector<float>& color_avg(itn->second);
      
      
      //Do same for writing to activation times
      stringstream ssAT;
      ssAT << setw(12)<< right << coarsening_.getCenterGid(color) <<" ";
      ssAT << setw(7)<< right << avg_valcolors_.nValues(color);
      
      ssAT << setprecision(5) << scientific;
      
      for (unsigned jj=0; jj<activationTime_[ii].size(); ++jj)		//AT-HACK, while this jj loop implies that multiple ATs may be recorded for each coarsened_anatomy gid (ii), this in fact will never happen, accoring to the way the eval() function works below, so this for loop is kind of just an artificat
        ssAT<< " " << setw(12)<< activationTime_[ii][jj];
	   
     
      ssAT << endl;
      string lineAT(ssAT.str());
      Pwrite(lineAT.c_str(), lineAT.size(), 1, fileAT);
      ii++;
   }
   
   Pclose(fileAT);
}


void DataVoronoiCoarsening::eval(double time, int loop)			// AT-HACK, evaluates/computes values you want printed out each time step, added some code here to evaluate AT
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);
   
   startTimer(sensorEvalTimer);
   
   times_.push_back(time);

   ro_array_ptr<double> VmArray = vdata_.VmTransport_.useOn(CPU);
   computeColorAverages(VmArray);
   
   
   const std::set<int>& owned_colors(coarsening_.getOwnedColors());
   int ii = 0;	//AT-HACK, a counter to index into activationTime array, follows the iterator "it" as it iterates through owned_colors, i.e., the colors that are averaged together to get final averaged color in Voronoi Coarsening
   for(set<int>::const_iterator it = owned_colors.begin();
                                it!= owned_colors.end();
                              ++it)
   {
      const int color=(*it);
      
      vector<float>& color_avg(averages_[color]);
      color_avg.push_back( avg_valcolors_.averageValue(color) );
      
      // AT-HACK, lines 249 to 266,aA hack to get out coarsened activation times fron Voronoi coarsened Vm values, when coarsened Vm value reaches certain threshold (-40 mV in this case)
      if (active_[ii] && avg_valcolors_.averageValue(color) < -40 )	//AT-HACK, once coarsened Vm falls below -40 mV and this coarsened_anatomy gid is currently active, reset it to be inactive and default AT of -1000 mS, in preparation for next cardiac cycle heartbeat
      {
        active_[ii] = false;
	activationTime_[ii].clear();
        activationTime_[ii].reserve(10);
        activationTime_[ii].push_back(-1000);
      }  
      if (! active_[ii] && avg_valcolors_.averageValue(color) > -40 )	//AT-HACK, once coarsened Vm goes above -40 mV and this coarsened_anatomy gid is currently inactive, set it to be active and record AT!
      {
         active_[ii] = true;
	 activationTime_[ii].clear();
         activationTime_[ii].reserve(10);
         activationTime_[ii].push_back(time);
      }
      ii++;    
   }
   
   stopTimer(sensorEvalTimer);
   
}

void DataVoronoiCoarsening::print(double time, int loop)				//calls the writeAverages functions and actually prints to files in folder snapshot.0xxxxx
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
   string fullnameAT = name.str() + "/" + filename_ + "AT";		//AT-HACK, filenames for AT will just be the same as for regular files, with "AT" appended
    
   writeAverages(fullname,time, loop); 		//Print coarsened Vm values to appropriate files
   writeAveragesAT(fullnameAT,time,loop);  	//AT-HACK, print AT for coarsened anatomy gids to appropriate files

   times_.clear();
   for(map<int,std::vector<float> >::iterator itg =averages_.begin();
                                              itg!=averages_.end();
                                            ++itg)
   {
      (itg->second).clear();
   }
   
   stopTimer(sensorPrintTimer);
}

void DataVoronoiCoarsening::clear()					//AT-HACK, this whole function was added as part of the hack, it resets all select gid ATs to -1000 ms and active status to false
{
   const std::set<int>& owned_colors(coarsening_.getOwnedColors());
   for (unsigned ii=0; ii<owned_colors.size(); ++ii)
   {
      active_[ii] = false;
      activationTime_[ii].clear();
      activationTime_[ii].reserve(10);
      activationTime_[ii].push_back(-1000);
   }
}
