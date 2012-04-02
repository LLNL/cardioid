#ifndef VORONOICOARSENING_HH
#define VORONOICOARSENING_HH

#include "Vector.hh"
#include "AnatomyCell.hh"
#include "Anatomy.hh"
#include "IndexToVector.hh"
#include "Long64.hh"

#include <mpi.h>

#include <map>
#include <set>
#include <vector>
#include <cassert>
#include <iostream>

// packed format for communications
struct PackedData
{
   int color;
   int n;
   double value;
};

// local sums of values for each color
class LocalSums
{
 private:
   std::map<int,int>    nval_; // number of values summed up
   std::map<int,double> sum_;
   
 public:
   LocalSums()
   {
   }
   
   void add1value(const int color, const double value)
   {
      assert( value==value );
      nval_[color]++;
      sum_[color]+=value;
   }
   void addnvalues(const int color, const double value, const int n)
   {
      assert( value==value );
      nval_[color]+=n;
      sum_[color]+=value;
   }
   
   void packData(PackedData* packeddata, const int ndata)
   {
      int i=0;
      std::map<int,double>::iterator itv=sum_.begin();
      for(std::map<int,int>::iterator itn = nval_.begin();
                                      itn != nval_.end();
                                      ++itn)
      {
         assert( itn->first==itv->first );
         
         packeddata[i].color=itn->first;
         packeddata[i].n    =itn->second;
         packeddata[i].value=itv->second;
         
         itv++;
         i++;
      }
      assert( i==ndata );
   }
   
   void clear()
   {
      nval_.clear();
      sum_.clear();
   }
   
   double averageValue(const int color)const
   {
      assert( nval_.find(color)!=nval_.end() );
      std::map<int,double>::const_iterator is=sum_.find(color);
      std::map<int,int>::const_iterator in=nval_.find(color);
      assert( is!=sum_.end() );
      assert( in!=nval_.end() );
      assert( in->second>0 );
      if( is->second!=is->second )
         std::cerr<<"ERROR!!! is->second="<<is->second<<std::endl;
      
      return is->second/(double)(in->second);
   }
   
   int nValues(const int color)const
   {
      std::map<int,int>::const_iterator in=nval_.find(color);
      return in->second;
   }
};


class VoronoiCoarsening
{
 private:
   MPI_Comm comm_;

   const Anatomy& anatomy_;

   // shared globally
   std::vector<Vector> centers_;

   // local only
   const std::vector<AnatomyCell>& cells_; // local cells
   std::vector<int> colors_; // colors of local cells
   
   std::map<int,int> ncolors_; // color -> number of local cells of that color
   
   std::map<int,int> ncolors_to_send_; // remote pe -> number of colors
   std::map<int,int> ncolors_to_recv_; // remote pe -> number of colors

   // set of colors local task is responsible for
   std::set<int> local_colors_;
   
   std::set<int> remote_tasks_;
   std::set<int> dst_tasks_;
   std::set<int> src_tasks_;

   IndexToVector indexToVector_;

   LocalSums valcolors_;

 public:
   VoronoiCoarsening(const Anatomy& anatomy,
                     const std::vector<Long64>& gid,
                     MPI_Comm comm);
   int bruteForceColoring();
   void computeRemoteTasks();
   void computeColorAverages(const std::vector<double>& val);
   void writeAverages(const std::string& filename,
                      const double current_time,
                      const int current_loop)const;
};

#endif
