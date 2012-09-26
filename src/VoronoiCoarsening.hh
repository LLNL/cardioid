#ifndef VORONOICOARSENING_HH
#define VORONOICOARSENING_HH

#include "Vector.hh"
#include "Anatomy.hh"
#include "IndexToVector.hh"
#include "Long64.hh"
#include "VectorDouble32.hh"

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
   std::map<int,int>    nval_; // number of values summed up for each color
   std::map<int,double> sum_;
   
 public:
   LocalSums()
   {
   }
   
   size_t size()
   {
      return nval_.size();
   }
   
   void add1value(const int color, const double value)
   {
      assert( value==value );
      nval_[color]++;
      sum_[color]+=value;
   }
   void setSum(const int color, const int nval, const double sum)
   {
      assert( sum==sum );
      
      nval_[color]=nval;
      sum_[color] =sum;
   }
   void addnvalues(const int color, const double value, const int n)
   {
      assert( value==value );
      nval_[color]+=n;
      sum_[color]+=value;
   }
   
   void packData(PackedData* packeddata, const int ndata)
   {
      assert( nval_.size()==sum_.size() );
      assert( ndata==nval_.size() );
      
      int i=0;
      std::map<int,double>::const_iterator itv=sum_.begin();
      for(std::map<int,int>::const_iterator itn  = nval_.begin();
                                            itn != nval_.end();
                                          ++itn)
      {
         assert( itn->first==itv->first );
         
         packeddata[i].color=itn->first;
         packeddata[i].n    =itn->second;
         packeddata[i].value=itv->second;
         
         ++itv;
         ++i;
      }
      if( i!=ndata )
      {
         std::cerr<<"i="<<i<<", ndata="<<ndata<<std::endl;
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
   
   double value(const int color)const
   {
      std::map<int,double>::const_iterator is=sum_.find(color);
      assert( is!=sum_.end() );
      if( is->second!=is->second )
         std::cerr<<"ERROR!!! is->second="<<is->second<<std::endl;
      
      return is->second;
   }
   
   int nValues(const int color)const
   {
      std::map<int,int>::const_iterator in=nval_.find(color);
      assert( in!=nval_.end() );
      assert( in->second>0 );
      return in->second;
   }
};


class VoronoiCoarsening
{
 private:
   std::map<int,int> ncolors_to_send_; // remote pe -> number of colors
   std::map<int,int> ncolors_to_recv_; // remote pe -> number of colors

   std::set<int> remote_tasks_;
   std::set<int> dst_tasks_;
   std::set<int> src_tasks_;

   std::set<int> local_colors_;

   void computeColorAverages(const std::vector<double>& val);
   void computeColorCenterValues(const std::vector<double>& val);

   std::map<int,int> ncolors_; // color -> number of local cells of that color
   
   MPI_Comm comm_;

   const Anatomy& anatomy_;

   // shared globally
   std::vector<Vector> centers_;
   std::vector<short> multiplicities_;

   // local only
   const std::vector<AnatomyCell>& cells_; // local cells
   std::vector<int> colors_; // colors of local cells

   // set of colors local task is responsible for
   std::set<int> owned_colors_;
   
   IndexToVector indexToVector_;

   void setupComm(const std::map< int, int* >& nremote_colors_for_task, 
                const int size_nremote_colors_for_task);
   void setOwnedColors(const std::map< int, int* >& nremote_colors_for_task, 
                       const int size_nremote_colors_for_task);

 public:
   VoronoiCoarsening(const Anatomy& anatomy,
                     const std::vector<Long64>& gid,
                     MPI_Comm comm);
   void computeRemoteTasks();
   void exchangeAndSum(LocalSums& valcolors);
   void exchangeAndSum(std::vector<LocalSums*> valcolors);
   int bruteForceColoring();
   void colorDisplacements(std::vector<double>& dx,
                           std::vector<double>& dy,
                           std::vector<double>& dz);
   void accumulateValues(const VectorDouble32& val, LocalSums& valcolors);
   
   const Vector& center(const int color)const
   {
      return centers_[color];
   }
   
   short multiplicity(const int color)const
   {
      return multiplicities_[color];
   }
   
   const std::set<int>& getOwnedColors()const
   {
      return owned_colors_;
   }
   
   int getColor(const int ic)const
   {
      return colors_[ic];
   }
};

#endif
