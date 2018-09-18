#ifndef VORONOICOARSENING_HH
#define VORONOICOARSENING_HH

#include "Vector.hh"
#include "Anatomy.hh"
#include "IndexToVector.hh"
#include "Long64.hh"
#include "lazy_array.hh"

#include <mpi.h>

#include <map>
#include <set>
#include <vector>
#include <cassert>
#include <iostream>
class CommTable;
class LocalSums;

class VoronoiCoarsening
{
 public:
   VoronoiCoarsening(const Anatomy& anatomy,
                     std::vector<Long64>& sensorPoint,
                     const double maxDistance,
                     const CommTable* commtable);
   void exchangeAndSum(LocalSums& valcolors);
   void exchangeAndSum(std::vector<LocalSums*> valcolors);
   void colorDisplacements(std::vector<double>& dx,
                           std::vector<double>& dy,
                           std::vector<double>& dz);
   void accumulateValues(ro_array_ptr<double> val, LocalSums& valcolors);
   
   const Long64& getCenterGid(int color) const
   {
      std::map<int, Long64>::const_iterator here = colorToGidMap_.find(color);
      assert( here!=colorToGidMap_.end() );
      return here->second;
   }
   
   const std::set<int>& getOwnedColors()const {return ownedColors_;}
   const std::set<int>& getLocalColors()const {return localColors_;}
   int getColor(const int ic)const {return cell_colors_[ic]; }

 private:

   int gaoColoring(const double maxDistance,
                   const std::vector<Long64>& sensorPoint);   
   void computeRemoteTasks();
   void computeColorAverages(const std::vector<double>& val);
   void computeColorCenterValues(const std::vector<double>& val);
   void setupComm(const std::map< int, int* >& nremote_colors_for_task, 
                  const int size_nremote_colors_for_task);

   Vector getDomaincenter()const;
   double getDomainRadius(const Vector& domain_center)const;
   std::map<int,Vector> getCloseCenters(const Vector& domain_center, 
                                        const double)const;
   


   std::map<int,int> ncolors_to_send_; // remote pe -> number of colors
   std::map<int,int> ncolors_to_recv_; // remote pe -> number of colors

   std::set<int> remote_tasks_;
   std::set<int> dst_tasks_;
   std::set<int> src_tasks_;

   std::set<int> localColors_;

   std::map<int,int> ncolors_; // color -> number of local cells of that color
   
   MPI_Comm comm_;
   const CommTable* commTable_;

   const Anatomy& anatomy_;

   // shared globally
   std::map<int, Long64> colorToGidMap_; // color->cell gid

   // local only
   std::vector<int> cell_colors_; // colors of local cells

   // set of colors local task is responsible for
   std::set<int> ownedColors_;
   
   IndexToVector indexToVector_;

};

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
   std::map<int,double> sum_;  // sum of values for each color
   
 public:
   LocalSums()
   {
   }
   
   size_t size()const
   {
      return nval_.size();
   }
   
   void increaseValue(const int color, const double value)
   {
      //assert( value==value );
      sum_[color]+=value;
   }
   void add1value(const int color, const double value)
   {
      //assert( value==value );
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
      //assert( sum_.size()>0 );
      std::map<int,double>::const_iterator is=sum_.find(color);
      //if( is==sum_.end() )
      //{
      //   int myRank;
      //   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      //   std::cerr<<"ERROR: rank="<<myRank<<", value(color="<<color<<") not in map of size "<<sum_.size()<<std::endl;
      //   usleep(1);
      //}
      assert( is!=sum_.end() );
      if( is->second!=is->second )
         std::cerr<<"ERROR!!! is->second="<<is->second<<std::endl;
      
      return is->second;
   }
   
   int nValues(const int color)const
   {
      std::map<int,int>::const_iterator in=nval_.find(color);
      if( in==nval_.end() )return 0;
      assert( in->second>0 );
      return in->second;
   }
};


#endif
