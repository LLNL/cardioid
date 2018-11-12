#ifndef GRADIENTVORONOICOARSENING_HH
#define GRADIENTVORONOICOARSENING_HH

#include "VoronoiCoarsening.hh"
#include "Sensor.hh"

#include <string.h>

class Anatomy;
class PotentialData;


class ColoredCell
{
   static double tol_;
   
   int index_;     // index of cell in storage array
   int color_;
   double norm2_;  // distance from color center
   double inorm2_; // 1./norm2_
   
public:
   ColoredCell(const int index, const int color, const double norm2)
      :index_(index),
       color_(color),
       norm2_(norm2)
      
   {
      if( norm2_>tol_ )inorm2_=1./norm2_;
      else             inorm2_=-1.;
   }
   
   bool normLargerThanTol()const
   {
      return (norm2_>tol_);
   }
   
   int index()const
   {
      return index_;
   }
   int color()const
   {
      return color_;
   }
   double norm2()const
   {
      return norm2_;
   }
   double inorm2()const
   {
      return inorm2_;
   }
};

class GradientVoronoiCoarsening : public Sensor
{
 private:
   VoronoiCoarsening coarsening_;
 
   const std::string filename_;
   const int nFiles_; // number of files for IO
   const Anatomy& anatomy_;
   const PotentialData& vdata_;

   const MPI_Comm comm_;
   const std::string format_;
   
   // flag telling us if values dumped for gradient are partial contribution and need
   // to be consolidated with other tasks output
   const bool use_communication_avoiding_algorithm_;
   
   std::vector<double> dx_;
   std::vector<double> dy_;
   std::vector<double> dz_;
   
   std::vector<ColoredCell>  colored_cells_;
   
   std::set<int> included_eval_colors_;
   
   Long64 nb_sampling_pts_;
   int nb_excluded_pts_;
   
   LocalSums valcolors_;

   LocalSums valMat00_;
   LocalSums valMat01_;
   LocalSums valMat02_;
   LocalSums valMat11_;
   LocalSums valMat12_;
   LocalSums valMat22_;
   std::map<int,double> invDetMat_;

   LocalSums valRHS0_;
   LocalSums valRHS1_;
   LocalSums valRHS2_;
   
   // matrices for LS systems for each local color
   std::map<int,double*> matLS_;
   
   // number of eval times
   int eval_count_;
   double time0_;
   double dt_;
   
   // gradient for each local color
   std::map<int,std::vector<float> > gradients_;
   
   void computeLeastSquareGradients(const double current_time,
                                    const int current_loop);
   void writeGradients(const std::string& filename,
                       const double current_time,
                       const int current_loop)const;
   void setupLSsystem(ro_array_ptr<double> val);
   void computeColorCenterValues(ro_array_ptr<double> val);
   void setupLSmatrix();
   void prologComputeLeastSquareGradients();
   int solve3x3(const double s[6], const double r[3], double x[3], const int color);
   
 public:
   GradientVoronoiCoarsening(const SensorParms& sp,
                             std::string filename,
                             unsigned nFiles, 
                             const Anatomy& anatomy,
                             std::vector<Long64>& gid,
                             const PotentialData& vdata_,
                             const CommTable* commtable,
                             const std::string format,
                             const double max_distance,
                             const bool use_communication_avoiding_algorithm=false);
   
   ~GradientVoronoiCoarsening()
   {
      for(std::map<int,double*>::const_iterator it = matLS_.begin();
                                                it!= matLS_.end();
                                              ++it)
      {
         delete[] it->second;
      }
   }
   
   void eval(double time, int loop);
   void print(double time, int loop);
};

#endif
