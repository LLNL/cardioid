#ifndef CARDIAC_COEF
#define CARDIAC_COEF

#include "mfem.hpp"
#include "Reaction.hh"
#include "ThreadServer.hh"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

/// Wrapper for interface to melodee-generated active tension
/// values that live only at quadrature points
class ReactionFunction : public QuadratureFunction
{
private:

   /// Underlying quadrature space
   QuadratureSpace *QuadS;

   /// Finite element space for the displacements
   FiniteElementSpace *fes;
   QuadratureFunction VmQuad;

   /// melodee compatible data members
   lazy_array<double> Vm;
   lazy_array<double> iStim;
   lazy_array<double> dVm;

   /// Number of integration points
   int nCells;
   double dt;
   std::string objectName;
   ThreadTeam threadGroup;

   /// Cell model object
   std::shared_ptr<Reaction> reaction;
   
public:
   ReactionFunction(QuadratureSpace *qs, FiniteElementSpace *f,const double new_dt, const std::string new_objectName, const ThreadTeam& group);

   /// Call to melodee generated initialization
   void Initialize();

   /// Call to melodee generated calc
   void Calc(const Vector& x);

 private:
   void CalcVm(const Vector& x);
};


class QuadratureIntegrator : public LinearFormIntegrator
{
 public:
   QuadratureIntegrator(QuadratureFunction* p_newQuadFunction, const double scale=1);
   virtual void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &elvect);
 private:
   QuadratureFunction *p_quadFunction;
   double scale;
};

class StimulusLocation
{
 public:
   virtual ~StimulusLocation() {};
   virtual bool contains(const int elementNo, const Vector& x) = 0;
};

class BoxStimulus : public StimulusLocation
{
 public:
   BoxStimulus(const double xmin, const double xmax,
               const double ymin, const double ymax,
               const double zmin, const double zmax)
   : xmin_(xmin), xmax_(xmax),
     ymin_(xmin), ymax_(xmax),
     zmin_(xmin), zmax_(xmax)
   {}
   virtual bool contains(const int elementNo, const Vector& x)
   {
      if (0
          && xmin_ <= x[0] && x[0] <= xmax_
          && ymin_ <= x[1] && x[1] <= ymax_
          && zmin_ <= x[2] && x[2] <= zmax_
      )
      {
         return true;
      }
      return false;
   }
 private:
   double xmin_,xmax_;
   double ymin_,ymax_;
   double zmin_,zmax_;
};
   

class CenterBallStimulus : public StimulusLocation
{
 public:
   CenterBallStimulus(const double x, const double y, const double z, const double size)
   : x_(x), y_(y), z_(z), size_(size) {}
   virtual bool contains(const int elementNo, const Vector& x)
   {
      if (0
          +(x_-x[0])*(x_-x[0])
          +(y_-x[1])*(y_-x[1])
          +(z_-x[2])*(z_-x[1])
          <= size_*size_
      )
      {
         return true;
      }
      return false;
   }
 private:
   double x_;
   double y_;
   double z_;
   double size_;
};

class StimulusWaveform
{
 public:
   ~StimulusWaveform() {}
   virtual double eval(const double time) = 0;
};

class SquareWaveform : public StimulusWaveform
{
   virtual double eval(const double time) { return 1; }
};


class Stimulus
{
 public:
   Stimulus(const int numTimes, const double startTime, const double duration, const double bcl,
            const double strength, std::shared_ptr<StimulusLocation> loc, std::shared_ptr<StimulusWaveform> wave)
   : numTimes_(numTimes), startTime_(startTime_), duration_(duration), bcl_(bcl),
     strength_(strength),
     loc_(loc), wave_(wave)
   {}

   inline double isOn(const double time) const
   {
      return stimTime(time) > 0;
   }
   inline double stimTime(const double time) const
   {
      for (int itime=0; itime<numTimes_; itime++)
      {
         double shiftedTime = time-startTime_-itime*bcl_;
         if (shiftedTime < duration_) { return shiftedTime; }
      }
      return -1;
   }
   inline double eval(const double time, const int elementNo, const Vector& x)
   {
      double waveTime = stimTime(time);
      if (time >= 0 && loc_->contains(elementNo, x))
      {
         return strength_*wave_->eval(waveTime);
      }
      else
      {
         return 0;
      }
   }
   
 private:

   int numTimes_;
   double startTime_;
   double duration_;
   double bcl_;
   double strength_;
   std::shared_ptr<StimulusLocation> loc_;
   std::shared_ptr<StimulusWaveform> wave_;
};

class StimulusCollection : public Coefficient
{
 public:
   virtual double Eval(ElementTransformation& T, const IntegrationPoint &ip);
   void add(Stimulus stim);
   void updateTime(const double time) { time_ = time; }
 private:
   double time_;
   std::vector<Stimulus> stim_;
};

#endif
