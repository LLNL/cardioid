#ifndef CARDIAC_COEF
#define CARDIAC_COEF

#include "mfem.hpp"
#include "ReactionManager.hh"
#include "ThreadServer.hh"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

class ReactionWrapper
{
 public:
   ReactionWrapper(const double new_dt, const std::vector<std::string> new_objectNames, const ThreadTeam& group, const std::vector<int>& cellTypes);
   void Initialize();
   void Calc();
   Vector& getVmReadwrite();
   Vector& getIionReadwrite();
   const Vector& getVmReadonly() const;
   const Vector& getIionReadonly() const;

   //private:
   lazy_array<double> Vm;
   lazy_array<double> iStim;
   lazy_array<double> dVm;

   int nCells;
   double dt;
   std::string objectName;
   ThreadTeam threadGroup;

   ReactionManager reaction;

   Vector Vm_vector;
   Vector Iion_vector;
};



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

   /// Cell model object
   ReactionWrapper* reactionWrapper;
   
public:
   ReactionFunction(QuadratureSpace *qs, FiniteElementSpace *f,ReactionWrapper* rw);

   /// Call to melodee generated calc
   void Calc(const Vector& x);

 private:
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
      if (1
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
          +(z_-x[2])*(z_-x[2])
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
         if (shiftedTime <= duration_) { return shiftedTime; }
      }
      return -1;
   }
   inline double eval(const double time, const int elementNo, const Vector& x)
   {
      double waveTime = stimTime(time);
      if (waveTime >= 0 && loc_->contains(elementNo, x))
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
   StimulusCollection(const double new_dt) : dt_(new_dt) {}
   virtual double Eval(ElementTransformation& T, const IntegrationPoint &ip);
   void add(Stimulus stim);
   void updateTime(const double time) { time_ = time; }
 private:
   double time_;
   double dt_;
   std::vector<Stimulus> stim_;
};

#endif
