#ifndef SENSOR_HH
#define SENSOR_HH
#include <vector>
#include <iostream>

struct SensorParms
{
   int evalRate;
   int printRate;
   double startTime;
   double endTime;
   double value;
};


class Sensor
{
 public:
   Sensor(const SensorParms& p)
   : evalRate_(p.evalRate),
     printRate_(p.printRate),
     startTime_(p.startTime),
     endTime_(p.endTime),
     value_(p.value)
   {}
   virtual ~Sensor() {};

   virtual void run(double time, int loop)
   {
      if (time < startTime_ || time > endTime_)
         return;
      if (loop % evalRate_ == 0)
      {
         eval( time, loop);
      }
      if (loop % printRate_ == 0)
      {
         print(time, loop);
      }
   }

   int printRate()const{return printRate_;}
   int evalRate()const{return evalRate_;}
   
   // to be implemented if sensor needs to know 
   // about reaction data
   virtual void bufferReactionData(const int loop)
   { return; };
   // This version called in threaded regions where each thread buffers
   // a portion of the cells from begin to end.
   virtual void bufferReactionData(const int begin, const int end, const int loop)
   { return; };
   

   bool checkPrintAtStep(const int loop)const
   {
      return (loop % printRate_ == 0);
   }

   bool checkEvalAtStep(const int loop)const
   {
      return (loop % evalRate_ == 0);
   }

 private:
   int evalRate_;
   int printRate_;
   double startTime_;
   double endTime_;
   double value_;
    
   virtual void print(double time, int loop) = 0;
   virtual void eval(double time, int loop) = 0;
};

#endif
