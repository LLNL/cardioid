#ifndef SENSOR_HH
#define SENSOR_HH
#include <vector>

struct SensorParms
{
   int evalRate;
   int printRate;
};


class Sensor
{
 public:
   Sensor(const SensorParms& p)
   : evalRate_(p.evalRate),
     printRate_(p.printRate)
   {}
   virtual ~Sensor() {};

   virtual void run(double time, int loop)
   {
      if (loop % printRate_ == 0)
      {
         eval( time, loop);
         print(time, loop);
      }
   }

   int printRate()const{return printRate_;}
   int evalRate()const{return evalRate_;}
   
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

   virtual void print(double time, int loop) = 0;
   virtual void eval(double time, int loop) = 0;
};

#endif
