#include "RandomPulse.hh"
#include <cmath>
#include <cassert>
#include <stdlib.h>
#include <iostream>

using namespace std;


// random state for initial run
unsigned short _init_rand_state[3]={2691,9186,2737};

RandomPulse::RandomPulse(const double minperiod, const double maxperiod, 
                         const double duration,
                         const double amplitude, const double phase,
                         const double begintime,
                         const unsigned short rand_state[3])
: minperiod_(minperiod),
  maxperiod_(maxperiod),
  duration_(duration),
  amplitude_(amplitude),
  phase_(phase),
  begintime_(begintime)
{
   assert(minperiod_ > 0);
   assert(maxperiod_ > 0);
   assert(phase_ < minperiod_);
   assert(duration_ > 0);
   assert(phase_+duration_ <= minperiod_);
   
   endtime_=begintime_;

   if( rand_state==0 ){
      for(short d=0;d<3;d++)
         rand_state_[d]=_init_rand_state[d];  // default
   }else{
      for(short d=0;d<3;d++)
         rand_state_[d]=rand_state[d];
   
   }
}

double RandomPulse::eval(double time)
{
   assert( time>=begintime_ );
   if( time>endtime_ )
   {
      begintime_=endtime_;
      endtime_+=(minperiod_+erand48(rand_state_)*(maxperiod_-minperiod_));
      //std::cout<<"New end time for period: "<<endtime_<<std::endl;
   }
   double t = time - begintime_;
   
   if (t >= phase_ && t < phase_+duration_)
      return amplitude_;
   return 0.0;
}
