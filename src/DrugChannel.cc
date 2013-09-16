#include "DrugChannel.hh"
#include <cmath>

DrugChannel::DrugChannel(string current, double low, double high, double nh, double xc50) :
    current_(current),low_(low),high_(high),nh_(nh),xc50_(xc50)
{
}

DrugChannel::~DrugChannel()
{
}

double DrugChannel::scalingFactor(double concentration)
{
   double rescale = 1.0;
   if (concentration > 0.0)
   {
      double tpow = pow(xc50_/concentration,nh_);
      rescale = low_ + (high_-low_)/( 1. + tpow);
   }
   return rescale;
}
