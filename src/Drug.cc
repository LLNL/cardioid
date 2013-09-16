#include "Drug.hh"
#include "DrugChannel.hh"
#include <cassert>

Drug::Drug(const string name, const double concentration) : name_(name), concentration_(concentration)
{
}

Drug::~Drug()
{
   for (int ii=0; ii<channels_.size(); ii++)
      if (channels_[ii] != 0)
         delete channels_[ii];
}

void Drug::addChannel(const string current, const double low, const double high,
                      const double nh, const double xc50)    
{

   DrugChannel* dchannel = new DrugChannel(current,low,high,nh,xc50);
   channels_.push_back(dchannel);
   return;
}
double Drug::scaleFactor(const string current)
{
   return scaleFactor(current,concentration_);
}
double Drug::scaleFactor(const int ind)
{
   string current = channels_[ind]->current();
   return scaleFactor(current,concentration_);
}
double Drug::scaleFactor(const string current, const double concentration)
{
   double factor = -1.0;
   for (int ii=0; ii<channels_.size(); ii++)
      if (channels_[ii]->current() == current)
         factor = channels_[ii]->scalingFactor(concentration);

   assert(factor >= 0.0);
   return factor;
}
