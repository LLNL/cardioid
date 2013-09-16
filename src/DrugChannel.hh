#ifndef DRUGCHANNEL_HH
#define DRUGCHANNEL_HH
#include <string>
using std::string;

class DrugChannel
{
  public:
    DrugChannel(string current, double low, double high, double nh, double xc50);
    ~DrugChannel();
    double scalingFactor(double concentration);
    string current(void) { return current_; }
    
  private:
    string current_;
    double low_;
    double high_;
    double nh_;
    double xc50_;
};
#endif
