#ifndef DRUG_HH
#define DRUG_HH
#include "DrugChannel.hh"
#include <vector>
#include <string>
using std::vector;
using std::string;

class Drug
{
 public:
    Drug(const string name, const double concentration);
    ~Drug();
    void addChannel(const string current, const double low,
                    const double high, const double nh,
                    const double xc50);

    // return rescaling factor using stored concentration
    double scaleFactor(const string current);                        
    double scaleFactor(const int ind);                        
    // return rescaling factor for given concentration
    double scaleFactor(const string current, double concentration);  

    int nChannels(void) { return channels_.size(); };
    string name(void) { return name_; };
    double concentration(void) { return concentration_; };
    string current(int ind) { return channels_[ind]->current(); };
    
 private:
    string name_;                        // name of drug
    double concentration_;               // concentration of drug in micromol
    vector<DrugChannel*> channels_;      // channel currents modified by drug
};
#endif
