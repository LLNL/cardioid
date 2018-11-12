#include "drugFactory.hh"
#include <iostream>
#include <cassert>
#include <string>
#include <cmath>

#include "object_cc.hh"
#include "Drug.hh"
#include "Simulate.hh"

using namespace std;

/*!
  @page obj_DOSE DOSE object

  Used to model drug effects by rescaling currents.

  The following keywords are used to define the dosage of a given drug.  Details on how a
  given drug rescales the currents are given in the drug section.

  @beginkeywords

    @kw{drug, Name of DRUG block in object.data that specifies detailed information on
      how the drug changes the channel currents., No default}
    @kw{concentration, Effective free therapeutic plasma concentration (EFTPC) in micromols.,0.0}
    @endkeywords

    @subpage DOSE_drug

    @page DOSE_drug DRUG object

    Sets the four parameters that define how the given drug changes a given channel
    current as a function of concentration, for some or all of the currents.  The syntax is

    [current] = [low] [high] [NH] [XC50 (micromol)]

    Low sets the rescaling factor when the drug concentration is zero,
    high sets the rescaling factor when the drug concentration is
    effectively infinite.  NH is the Hill coefficient and XC50 sets
    the concentration of the rescaling midpoint (50\% rescaling).  The formula
    is:

    \f$\mbox{rescaling} = \mbox{low} + \frac{high-low}{1+(\frac{XC50}{concentration})^{NH}} \f$
    
    The most common cases will be inhibitory compounds (IC50):
    
      low = 1.0 and high = 0.0.

      and drugs that enhance current function (EC50):

      low = 1.0 and high = max rescaling.

    The full list of currents that can be rescaled is listed here:
      
    @beginkeywords
    @kw{I_K1, Parameters for K1 current,No default}
    @kw{I_Kr, Parameters for Kr current,No default}
    @kw{I_Ks, Parameters for Ks current,No default}
    @kw{I_Na, Parameters for Na current,No default}
    @kw{I_bNa, Parameters for bNa current,No default}
    @kw{I_CaL, Parameters for CaL current,No default}
    @kw{I_bCa, Parameters for bCa current,No default}
    @kw{I_to, Parameters for to current,No default}
    @kw{I_NaK, Parameters for NaK current,No default}
    @kw{I_NaCa, Parameters for NaCa current,No default}
    @kw{I_pCa, Parameters for pCa current,No default}
    @kw{I_pK, Parameters for pK current,No default}
    @kw{I_NaL, Parameters for NaL current,No default}
    @endkeywords
*/

Drug* drugFactory(const std::string& dosename, const Simulate& sim)
{
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  OBJECT* obj = objectFind(dosename, "DOSE");
  string drugobj;
  objectGet(obj, "drug", drugobj, "undefined");
  if (drugobj == "undefined")
    assert(false);
  double concentration;
  objectGet(obj, "concentration", concentration, "-1.0");
  if (concentration < 0.0)
     assert(false);

  Drug* drug = new Drug(dosename,concentration);

  // vector currentNames contains all possible currents allowed by
  // any reaction model.  Not all currents may be supported by
  // a given reaction model:  the reaction model is expected to throw an
  // error if called with a current name it doesn't recognize.
  vector<string> currentNames;
  currentNames.push_back("I_K1"); currentNames.push_back("I_Kr"); currentNames.push_back("I_Ks");
  currentNames.push_back("I_Na"); currentNames.push_back("I_bNa"); currentNames.push_back("I_CaL");
  currentNames.push_back("I_bCa"); currentNames.push_back("I_to"); currentNames.push_back("I_NaK");
  currentNames.push_back("I_NaCa"); currentNames.push_back("I_pCa"); currentNames.push_back("I_pK");
  currentNames.push_back("I_NaL"); currentNames.push_back("I_leak"); currentNames.push_back("I_up"); 
  currentNames.push_back("I_rel"); currentNames.push_back("I_xfer");
  obj = objectFind(drugobj, "DRUG");
  for (int ii=0; ii<currentNames.size(); ++ii)
  {
     vector<string> currentArgs;
     objectGet(obj, currentNames[ii], currentArgs);

     // read low, hi, nh, xc50 values
     if (currentArgs.size() > 3)
     {
        double low = atof(currentArgs[0].c_str());
        double high = atof(currentArgs[1].c_str());
        double nh = atof(currentArgs[2].c_str());
        double xc50 = atof(currentArgs[3].c_str());
        drug->addChannel(currentNames[ii],low,high,nh,xc50);
     }
  }
  return drug;
}
