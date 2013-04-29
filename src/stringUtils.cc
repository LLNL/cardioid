#include "stringUtils.hh"

using namespace std;

string concat(const vector<string> vv)
{
   if (vv.size() == 0)
      return "";
   string tmp = vv[0];
   for (unsigned ii=1; ii<vv.size(); ++ii)
      tmp += " " + vv[ii];
   return tmp;
}
