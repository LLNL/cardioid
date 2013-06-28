#include "getUserInfo.hh"
#include <unistd.h>
#include <pwd.h>

using namespace std;


#if defined (BGP)
#define NO_GETPWUID
#endif
// getpwuid Not supported on all platforms
// Make the default behavior to be has getpwuid independent of system arch.
// Always do not have getpwuid on BGP

#if not defined (NO_GETPWUID)
string getUserName()
{
   return getpwuid(getuid())->pw_name;
}
#else
string getUserName()
{
   return "unavailable";
}
#endif

string getHostName()
{
   char tmp[256];
   gethostname(tmp, 256);
   return tmp;
}
