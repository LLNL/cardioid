#include "getUserInfo.hh"
#include <unistd.h>
#include <pwd.h>

using namespace std;

#define HAS_GETPWUID
#if defined (BGP)
#undef HAS_GETPWUID 
#endif


// getpwuid Not supported on all platforms
#if defined(HAS_GETPWUID)
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
