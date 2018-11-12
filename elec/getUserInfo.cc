#include "getUserInfo.hh"
#include <unistd.h>
#include <pwd.h>

using namespace std;


static char unknown[] = {'u','n','k','n','o','w','n','\0'} ; 
#if defined (BGP)
#define NO_GETPWUID
#endif
// getpwuid Not supported on all platforms
// Make the default behavior to be has getpwuid independent of system arch.
// Always do not have getpwuid on BGP

#if not defined (NO_GETPWUID)
string getUserName()
{
   struct passwd *pwd = getpwuid(getuid());
   string name("unknown"); 
   if (pwd != NULL) name = pwd->pw_name; 
   return name; 
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
