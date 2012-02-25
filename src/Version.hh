#ifndef VERSION_HH
#define VERSION_HH

#include <string>
#include <iosfwd>

class Version
{
 public:
   static const Version& getInstance();

   std::ostream& versionPrint(std::ostream& out) const;
   
 private:
   Version();
   
   std::string svnVersion_;
   std::string compileTime_;
   std::string compileDate_;
   std::string srcPath_;
   std::string cxxFlags_;
   std::string cFlags_;
   std::string ldFlags_;
   std::string buildTarget_;
   std::string buildArch_;
   std::string host_;
   std::string user_;
};

#endif
