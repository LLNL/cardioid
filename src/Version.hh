#ifndef VERSION_HH
#define VERSION_HH

#include <string>
#include <iosfwd>

class Version
{
 public:
   static const Version& getInstance();

   std::ostream& versionPrint(std::ostream& out) const;
   std::string version() const {return buildTarget_ + " r" + svnVersion_;}
   std::string compileTime() const {return compileDate_ + " " + compileTime_;}
   const std::string& srcPath() const {return srcPath_;}
   const std::string& host()    const {return host_;}
   const std::string& user()    const {return user_;}
   
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
