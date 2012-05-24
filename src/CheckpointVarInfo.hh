#ifndef CHECKPOINT_VAR_INFO_HH
#define CHECKPOINT_VAR_INFO_HH

#include <map>
#include <string>

class CheckpointVarInfo
{
 public:
   /** -1 is the handle value for an undefined name */
   CheckpointVarInfo()
   : handle_(-1), checkpoint_(false), unit_("1")
   {};
   CheckpointVarInfo(int handle, bool checkpoint, std::string unit)
   : handle_(handle), checkpoint_(checkpoint), unit_(unit)
   {};
   
   int         handle_;
   bool        checkpoint_;
   std::string unit_; // output unit
};

typedef std::map<std::string, CheckpointVarInfo> HandleMap; 

#endif
