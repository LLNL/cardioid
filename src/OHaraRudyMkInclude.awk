BEGIN { 
setValuehh     = "OHaraRudySetValue.hh"
getValuehh     = "OHaraRudyGetValue.hh"
enumhh         = "OHaraRudyEnum.hh"
nameArrayh        = "OHaraRudyNameArray.h"
getHandleMaphh = "OHaraRudyGetHandleMap.hh"
flag=0;
units = 2

# preamble for setValue
print "void OHaraRudy::setValue(int varHandle, double value)" > setValuehh
print "{" >> setValuehh
print " " >> setValuehh
print "   switch (varHandle)" >>  setValuehh
print "   {" >> setValuehh 
print "      case undefinedName   : assert(false)                   ; break;" >> setValuehh 

# preamble for getValue
print "double OHaraRudy::getValue(int varHandle) const" > getValuehh
print "{" >> getValuehh
print " " >> getValuehh
print "   switch (varHandle)" >>  getValuehh
print "   {" >> getValuehh 
print "      case undefinedName   : assert(false)                 ; break;" >> getValuehh 

# preamble for getHandleMap
print "/* Remember that down in the cell models the units don't necessarily" > getHandleMaphh
print " *  correspond to the internal units of Cardioid.  The units in this map" >> getHandleMaphh
print " *  are the units the cell model expects the variables to have." >> getHandleMaphh
print "*/" >> getHandleMaphh
print "HandleMap& OHaraRudy::getHandleMap()" > getHandleMaphh
print "{" > getHandleMaphh
print "  static HandleMap handleMap;" >> getHandleMaphh
print "   if (handleMap.size() == 0)" >> getHandleMaphh
print "   {" >> getHandleMaphh

#preamble for enum
printf "enum VarHandle\n{\n" > enumhh
printf "   %s\n","undefinedName = -1," >> enumhh
}

{
sub(";" ," ",$0)
if ( flag > 0  ) 
{
if (index($0,"STATE")>0) flag=0;  
if (index($0,"CELLPARMS")>0) flag=0;  
}
if (flag>0 && $1 == "double") 
{
   printf "      case %-16s: %s%-16s = value; break;\n", $2, variable, $2  >> setValuehh
   printf "      case %-16s: return %s%-16s; break;\n", $2, variable, $2  >> getValuehh
   printf "      handleMap[%-16s] = CheckpointVarInfo(%-16s, %-5s, %-8s);\n", "\""$2"\"",$2, $5, "\""$4"\"" >> getHandleMaphh
   printf "   %-13s\n",$2"," >> enumhh
   printf "   %s\n", "\""$2"\"," >> nameArrayh
}
if ( flag == 0 && $3  == "state_st") 
{ 
   printf "char *stateVarNames[] = \n{\n" >  nameArrayh 
   flag=1;  
   variable = "state_."
}
if ( flag == 0 && $3  == "cellParms_st") 
{
   printf "// Start of cell parmeters \n" >> enumhh
   printf "};\n" >>  nameArrayh 
   printf "char *cellParmNames[] = \n{\n" >  nameArrayh 
   flag=1; 
   variable = "cellParms_->"; 
}
}
END {
# postamble for setValue
print "      case nVars           : assert(false)                  ; break;"  >> setValuehh
print "   }"  >> setValuehh
print "}"  >> setValuehh

# postamble for getValue
print "      case nVars           : assert(false)                  ; break;"  >> getValuehh
print "   }"  >> getValuehh
print "   return 0.;"  >> getValuehh
print "}"  >> getValuehh

# postamble for getHandleMap
print "      assert(handleMap.size() == nVars);" >> getHandleMaphh
print "   }" >> getHandleMaphh
print "   return handleMap;" >> getHandleMaphh
print "}" >> getHandleMaphh

# postamble for enum 
printf "   nVars\n};\n" >> enumhh

# postable for nameArray
printf "};\n" >>  nameArrayh 
}


