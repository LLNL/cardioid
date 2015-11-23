BEGIN { 
typeMap["PARAMETER"] = "PARAMETER_TYPE"; 
typeMap["PSTATE"] = "PSTATE_TYPE"; 
nVar = 0; ns=0; np=0;  nt=0 ; 
}
{
   if ( !($1 == "\/\/" || NF == 0))  
   {
       nVar++;
       model[nVar] = $1 
       type[nVar]  = $2
       root[nVar]  = $1"_"$3
       var[nVar]   = $4
       default[nVar] = defaultString($5); 
       units[nVar] = $6
       if (root[nVar] != root[nVar-1]) 
       { 
          rootList[++nt] = root[nVar]; 
       }
   }
}
END {
printf "/* ******   THIS FILE IS AUTOMATICLY GENERATED.  DO NOT EDIT!! ******* */\n"
printf "#include <assert.h>\n"; 
printf "enum enumIndex{"  ;
for (i=1;i<=nVar;i++) printf " %sIndex,",var[i];
printf " %s","nVar";
printf "};\n"

if (nVar > 0) {
printf "static VARINFO varInfo[] =\n";
printf "{\n"; 
for (i=1;i<nVar;i++) {printf "   {\"%s\",%s,%sIndex,%s,\"%s\"},\n",var[i],typeMap[type[i]],var[i],default[i],units[i];}
printf "   {\"%s\",%s,%sIndex,%s,\"%s\"}\n",var[i],typeMap[type[i]],var[i],default[i],units[i];
printf "};\n";
}


np = 0; 
for (i=1;i<=nVar;i++) if ( type[i] == "PARAMETER") { parms[++np] = var[i]; parmsDefault[np] = default[i]; }
if (np > 0) 
{
  printf "typedef struct parameters_str { double ";
  for (i=1;i<np;i++) printf " %s,", parms[i];
  printf " %s", parms[i];
  printf ";} PARAMETERS;\n"; 

#  printf "static PARAMETERS parmsDefault = { ";
#  for (i=1;i<np;i++) printf " %s,", parmsDefault[i];
#  printf " %s", parmsDefault[i];
#  printf "};\n"; 
}
ns = 0; 
for (i=1;i<=nVar;i++) if ( type[i] == "PSTATE") { pstate[++ns] = var[i]; pstateDefault[ns] = default[i]; }
if (ns > 0) 
{
  printf "typedef struct pstate_str { double ";
  for (i=1;i<ns;i++) printf " %s,", pstate[i];
  printf " %s", pstate[i];
  printf ";} PSTATE;\n"; 

#  printf "static PSTATE  pstateDefault = { ";
#  for (i=1;i<ns;i++) printf " %s,", pstateDefault[i];
#  printf " %s", pstateDefault[i];
#  printf "};\n"; 
}

printf "void %sFunc(CELLPARMS *parmsPtr, double *state, int pOffset, DERIVED *derived, double dt );\n",root[1];

printf "void %sAccess(int type,int index,double *value, double  *parmsPtr, double *statePtr)\n",root[1]; 
print "{\n"; 
if (ns  > 0) printf "   PSTATE *state = (PSTATE *)statePtr;\n" 
if (np  > 0) printf "   PARAMETERS *parms = (PARAMETERS *)parmsPtr;\n" 
printf "   if (type == READ)\n";
printf "   {\n";
printf "      switch (index)\n";
printf "      {\n";
            for (i=1;i<=nVar;i++) 
            {
                 typeName = "parms";
                 if (type[i] == "PSTATE") typeName = "state"; 
printf "         case %sIndex:\n",var[i];
printf "            *value = %s->%s; \n",typeName,var[i];
printf "            break;\n";
            }
printf "         default:\n";
printf "            assert(0); \n";
printf "      }\n";
printf "   }\n";
printf "   if (type == WRITE)\n";
printf "   {\n";
printf "      switch (index)\n";
printf "      {\n";
            for (i=1;i<=nVar;i++) 
            {
                 typeName = "parms";
                 if (type[i] == "PSTATE") typeName = "state"; 
printf "         case %sIndex:\n",var[i];
printf "            %s->%s = *value;\n",typeName,var[i];
printf "            break;\n";
            }
printf "            assert(0); \n";
printf "      }\n";
printf "   }\n";
printf "}\n";
constFunc = root[1]"Constants()";
cmd = sprintf( "fgrep '%s' %s.c ", constFunc, root[1]) ;
cmd | getline line 
constExist=0; 
if (length(line) > 0) constExist=1; 

if (constExist) printf "void   %s;\n",constFunc
printf "COMPONENTINFO %sInit()\n{\n", root[1];
printf "   COMPONENTINFO info;\n"
if (constExist) printf "   %s;\n",constFunc
printf "   if (FRT  < 0) FRT = F/(R*T);\n"
if (ns > 0) printf "   info.pStateSize = sizeof(PSTATE);\n"; 
if (np > 0) printf "   info.parmsSize = sizeof(PARAMETERS);\n"; 
printf "   info.nVar = %d;\n",nVar; 
printf "   info.varInfo = varInfo;\n"
printf "   info.func = %sFunc;\n",root[1]; 
printf "   info.access = %sAccess;\n",root[1]; 
printf "   return info;\n}\n"; 

print "#ifdef doEnd"
print "#define ENDCODE() {\\" 
k=0; 
for (i=1;i<=nVar;i++) 
{
     if (type[i] != "PSTATE") continue; 
     line = sprintf( "     printf(\"%-16s: %%2d %%22.15e %%22.15e\\n\",pOffset+%d,pState->%s,d%s);", var[i],k,var[i],var[i]) ;
     k++; 
     print line,"\\"; 
}
print "}"; 
print "#else"
print "#define ENDCODE()"
print "#endif"
}
function defaultString(string)
{
 n=split(string,field,"[:\*]")
 s0=length(field[1])+1;
 fs0 = substr(string,s0,1)
 s1=s0+length(field[2])+1;
 fs1 = substr(string,s1,1)
 ENDO = field[1]; 
 EPI=M=ENDO; 
 if (n == 2) 
 {
   EPI = field[2]; 
   if (fs0 == "*") EPI = ENDO*field[2]; 
 }
 if (n == 3) 
 {
   EPI = field[2]; 
   M   = field[3]; 

   if (fs0 == "*") EPI  = ENDO*field[2]; 
   if (fs1 == "*") M    = ENDO*field[3]; 
 }
 return sprintf("%.14g,%.14g,%.14g",ENDO,EPI,M); 
}
