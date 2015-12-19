#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "OHaraRudy.h"
#include "OHaraRudyInit.h"

static char *currentNames_[]={"INaCai","INaCass","INaK","INab","ICab","IKb","IpCa","INaFast","INaL","Ito","IKr","IKs","IK1","ICa",""};
static char *currentModels_[]={"OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy","OHaraRudy",""};

static int nCurrents_ = 14; 
static int nComp_ ; 
static int * privateStateOffset_=0;; 
static int * privateParmsOffset_=0;; 
static COMPONENTINFO* info_; 

static int nTypes_; 

static int nParms_=0; 
static double *typeParmsBuffer_=0; 
static double **typeParms_=0;
static char   **typeName_=0;

static int nState_=0; 
static double *initialStateBuffer_=0; 
static double **initialState_=0;
static int *cellType_; 
static char **stateName_=0; 
static int  *stateCompMap_=0; 
static int  *parmCompMap_=0; 


static double time_=0.0;
static int loop_=0;
static int nCells_ =1;
static double dt_=0.0; 
static double **cellParms_=0;

static double *state_=0; 
static double *dState_=0; 
static DERIVED *derived_=0;

void OHaraRudyDefineComps();
void OHaraRudyDefineDefaultTypes();

void OHaraRudyInit(double dt, int nCells, int *cellType)
{
   dt_ = dt; 
   nCells_ = nCells; 
   reversalPotentialsInit();
   OHaraRudyDefineComps(); 
   OHaraRudyDefineDefaultTypes();
   state_  = (double *)malloc(nCells_*nState_*sizeof(double)); 
   dState_ = (double *)malloc(nCells_*nState_*sizeof(double)); 
   cellType_ = (int *)malloc(nCells_*sizeof(int)); 
   derived_ = (DERIVED *)malloc(nCells_*sizeof(DERIVED)); 
   cellParms_ = (double **)malloc(nCells_*sizeof(char*)); 
   for (int i=0;i<nCells_;i++) 
   {
      cellType_[i] = cellType[i]; 
      cellParms_[i] = typeParms_[cellType[i]]; 
      double *state = state_ +i*nState_ ; 
      double *dState = dState_ +i*nState_ ; 
      derived_[i].dState = dState; 
      for (int j=0;j<nState_;j++) state[j] = initialState_[cellType[i]][j]; 
   }
}
void OHaraRudyDefineComps()
{
   nComp_     = nCurrents_ + 4; 
   typedef COMPONENTINFO (*CINIT)(); 
   CINIT *cInit = (CINIT *)calloc(nComp_,sizeof(CINIT)); 
   int k=0; 
   cInit[k++] = OHaraRudy_VoltageInit;   // Concentration must be second and CaMK trap third. 
   cInit[k++] = OHaraRudy_ConcentInit;   
   cInit[k++] = OHaraRudy_CaMKtrapInit; 
   for (int i=0;i<nCurrents_;i++)   // Order of currents initialization does not matter. 
   {
      char name[256]; 
      sprintf(name,"%s_%s",currentModels_[i],currentNames_[i]) ;
      if ( strcmp(name,"OHaraRudy_INaCai")     ==0) cInit[k++] = OHaraRudy_INaCaiInit; 
      if ( strcmp(name,"OHaraRudy_INaCass")    ==0) cInit[k++] = OHaraRudy_INaCassInit; 
      if ( strcmp(name,"OHaraRudy_INaK")       ==0) cInit[k++] = OHaraRudy_INaKInit; 
      if ( strcmp(name,"OHaraRudy_INab")       ==0) cInit[k++] = OHaraRudy_INabInit; 
      if ( strcmp(name,"OHaraRudy_ICab")       ==0) cInit[k++] = OHaraRudy_ICabInit; 
      if ( strcmp(name,"OHaraRudy_IKb")        ==0) cInit[k++] = OHaraRudy_IKbInit; 
      if ( strcmp(name,"OHaraRudy_IpCa")       ==0) cInit[k++] = OHaraRudy_IpCaInit; 
      if ( strcmp(name,"OHaraRudy_INaFast")    ==0) cInit[k++] = OHaraRudy_INaFastInit; 
      if ( strcmp(name,"OHaraRudy_INaL")       ==0) cInit[k++] = OHaraRudy_INaLInit; 
      if ( strcmp(name,"OHaraRudy_Ito")        ==0) cInit[k++] = OHaraRudy_ItoInit; 
      if ( strcmp(name,"OHaraRudy_IKr")        ==0) cInit[k++] = OHaraRudy_IKrInit; 
      if ( strcmp(name,"OHaraRudy_IKs")        ==0) cInit[k++] = OHaraRudy_IKsInit; 
      if ( strcmp(name,"OHaraRudy_IK1")        ==0) cInit[k++] = OHaraRudy_IK1Init; 
      if ( strcmp(name,"OHaraRudy_ICa")        ==0) cInit[k++] = OHaraRudy_ICaInit; 

      if ( strcmp(name,"OHaraRudyMod_INaFast") ==0) cInit[k++] = OHaraRudyMod_INaFastInit; 
      if ( strcmp(name,"RTYSC14A_IKr")         ==0) cInit[k++] = RTYSC14A_IKrInit; 
      if ( strcmp(name,"MYBGBKC_INaL")         ==0) cInit[k++] = MYBGBKC_INaInit; 
      if ( strcmp(name,"MYBGBKC_INaFast")         ==0) cInit[k++] = null_INullInit; 
   }
   cInit[k++] = OHaraRudy_FluxesInit;      // must be initialize after ICa
   assert(k == nComp_); 

   info_=(COMPONENTINFO *)malloc(sizeof(COMPONENTINFO)*nComp_); 
   int stateOffset[nComp_]; 
   int parmsOffset[nComp_]; 

   nState_ = 0; 
   nParms_ = 0; 
   for (int i=0;i<nComp_;i++)
   {
      stateOffset[i] =  nState_;
      parmsOffset[i] =  nParms_;
      info_[i] = cInit[i](); 
      VARINFO *varInfo = info_[i].varInfo; 
      for (int j=0;j<info_[i].nVar;j++) 
      {
         unsigned int index = (i << 16) + varInfo[j].index;
         if (varInfo[j].type == PSTATE_TYPE)    nState_++; 
         if (varInfo[j].type == PARAMETER_TYPE) nParms_++;
      }
   } 
   privateStateOffset_=(int *)malloc(sizeof(int)*nState_); 
   stateName_ = (char **)malloc(nState_*sizeof(char *)); 
   for (int j=0;j<nState_;j++) privateStateOffset_[j] = stateOffset[j]; 

   privateParmsOffset_=(int *)malloc(sizeof(int)*nParms_); 
   typeName_ = (char **)malloc(nParms_*sizeof(char *)); 
   for (int j=0;j<nParms_;j++) privateParmsOffset_[j] = parmsOffset[j]; 
   stateCompMap_ = (int *)malloc(nState_*sizeof(int)); 
   parmCompMap_ = (int *)malloc(nParms_*sizeof(int)); 

   int nState=0; 
   int nParms=0; 
   for (int i=0;i<nComp_;i++)
   {
      VARINFO *varInfo = info_[i].varInfo; 
      for (int j=0;j<info_[i].nVar;j++) 
      {
         if (varInfo[j].type == PSTATE_TYPE)    {stateName_[nState] = strdup(varInfo[j].name); stateCompMap_[nState++] = j;}
         if (varInfo[j].type == PARAMETER_TYPE) {typeName_[nParms] = strdup(varInfo[j].name);  parmCompMap_[nParms++] = j; }
      }
   } 
}
void OHaraRudyDefineDefaultTypes()
{
   nTypes_= 3; 
   initialStateBuffer_ = (double *)malloc(nTypes_*nState_*sizeof(double)); 
   typeParmsBuffer_    = (double *)malloc(nTypes_*nParms_*sizeof(double)); 
   initialState_       = (double **)malloc(nTypes_*sizeof(double *)); 
   typeParms_          = (double **)malloc(nTypes_*sizeof(double *)); 
   for (int i = 0; i < nTypes_; i++) 
   {
      double *initialState = initialState_[i] = initialStateBuffer_+nState_*i; 
      double *typeParms    = typeParms_[i]    = typeParmsBuffer_+nParms_*i; 
      int nState=0; 
      int nParms=0; 
      for (int k = 0; k < nComp_; k++) 
      {
         VARINFO *varInfo = info_[k].varInfo; 
         for (int j = 0; j < info_[k].nVar; j++) 
         {
            double value; 
            switch (i) 
            {
               case ENDO_CELL: 
                  value  = varInfo[j].defaultValueENDO;
                  break; 
               case EPI_CELL: 
                  value  = varInfo[j].defaultValueEPI;
                  break; 
               case M_CELL: 
                  value  = varInfo[j].defaultValueM;
                  break; 
               default: 
                  value = 0.0; 
            }
            if (varInfo[j].type == PSTATE_TYPE)    initialState_[i][nState++] = value ;
            if (varInfo[j].type == PARAMETER_TYPE) typeParms_[i][nParms++] =   value ;
         }
      }
   }
}
int OHaraRudyGet_nComp() {return nComp_;} 
COMPONENTINFO* OHaraRudyGet_compInfo() {return info_;} 
void OHaraRudySetValue(int iCell, int varHandle, double value) 
{
   int iComp = varHandle >> 16; 
   int index = varHandle &  0xffff; 
   double *pstate = state_+ nState_*iCell  + privateStateOffset_[iComp]; 
   double *cellParms = cellParms_[iCell] + privateParmsOffset_[iComp]; 
   info_[iComp].access(WRITE,index,&value,cellParms,pstate); 
}
double OHaraRudyGetValue(int iCell, int varHandle) 
{
   double value; 
   int iComp = varHandle >> 16; 
   int index = varHandle &  0xffff; 
   double *pstate = state_+ nState_*iCell  + privateStateOffset_[iComp]; 
   double *cellParms = cellParms_[iCell] + privateParmsOffset_[iComp]; 
   info_[iComp].access(READ,index,&value,cellParms,pstate); 
   return value; 
}
void  OHaraRudyPut(double dt, int nCells, const double *Vm, const double *iStim)
{
      dt_= dt; 
      nCells_ = nCells; 
      double *vv= (double *)Vm; 
      for (int ii=0;ii<nCells_;ii++) 
      {
         double *state     = state_+ii*nState_; 
         VOLTAGE *voltage = (VOLTAGE*)state; 
         voltage->Vm = Vm[ii]; 
         voltage->iStim = iStim[ii]; 
      }
}
void  OHaraRudyGet(double *dVm)
{
      for (int ii=0;ii<nCells_;ii++) 
      {
         double *state     = state_+ii*nState_; 
         VOLTAGE *voltage = (VOLTAGE*)state; 
         dVm[ii] = voltage->dVm; 
      }
}
void  OHaraRudyCalc()
{
   for (int ii=0;ii<nCells_;ii++) 
   {
      CONCENTRATIONS *concentrations     = (CONCENTRATIONS *)(state_+ii*nState_+privateStateOffset_[1]); 
      double *cellParms = cellParms_[ii]; 
      double Nai = concentrations->Nai; 
      double Ki  = concentrations->Ki; 

      reversalPotentials(Nai,Ki,derived_+ii); 

      for (int i=2;i<=17;i++) 
      {
         info_[i].func(cellParms+privateParmsOffset_[i], state_, privateStateOffset_[ i], derived_+ii, dt_);
      }
      info_[1].func(cellParms+privateParmsOffset_[0], state_, privateStateOffset_[1], derived_+ii, dt_);
/*
      double *I =  (double *)&(derived_[ii].I.NaCai); 
      for (int i=0;i<nCurrents_;i++)
         printf("%2d %10s %20.12f\n",i,currentNames_[i],I[i]); 
         
*/
   }
}
#ifdef SA
#include <unistd.h>
typedef struct parms_st 
{
   int loopMax; 
   int printInterval;
   int fixVmFlag; 
   double fixVm;
   double dt; 
   char * initialConditions; 
   char *currents; 
} PARMS;
PARMS commandLineArgs(int iargc, char *argv[])
{
   PARMS parms={100000000,100000,0,0.0,0.0,NULL,NULL}; 
   char opt; 
   while ((opt = getopt(iargc, argv, "l:p:i:v:I:h:")) != -1) 
   {
      switch (opt) 
      {
         case 'h':
            parms.dt = atof(optarg); 
            break;
         case 'I':
            parms.currents = strdup(optarg); 
            break;
         case 'i':
            parms.initialConditions = strdup(optarg); 
            break;
         case 'l':
            parms.loopMax = atoi(optarg); 
            break;
         case 'p':
            parms.printInterval = atoi(optarg); 
            break;
         case 'v':
            parms.fixVmFlag =1; 
            parms.fixVm = atof(optarg); 
            break;
         case '?':
            break;
         default:
            {

            }
      }
   }


   iargc -= optind;
   argv += optind;
   return parms; 
}
void readCurrents(char *filename)
{
   FILE *file=fopen(filename,"r"); 
   char name[32];
   char model[32];
   for( int i=0;i<14;i++)
   {
      
      fscanf(file,"%s %s",name,model); 
      currentNames_[i]=strdup(name);
      currentModels_[i]=strdup(model);
   }
   fclose(file); 
   
}
void doIO(int loop, double t, double *state, double *dState)
{
         printf("%10.3f %24.15f ",t,state[0]); 
         for (int i=7;i<7;i++) 
         printf("%24.15f ",state[i]); 
         printf("\n"); 
         fflush(stdout); 
}
void checkpoint(int loop, double t, int nCell, int nState, double *state, double *dState)
{
   static FILE *file = 0;
   if (file == 0) file = fopen("init.data","w"); 
   rewind(file); 
   fprintf(file,"%d  %12.6f\n",loop,t); 
   for (int j=0;j<nCell;j++)
   {
      int k=0; 
      for (int i=0;i<nComp_;i++)
      {
         char *compName = info_[i].compName; 
         VARINFO *varInfo = info_[i].varInfo; 
         for (int j=0;j<info_[i].nVar;j++) 
         {
            if (varInfo[j].type == PSTATE_TYPE) 
            {
               char *stateName= varInfo[j].name; 
               if (i> 0 | j == 0 ) fprintf(file,"%10s %16s  %22.15e %22.15e\n",compName,stateName,state[k],dState[k]); 
               k++; 
            }
         }
      } 
   }
   fflush(file); 
}
void readInitialConditions(char *filename) 
{
   FILE *file = fopen(filename,"r");
   int iS=0; 
   int iP=0; 
   char line[256]; 
   fgets(line,256,file); 
   sscanf(line,"%d %g",&loop_,&time_); 
   for (int i=0;i<nState_;i++) 
   {
      double value; 
      fgets(line,256,file); 
      sscanf(line,"%*s %*s %lg",&value); 
      if (i < nState_) state_[iS++] = value;
      else cellParms_[0][iP++] = value;
      if (i==0) {i+=2; iS+=2;}
   }
   fclose(file); 
}

int main(int iargc, char *argv[])
{
   PARMS parms = commandLineArgs(iargc, argv);
   double dt=parms.dt ;
   int nCells = 1; 
   if (parms.currents != NULL)    readCurrents(parms.currents); 
   int cellType[nCells]; 
   for (int ii=0;ii<nCells;ii++) cellType[ii] = ENDO_CELL; 
   OHaraRudyInit(dt,nCells,cellType);
   if (parms.initialConditions != NULL)    readInitialConditions(parms.initialConditions); 
   VOLTAGE *voltage0 = (VOLTAGE *)(state_+0*nState_) ;
   double amp = -80; 
   if (parms.fixVmFlag) 
   {
      voltage0->Vm=parms.fixVm; 
      amp = 0.0 ; 
   }
   double state[nState_*nCells_]; 
   for ( ;loop_< parms.loopMax;loop_++) 
   {
      double iStim =0; 
      if (0<=fmod(time_,1000.0) && fmod(time_,1000.0)<0.5 )iStim = amp; 

      for (int i=0;i<nState_*nCells_;i++) state[i] = state_[i]; 
      for (int i=0;i<nCells_;i++)
      {
         VOLTAGE *voltage = (VOLTAGE *)(state_+i*nState_) ;
         voltage->iStim=iStim;
         voltage->Vm += (voltage->dVm-iStim)*dt;
         if (parms.fixVmFlag) voltage->Vm = parms.fixVm; 
      }
      OHaraRudyCalc();
      double time = time_-parms.loopMax*dt+1000; 
      if (loop_%parms.printInterval==0 && time >=-0.0000000001) 
      {
         doIO(loop_,time,state,dState_); 
         checkpoint(loop_,time_,nCells_,nState_,state,dState_); 
      }
      time_ += dt; 
   }
   dt = 0.0; 
   OHaraRudyCalc();
   double   time = time_-parms.loopMax*dt+1000; 
   if (loop_%parms.printInterval==0) doIO(loop_,time_,state,dState_); 
   //   if (loop > 1000) checkpoint(loop_,time_,nCells_,nState_,state_,dState_); 
   //   for (int i=0;i<nParms_;i++) fprintf(initData,"%2d %16s  %22.15e\n",i,typeName_[i],cellParms_[0][i]); 
}
#endif
