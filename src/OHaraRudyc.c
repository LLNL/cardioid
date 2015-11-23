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
static char **stateName_=0; 

static int nCells_ =1;
static double dt_=0.0; 
static double **cellParms_=0;

static double *state_; 


void OHaraRudyDefineComps();
void OHaraRudyDefineDefaultTypes();

void OHaraRudyInit(double dt, int nCells)
{
   dt_ = dt; 
   nCells_ = nCells; 
   reversalPotentialsInit();
   OHaraRudyDefineComps(); 
   OHaraRudyDefineDefaultTypes();
   state_ = (double *)malloc(nCells_*nState_*sizeof(double)); 
   cellParms_ = (double **)malloc(nCells_*sizeof(char*)); 
   for (int i=0;i<nCells_;i++) 
   {
      cellParms_[i] = typeParms_[ENDO_CELL]; 
      double *state = state_ +i*nState_ ; 
      for (int j=0;j<nState_;j++) state[j] = initialState_[ENDO_CELL][j]; 
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
      if ( strcmp(name," MYBGBKC_INa")         ==0) cInit[k++] = MYBGBKC_INaInit; 
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

   int nState=0; 
   int nParms=0; 
   for (int i=0;i<nComp_;i++)
   {
      VARINFO *varInfo = info_[i].varInfo; 
      for (int j=0;j<info_[i].nVar;j++) 
      {
         if (varInfo[j].type == PSTATE_TYPE)    stateName_[nState++] = strdup(varInfo[j].name); 
         if (varInfo[j].type == PARAMETER_TYPE) typeName_[nParms++] = strdup(varInfo[j].name); 
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
void  OHaraRudyPut(double dt, int nCells, const double *Vm, const double *iStim)
{
      dt_= dt; 
      nCells_ = nCells; 
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
   DERIVED derived; 
   for (int ii=0;ii<nCells_;ii++) 
   {
      CONCENTRATIONS *concentrations     = (CONCENTRATIONS *)(state_+ii*nState_+privateStateOffset_[1]); 
      double *cellParms = cellParms_[ii]; 
      double Nai = concentrations->Nai; 
      double Ki  = concentrations->Ki; 

      reversalPotentials(Nai,Ki,&derived); 

      for (int i=2;i<=17;i++) 
      {
         info_[i].func(cellParms+privateParmsOffset_[i], state_, privateStateOffset_[ i], &derived, dt_);
      }
      info_[1].func(cellParms+privateParmsOffset_[0], state_, privateStateOffset_[1], &derived, dt_);
/*
      double *I =  (double *)&(derived.I.NaCai); 
      for (int i=0;i<nCurrents_;i++)
         printf("%2d %10s %20.12f\n",i,currentNames_[i],I[i]); 
         
*/
   }
}
#ifdef SA
#include <unistd.h>
typedef struct parms_st 
{int loopMax, 
printInterval;
char *filename; 
} PARMS;
PARMS comandLineArgs(int iargc, char *argv[])
{
   PARMS parms={10000000,100000," "}; 
   char opt; 
   while ((opt = getopt(iargc, argv, "l:p:")) != -1) 
   {
      switch (opt) 
     {
         case 'l':
            parms.loopMax = atoi(optarg); 
            break;
         case 'p':
            parms.printInterval = atoi(optarg); 
            break;
         case '?':
         break;
         default:
         {

         }
      }
   if (optind == iargc-1) parms.filename = strdup(argv[optind]); 
   printf("opt=%s %d %d\n",optarg,optind,iargc); 
   }
   
   
   iargc -= optind;
   argv += optind;
   return parms; 
}

int main(int iargc, char *argv[])
{
   PARMS parms = comandLineArgs(iargc, argv);
   double dt = 0.01; 
   double nCells = 1; 
   OHaraRudyInit(dt,nCells);
   if (iargc == 2 )
   {
      char *filename = argv[1]; 
      FILE *file = fopen(filename,"r"); 

      int iS=0; 
      int iP=0; 
      for (int i=0;i<nState_+nParms_;i++) 
      {
         char line[256]; 
         double value; 
         fgets(line,256,file); 
         sscanf(line,"%*s %*s %lg",&value); 
         if (i < nState_) state_[iS++] = value;
         else cellParms_[0][iP++] = value;
      }
      fclose(file); 
   }
   double t=0.0; 
   OHaraRudyCalc();
   VOLTAGE *voltage0 = (VOLTAGE *)(state_+0*nState_) ;
#ifdef doEnd
   printf("Vm              : %2d %22.15e %22.15e\n",0,voltage0->Vm,voltage0->dVm); 
   exit(0); 
#endif
   int loop = 0; 
   for (loop ;loop< parms.loopMax;loop++) 
   {
      if (loop%parms.printInterval==0)
      {
         printf("%10d %10.3f %24.15f %24.15f\n",loop,loop*dt,voltage0->Vm,voltage0->dVm); 
         fflush(stdout); 
      }
      for (int i=0;i<nCells_;i++)
      {
         VOLTAGE *voltage = (VOLTAGE *)(state_+i*nState_) ;
         voltage->Vm += voltage->dVm*dt;
      }
      OHaraRudyCalc();
   }
   for (int i=0;i<nCells_;i++) printf("%10d %10.3f %24.15f %24.15f\n",loop,loop*dt,voltage0->Vm,voltage0->dVm); 
   FILE *initData = fopen("init.data","w"); 
   for (int i=0;i<nState_;i++) fprintf(initData,"%2d %16s  %22.15e\n",i,stateName_[i],state_[i]); 
   for (int i=0;i<nParms_;i++) fprintf(initData,"%2d %16s  %22.15e\n",i,typeName_[i],cellParms_[0][i]); 
}
#endif
