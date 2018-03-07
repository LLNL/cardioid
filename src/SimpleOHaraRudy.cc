/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "SimpleOHaraRudy.hh"
#include "object_cc.hh"
#include "Anatomy.hh"
#include "reactionFactory.hh"
#include <cmath>

#define pi 3.141592653589793238462643383279502884197169399375105820974944592307816406286

using namespace std;


#define setDefault(name, value) objectGet(obj, #name, reaction->name, #value)
   
REACTION_FACTORY(SimpleOHaraRudy)(OBJECT* obj, const double, const int numPoints, const ThreadTeam&)
{
   SimpleOHaraRudy::ThisReaction* reaction = new SimpleOHaraRudy::ThisReaction(numPoints);

   //override the defaults
   //EDIT_FLAGS
   setDefault(celltype, 0);
   setDefault(useINaFromTT06, 0);

   //EDIT_PARAMETERS
   setDefault(GCaB, 6.0643e-4);    // [uA/uF] 3
   setDefault(JrelStiffConst, 0.005);  // ms
      
   return reaction;
}
#undef setDefault


namespace SimpleOHaraRudy 
{

inline double pow(const double x, const int p)
{
   double ret=1;
   if (p > 0) 
   {
      for (int ii=0; ii<p; ii++) 
      {
         ret *= x;
      }
   }
   else
   {
      for (int ii=0; ii<-p; ii++) 
      {
         ret /= x;
      }
   }
   return ret;
}
   
string ThisReaction::methodName() const
{
   return "SimpleOHaraRudy";
}
const char* varNames[] = 
{
   //EDIT_STATE
   "nai",
   "nass",
   "ki",
   "kss",
   "cai",
   "cass",
   "cansr",
   "cajsr",
   "m",
   "hf",
   "hs",
   "j",
   "hsp",
   "jp",
   "mL",
   "hL",
   "hLp",
   "a",
   "iF",
   "iS",
   "ap",
   "iFp",
   "iSp",
   "d",
   "ff",
   "fs",
   "fcaf",
   "fcas",
   "jca",
   "nca",
   "ffp",
   "fcafp",
   "xrf",
   "xrs",
   "xs1",
   "xs2",
   "xk1",
   "Jrelnp",
   "Jrelp",
   "CaMKt"
};
#define NUMVARS (sizeof(varNames)/sizeof(char*))

int getVarOffset(const std::string& varName)
{
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      if (varNames[ivar] == varName) 
      {
         return ivar;
      }
   }
   return -1;
}

void assertStateOrderAndVarNamesAgree(void)
{
   State s;
#define checkVarOrder(x) assert(reinterpret_cast<double*>(&s)+getVarOffset(#x) == &s . x)

   int STATIC_ASSERT_checkAllDouble[(NUMVARS == sizeof(s)/sizeof(double))? 1: 0];

   //EDIT_STATE
   checkVarOrder(nai);
   checkVarOrder(nass);
   checkVarOrder(ki);
   checkVarOrder(kss);
   checkVarOrder(cai);
   checkVarOrder(cass);
   checkVarOrder(cansr);
   checkVarOrder(cajsr);
   checkVarOrder(m);
   checkVarOrder(hf);
   checkVarOrder(hs);
   checkVarOrder(j);
   checkVarOrder(hsp);
   checkVarOrder(jp);
   checkVarOrder(mL);
   checkVarOrder(hL);
   checkVarOrder(hLp);
   checkVarOrder(a);
   checkVarOrder(iF);
   checkVarOrder(iS);
   checkVarOrder(ap);
   checkVarOrder(iFp);
   checkVarOrder(iSp);
   checkVarOrder(d);
   checkVarOrder(ff);
   checkVarOrder(fs);
   checkVarOrder(fcaf);
   checkVarOrder(fcas);
   checkVarOrder(jca);
   checkVarOrder(nca);
   checkVarOrder(ffp);
   checkVarOrder(fcafp);
   checkVarOrder(xrf);
   checkVarOrder(xrs);
   checkVarOrder(xs1);
   checkVarOrder(xs2);
   checkVarOrder(xk1);
   checkVarOrder(Jrelnp);
   checkVarOrder(Jrelp);
   checkVarOrder(CaMKt);
}

   
ThisReaction::ThisReaction(const int numPoints)
: nCells_(numPoints)
{
   assertStateOrderAndVarNamesAgree();
   state_.resize(nCells_);
   perCellFlags_.resize(nCells_);
   perCellParameters_.resize(nCells_);
}

void ThisReaction::calc(double dt, const VectorDouble32& Vm,
                       const vector<double>& iStim , VectorDouble32& dVm)
{
   for (unsigned ii=0; ii<nCells_; ++ii)
   {

      //set Vm
      const double v = Vm[ii];
      const double istim = iStim[ii];

      //set all state variables
      //EDIT_STATE
      const double nai=state_[ii].nai;
      const double nass=state_[ii].nass;
      const double ki=state_[ii].ki;
      const double kss=state_[ii].kss;
      const double cai=state_[ii].cai;
      const double cass=state_[ii].cass;
      const double cansr=state_[ii].cansr;
      const double cajsr=state_[ii].cajsr;
      const double m=state_[ii].m;
      const double hf=state_[ii].hf;
      const double hs=state_[ii].hs;
      const double j=state_[ii].j;
      const double hsp=state_[ii].hsp;
      const double jp=state_[ii].jp;
      const double mL=state_[ii].mL;
      const double hL=state_[ii].hL;
      const double hLp=state_[ii].hLp;
      const double a=state_[ii].a;
      const double iF=state_[ii].iF;
      const double iS=state_[ii].iS;
      const double ap=state_[ii].ap;
      const double iFp=state_[ii].iFp;
      const double iSp=state_[ii].iSp;
      const double d=state_[ii].d;
      const double ff=state_[ii].ff;
      const double fs=state_[ii].fs;
      const double fcaf=state_[ii].fcaf;
      const double fcas=state_[ii].fcas;
      const double jca=state_[ii].jca;
      const double nca=state_[ii].nca;
      const double ffp=state_[ii].ffp;
      const double fcafp=state_[ii].fcafp;
      const double xrf=state_[ii].xrf;
      const double xrs=state_[ii].xrs;
      const double xs1=state_[ii].xs1;
      const double xs2=state_[ii].xs2;
      const double xk1=state_[ii].xk1;
      const double Jrelnp=state_[ii].Jrelnp;
      const double Jrelp=state_[ii].Jrelp;
      const double CaMKt=state_[ii].CaMKt;

      //set per-cell flags
      //EDIT_PERCELL_FLAGS
      
      //set per-cell parameters
      //EDIT_PERCELL_PARAMETERS

      // Constants
      //extracellular ionic concentrations
      const double nao=140.0;
      const double cao=1.8;
      const double ko=5.4;

      //physical constants
      const double R=8314.0;
      const double T=310.0;
      const double F=96485.0;

      //cell geometry
      const double L=0.01;
      const double rad=0.0011;
      const double vcell=1000*3.14*rad*rad*L;
      const double Ageo=2*3.14*rad*rad+2*3.14*rad*L;
      const double Acap=2*Ageo;
      const double vmyo=0.68*vcell;
      const double vnsr=0.0552*vcell;
      const double vjsr=0.0048*vcell;
      const double vss=0.02*vcell;

      //reversal potentials
      const double ENa=(R*T/F)*log(nao/nai);
      const double EK=(R*T/F)*log(ko/ki);
      const double PKNa=0.01833;
      const double EKs=(R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai));

      //convenient shorthand calculations
      const double vffrt=v*F*F/(R*T);
      const double vfrt=v*F/(R*T);

      //CaMK constants
      const double KmCaMK=0.15;
      const double aCaMK=0.05;
      const double bCaMK=0.00068;
      const double CaMKo=0.05;
      const double KmCaM=0.0015;

      //update CaMK
      const double CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
      const double CaMKa=CaMKb+CaMKt;
      const double diff_CaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt;

      //// Membrane Currents
      double mss;
      double tm;
      double hss;
      double thf;
      double jss;
      double tj;
      if (useINaFromTT06)
      {
         const double aa_mTT2=1.0/(1.0+exp((-60.0-v)/5.0));
         const double bb_mTT2=0.1/(1.0+exp((v+35.0)/5.0))+0.1/(1.0+exp((v-50.0)/200.0));
         const double tau_mTT2=aa_mTT2*bb_mTT2;
         const double mTT2_inf=1.0/((1.0+exp((-56.86-v)/9.03))*(1.0+exp((-56.86-v)/9.03)));
         double aa_hTT2;
         double bb_hTT2;
         if (v>=-40.0)
         {
            aa_hTT2=0.0;
            bb_hTT2=(0.77/(0.13*(1.0+exp(-(v+10.66)/11.1))));
         }
         else
         {
            aa_hTT2=(0.057*exp(-(v+80.0)/6.8));
            bb_hTT2=(2.7*exp(0.079*v)+(3.1e5)*exp(0.3485*v));
         }
         const double tau_hTT2=1.0/(aa_hTT2+bb_hTT2);
         const double hTT2_inf=1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
         double aa_jTT2;
         double bb_jTT2;
         if (v>=-40)
         {
            aa_jTT2=0.0;
            bb_jTT2 =(0.6*exp((0.057)*v)/(1.0+exp(-0.1*(v+32.0))));
         }
         else
         {
            aa_jTT2=(((-2.5428e4)*exp(0.2444*v)-(6.948e-6)*exp(-0.04391*v))*(v+37.78)/(1.0+exp(0.311*(v+79.23))));
            bb_jTT2=(0.02424*exp(-0.01052*v)/(1.+exp(-0.1378*(v+40.14))));
         }
         const double tau_jTT2=1.0/(aa_jTT2+bb_jTT2);
         const double jTT2_inf=hTT2_inf;

         mss = mTT2_inf;
         tm = tau_mTT2;

         hss = hTT2_inf;
         thf = tau_hTT2;

         jss = jTT2_inf;
         tj = tau_jTT2;
      }
      else
      {
         mss=1.0/(1.0+exp((-(v+39.57))/9.871));
         tm=1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
         hss=1.0/(1+exp((v+82.90)/6.086));
         thf=1.0/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
         jss=hss;
         tj=2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
      }
      const double diff_m=(mss-m)/tm;
      const double diff_hf=(hss-hf)/thf;
      const double diff_j=(jss-j)/tj;

      const double ths=1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
      const double diff_hs=(hss-hs)/ths;
      const double hspss=1.0/(1+exp((v+89.1)/6.086));
      const double thsp=3.0*ths;
      const double diff_hsp=(hspss-hsp)/thsp;
      const double tjp=1.46*tj;
      const double diff_jp=(jss-jp)/tjp;

      double GNa;
      double INa;
      if (useINaFromTT06) {
         GNa=14.838;
         INa=GNa*(v-ENa)*m*m*m*hf*j;
      }
      else
      {
         const double Ahf=0.99;
         const double Ahs=1.0-Ahf;
         const double h=Ahf*hf+Ahs*hs;
         const double hp=Ahf*hf+Ahs*hsp;
         const double fINap=(1.0/(1.0+KmCaMK/CaMKa));
         GNa=75;
         INa=GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);
      }

      const double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
      const double tmL=tm;
      const double diff_mL=(mLss-mL)/tmL;
      const double hLss=1.0/(1.0+exp((v+87.61)/7.488));
      const double thL=200.0;
      const double diff_hL=(hLss-hL)/thL;
      const double hLpss=1.0/(1.0+exp((v+93.81)/7.488));
      const double thLp=3.0*thL;
      const double diff_hLp=(hLpss-hLp)/thLp;
      double GNaL=0.0075;
      if (celltype==1)
	{
	  GNaL*=0.6;
	}
      const double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
      const double INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

      const double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
      const double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
      const double diff_a=(ass-a)/ta;
      const double iss=1.0/(1.0+exp((v+43.94)/5.711));
      double delta_epi;
      if (celltype==1)
	{
	  delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
	}
      else
	{
	  delta_epi=1.0;
	}
      double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
      double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
      tiF*=delta_epi;
      tiS*=delta_epi;
      const double AiF=1.0/(1.0+exp((v-213.6)/151.2));
      const double AiS=1.0-AiF;
      const double diff_iF=(iss-iF)/tiF;
      const double diff_iS=(iss-iS)/tiS;
      const double i=AiF*iF+AiS*iS;
      const double apss=1.0/(1.0+exp((-(v-24.34))/14.82));
      const double diff_ap=(apss-ap)/ta;
      const double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
      const double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
      const double tiFp=dti_develop*dti_recover*tiF;
      const double tiSp=dti_develop*dti_recover*tiS;
      const double diff_iFp=(iss-iFp)/tiFp;
      const double diff_iSp=(iss-iSp)/tiSp;
      const double ip=AiF*iFp+AiS*iSp;
      double Gto=0.02;
      if (celltype==1)
	{
	  Gto*=4.0;
	}
      else if (celltype==2)
	{
	  Gto*=4.0;
	}
      const double fItop=(1.0/(1.0+KmCaMK/CaMKa));
      const double Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

      const double dss=1.0/(1.0+exp((-(v+3.940))/4.230));
      const double td=0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
      const double diff_d=(dss-d)/td;
      const double fss=1.0/(1.0+exp((v+19.58)/3.696));
      const double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
      const double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
      const double Aff=0.6;
      const double Afs=1.0-Aff;
      const double diff_ff=(fss-ff)/tff;
      const double diff_fs=(fss-fs)/tfs;
      const double f=Aff*ff+Afs*fs;
      const double fcass=fss;
      const double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
      const double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
      const double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
      const double Afcas=1.0-Afcaf;
      const double diff_fcaf=(fcass-fcaf)/tfcaf;
      const double diff_fcas=(fcass-fcas)/tfcas;
      const double fca=Afcaf*fcaf+Afcas*fcas;
      const double tjca=75.0;
      const double diff_jca=(fcass-jca)/tjca;
      const double tffp=2.5*tff;
      const double diff_ffp=(fss-ffp)/tffp;
      const double fp=Aff*ffp+Afs*fs;
      const double tfcafp=2.5*tfcaf;
      const double diff_fcafp=(fcass-fcafp)/tfcafp;
      const double fcap=Afcaf*fcafp+Afcas*fcas;
      const double Kmn=0.002;
      const double k2n=1000.0;
      const double km2n=jca*1.0;
      const double anca=1.0/(k2n/km2n+pow(1.0+Kmn/cass,4.0));
      const double diff_nca=anca*k2n-nca*km2n;
      const double PhiCaL=4.0*vffrt*(cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
      const double PhiCaNa=1.0*vffrt*(0.75*nass*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
      const double PhiCaK=1.0*vffrt*(0.75*kss*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
      const double zca=2.0;
      double PCa=0.0001;
      if (celltype==1)
	{
	  PCa*=1.2;
	}
      else if (celltype==2)
	{
	  PCa*=2.5;
	}
      const double PCap=1.1*PCa;
      const double PCaNa=0.00125*PCa;
      const double PCaK=3.574e-4*PCa;
      const double PCaNap=0.00125*PCap;
      const double PCaKp=3.574e-4*PCap;
      const double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
      const double ICaL=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
      const double ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
      const double ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);

      const double xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
      const double txrf=12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
      const double txrs=1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
      const double Axrf=1.0/(1.0+exp((v+54.81)/38.21));
      const double Axrs=1.0-Axrf;
      const double diff_xrf=(xrss-xrf)/txrf;
      const double diff_xrs=(xrss-xrs)/txrs;
      const double xr=Axrf*xrf+Axrs*xrs;
      const double rkr=1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
      double GKr=0.046;
      if (celltype==1)
	{
	  GKr*=1.3;
	}
      else if (celltype==2)
	{
	  GKr*=0.8;
	}
      const double IKr=GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);

      const double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
      const double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
      const double diff_xs1=(xs1ss-xs1)/txs1;
      const double xs2ss=xs1ss;
      const double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
      const double diff_xs2=(xs2ss-xs2)/txs2;
      const double KsCa=1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
      double GKs=0.0034;
      if (celltype==1)
	{
	  GKs*=1.4;
	}
      const double IKs=GKs*KsCa*xs1*xs2*(v-EKs);

      const double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
      const double txk1=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
      const double diff_xk1=(xk1ss-xk1)/txk1;
      const double rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
      double GK1=0.1908;
      if (celltype==1)
	{
	  GK1*=1.2;
	}
      else if (celltype==2)
	{
	  GK1*=1.3;
	}
      const double IK1=GK1*sqrt(ko)*rk1*xk1*(v-EK);

      double kna1=15.0;
      double kna2=5.0;
      double kna3=88.12;
      double kasymm=12.5;
      double wna=6.0e4;
      double wca=6.0e4;
      double wnaca=5.0e3;
      double kcaon=1.5e6;
      double kcaoff=5.0e3;
      double qna=0.5224;
      double qca=0.1670;
      double hca=exp((qca*v*F)/(R*T));
      double hna=exp((qna*v*F)/(R*T));
      double h1=1+nai/kna3*(1+hna);
      double h2=(nai*hna)/(kna3*h1);
      double h3=1.0/h1;
      double h4=1.0+nai/kna1*(1+nai/kna2);
      double h5=nai*nai/(h4*kna1*kna2);
      double h6=1.0/h4;
      double h7=1.0+nao/kna3*(1.0+1.0/hna);
      double h8=nao/(kna3*hna*h7);
      double h9=1.0/h7;
      double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
      double h11=nao*nao/(h10*kna1*kna2);
      double h12=1.0/h10;
      double k1=h12*cao*kcaon;
      double k2=kcaoff;
      double k3p=h9*wca;
      double k3pp=h8*wnaca;
      double k3=k3p+k3pp;
      double k4p=h3*wca/hca;
      double k4pp=h2*wnaca;
      double k4=k4p+k4pp;
      double k5=kcaoff;
      double k6=h6*cai*kcaon;
      double k7=h5*h2*wna;
      double k8=h8*h11*wna;
      double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
      double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
      double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
      double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
      double E1=x1/(x1+x2+x3+x4);
      double E2=x2/(x1+x2+x3+x4);
      double E3=x3/(x1+x2+x3+x4);
      double E4=x4/(x1+x2+x3+x4);
      double KmCaAct=150.0e-6;
      double allo=1.0/(1.0+(KmCaAct/cai)*(KmCaAct/cai));
      double zna=1.0;
      double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
      double JncxCa=E2*k2-E1*k1;
      double Gncx=0.0008;
      if (celltype==1)
	{
	  Gncx*=1.1;
	}
      else if (celltype==2)
	{
	  Gncx*=1.4;
	}
      const double INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);

      h1=1+nass/kna3*(1+hna);
      h2=(nass*hna)/(kna3*h1);
      h3=1.0/h1;
      h4=1.0+nass/kna1*(1+nass/kna2);
      h5=nass*nass/(h4*kna1*kna2);
      h6=1.0/h4;
      h7=1.0+nao/kna3*(1.0+1.0/hna);
      h8=nao/(kna3*hna*h7);
      h9=1.0/h7;
      h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
      h11=nao*nao/(h10*kna1*kna2);
      h12=1.0/h10;
      k1=h12*cao*kcaon;
      k2=kcaoff;
      k3p=h9*wca;
      k3pp=h8*wnaca;
      k3=k3p+k3pp;
      k4p=h3*wca/hca;
      k4pp=h2*wnaca;
      k4=k4p+k4pp;
      k5=kcaoff;
      k6=h6*cass*kcaon;
      k7=h5*h2*wna;
      k8=h8*h11*wna;
      x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
      x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
      x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
      x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
      E1=x1/(x1+x2+x3+x4);
      E2=x2/(x1+x2+x3+x4);
      E3=x3/(x1+x2+x3+x4);
      E4=x4/(x1+x2+x3+x4);
      KmCaAct=150.0e-6;
      allo=1.0/(1.0+(KmCaAct/cass)*(KmCaAct/cass));
      JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
      JncxCa=E2*k2-E1*k1;
      const double INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);

      const double k1p=949.5;
      const double k1m=182.4;
      const double k2p=687.2;
      const double k2m=39.4;
      k3p=1899.0;
      const double k3m=79300.0;
      k4p=639.0;
      const double k4m=40.0;
      const double Knai0=9.073;
      const double Knao0=27.78;
      const double delta=-0.1550;
      const double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
      const double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
      const double Kki=0.5;
      const double Kko=0.3582;
      const double MgADP=0.05;
      const double MgATP=9.8;
      const double Kmgatp=1.698e-7;
      const double H=1.0e-7;
      const double eP=4.2;
      const double Khp=1.698e-7;
      const double Knap=224.0;
      const double Kxkur=292.0;
      const double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
      const double a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
      const double b1=k1m*MgADP;
      const double a2=k2p;
      const double b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
      const double a3=(k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
      const double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
      const double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
      const double b4=(k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
      x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
      x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
      x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
      x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
      E1=x1/(x1+x2+x3+x4);
      E2=x2/(x1+x2+x3+x4);
      E3=x3/(x1+x2+x3+x4);
      E4=x4/(x1+x2+x3+x4);
      const double zk=1.0;
      const double JnakNa=3.0*(E1*a3-E2*b3);
      const double JnakK=2.0*(E4*b1-E3*a1);
      double Pnak=30;
      if (celltype==1)
	{
	  Pnak*=0.9;
	}
      else if (celltype==2)
	{
	  Pnak*=0.7;
	}
      const double INaK=Pnak*(zna*JnakNa+zk*JnakK);

      const double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
      double GKb=0.003;
      if (celltype==1)
	{
	  GKb*=0.6;
	}
      const double IKb=GKb*xkb*(v-EK);

      const double PNab=3.75e-10;
      const double INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

      const double PCab=2.5e-8;
      const double ICab=PCab*4.0*vffrt*(cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);

      const double GpCa=0.0005;
      const double IpCa=GpCa*cai/(0.0005+cai);

      //Fluxes and Buffers

      const double JdiffNa=(nass-nai)/2.0;
      const double JdiffK=(kss-ki)/2.0;
      const double Jdiff=(cass-cai)/0.2;

      const double bt=4.75;
      const double a_rel=0.5*bt;
      double Jrel_inf=a_rel*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
      if (celltype==2)
	{
	  Jrel_inf*=1.7;
	}
      const double tau_rel=bt/(1.0+0.0123/cajsr);
      double tau_rel_use;
      if (tau_rel<0.001)
	{
	  tau_rel_use=0.001; 
	}
      else
	{
	  tau_rel_use=tau_rel;
	}
      const double diff_Jrelnp=(Jrel_inf-Jrelnp)/tau_rel_use;

      const double btp=1.25*bt;
      const double a_relp=0.5*btp;
      double Jrelp_inf=a_relp*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
      if (celltype==2)
	{
	  Jrelp_inf=Jrelp_inf*1.7;
	}
      const double tau_relp=btp/(1.0+0.0123/cajsr);
      double tau_relp_use;
      if (tau_relp<0.001)
	{
	  tau_relp_use=0.001; 
	}
      else
	{
	  tau_relp_use=tau_relp;
	}
      const double diff_Jrelp=(Jrelp_inf-Jrelp)/tau_relp_use;
      const double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
      /* May 2011:
         RCB:
         Ok, I came up with a hack for Jrel that doesn't involve clamping if
         cajsr goes negative.  This gives the same results as an adaptive
         integration scheme.

         jrel_stiff_const should be 0.1*dt, roughly.


         -Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;
         +Jrel_canidate=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;
         +Jrel=Jrel_canidate;
         +
         +jrel_stiff_const = 0.005; .param();
         +if (Jrel_canidate*jrel_stiff_const > cajsr) {
         +   Jrel = cajsr/jrel_stiff_const;
         +}

         Tom:
         
         This is terrific news.  My evening beers are dedicated to your
         efforts!  She runs fast and smooth now?  I want very much to make it
         so.  Please let me know what's next.

         Rob:
         I think she runs fast and smooth, at least in a single cell.  Simulation
         results are going to take a few days to run and process.

         The Jrel hack should maintain conservation of calcium.  I spent a couple
         of hours writing equations to figure out what was happening.  Basically,
         near zero the differential equation for cajsr turns into something like
         this (I forget the exact details):
         
         y'' = -(1+c/y)*y'-b*y - d
         
         All my fix does is damp that (1+c/y) term when y is near zero.  cajsr
         might take an extra timestep or two to reach zero, but it won't go negative.
      */
      double Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;
      if (Jrel*JrelStiffConst > cajsr)
      {
         Jrel = cajsr/JrelStiffConst;
      }
      
      double Jupnp=0.004375*cai/(cai+0.00092);
      double Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);
      if (celltype==1)
	{
	  Jupnp*=1.3;
	  Jupp*=1.3;
	}
      const double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
      const double Jleak=0.0039375*cansr/15.0;
      const double Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

      const double Jtr=(cansr-cajsr)/100.0;

      double cmdnmax=0.05;
      if (celltype==1)
	{
	  cmdnmax*=1.3;
	}
      const double kmcmdn=0.00238;
      const double trpnmax=0.07;
      const double kmtrpn=0.0005;
      const double BSRmax=0.047;
      const double KmBSR=0.00087;
      const double BSLmax=1.124;
      const double KmBSL=0.0087;
      const double csqnmax=10.0;
      const double kmcsqn=0.8;

      const double diff_nai=-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;
      const double diff_nass=-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;

      const double diff_ki=-(Ito+IKr+IKs+IK1+IKb-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
      const double diff_kss=-(ICaK)*Acap/(F*vss)-JdiffK;

      const double Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
      const double diff_cai=Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);

      const double Bcass=1.0/(1.0+BSRmax*KmBSR/pow(KmBSR+cass,2.0)+BSLmax*KmBSL/pow(KmBSL+cass,2.0));
      const double diff_cass=Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);

      const double diff_cansr=Jup-Jtr*vjsr/vnsr;

      const double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow(kmcsqn+cajsr,2.0));
      const double diff_cajsr=Bcajsr*(Jtr-Jrel);


      //// Membrane Potential
      ////
      const double diff_v=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab);
      dVm[ii] = diff_v;

      if (1) 
      {
         bool foundError=false;
#define CHECK_BLOWUP(x) do { if (!isfinite(x)) { fprintf(stderr, "Error in node %d, variable " #x " = %g\n", ii, (x)); foundError=true; } } while(0)
         CHECK_BLOWUP(diff_v);
            
         //EDIT_STATE
         CHECK_BLOWUP(diff_nai);
         CHECK_BLOWUP(diff_nass);
         CHECK_BLOWUP(diff_ki);
         CHECK_BLOWUP(diff_kss);
         CHECK_BLOWUP(diff_cai);
         CHECK_BLOWUP(diff_cass);
         CHECK_BLOWUP(diff_cansr);
         CHECK_BLOWUP(diff_cajsr);
         CHECK_BLOWUP(diff_m);
         CHECK_BLOWUP(diff_hf);
         CHECK_BLOWUP(diff_hs);
         CHECK_BLOWUP(diff_j);
         CHECK_BLOWUP(diff_hsp);
         CHECK_BLOWUP(diff_jp);
         CHECK_BLOWUP(diff_mL);
         CHECK_BLOWUP(diff_hL);
         CHECK_BLOWUP(diff_hLp);
         CHECK_BLOWUP(diff_a);
         CHECK_BLOWUP(diff_iF);
         CHECK_BLOWUP(diff_iS);
         CHECK_BLOWUP(diff_ap);
         CHECK_BLOWUP(diff_iFp);
         CHECK_BLOWUP(diff_iSp);
         CHECK_BLOWUP(diff_d);
         CHECK_BLOWUP(diff_ff);
         CHECK_BLOWUP(diff_fs);
         CHECK_BLOWUP(diff_fcaf);
         CHECK_BLOWUP(diff_fcas);
         CHECK_BLOWUP(diff_jca);
         CHECK_BLOWUP(diff_nca);
         CHECK_BLOWUP(diff_ffp);
         CHECK_BLOWUP(diff_fcafp);
         CHECK_BLOWUP(diff_xrf);
         CHECK_BLOWUP(diff_xrs);
         CHECK_BLOWUP(diff_xs1);
         CHECK_BLOWUP(diff_xs2);
         CHECK_BLOWUP(diff_xk1);
         CHECK_BLOWUP(diff_Jrelnp);
         CHECK_BLOWUP(diff_Jrelp);
         CHECK_BLOWUP(diff_CaMKt);

#undef CHECK_BLOWUP
         
         if (foundError) 
         {
#define PRINT_STATE(x) do { fprintf(stderr, "node %d: " #x " = %g\n", ii, (x)); } while(0)
            //EDIT_STATE
            PRINT_STATE(nai);
            PRINT_STATE(nass);
            PRINT_STATE(ki);
            PRINT_STATE(kss);
            PRINT_STATE(cai);
            PRINT_STATE(cass);
            PRINT_STATE(cansr);
            PRINT_STATE(cajsr);
            PRINT_STATE(m);
            PRINT_STATE(hf);
            PRINT_STATE(hs);
            PRINT_STATE(j);
            PRINT_STATE(hsp);
            PRINT_STATE(jp);
            PRINT_STATE(mL);
            PRINT_STATE(hL);
            PRINT_STATE(hLp);
            PRINT_STATE(a);
            PRINT_STATE(iF);
            PRINT_STATE(iS);
            PRINT_STATE(ap);
            PRINT_STATE(iFp);
            PRINT_STATE(iSp);
            PRINT_STATE(d);
            PRINT_STATE(ff);
            PRINT_STATE(fs);
            PRINT_STATE(fcaf);
            PRINT_STATE(fcas);
            PRINT_STATE(jca);
            PRINT_STATE(nca);
            PRINT_STATE(ffp);
            PRINT_STATE(fcafp);
            PRINT_STATE(xrf);
            PRINT_STATE(xrs);
            PRINT_STATE(xs1);
            PRINT_STATE(xs2);
            PRINT_STATE(xk1);
            PRINT_STATE(Jrelnp);
            PRINT_STATE(Jrelp);
            PRINT_STATE(CaMKt);
            
#undef PRINT_STATE
            
            exit(255);
         }
      }
      
      
      //EDIT_STATE
      state_[ii].nai += diff_nai*dt;
      state_[ii].nass += diff_nass*dt;
      state_[ii].ki += diff_ki*dt;
      state_[ii].kss += diff_kss*dt;
      state_[ii].cai += diff_cai*dt;
      state_[ii].cass += diff_cass*dt;
      state_[ii].cansr += diff_cansr*dt;
      state_[ii].cajsr += diff_cajsr*dt;
      state_[ii].m = mss-(mss-m)*exp(-dt/tm);
      state_[ii].hf = hss-(hss-hf)*exp(-dt/thf);
      state_[ii].hs = hss-(hss-hs)*exp(-dt/ths);
      state_[ii].j = jss-(jss-j)*exp(-dt/tj);
      state_[ii].hsp = hspss-(hspss-hsp)*exp(-dt/thsp);
      state_[ii].jp = jss-(jss-jp)*exp(-dt/tjp);
      state_[ii].mL = mLss-(mLss-mL)*exp(-dt/tmL);
      state_[ii].hL = hLss-(hLss-hL)*exp(-dt/thL);
      state_[ii].hLp = hLpss-(hLpss-hLp)*exp(-dt/thLp);
      state_[ii].a = ass-(ass-a)*exp(-dt/ta);
      state_[ii].iF = iss-(iss-iF)*exp(-dt/tiF);
      state_[ii].iS = iss-(iss-iS)*exp(-dt/tiS);
      state_[ii].ap = apss-(apss-ap)*exp(-dt/ta);
      state_[ii].iFp = iss-(iss-iFp)*exp(-dt/tiFp);
      state_[ii].iSp = iss-(iss-iSp)*exp(-dt/tiSp);
      state_[ii].d = dss-(dss-d)*exp(-dt/td);
      state_[ii].ff = fss-(fss-ff)*exp(-dt/tff);
      state_[ii].fs = fss-(fss-fs)*exp(-dt/tfs);
      state_[ii].fcaf = fcass-(fcass-fcaf)*exp(-dt/tfcaf);
      state_[ii].fcas = fcass-(fcass-fcas)*exp(-dt/tfcas);
      state_[ii].jca = fcass-(fcass-jca)*exp(-dt/tjca);
      state_[ii].ffp = fcass-(fcass-ffp)*exp(-dt/tffp);
      state_[ii].fcafp = fcass-(fcass-fcafp)*exp(-dt/tfcafp);
      state_[ii].nca += diff_nca*dt;
      state_[ii].xrf = xrss-(xrss-xrf)*exp(-dt/txrf);
      state_[ii].xrs = xrss-(xrss-xrs)*exp(-dt/txrs);
      state_[ii].xs1 = xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
      state_[ii].xs2 = xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
      state_[ii].xk1 = xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
      state_[ii].Jrelnp = Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
      state_[ii].Jrelp = Jrelp_inf-(Jrelp_inf-Jrelp)*exp(-dt/tau_relp);
      state_[ii].CaMKt += diff_CaMKt*dt;
      
   }
}

void ThisReaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= nCells_);
   Vm.assign(Vm.size(), -85.0);
   State initState;
   //EDIT_STATE

   initState.nai=7;
   initState.nass=7;
   initState.ki=145;
   initState.kss=145;
   initState.cai=1.0e-4;
   initState.cass=1.0e-4;
   initState.cansr=1.2;
   initState.cajsr=1.2;
   initState.m=0;
   initState.hf=1;
   initState.hs=1;
   initState.j=1;
   initState.hsp=1;
   initState.jp=1;
   initState.mL=0;
   initState.hL=1;
   initState.hLp=1;
   initState.a=0;
   initState.iF=1;
   initState.iS=1;
   initState.ap=0;
   initState.iFp=1;
   initState.iSp=1;
   initState.d=0;
   initState.ff=1;
   initState.fs=1;
   initState.fcaf=1;
   initState.fcas=1;
   initState.jca=1;
   initState.nca=0;
   initState.ffp=1;
   initState.fcafp=1;
   initState.xrf=0;
   initState.xrs=0;
   initState.xs1=0;
   initState.xs2=0;
   initState.xk1=1;
   initState.Jrelnp=0;
   initState.Jrelp=0;
   initState.CaMKt=0;

   state_.resize(nCells_);
   state_.assign(state_.size(), initState);

}

const string ThisReaction::getUnit(const std::string& varName) const
{
   //deliberatly broken for now, if this code still is being used past 2016-11-01 something has gone wrong.
   return "1";
}

#define HANDLE_OFFSET 1000
int ThisReaction::getVarHandle(const std::string& varName) const
{
   int retVal = getVarOffset(varName);
   if (retVal >= 0) {
      retVal += HANDLE_OFFSET;
   }
   return retVal;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
   reinterpret_cast<double*>(&state_[iCell])[varHandle-HANDLE_OFFSET] = value;
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
   return reinterpret_cast<const double*>(&state_[iCell])[varHandle-HANDLE_OFFSET];
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.resize(NUMVARS);
   fieldUnits.resize(NUMVARS);
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      fieldNames[ivar] = varNames[ivar];
      fieldUnits[ivar] = getUnit(fieldNames[ivar]);
   }
}

}
