/* $Id$ */
#include "units.h"
#include <stdlib.h>
#include <string.h>
#include "codata.h"
#include "error.h"
#include "ddcMath.h"
#if 1
#include "ddcMalloc.h"
#define Free ddcFree 
#define Malloc ddcMalloc
#else
#define Free free
#define Malloc malloc
#endif 

enum UNITS_ENUM { LENGTH,MASS,TIME,CURRENT,TEMPERATURE,AMOUNT,LUMINOUS_INTENSITY,ENERGY,PRESSURE,VELOCITY,EXTERNAL}; 
typedef struct symbol_struc { char *name; double length,mass,time,current,temperature,amount,luminous_intensity,mks_value;} SYMBOL;
double kB; 
double hbar;
double mu0;
double ke; 
double hc;
double mec2;
double re;
double qelectron;
double M_e;
double M_p;
double M_n;
double M_D;
double M_T;
double M_He4;
double M_He3;
static SYMBOL symbol_table[] = 
{
	{"length_internal",                    	1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"mass_internal",                      	0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"time_internal",                      	0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"current_internal",                   	0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
	{"temperature_internal",               	0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
	{"amount_internal",                    	0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0},
	{"luminous_intensity_internal",        	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0},
	{"energy_internal",						2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"pressure_internal",			   	-1.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"velocity_internal",			   	1.0, 1.0,1.0, 0.0, 0.0, 0.0, 0.0, 1.0},

	{"l",        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"m",        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"t",        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"i",        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
	{"T",        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
	{"n",        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0},
	{"I",        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0},
	{"energy",	 2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"pressure",-1.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0},
	{"velocity", 1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0},

	{"Ang",      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-10},
	{"Angstrom", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-10},
	{"Bohr",     1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a0_MKS},
	{"a0",       1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a0_MKS},
	{"eV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0},
	{"Hartree",  2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, Eh},
	{"keV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0e+3},
	{"MeV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0e+6},
	{"meV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0e-3},
	{"ueV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0e-6},
	{"kB",		 2.0, 1.0,-2.0, 0.0,-1.0, 0.0, 0.0, kB_MKS},
	{"amu",		 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, u_MKS},
	{"fs",	 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1e-15},
	{"GPa",-1.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1e9},
	{"K",        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
	{"e",        0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, e_MKS},
	{"M_e",      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, me_MKS},
	{"M_p",      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mp_MKS},
	{"M_n",      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mn_MKS},
	{"M_D",      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mD_MKS},
	{"M_T",      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mT_MKS},
	{"M_He3",  0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mh_MKS},
	{"M_He4",  0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mAlpha_MKS},
	{"END",      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
}; 
double units_parse(char *unit_string, double *expSum)
{
	//l m t I T A L
	char *ptr, *next, *string,op;
	double  Sexp, *exp,mks_value; 
	SYMBOL *symbol; 
	int i,j,l; 
	double post_op = 1.0; 
	for (i=0;i<7;i++) expSum[i]=0.0; 
	l = strlen(unit_string)+2; 
	string = (char*) Malloc(l); 
	j=0; 
	for (i=0;i<l-2;i++) if ( unit_string[i] != ' ' ) string[j++] = unit_string[i]; 
	string[j] = '*'; 
	string[j+1] = (char)0; 
	op = '*'; 
	mks_value=1.0; 
	char *tok = string; 
	while (*tok != (char)0)
	{
		if (op  == '*' ) post_op= 1.0; 
		if (op  == '/' ) post_op=-1.0; 
		next   = strpbrk(tok, "/*");
		op = *next; 
		*next = (char)0;
		ptr  = strpbrk(tok, "^");
		if (ptr == NULL) Sexp = 1.0; 
		else {Sexp= strtod(ptr+1,NULL); *ptr = (char)0; }
		Sexp *= post_op;
		symbol = symbol_table;
		char *start = tok; 
		double value = strtod(start,&tok);
		if (start != tok) mks_value *= value;
		if (*tok != (char)0)
		{
			while(strcmp(symbol->name,"END")!=0)
			{
				if (strcmp(symbol->name,tok)==0)
				{
						exp = &(symbol->length);
						mks_value *= pow(symbol->mks_value,Sexp);
						for (i=0;i<7;i++) expSum[i] += Sexp * exp[i];
						break;
				}
				symbol++;
			}
		}
		tok = next + 1; 
	}
	Free(string); 
	return mks_value; 
}
void units_internal(double length, double mass, double time, double current, double  temperature, double amount, double luminous_intensity)
{
	double energy; 
	symbol_table[LENGTH].mks_value = length; 
	symbol_table[MASS].mks_value = mass; 
	symbol_table[TIME].mks_value = time; 
	symbol_table[CURRENT].mks_value = current; 
	symbol_table[TEMPERATURE].mks_value = temperature; 
	symbol_table[AMOUNT].mks_value = amount; 
	symbol_table[LUMINOUS_INTENSITY].mks_value = luminous_intensity; 
	symbol_table[ENERGY].mks_value= energy = mass*length*length/(time*time);
	double charge = current*time; 
	kB = kB_MKS * temperature/energy; 
	hbar = hbar_MKS/(energy*time); 
	mu0 = mu0_MKS/((energy/length) * current*current); 
	ke = ke_MKS*charge*charge/(energy*length); 
	hc =  h_MKS*c_MKS/(energy*length);
	mec2 = mec2_MKS/energy; 
	re = re_MKS/length; 
	M_e = me_MKS/mass; 
	M_p = mp_MKS/mass; 
	M_n = mn_MKS/mass; 
	M_D = mD_MKS/mass; 
	M_T = mT_MKS/mass; 
	M_He3 = mh_MKS/mass; 
	M_He4 = mAlpha_MKS/mass; 
	qelectron = e_MKS/(current*time); 
	
}
void units_external(double length, double mass, double time, double current, double  temperature, double amount, double luminous_intensity)
{
	SYMBOL *symbol=symbol_table+EXTERNAL; 
	symbol[LENGTH].mks_value = length; 
	symbol[MASS].mks_value = mass; 
	symbol[TIME].mks_value = time; 
	symbol[CURRENT].mks_value = current; 
	symbol[TEMPERATURE].mks_value = temperature; 
	symbol[AMOUNT].mks_value = amount; 
	symbol[LUMINOUS_INTENSITY].mks_value = luminous_intensity; 
	symbol[ENERGY].mks_value= mass*length*length/(time*time);
	symbol[PRESSURE].mks_value= mass/(time*time*length);
	symbol[VELOCITY].mks_value= length/time;
}
int units_check(char *u1, char*u2) 
{
	double mks,exp1[7],exp2[7]; 
	if ((u1 == NULL) || (u2 == NULL)) return 0; 
	mks = units_parse(u1,exp1); 
	mks = units_parse(u2,exp2); 
	double sum =0.0; 
	for (int i=0;i<7;i++) sum+= fabs(exp1[i]-exp2[i]); 
	if ( sum > 1e-8)  return 1 ;
	return 0; 
}
double units_convert(double value,char *from , char *to )
{
	
	double to_mks=0.0,from_mks=0.0;;
	double fexp[7],texp[7];  
	int i; 
	if (from == NULL && to == NULL ) error_action("Can't determine conversion both, from and to are NULL", ERROR_IN("unit_convert", ABORT));
	if (from == NULL) 
	{
		to_mks=units_parse(to,texp); 
		from_mks = 1.0; 
		for (i=0;i<7;i++) 
		{
			from_mks *= pow(symbol_table[i].mks_value,texp[i]); 
		}
		return value*from_mks/to_mks;
	}
	if (to == NULL) 
	{
		from_mks = units_parse(from,fexp); 
		to_mks=1.0; 
		for (i=0;i<7;i++) to_mks *= pow(symbol_table[i].mks_value,fexp[i]); 
		return value*from_mks/to_mks;
	}
	to_mks=units_parse(to,texp); 
	from_mks = units_parse(from,fexp); 
	double sum =0.0; 
	for (i=0;i<7;i++) sum+= fabs(fexp[i]-texp[i]); 
	if ( sum > 1e-8) error_action("unit mismatch", ERROR_IN("unit_convert", ABORT));
	return  value*from_mks/to_mks; 
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
