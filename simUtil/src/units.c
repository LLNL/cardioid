/* $Id$ */
#include "units.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "codata.h"
#include "units.h"
#include "error.h"
#include "assert.h"

enum UNITS_ENUM
{
   LENGTH, MASS, TIME, CURRENT, TEMPERATURE, AMOUNT, LUMINOUS_INTENSITY,
   ENERGY, PRESSURE, VELOCITY, CAPACITANCE, RESISTIVITY, VOLTAGE, EXTERNAL
}; 

typedef struct symbol_struc { char *name; double length,mass,time,current,temperature,amount,luminous_intensity,mks_value;} SYMBOL;

double kB; 
double hbar;
double mu0;
double ke; 
double hc;
double mec2;
double bohr;
double re;
double qelectron;
double M_e;
double M_p;
double M_n;
double M_D;
double M_T;
double M_He4;
double M_He3;
double amu;
double Eh;

/** A note on capitalization:
 *  -------------------------
 *  Officially, the names of SI units are not capitalized (unless
 *  grammar requires such as at the start of a sentence.  There are a
 *  few unit names in the table below that are capitalized.  They were
 *  created before we began to observe this convention and we leave them
 *  in place to avoid breaking things.
 *
 *  As for abbreviations, SI units do not have abbreviations, they have
 *  symbols.  For example, V is not the abbreviation for volt, it is the
 *  symbol for volt.  Symbols are like mathematical symbols that are
 *  case sensitive, (i.e., s = seconds, S = siemens).  Hence, there is
 *  no systematic capitalization rule to apply to symbols.
 *
 *
 *  Adding internal or external units:
 *  ----------------------------------
 *
 *  If you want to add internal or external units, you must do three
 *  things:
 *  1.  You must add them in pairs.  Do not add an internal unit without
 *      adding the corresponding external unit and vice-versa.
 *  2.  You must add an entry to the UNITS_ENUM that corresponds to
 *      the unit you have added.  Keep everything in the right order and
 *      EXTERNAL must be last.
 *  3.  You must add code to the units_internal and units_external
 *      functions to properly convert your new units from the SI
 *      definitions that you enetered in the table to whatevery (crazy)
 *      system the caller is requesting.
 *  If you fail to do any of these things the unit conversion system may
 *  fail in all sorts of strange and interesting ways.  You have been
 *  warned. 
 */
static SYMBOL symbol_table[] = 
{
   // Read comments above before adding internal units.
   //                               l    m    t    i    T    n    I    const
   {"length_internal",              1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"mass_internal",                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"time_internal",                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"current_internal",             0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
   {"temperature_internal",         0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
   {"amount_internal",              0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0},
   {"luminous_intensity_internal",  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0},

   {"energy_internal",              2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"pressure_internal",           -1.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"velocity_internal",            1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"capacitance_internal",        -2.0,-1.0, 4.0, 2.0, 0.0, 0.0, 0.0, 1.0}, 
   {"resistivity_internal",        -2.0,-1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0},
   {"voltage_internal",             2.0, 1.0,-3.0,-1.0, 0.0, 0.0, 0.0, 1.0},

   // read comments above before adding external units.
   // external units
   //                   l    m    t    i    T    n    I    const
   {"l",                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"m",                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"t",                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"i",                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
   {"T",                0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
   {"n",                0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0},
   {"I",                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0},
   {"energy",           2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"pressure",        -1.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"velocity",         1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"capacitance",     -2.0,-1.0, 4.0, 2.0, 0.0, 0.0, 0.0, 1.0}, 
   {"resistivity",     -2.0,-1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0},
   {"voltage",          2.0, 1.0,-3.0,-1.0, 0.0, 0.0, 0.0, 1.0},

   // named units
   //            l    m    t    i    T    n    I    const
   {"Ang",       1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-10},
   {"Angstrom",  1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-10},
   {"Bohr",      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a0_MKS},
   {"a0",        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a0_MKS},
   {"meter",     1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"mm",        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-3},
   {"um",        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-6},
   {"nm",        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-9},
   {"gram",      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-3},
   {"g",         0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-3},
   {"kg",        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
   {"mg",        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-6},
   {"ug",        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-9},
   {"eV",        2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0},
   {"keV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0e+3},
   {"MeV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0e+6},
   {"meV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0e-3},
   {"ueV",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, e_MKS*1.0e-6},
   {"Hartree",   2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, Eh_MKS},
   {"Ry",        2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, Eh_MKS*0.5},
   {"J",         2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1},
   {"kJ",        2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0e+3},
   {"cal",       2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, cal_MKS},
   {"kcal",      2.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, cal_MKS*1.0e+3},
   {"amu",       0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, u_MKS},
   {"second",    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1},
   {"s",         0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1},
   {"ms",        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1e-3},
   {"us",        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1e-6},
   {"ns",        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1e-9},
   {"ps",        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1e-12},
   {"fs",        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1e-15},
   {"amp",       0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
   {"A",         0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
   {"mA",        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0e-3},
   {"uA",        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0e-6},
   {"nA",        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0e-9},
   {"pA",        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0e-12},
   {"cc",        3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-6},
   {"mol",       0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, NA},  // mol avogadro
   {"molar",    -3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0e3}, // 1 molar = 1mol/liter
   {"M",        -3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0e3}, // 1 M = 1 molar
   {"mM",       -3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0},   // 1 mM = 1 mol/meter^3
   {"uM",       -3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0e-3},   
   {"GPa",      -1.0, 1.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1e9},
   {"K",         0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
   {"volt",      2.0, 1.0,-3.0,-1.0, 0.0, 0.0, 0.0, 1.0},
   {"V",         2.0, 1.0,-3.0,-1.0, 0.0, 0.0, 0.0, 1.0},
   {"mV",        2.0, 1.0,-3.0,-1.0, 0.0, 0.0, 0.0, 1.0e-3},
   {"uV",        2.0, 1.0,-3.0,-1.0, 0.0, 0.0, 0.0, 1.0e-6},
   {"siemens",  -2.0,-1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0},
   {"S",        -2.0,-1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0},
   {"mS",       -2.0,-1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0e-3},
   {"uS",       -2.0,-1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0e-6},
   {"nS",       -2.0,-1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0e-9},
   {"pS",       -2.0,-1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0e-12},
   {"farad",    -2.0,-1.0, 4.0, 2.0, 0.0, 0.0, 0.0, 1.0}, 
   {"F",        -2.0,-1.0, 4.0, 2.0, 0.0, 0.0, 0.0, 1.0}, 
   {"mF",       -2.0,-1.0, 4.0, 2.0, 0.0, 0.0, 0.0, 1.0e-3}, 
   {"uF",       -2.0,-1.0, 4.0, 2.0, 0.0, 0.0, 0.0, 1.0e-6}, 
   {"nF",       -2.0,-1.0, 4.0, 2.0, 0.0, 0.0, 0.0, 1.0e-9}, 
   {"pF",       -2.0,-1.0, 4.0, 2.0, 0.0, 0.0, 0.0, 1.0e-12}, 
   {"coulomb",   0.0, 0.0,-1.0, 1.0, 0.0, 0.0, 0.0, 1.0}, 
   {"C",         0.0, 0.0,-1.0, 1.0, 0.0, 0.0, 0.0, 1.0}, 
   {"mC",        0.0, 0.0,-1.0, 1.0, 0.0, 0.0, 0.0, 1.0e-3}, 
   {"uC",        0.0, 0.0,-1.0, 1.0, 0.0, 0.0, 0.0, 1.0e-6}, 

   // constants
   //            l    m    t    i    T    n    I    const
   {"kB",        2.0, 1.0,-2.0, 0.0,-1.0, 0.0, 0.0, kB_MKS},
   {"e",         0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, e_MKS},
   {"M_e",       0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, me_MKS},
   {"M_p",       0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mp_MKS},
   {"M_n",       0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mn_MKS},
   {"M_D",       0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mD_MKS},
   {"M_T",       0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mT_MKS},
   {"M_He3",     0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mh_MKS},
   {"M_He4",     0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mAlpha_MKS},
   {"END",       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
}; 

/* 
expression : factor * expression
           | factor / expression
           | factor

factor : paren ^ number
       | paren

paren : ( expression )
      | var

var : name

*/

static int match(char *s[], const char* token) {
   if (strncmp(*s, token, strlen(token)) == 0) {
      *s += strlen(token);
      return 1;
   }

   return 0;
}

static int parse_simple(char *s[], double expSum[], double* expVal) {
   char *original=*s;

   if (**s == '1') {
      ++*s;
      for (int i=0;i<7;i++) expSum[i] = 0;
      *expVal = 1;
      return 1;
   }
   
   int length=0;
   while (('a' <= **s && **s <= 'z')
          ||
          ('A' <= **s && **s <= 'Z')
          ||
          (**s == '_')
   ) {
      ++*s;
      length++;
   }
   
   for (int ii=0; ii<sizeof(symbol_table)/sizeof(symbol_table[0]); ++ii) {
      char* token = symbol_table[ii].name;
      if (strncmp(original, token, length) == 0) {
         double *exp = &(symbol_table[ii].length);
         for (int i=0;i<7;i++) expSum[i] = exp[i];
         *expVal = symbol_table[ii].mks_value;
         return 1;
      }
   }

   *s = original;
   return 0;
}

static int parse_number(char *s[], double *num) {
   char *original=*s;

   if (**s == '+') { ++*s; }
   else if (**s == '-') {++*s; }
   int beforeDot = 0;
   int afterDot = 0;
   while ('0' < **s && **s < '9') {
      ++*s;
      beforeDot = 1;
   }
   if (**s == '.') {
      while ('0' < **s && **s < '9') {
         ++*s;
         afterDot = 1;
      }
   }
   if (!beforeDot && !afterDot) {
      *s = original;
      return 0;
   }
   if (**s == 'e' || **s == 'E') {
      ++*s;
      if (**s == '+') { ++*s; }
      else if (**s == '-') {++*s; }
      int afterE = 0;
      while ('0' < **s && **s < '9') {
         ++*s;
         afterE = 1;
      }
      if (!afterE) {
         *s = original;
         return 0;
      }
   }
   char tmp = **s;
   **s = '\0';
   *num = atof(original);
   **s = tmp;

   return 1;
}

static int parse_expression(char *s[], double expSum[], double* expVal);

static int parse_paren(char *s[], double expSum[], double* expVal) {
   char *original=*s;
   if (match(s, "(") && parse_expression(s,expSum,expVal) && match(s,")")) {
      return 1;
   }
   *s = original;

   if (parse_simple(s,expSum,expVal)) {
      return 1;
   }
   *s = original;

   return 0;
}

static int parse_factor(char *s[], double expSum[], double* expVal) {
   char *original=*s;

   if (parse_paren(s,expSum,expVal)) {
      double exponent = 1;
      char* secondOrig = *s;
      if (match(s, "^") && parse_number(s, &exponent)) {
         for (int i=0;i<7;i++) expSum[i] *= exponent;
         *expVal = pow(*expVal,exponent);
         return 1;
      }
      *s = secondOrig;

      return 1;
   }
   *s = original;

   return 0;
}

static int parse_expression(char *s[], double expSum[], double* expVal) {
   char *original=*s;

   if (parse_factor(s,expSum,expVal)) {
      double secondExpSum[7];
      double secondExpVal;
      char* secondOrig = *s;
      if (match(s, "/") && parse_expression(s,secondExpSum,&secondExpVal)) {
         //do the division
         for (int i=0;i<7;i++) expSum[i] -= secondExpSum[i];
         *expVal /= secondExpVal;
         return 1;
      }
      *s = secondOrig;
      
      if (match(s, "*") && parse_expression(s,secondExpSum,&secondExpVal)) {
         //do the multiplication
         for (int i=0;i<7;i++) expSum[i] += secondExpSum[i];
         *expVal *= secondExpVal;
         return 1;
      }
      *s = secondOrig;

      return 1;
   }
   
   *s = original;
   return 0;

}


double units_parse(const char* unit_string, double* expSum)
{
   double retVal;
   char* orig_string = strdup(unit_string);
   char* mod_string = orig_string;
   if (!parse_expression(&mod_string, expSum, &retVal) || mod_string[0] != '\0') {
      printf("Couldn't understand unit %s\n", orig_string);
      assert(0);
   }
   free(orig_string);
   return retVal;
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
   symbol_table[PRESSURE].mks_value = mass/(time*time*length);
   symbol_table[VELOCITY].mks_value = length/time;
   symbol_table[CAPACITANCE].mks_value = time*time*time*time*current*current/(length*length*mass);
   symbol_table[RESISTIVITY].mks_value = current*current*time*time*time/(mass*length*length);
   symbol_table[VOLTAGE].mks_value = mass*length*length/(current*time*time*time);
   double charge = current*time; 
   kB = kB_MKS * temperature/energy; 
   hbar = hbar_MKS/(energy*time); 
   mu0 = mu0_MKS/((energy/length) * current*current); 
   ke = ke_MKS*charge*charge/(energy*length); 
   hc =  h_MKS*c_MKS/(energy*length);
   mec2 = mec2_MKS/energy; 
   bohr = a0_MKS/length; 
   re = re_MKS/length; 
   M_e = me_MKS/mass; 
   M_p = mp_MKS/mass; 
   M_n = mn_MKS/mass; 
   M_D = mD_MKS/mass; 
   M_T = mT_MKS/mass; 
   M_He3 = mh_MKS/mass; 
   M_He4 = mAlpha_MKS/mass; 
   amu   = u_MKS/mass; 
   Eh =   Eh_MKS/energy; 
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
   symbol[CAPACITANCE].mks_value = time*time*time*time*current*current/(length*length*mass);
   symbol[RESISTIVITY].mks_value = current*current*time*time*time/(mass*length*length);
   symbol[VOLTAGE].mks_value = mass*length*length/(current*time*time*time);
}
int units_check(const char* u1, const char* u2) 
{
   double exp1[7],exp2[7]; 
   if ((u1 == NULL) || (u2 == NULL)) return 0; 
   units_parse(u1,exp1); 
   units_parse(u2,exp2); 
   double sum =0.0; 
   for (int i=0;i<7;i++) sum+= fabs(exp1[i]-exp2[i]); 
   if ( sum > 1e-8)  return 1 ;
   return 0; 
}
double units_convert(double value, const char* from, const char* to)
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
