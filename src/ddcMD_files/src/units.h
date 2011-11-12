/* $Id$ */
#ifndef UNITS_H
#define UNITS_H

#ifdef __cplusplus
extern "C" {
#endif

/** When NULL is passed to either to or from argument it stands for
 *  internal units.
 *
 *  Examples:
 *  To obtain the conversion factor from internal units to Angstroms:
 *  length_convert = units_convert(1.0, NULL, "Angstrom"); 
 */
int units_check(char *u1, char*u2) ;
double units_convert(double value,char *from , char *to );
int units_check(char* from, char* to );
void units_internal(double length, double mass, double time,
		    double current, double  temperature,
		    double amount, double luminous_intensity);
void units_external(double length, double mass, double time,
		    double current, double  temperature,
		    double amount, double luminous_intensity);

#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
