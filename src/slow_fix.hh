#ifndef SLOW_FIX
//#define SLOW_FIX
#endif


// The LEGACY_SORTING macro controls minor sorting
// order of cells in the TT06Dev_Reaction constructor.
// Setting this to 1 (one) uses the old and slower
// sorting scheme, and setting it to 0 (zero) enables
// an improved scheme.
#ifndef LEGACY_SORTING
#define LEGACY_SORTING 0
#endif

/* The LEGACY_NG_WORKPARTITION macro controls workpartitioning
   for the non-gate code and the integrator. Not defining it uses
   new code, which also enables per core barriers. It runs about
   1.3us faster (2%) than the old code for a 16x16x14 block
   performance run, and 2us faster (0.5%) for a 48x48x14 block
   performance run.
*/
#ifndef LEGACY_NG_WORKPARTITION
//#define LEGACY_NG_WORKPARTITION
#endif


/* Defining PRINT_WP to something ther than 0 prints the
   work partitioning information for the gate, non-gate,
   and integrator subroutines. This incurs a slowdown on
   the first iteration, and generates for each node, three
   lines of output for each thread.
*/
#ifndef PRINT_WP
#define PRINT_WP 0
#endif
