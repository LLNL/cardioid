#ifndef SLOW_FIX_HH
#define SLOW_FIX_HH

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

inline int
convertActualSizeToBufferSize(int localSize) {
   //FIXME FIXME FIXME
   /* This code is to fix the buffer sizes to match the sizes created
      by both TT06 and Simulate classes.  The code was first
      introduced in 0b947032b9876b89fdf124c4ef0259c96e98ebaf, r1100
      with this commit.
      
    Added code in Simulate.hh and TT06Dev_Reaction.cc to change the
    memory allocation sizes so that the problem of nodes with
    particular numbers of cells running slowly.

    The new code is guarded by a macro defined in the new file slow_fix.hh.
    As of checkin, the macro definition is commented out, so that the new
    code is turned off. Enable by uncommenting the SLOW_FIX macro definition.

    The new code, when turned on, is expected to improve nonGate-times for
    nodes with unfavourable numbers of cells by 4us, assuming about 5300
    cells per node.

      Author: Tomas Oppelstrup

      blake14: this code sends 0-8 -> 8, 9-40 -> 40, 41-72 -> 72 and
      so on.  This seems awfully strange, and I don't know why it's
      needed.
   */
      int paddedSize = 4*((localSize+3)/4);
      {
	int nFourVecs = paddedSize>>2;
	if(0) paddedSize += 4*((10 - (nFourVecs % 8)) % 8);
	else  paddedSize += ((10 - (nFourVecs & 7)) & 7) << 2;
       
      }
      return paddedSize;
}

#endif
