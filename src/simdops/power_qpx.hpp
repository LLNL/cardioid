#pragma once

#if defined(SIMDOPS_ARCH_POWER_QPX)
#include <builtins.h>


namespace simdops {
   
typedef vector4double native_vector_type;

#define SIMDOPS_SIZE 4

#define load(x) vec_ld(&x)
#define store(x,y) vec_store(y,0,x)
#define splat(x) vec_spats(x)
#define make_float(x) splat(x)

#define add vec_add
#define sub vec_sub
#define mul vec_mul
#define div vec_swdiv_nochk

}

#endif
