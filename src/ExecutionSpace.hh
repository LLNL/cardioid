#pragma once

typedef int ExecutionSpace;

#define NONE (-1)
#define CPU (0)
#ifdef USE_CUDA
# define GPU (1)
# define NUMSPACES (2)
#else
# define NUMSPACES (1)
#endif
