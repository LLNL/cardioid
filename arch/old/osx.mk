#------------------------------------------------------------------------

CXX=mpicxx -fopenmp
CC =mpicc --std=gnu99
LD=$(CXX)
NVCC = nvcc

DFLAGS = -DOSX -DWITH_PIO -DWITH_MPI \
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

INCLUDE = -I/opt/local/include/ -I/Developer/NVIDIA/CUDA-8.0/include

CFLAGS_BASE =   $(INCLUDE) $(DFLAGS)
CXXFLAGS_BASE = $(INCLUDE) $(DFLAGS)

DFLAGS += -DHAVE_LAPACK -DHAVE_CUDA
LDFLAGS_BASE = -L/Developer/NVIDIA/CUDA-8.0/lib -framework CUDA -ldl -lcudart -lnvrtc -llapack -lblas

CFLAGS_OPT =   $(CFLAGS_BASE) -g -O3
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -ggdb -O0 
CFLAGS_PROF =  $(CFLAGS_BASE) -g -pg -O3 -DPROFILE

CXXFLAGS_OPT =   $(CXXFLAGS_BASE) -g -O3
CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g -ggdb -O0 
CXXFLAGS_PROF =  $(CXXFLAGS_BASE) -g -pg -O3 -DPROFILE

LDFLAGS_OPT   = $(LDFLAGS_BASE) $(CFLAGS_OPT) $(CXXFLAGS_OPT) 
LDFLAGS_DEBUG = $(LDFLAGS_BASE) $(CFLAGS_DEBUG) $(CXXFLAGS_DEBUG)
LDFLAGS_PROF  = $(LDFLAGS_BASE) $(CFLAGS_PROF) $(CXXFLAGS_PROF)
