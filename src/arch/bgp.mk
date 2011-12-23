#-------------------------------------------------------------------------------

CXX=mpixlcxx_r
CC=mpixlc_r
LD=$(CXX)

DFLAGS = -DWITH_PIO -DWITH_MPI -DBGP \
	-DADD_ -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK

INCLUDE =  -I/bgsys/drivers/ppcfloor/arch/include \
	-I/usr/local/tools/gsl/include/

CFLAGS_BASE =   -qsmp=omp -qarch=450d $(INCLUDE) $(DFLAGS)
CXXFLAGS_BASE = -qsmp=omp -qarch=450d $(INCLUDE) $(DFLAGS)


HAVE_GSL = 1
ifeq ($(HAVE_GSL),1) 
   CFLAGS_BASE  += -DHAVE_GSL
   CXXFLAGS_BASE  += -DHAVE_GSL
   LDFLAGS_BASE += -L/usr/local/tools/gsl/lib -lgsl -lgslcblas
endif

CFLAGS_OPT =   $(CFLAGS_BASE) -g -O3 -qtune=450 
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0
CFLAGS_PROF =  $(CFLAGS_BASE) -g -pg -O3 -DPROFILE

CXXFLAGS_OPT =   $(CXXFLAGS_BASE) -g -O3 -qtune=450 
CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g -O0
CXXFLAGS_PROF =  $(CXXFLAGS_BASE) -g -pg -O3 -DPROFILE


LDFLAGS_OPT   = $(LDFLAGS_BASE) $(CFLAGS_OPT) $(CXXFLAGS_OPT)
LDFLAGS_DEBUG = $(LDFLAGS_BASE) $(CFLAGS_DEBUG) $(CXXFLAGS_DEBUG)
LDFLAGS_PROF  = $(LDFLAGS_BASE) $(CFLAGS_PROF) $(CXXFLAGS_PROF)
