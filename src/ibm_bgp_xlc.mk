#
 PLAT=IBM_BGP
#-------------------------------------------------------------------------------

 CXX=mpixlcxx_r
 CC=mpixlc_r

 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = -I/usr/local/tools/gsl/include/
 OPTFLAGS = -g -O3 -qarch=450 -qtune=450 
 CXXFLAGS = -qsmp=omp $(OPTFLAGS) $(INCLUDE) $(DFLAGS)
 CFLAGS =   -qsmp=omp $(OPTFLAGS) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L/usr/local/tools/gsl/lib
 LIBS =    -lgsl -lgslcblas

 LDFLAGS = -qsmp=omp $(LIBPATH) $(OPTFLAGS) $(LIBS)
