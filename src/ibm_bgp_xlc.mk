#
 PLAT=IBM_BGP
#-------------------------------------------------------------------------------

 CXX=mpixlcxx_r
 CC=mpixlc_r

 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = 
 OPTFLAGS = -g -O3 -qarch=450 -qtune=450 
 CXXFLAGS = -qsmp=omp $(OPTFLAGS) $(INCLUDE) $(DFLAGS)
 CFLAGS =   -qsmp=omp $(OPTFLAGS) $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = -qsmp=omp $(LIBPATH) $(OPTFLAGS) $(LIBS)
