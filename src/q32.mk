#
 PLAT=IBM_BGQ
#-------------------------------------------------------------------------------

 CXX=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r
 CC=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc_r

 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = 
 OPTFLAGS = -g -O3
 CXXFLAGS = -qsmp=omp $(OPTFLAGS) $(INCLUDE) $(DFLAGS)
 CFLAGS =   -qsmp=omp $(OPTFLAGS) $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = -qsmp=omp $(LIBPATH) $(OPTFLAGS) $(LIBS)
