#
 PLAT=IBM_BGQ
#-------------------------------------------------------------------------------

 CXX=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx
 CC=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc

 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = 
 OPTFLAGS = -g -O3
 CXXFLAGS = $(OPTFLAGS) $(INCLUDE) $(DFLAGS)
 CFLAGS =   $(OPTFLAGS) $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(OPTFLAGS) $(LIBS)
