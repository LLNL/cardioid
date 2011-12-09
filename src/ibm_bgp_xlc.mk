#
 PLAT=IBM_BGP
#-------------------------------------------------------------------------------

 CXX=mpixlcxx
 CC=mpixlc

 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = 
 OPTFLAGS = -g -O3 -qarch=450 -qtune=450 
 CXXFLAGS = $(OPTFLAGS) $(INCLUDE) $(DFLAGS)
 CFLAGS =   $(OPTFLAGS) $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(OPTFLAGS) $(LIBS)
