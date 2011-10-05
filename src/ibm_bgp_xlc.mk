#
 PLAT=IBM_BGP
#-------------------------------------------------------------------------------

 CXX=mpixlcxx
 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = 
 CXXFLAGS= -g -O3 -qarch=450 -qtune=450 $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(LIBS)
