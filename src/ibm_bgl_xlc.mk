#
 PLAT=IBM_BGL
#-------------------------------------------------------------------------------

 CXX=mpxlC
 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = 
 CXXFLAGS= -g -O3 -qarch=440 -qtune=440 $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(LIBS)
