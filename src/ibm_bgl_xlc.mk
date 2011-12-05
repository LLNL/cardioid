#
 PLAT=IBM_BGL
#-------------------------------------------------------------------------------

 CXX=mpxlC
 CC=mpxlc
 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = 
 OPTFLAGS = -g -O3 -qarch=440 -qtune=440 
 CXXFLAGS = $(OPTFLAGS) $(INCLUDE) $(DFLAGS)
 CFLAGS = $(OPTFLAGS) $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(LIBS)
