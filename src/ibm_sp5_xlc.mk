#
 PLAT=IBM_SP5
#-------------------------------------------------------------------------------

 CXX=mpxlC
 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS 
 INCLUDE = 
 CXXFLAGS= -g -O3 -qarch=pwr5 -qtune=pwr5 -qmaxmem=-1  $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(LIBS)
