#
 PLAT=OSX
#-------------------------------------------------------------------------------

 CXX=mpic++
 CC =mpicc
 LD=$(CXX)

 DFLAGS += -D$(PLAT) \
	-DWITH_MPI -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 INCLUDE = 
 CXXFLAGS= -g -ggdb $(INCLUDE) $(DFLAGS)
 CFLAGS= --std=gnu99 -g -ggdb  $(INCLUDE) $(DFLAGS)
# CXXFLAGS= -g -O3 $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(LIBS)
