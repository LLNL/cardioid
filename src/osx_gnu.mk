#
 PLAT=OSX
#-------------------------------------------------------------------------------

 CXX=mpic++
 CC =mpicc
 LD=$(CXX)

 DFLAGS += -D$(PLAT) \
	-DWITH_MPI -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 INCLUDE = 
 OPTFLAGS = -g -O3 -I/opt/local/include
 CXXFLAGS= $(OPTFLAGS) $(INCLUDE) $(DFLAGS)
 CFLAGS= --std=gnu99 $(OPTFLAGS) $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(LIBS) -L/opt/local/lib -lgsl -lgslcblas

