#
 PLAT=LINUX_X86_64
#-------------------------------------------------------------------------------

 CC=/usr/local/bin/mpiicc
 CXX=/usr/local/bin/mpiicpc
 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS \
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 INCLUDE = 
 CXXFLAGS= -g -O3 -xW $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS = -openmp

 LDFLAGS = $(LIBPATH) $(LIBS)
