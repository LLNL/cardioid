#
 PLAT=LINUX_X86_64
#-------------------------------------------------------------------------------

 CC=/usr/local/bin/mpiicc
 CXX=/usr/local/bin/mpiicpc
 CC=/usr/local/bin/mpiicc
 CXX=mpig++
 CC=mpicc

 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS \
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 INCLUDE = 
 CXXFLAGS= -g  -O3 $(INCLUDE) $(DFLAGS)
 CXXFLAGS= -g   $(INCLUDE) $(DFLAGS)
 CFLAGS= -std=c99 -g -xW $(INCLUDE) $(DFLAGS)
 CFLAGS= --std=gnu99 -g -O3 $(INCLUDE) $(DFLAGS)
 LIBPATH = 
 LIBS = 

 LDFLAGS = $(LIBPATH) $(LIBS)
