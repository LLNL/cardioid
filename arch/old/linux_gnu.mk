#
 PLAT=LINUX
#-------------------------------------------------------------------------------

 CXX=g++
 LD=$(CXX)

 DFLAGS += -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS \
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 INCLUDE = 
 CXXFLAGS= -g -ggdb -O3 $(INCLUDE) $(DFLAGS)
# CXXFLAGS= -g -O3 $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS =

 LDFLAGS = $(LIBPATH) $(LIBS)
