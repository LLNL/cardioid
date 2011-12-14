#
 PLAT=IBM_BGL
#-------------------------------------------------------------------------------

 CXX=mpxlC
 CC=mpxlc
 LD=$(CXX)

 DFLAGS += -DWITH_MPI -DADD_ -D$(PLAT) -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK
 INCLUDE = -I/bgl/BlueLight/ppcfloor/bglsys/include
 OPTFLAGS = -g -O3 -qarch=440 -qtune=440 
 CXXFLAGS = $(OPTFLAGS) $(INCLUDE) $(DFLAGS)
 CFLAGS = -qlanglvl=stdc99 $(OPTFLAGS) $(INCLUDE) $(DFLAGS)

 LIBPATH = 
 LIBS = -lc -lnss_files -lnss_dns -lresolv

 LDFLAGS = $(LIBPATH) $(LIBS)
