#-------------------------------------------------------------------------------

CXX=mpxlC
CC=mpxlc
LD=$(CXX)

DFLAGS = -DWITH_PIO -DWITH_MPI -DBGL \
	-DADD_ -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK

INCLUDE = -I/bgl/BlueLight/ppcfloor/bglsys/include

CFLAGS_BASE =   -qlanglvl=stdc99 -qarch=440d $(INCLUDE) $(DFLAGS)
CXXFLAGS_BASE = -qarch=440 $(INCLUDE) $(DFLAGS)
LDFLAGS_BASE = -lc -lnss_files -lnss_dns -lresolv


CFLAGS_OPT =   $(CFLAGS_BASE) -g -O3 -qtune=440 
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0
CFLAGS_PROF =  $(CFLAGS_BASE) -g -pg -O3 -DPROFILE

CXXFLAGS_OPT =   $(CXXFLAGS_BASE) -g -O3 -qtune=440 
CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g -O0
CXXFLAGS_PROF =  $(CXXFLAGS_BASE) -g -pg -O3 -DPROFILE


LDFLAGS_OPT   = $(LDFLAGS_BASE) $(CFLAGS_OPT) $(CXXFLAGS_OPT)
LDFLAGS_DEBUG = $(LDFLAGS_BASE) $(CFLAGS_DEBUG) $(CXXFLAGS_DEBUG)
LDFLAGS_PROF  = $(LDFLAGS_BASE) $(CFLAGS_PROF) $(CXXFLAGS_PROF)


# LIBS = -lc -lnss_files -lnss_dns -lresolv

