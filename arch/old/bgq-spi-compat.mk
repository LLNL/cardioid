BGQ_SDK_PATH = /bgsys/drivers/ppcfloor
#BGQ_SDK_PATH=/bgsys/drivers/DRV2012_0229_1625/ppc64-rhel60
#### SPI ####
BGSYS_INC := -I$(BGQ_SDK_PATH)/comm/sys/include       \
        -I$(BGQ_SDK_PATH)/spi/include                 \
        -I$(BGQ_SDK_PATH)/spi/include/kernel/cnk      \
        -I$(BGQ_SDK_PATH)/spi/include/mu/default/

SPI_INC := -I$(SPI_PATH)/libutil/include
BGSYS_LIBS := -L$(BGQ_SDK_PATH)/lib -lrt -L$(BGQ_SDK_PATH)/spi/lib -lSPI -lSPI_cnk -L$(BGQ_SDK_PATH)/bgpm/lib -lbgpm
############

CXX=$(BGQ_SDK_PATH)/comm/xl/bin/mpixlcxx_r
CC=$(BGQ_SDK_PATH)/comm/xl/bin/mpixlc_r
LD=$(CXX)

DFLAGS = -DWITH_PIO -DWITH_MPI -DBGQ -DSPI \
	 -DADD_ -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK -DNO_GETPWUID
LDFLAGS_BASE = -lc -lnss_files -lnss_dns -lresolv -L/usr/local/tools/mpitrace/lib  -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm
#LDFLAGS_BASE = -lnss_files -lnss_dns -lresolv -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm

INCLUDE = -I$(SPI_PATH)/include
OPTFLAGS = -g -O3
CFLAGS_BASE =   -qsmp=omp -qhot=novector -qsimd=auto $(INCLUDE) $(DFLAGS) $(BGSYS_INC) $(SPI_INC)
CFLAGS_BASE =   -qsmp=omp $(INCLUDE) $(DFLAGS) $(BGSYS_INC) $(SPI_INC)
CXXFLAGS_BASE = -qsmp=omp $(INCLUDE) $(DFLAGS) $(BGSYS_INC) $(SPI_INC)


CFLAGS_OPT =   $(CFLAGS_BASE) -g -O3  
#CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0 -fnoinline
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0 
CFLAGS_PROF =  $(CFLAGS_BASE) -g -pg -O3 -DPROFILE

CXXFLAGS_OPT =   $(CXXFLAGS_BASE) -g  -O3
#CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g -O0 -fnoinline
CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g -O0 
CXXFLAGS_PROF =  $(CXXFLAGS_BASE) -g -pg -O3 -DPROFILE


LDFLAGS_OPT   = $(LDFLAGS_BASE) $(CFLAGS_OPT) $(CXXFLAGS_OPT)
LDFLAGS_DEBUG = $(LDFLAGS_BASE) $(CFLAGS_DEBUG) $(CXXFLAGS_DEBUG)
LDFLAGS_PROF  = $(LDFLAGS_BASE) $(CFLAGS_PROF) $(CXXFLAGS_PROF)
