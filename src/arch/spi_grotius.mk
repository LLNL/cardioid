BGQ_SDK_PATH = /bgsys/drivers/ppcfloor
BASE = /gsa/rchgsa-p2/16/xlcmpbld/run/vacpp/dev/bgq/daily/120126/opt/ibmcmp/vacpp/bg/12.1/bin
CXX = $(BASE)/bgxlC_r
CC = $(BASE)/bgxlc_r
MPI_PATH:= -L/bgsys/drivers/ppcfloor/comm/xl/lib -lmpich -lmpl -lopa \
           -L/bgsys/drivers/ppcfloor/comm/sys/lib -lpami \
           -L/bgsys/drivers/ppcfloor/spi/lib -lSPI_cnk -lrt

#MPI_TRACE := /bgusr/home2/walkup/mpi_trace/bgq/lib/libmpitrace.a
#MPI_TRACE := /bgusr/home1/hfwen/mpi_trace/bgq/src/libmpitrace.a

#### SPI ####
BGSYS_INC := -I$(BGQ_SDK_PATH)/comm/sys/include                 \
        -I$(BGQ_SDK_PATH) -I$(BGQ_SDK_PATH)/spi/include                 \
        -I$(BGQ_SDK_PATH)/spi/include/kernel/cnk                        \
        -I$(BGQ_SDK_PATH)/spi/include/mu/default/

SPI_INC := -I$(SPI_PATH)/libutil/include
BGSYS_LIBS := -L$(BGQ_SDK_PATH)/lib -lrt -L$(BGQ_SDK_PATH)/spi/lib -lSPI -lSPI_cnk
############


#CXX=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r
#CC=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc_r
LD=$(CXX)

DFLAGS = -DWITH_PIO -DWITH_MPI -DBGQ \
	 -DADD_ -DUSE_CSTDIO_LFS -DMPICH_IGNORE_CXX_SEEK

INCLUDE = -I/bgsys/drivers/ppcfloor -I/bgsys/drivers/ppcfloor/comm/xl/include
OPTFLAGS = -g -O3
CFLAGS_BASE =   -qsmp=omp $(INCLUDE) $(DFLAGS) -DBGQ $(BGSYS_INC) $(SPI_INC) -DSPI
CXXFLAGS_BASE = -qsmp=omp $(INCLUDE) $(DFLAGS) -DBGQ -DSPI


CFLAGS_OPT =   $(CFLAGS_BASE) -g -O3  
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0
CFLAGS_PROF =  $(CFLAGS_BASE) -g -pg -O3 -DPROFILE

CXXFLAGS_OPT =   $(CXXFLAGS_BASE) -g -O3  
CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g -O0
CXXFLAGS_PROF =  $(CXXFLAGS_BASE) -g -pg -O3 -DPROFILE

LDFLAGS_BASE = $(MPI_TRACE) $(MPI_PATH)
LDFLAGS_OPT   = $(LDFLAGS_BASE) $(CFLAGS_OPT) $(CXXFLAGS_OPT)
LDFLAGS_DEBUG = $(LDFLAGS_BASE) $(CFLAGS_DEBUG) $(CXXFLAGS_DEBUG)
LDFLAGS_PROF  = $(LDFLAGS_BASE) $(CFLAGS_PROF) $(CXXFLAGS_PROF)
