# $Id: Makefile 7 2010-03-11 21:50:18Z richards12 $

EXE = IOSpeed

HEART_FILES = \
	object.h object.c \
	pio.h pio.c \
	pioHelper.h pioHelper.c \
	pioFixedRecordHelper.h pioFixedRecordHelper.c \
	pioVariableRecordHelper.h pioVariableRecordHelper.c \
	utilities.h utilities.c \
	heap.h heap.c \
	tagServer.h tagServer.c \
	ddcMalloc.h ddcMalloc.c \
	mpiUtils.h mpiUtils.c \
	match.h match.c \
	units.h units.c \
	error.h error.c \
	hardwareInfo.h hardwareInfo.c \
	codata.h \
	ioUtils.h ioUtils.c \
	external.h \
	lessThan.h lessThan.c


CXX=mpixlcxx_r
CC=mpixlc_r

OPTFLAGS=-g -O3 -DBGQ -DWITH_MPI

CXXFLAGS += $(OPTFLAGS)
CFLAGS += $(OPTFLAGS)


HEART_CXX_SRCS = $(filter %.cc, $(HEART_FILES))
HEART_C_SRCS   = $(filter %.c, $(HEART_FILES))

C_SRCS = $(wildcard *.c)
CXX_SRCS = $(wildcard *.cc)

ALL_C_SRCS   = $(sort $(HEART_C_SRCS)   $(C_SRCS))
ALL_CXX_SRCS = $(sort $(HEART_CXX_SRCS) $(CXX_SRCS))

OBJECTS = $(patsubst  %.c,  %.o, $(ALL_C_SRCS))
OBJECTS += $(patsubst  %.cc,  %.o, $(ALL_CXX_SRCS))


$(EXE): $(HEART_FILES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS)

$(HEART_FILES):
	../tools/mkLinks_heart.sh $@

clean:
	rm -f *.o *.d
distclean:
	rm -f *.o *.d $(EXE) $(HEART_FILES)

links: $(HEART_FILES)
	@echo "links made"

-include $(C_SRCS:.c=.d)
-include $(CXX_SRCS:.cc=.d)

