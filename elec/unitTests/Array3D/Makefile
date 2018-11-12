# $Id: Makefile 7 2010-03-11 21:50:18Z richards12 $

EXE = Array3d_test

HEART_FILES = \
	Array3d.hh \




CXX=g++
CC=gcc

OPTFLAGS=-g -O3

CXXFLAGS += $(OPTFLAGS)
CXXFLAGS += -MMD
CXXFLAGS += -Wall -Wextra
CFLAGS += $(OPTFLAGS)
CFLAGS += -MMD
CFLAGS += -Wall -Wextra
CFLAGS += -std=gnu99


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

