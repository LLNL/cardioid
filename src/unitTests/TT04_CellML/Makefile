# $Id: Makefile 7 2010-03-11 21:50:18Z richards12 $

EXE = SingleCell

HEART_FILES = \
	TT04_CellML.hh \
	TT04_CellML_Endo.cc \
	TT04_CellML_Endo.hh \
	TT04_CellML_Epi.cc \
	TT04_CellML_Epi.hh \
	TT04_CellML_Mid.cc \
	TT04_CellML_Mid.hh \



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


C_SRCS = $(wildcard *.c)
CXX_SRCS = $(wildcard *.cc)



OBJECTS = $(patsubst  %.c,  %.o, $(C_SRCS))
OBJECTS += $(patsubst  %.cc,  %.o, $(CXX_SRCS))

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

