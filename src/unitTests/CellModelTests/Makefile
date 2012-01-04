# $Id: Makefile 7 2010-03-11 21:50:18Z richards12 $

.PHONY: DEFAULT clean distclean 

EXE = OneCell
DEFAULT: first

HEART_FILES = \
	TT04_CellML.hh \
	TT06Func.cc \
	TT06Func.hh \
	TT06Tau.cc \
	TT06Tau.hh \
	pade.hh \
	pade.cc \
	TT04_CellML_Endo.cc \
	TT04_CellML_Endo.hh \
	TT04_CellML_Epi.cc \
	TT04_CellML_Epi.hh \
	TT04_CellML_Mid.cc \
	TT04_CellML_Mid.hh \
	TT06_CellML.hh \
	TT06_CellML_Endo.cc \
	TT06_CellML_Endo.hh \
	TT06_CellML_Epi.cc \
	TT06_CellML_Epi.hh \
	TT06_CellML_Mid.cc \
	TT06_CellML_Mid.hh \
	IBM_tenTusscher04.cc \
	IBM_tenTusscher04.hh \
	IBM_tenTusscher04ParameterDefNLoad.hh \
	IBM_tenTusscher04_endoLUT.cc \
	IBM_tenTusscher04_endoLUT.hh \
	IBM_tenTusscher04_epiLUT.cc \
	IBM_tenTusscher04_epiLUT.hh \
	IBM_tenTusscher04_midLUT.cc \
	IBM_tenTusscher04_midLUT.hh \
	Reaction.hh \
	TT04Dev_Reaction.cc \
	TT04Dev_Reaction.hh \
	TT06Dev_Reaction.cc \
	TT06Dev_Reaction.hh \
	TT04_CellML_Reaction.cc \
	TT04_CellML_Reaction.hh \
	TT06_CellML_Reaction.cc \
	TT06_CellML_Reaction.hh \
	TT04_bbReaction.cc \
	TT04_bbReaction.hh \
	ReactionFHN.cc \
	ReactionFHN.hh \
	Anatomy.hh \
	AnatomyCell.hh \
	Tuple.hh \
	mpiUtils.h \
	mpiUtils.c \
	IndexToTuple.hh \
	Long64.hh \

CXX=mpig++
CC=mpigcc

LDFLAGS=  -lgsl -lgslcblas  
OPTFLAGS=-g -O3    
#OPTFLAGS=-g    

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

first: $(HEART_FILES)
	@$(MAKE) $(EXE)

$(EXE): $(HEART_FILES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

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

