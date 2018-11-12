# Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the MFEM library. For more information and source code
# availability see http://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

# Use the MFEM build directory
MFEM_DIR=deps/share/mfem
CONFIG_MK = $(MFEM_DIR)/config.mk
#TEST_MK = $(MFEM_DIR)/test.mk
# Use the MFEM install directory
# MFEM_DIR = ../mfem
# CONFIG_MK = $(MFEM_DIR)/config.mk

-include $(CONFIG_MK)

.PHONY: all clean clean-build clean-exec

TARGETS=ecg

all: $(SOURCE) $(TARGETS)

DDCMD_FILES = \
        ddcmdUtil/src/codata.h \
        ddcmdUtil/src/ddcMalloc.c \
        ddcmdUtil/src/ddcMalloc.h \
        ddcmdUtil/src/ddcMath.h \
        ddcmdUtil/src/error.c \
        ddcmdUtil/src/error.h \
        ddcmdUtil/src/external.h \
        ddcmdUtil/src/GridAssignmentObject.c \
        ddcmdUtil/src/GridAssignmentObject.h \
        ddcmdUtil/src/hardwareInfo.c \
        ddcmdUtil/src/hardwareInfo.h \
        ddcmdUtil/src/heap.c \
        ddcmdUtil/src/heap.h \
        ddcmdUtil/src/ioUtils.c \
        ddcmdUtil/src/ioUtils.h \
        ddcmdUtil/src/intQueue.c \
        ddcmdUtil/src/intQueue.h \
        ddcmdUtil/src/lessThan.c \
        ddcmdUtil/src/lessThan.h \
        ddcmdUtil/src/match.c \
        ddcmdUtil/src/match.h \
        ddcmdUtil/src/mpiUtils.c \
        ddcmdUtil/src/mpiUtils.h \
        ddcmdUtil/src/object.c \
        ddcmdUtil/src/object.h \
        ddcmdUtil/src/pio.c \
        ddcmdUtil/src/pio.h \
        ddcmdUtil/src/pioFixedRecordHelper.c \
        ddcmdUtil/src/pioFixedRecordHelper.h \
        ddcmdUtil/src/pioHelper.c \
        ddcmdUtil/src/pioHelper.h \
        ddcmdUtil/src/pioVariableRecordHelper.c \
        ddcmdUtil/src/pioVariableRecordHelper.h \
        ddcmdUtil/src/tagServer.c \
        ddcmdUtil/src/tagServer.h \
        ddcmdUtil/src/three_algebra.c \
        ddcmdUtil/src/three_algebra.h \
        ddcmdUtil/src/units.c \
        ddcmdUtil/src/units.h \
        ddcmdUtil/src/utilities.c \
        ddcmdUtil/src/utilities.h

DDCMDSRC = $(filter %.c, $(DDCMD_FILES))	
DDCMD_OBJECT = $(DDCMDSRC:.c=.o)
OBJECTS = $(DDCMD_OBJECT)

MY_FLAGS = -IddcmdUtil/include -D_USE_MATH_DEFINES -DDiff_Weight_Type_Single -DWITH_PIO -DWITH_MPI -IddcmdUtil/include -I/usr/tce/packages/impi/impi-2018.0-gcc-4.9.3/include -g -DM_PI=3.14159 -std=c++11 -DDEBUG

MFEM_CC = mpicxx #$(MFEM_CXX)
MFEM_CC := $(MFEM_CC:%c++=%cc)
MFEM_CC := $(MFEM_CC:%cxx=%cc)

ecg: $(OBJECTS) 
	$(MFEM_CXX) -std=c++11 $(MFEM_FLAGS) $(MY_FLAGS) $(OBJECTS) ecg.cpp -o $@ $(MFEM_LIBS) -L/usr/tce/packages/impi/impi-2018.0-gcc-4.9.3/lib

pecg: $(OBJECTS) 
	$(MFEM_CXX) -std=c++11 $(MFEM_FLAGS) $(MY_FLAGS) $(OBJECTS) pecg.cpp -o $@ $(MFEM_LIBS) -L/usr/tce/packages/impi/impi-2018.0-gcc-4.9.3/lib
	#$(MFEM_CXX) -std=c++11 $(MFEM_FLAGS) $(MY_FLAGS) $(OBJECTS) pecg.cpp -o $@ $(MFEM_LIBS) -L/usr/tce/packages/impi/impi-2018.0-gcc-4.9.3/lib -L../hypre/src/hypre/lib -I../hypre/src/hypre/include

%.o : %.cpp
	$(MFEM_CXX) -std=gnu++11 $(MFEM_FLAGS) $(MY_FLAGS) -c $(<) -o $(@)

ddcmdUtil/src/%.o : ddcmdUtil/src/%.c
	$(MFEM_CC) -std=c99 $(MFEM_FLAGS) $(MY_FLAGS) -c $< -o $@

# Testing: "test" target and mfem-test* variables are defined in config/test.mk

clean: clean-build clean-exec clean-test

clean-build:
	rm -f *.o *~ $(TARGETS)

clean-exec:
	@rm -f sphere_refined.* sol.* sol_u.* sol_p.*
	@rm -f deformed.* velocity.* elastic_energy.* mode_*

clean-test:

distclean:
	rm -f *.o *~ $(EXAMPLES)
	rm -rf *.dSYM *.TVD.*breakpoints
