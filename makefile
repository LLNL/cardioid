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
#MFEM_DIR=../mfem-3.3
CONFIG_MK = $(MFEM_DIR)/config/config.mk
TEST_MK = $(MFEM_DIR)/config/test.mk
# Use the MFEM install directory
# MFEM_DIR = ../mfem
# CONFIG_MK = $(MFEM_DIR)/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

SEQ_EXAMPLES = fiber 
#SEQ_EXAMPLES = mfemTest
PAR_EXAMPLES = fiberp

ifeq ($(MFEM_USE_MPI),NO)
   EXAMPLES = $(SEQ_EXAMPLES)
   SOURCE = io.cpp fiber.cpp solver.cpp utils.cpp triplet.cpp genfiber.cpp cardfiber.cpp cardgradients.cpp
   OBJECT = $(SOURCE:.cpp=.o)
else
   # MPI C Compiler for PIO.
   MPI_CC=mpicc

  ifeq ($(MFEM_DEBUG), NO)
       CC_FLAGS = -O3 --std=gnu99
  else
       CC_FLAGS = -g --std=gnu99
  endif

   DDCMD_FILES = \
        codata.h \
        ddcMalloc.c \
        ddcMalloc.h \
        ddcMath.h \
        error.c \
        error.h \
        external.h \
        GridAssignmentObject.c \
        GridAssignmentObject.h \
        hardwareInfo.c \
        hardwareInfo.h \
        heap.c \
        heap.h \
        ioUtils.c \
        ioUtils.h \
        intQueue.c \
        intQueue.h \
        lessThan.c \
        lessThan.h \
        match.c \
        match.h \
        mpiUtils.c \
        mpiUtils.h \
        object.c \
        object.h \
        pio.c \
        pio.h \
        pioFixedRecordHelper.c \
        pioFixedRecordHelper.h \
        pioHelper.c \
        pioHelper.h \
        pioVariableRecordHelper.c \
        pioVariableRecordHelper.h \
        tagServer.c \
        tagServer.h \
        three_algebra.c \
        three_algebra.h \
        units.c \
        units.h \
        utilities.c \
        utilities.h


   DDCMDSRC = $(filter %.c, $(DDCMD_FILES))	
   EXAMPLES = $(PAR_EXAMPLES) 
   FIBER_SOURCE = io.cpp fiberp.cpp solver.cpp utils.cpp triplet.cpp genfiber.cpp cardfiber.cpp cardgradientsp.cpp
   FIBER_OBJECT = $(FIBER_SOURCE:.cpp=.o)
   DDCMD_OBJECT = $(DDCMDSRC:.c=.o)
   SOURCE = $(FIBER_SOURCE) $(DDCMD_FILES)
   OBJECT = $(FIBER_OBJECT) $(DDCMD_OBJECT)

#print:
#	echo "DDCMDSRC", $(DDCMDSRC)
#	echo "OBJECT", $(OBJECT)
#	echo "SOURCE", $(SOURCE)


endif


.PHONY: all clean clean-build clean-exec

all: $(SOURCE) $(EXAMPLES)
	
$(EXAMPLES): $(OBJECT)
	$(MFEM_CXX) $(MFEM_FLAGS) $(OBJECT) -o $@ $(MFEM_LIBS)
	
%.o : %.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) -c $(<) -o $(@)

ifeq ($(MFEM_USE_MPI),YES)

$(DDCMD_FILES):
	./mkLinks_ddcMD.sh $@

%.o : %.c
	$(MPI_CC) $(CC_FLAGS) -c $(<) -o $(@)

endif

test: mfemTest kdtreeTest

mfemTest:
	$(MFEM_CXX) $(MFEM_FLAGS) -c mfemTest.cpp -o mfemTest.o
	$(MFEM_CXX) $(MFEM_FLAGS) mfemTest.o solver.o io.o -o $@ $(MFEM_LIBS)
	./mfemTest
	
kdtreeTest:
	$(MFEM_CXX) $(MFEM_FLAGS) kdtreeTest.cpp -o $@ 
	./kdtreeTest

surface:
	$(MFEM_CXX) $(MFEM_FLAGS) -c surface.cpp -o surface.o
	$(MFEM_CXX) $(MFEM_FLAGS) surface.o solver.o io.o -o $@ $(MFEM_LIBS)
	

# Testing: "test" target and mfem-test* variables are defined in config/test.mk

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec clean-test

clean-build:
	rm -f *.o *~ $(SEQ_EXAMPLES) $(PAR_EXAMPLES)
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -f refined.mesh displaced.mesh mesh.* ex5.mesh
	@rm -rf Example5* Example9* Example15*
	@rm -f sphere_refined.* sol.* sol_u.* sol_p.*
	@rm -f ex9.mesh ex9-mesh.* ex9-init.* ex9-final.*
	@rm -f deformed.* velocity.* elastic_energy.* mode_*

clean-test:
	rm -rf mfemTest.o mfemTest kdtreeTest.o kdtreeTest surface

distclean:
	rm -f *.o *~ $(SEQ_EXAMPLES) $(PAR_EXAMPLES)
	rm -rf *.dSYM *.TVD.*breakpoints
	rm -f $(DDCMD_FILES)
