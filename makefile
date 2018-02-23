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
TEST_MK = $(MFEM_DIR)/test.mk
# Use the MFEM install directory
# MFEM_DIR = ../mfem
# CONFIG_MK = $(MFEM_DIR)/config.mk

-include $(CONFIG_MK)

.PHONY: all clean clean-build clean-exec

TARGETS=ecg

all: $(SOURCE) $(TARGETS)


MPI_FLAGS=-I/usr/local/Cellar/mpich/3.2.1_1/include -L/usr/local/Cellar/mpich/3.2.1_1/lib -lmpicxx -lmpi -lpmpi	
ecg: $(OBJECTS)
	clang++ -std=c++11 -stdlib=libc++ $(MFEM_FLAGS) $(OBJECTS) ecg.cpp -o $@ $(MFEM_LIBS) $(MPI_FLAGS)
	
%.o : %.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) -c $(<) -o $(@)

# Testing: "test" target and mfem-test* variables are defined in config/test.mk

clean: clean-build clean-exec clean-test

clean-build:
	rm -f *.o *~ $(EXAMPLES)

clean-exec:
	@rm -f sphere_refined.* sol.* sol_u.* sol_p.*
	@rm -f deformed.* velocity.* elastic_energy.* mode_*

clean-test:

distclean:
	rm -f *.o *~ $(EXAMPLES)
	rm -rf *.dSYM *.TVD.*breakpoints
