###############################################################
# Automatically-generated sample makefile to illustrate how to  
# link against oomph-lib from outside the automake/autoconf
# framework. Do not edit this -- make a copy first
# 
# When customising this makefile, you should only have to change
# 
# - the variable OOMPH-LIB_LIBS:
#         Add any additional oomph-lib sub-libraries that 
#         you may wish to use in your code. 
# 
# - the specific dependencies for your driver code: 
#         Include any additional local dependencies such as 
#         your own header files etc.
# 
###############################################################
 
OOMPH_DIR = ../../../base/oomph
 
# Installation-specific information -- don't change any of this! 
#-------------------------------------------------------------- 
 
# Flags for C pre-processor 
AM_CPPFLAGS=-DHAVE_CONFIG_H -I. -I$(OOMPH_DIR)  -DOOMPH_HAS_MPI -DOOMPH_HAS_TRIANGLE_LIB -DOOMPH_HAS_TETGEN_LIB -DUSING_OOMPH_SUPERLU -DUSING_OOMPH_SUPERLU_DIST -I$(OOMPH_DIR)/build/include
 
# Library include directory: This is where all the header files live
OOMPH-LIB_INCLUDE_DIR=$(OOMPH_DIR)/build/include
 
# Library directory: This is where all of oomph-lib's sub-libraries live
OOMPH-LIB_LIB_DIR=$(OOMPH_DIR)/build/lib
 
# These are the external (3rd party) libraries that are distributed
# with oomph-lib and that we always link against
OOMPH-LIB_EXTERNAL_LIBS=-loomph_hsl -loomph_arpack -loomph_triangle -loomph_tetgen -loomph_superlu_3.0 -loomph_parmetis_3.1.1 -loomph_superlu_dist_2.3 -loomph_metis_from_parmetis_3.1.1 -loomph_lapack -loomph_flapack -loomph_blas
 
# This specifies where libraries built from third party 
# distributions can be found
EXTERNAL_DIST_LIBRARIES=
 
# This is additional machine-specific linking information that 
# allows mixed-language compilation/linking
FLIBS=-L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/../../../../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/../../.. -L/lib/../lib64 -L/usr/lib/../lib64 -lfrtbegin -lg2c -lm -lgcc_s
 
# Flags required for the use of shared libraries 
SHARED_LIBRARY_FLAGS=-Wl,--rpath -Wl,$(OOMPH_DIR)/build/lib
 
#Mac OSX: Replace the above line with the following
#SHARED_LIBRARY_FLAGS= --rpath=$(OOMPH_DIR)/build/lib
# Problem-specific information -- edit this for your driver code 
 
#---------------------------------------------------------------
# These are the specific oomph-lib sub-libraries that we have to link against
# for this driver code -- edit this according to your requirements
# but remember that the order of the libraries matters: List the
# the more specific ones before the more general ones!
OOMPH-LIB_LIBS=-lsolid -lconstitutive -lpoisson -lgeneric
 
OBJECTS =  main.o cardiac_constitutive_laws.o cardiac_elements.o cardiac_nodes.o cardiac_mesh.o cardiac_utilities.o cardiac_problem.o

# Dependencies for this driver code and compile instructions: 
# Which local source (usually *.cc or *.h) files does the
# driver code depend on?
%.o: %.cc
	 /usr/local/bin/mpig++ $(AM_CPPFLAGS) -c $< \
	       -I$(OOMPH-LIB_INCLUDE_DIR) 
  
# Linking instructions: Just declare the target (i.e. the name of the executable)
# and the dependencies (i.e. the object files created above). The rest 
# should not have to be changed. 
main: $(OBJECTS)
	 /usr/local/bin/mpig++ $(SHARED_LIBRARY_FLAGS) $(OBJECTS) -o $@ \
	       -L$(OOMPH-LIB_LIB_DIR) $(EXTERNAL_DIST_LIBRARIES) $(OOMPH-LIB_LIBS) \
	        $(OOMPH-LIB_EXTERNAL_LIBS) $(FLIBS) 
