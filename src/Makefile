# Get name of architecture-specific include file from hostname
SHELL = /bin/bash

.PHONY: DEFAULT COMMON_TARGETS makedirs clean distclean depend

include Makefile.arch
# If defined, use architecture file set with ARCH variable
ARCH ?= $(ARCHGUESS)
include arch/$(ARCH).mk
#------------------------------------------------------------------------------
#

BINDIR = ../bin

HEARTSRC=cardioid.cc \
	 Vector.cc PerformanceTimers.cc unionOfStrings.cc \
         GDLoadBalancer.cc CommTable.cc \
	 AnatomyReader.cc object_cc.cc \
	 Koradi.cc GridRouter.cc Grid3DStencil.cc writeCells.cc \
	 DomainInfo.cc \
	 initializeSimulate.cc \
	 initializeAnatomy.cc assignCellsToTasks.cc getRemoteCells.cc \
	 simulationLoop.cc \
	 diffusionFactory.cc reactionFactory.cc \
	 conductivityFactory.cc stimulusFactory.cc sensorFactory.cc \
	 FibreConductivity.cc UniformConductivity.cc JHUConductivity.cc\
	 Saleheen98Diffusion.cc SaleheenDev.cc \
	 ReactionFHN.cc \
	 TT04_CellML_Reaction.cc \
	 TT06_CellML_Reaction.cc  \
	 TT06Tau.cc \
	 TT06Func.cc \
         pade.cc \
	 TT06Dev_Reaction.cc \
	 TT04_CellML_Endo.cc TT04_CellML_Mid.cc TT04_CellML_Epi.cc \
	 TT06_CellML_Endo.cc TT06_CellML_Mid.cc TT06_CellML_Epi.cc \
	 TT06_JRG_Reaction.cc \
	 TestStimulus.cc PointStimulus.cc BoxStimulus.cc \
	 PointListSensor.cc ActivationTimeSensor.cc


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
	gsl.c \
	gsl.h \
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
	utilities.h \


DDCMDSRC = $(filter %.c, $(DDCMD_FILES))

HEARTOBJS = $(HEARTSRC:%.cc=$(OBJDIR)/%.o)
HEARTOBJS += $(DDCMDSRC:%.c=$(OBJDIR)/%.o)

BUILD_SUFFIX := $(ARCH)
OBJDIR := objs-$(BUILD_SUFFIX)
EXENAME = cardioid


DEFAULT: opt

COMMON_TARGETS: $(DDCMD_FILES) makedirs

opt:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/cardioid-$(ARCH) \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT)" CXXFLAGS="$(CXXFLAGS_OPT)" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=cardioid-$(ARCH)
debug:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/cardioid-$(ARCH)-debug \
	BUILD_SUFFIX=$(ARCH)-debug OBJDIR=objs-$(ARCH)-debug \
	CFLAGS="$(CFLAGS_DEBUG)" CXXFLAGS="$(CXXFLAGS_DEBUG)" \
	LDFLAGS="$(LDFLAGS_DEBUG)" EXENAME=cardioid-$(ARCH)-debug 
profile:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/cardioid-$(ARCH)-prof \
	BUILD_SUFFIX=$(ARCH)-prof OBJDIR=objs-$(ARCH)-prof \
	CFLAGS="$(CFLAGS_PROF)" CXXFLAGS="$(CXXFLAGS_PROF)" \
	LDFLAGS="$(LDFLAGS_PROF)"  EXENAME=cardioid-$(ARCH)-prof




$(OBJDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(BINDIR)/$(EXENAME): $(HEARTOBJS)
	$(LD) $(DFLAGS) -o $(BINDIR)/$(EXENAME) $(HEARTOBJS) $(LDFLAGS)

testLoadBalancer: makedirs $(HEARTOBJS) $(OBJDIR)/testLoadBalancer.o
	$(LD) $(DFLAGS) -o testLoadBalancer $(OBJDIR)/testLoadBalancer.o $(HEARTOBJS) $(LDFLAGS)

testGridRouter: makedirs $(HEARTOBJS) $(OBJDIR)/testGridRouter.o
	$(LD) $(DFLAGS) -o testGridRouter $(OBJDIR)/testGridRouter.o $(HEARTOBJS) $(LDFLAGS)

$(DDCMD_FILES):
	./mkLinks_ddcMD.sh $@

ddcMD_dist:
	./mkDist_ddcMD.sh $(DDCMD_FILES)

makedirs:
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR) ; fi
	@if [ ! -e .depend_$(ARCH) ]; then touch .depend_$(ARCH) ; fi
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR) ; fi

#------------------------------------------------------------------------------
# generate dependencies in makefile: use -Y to avoid library header files
# that are likely to be different on other platforms.
 depend :
	@makedepend -f .depend_$(BUILD_SUFFIX) -p'$(OBJDIR)/' -o '.o' -Y $(CFLAGS) $(CXXFLAGS) *.cc *.c *.hh *.h 2> /dev/null

.depend_$(BUILD_SUFFIX): $(HEARTSRC)
	@touch .depend_$(BUILD_SUFFIX)
	@$(MAKE) depend

#------------------------------------------------------------------------------
#  Cleanup object files
clean :
	rm -rf $(OBJDIR)-*/*.o
	rm -f $(BINDIR)/$(EXENAME)-$(ARCH)*
	rm -f .depend_$(ARCH)*
	rm -f testLoadBalancer
	rm -f testGridRouter

distclean: clean
	rm -f $(DDCMD_FILES)
	rm -rf ddcMD_files
	rm -rf $(OBJDIR)*
	rm -f .depend*
	rm -rf $(BINDIR)

-include .depend_$(ARCH)
