# Get name of architecture-specific include file from hostname
SHELL = /bin/bash

.PHONY: DEFAULT COMMON_TARGETS makedirs clean distclean depend

include Makefile.arch
# If defined, use architecture file set with ARCH variable
ARCH ?= $(ARCHGUESS)
include arch/$(ARCH).mk

SVNVERSION := $(shell svnversion ./ )

#------------------------------------------------------------------------------
#

BINDIR = ../bin

HEARTSRC= \
	ThreadServer.cc \
	ThreadUtils.cc \
	Vector.cc SymmetricTensor.cc ConstantModel.cc ConstantReaction.cc \
	PerformanceTimers.cc unionOfStrings.cc getUserInfo.cc \
        GDLoadBalancer.cc BlockLoadBalancer.cc workBoundBalancer.cc \
	pioBalancer.cc \
	CommTable.cc \
	BucketOfBits.cc \
	stateLoader.cc \
	readPioFile.cc \
	AnatomyReader.cc object_cc.cc \
	Koradi.cc GridRouter.cc Grid3DStencil.cc writeCells.cc \
	checkpointIO.cc \
	DomainInfo.cc \
	initializeSimulate.cc \
	Simulate.cc \
	initializeAnatomy.cc setConductivity.cc assignCellsToTasks.cc \
	getRemoteCells.cc \
	simulationLoop.cc \
	diffusionFactory.cc reactionFactory.cc \
	stimulusFactory.cc sensorFactory.cc \
	FibreConductivity.cc JHUConductivity.cc\
	DiffusionUtils.cc \
	FGRDiffusion.cc \
	FGRDiffusionOMP.cc \
	FGRDiffusionThreads.cc \
	FGRDiffusionStrip.cc \
	FGRUtils.cc \
	Reaction.cc \
	ReactionFHN.cc \
	TT04_CellML_Reaction.cc \
	TT06_CellML_Reaction.cc  \
	TT06Dev_Reaction.cc \
	TT06_RRG_Reaction.cc \
	TT06Tau.cc \
	TT06Func.cc \
	pade.cc \
	 TT04_CellML_Endo.cc TT04_CellML_Mid.cc TT04_CellML_Epi.cc \
	 TT06_CellML_Endo.cc TT06_CellML_Mid.cc TT06_CellML_Epi.cc \
	 TT06_RRG.cc \
	 TestStimulus.cc PointStimulus.cc BoxStimulus.cc \
	 PeriodicPulse.cc RandomPulse.cc \
	 PointListSensor.cc ActivationTimeSensor.cc readSnapshotCellList.cc \
	VoronoiCoarsening.cc DataVoronoiCoarsening.cc GradientVoronoiCoarsening.cc \
	CaAverageSensor.cc \
	MinMaxSensor.cc ProcBox.cc MaxDVSensor.cc \
   StateVariableSensor.cc

SPISRC=	 spi_impl.c
SIMDSRC=  clooper.c TT06GatesSimd.c TT06Gates.c TT06NonGates.c TT06NonGatesSimd.c workBound.c 

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
	utilities.h 

DDCMDSRC = $(filter %.c, $(DDCMD_FILES))

HEARTOBJS = $(HEARTSRC:%.cc=$(OBJDIR)/%.o)
HEARTOBJS += $(DDCMDSRC:%.c=$(OBJDIR)/%.o)
HEARTOBJS += $(SPISRC:%.c=$(OBJDIR)/%.o)
HEARTOBJS += $(SIMDSRC:%.c=$(OBJDIR)/%.o)

CARDIOIDSRC = $(HEARTSRC) cardioid.cc 
CARDIOIDOBJS = $(CARDIOIDSRC:%.cc=$(OBJDIR)/%.o)
CARDIOIDOBJS += $(DDCMDSRC:%.c=$(OBJDIR)/%.o)
CARDIOIDOBJS += $(SPISRC:%.c=$(OBJDIR)/%.o)
CARDIOIDOBJS += $(SIMDSRC:%.c=$(OBJDIR)/%.o)
ifdef TIMEOBJ
  CARDIOIDOBJS += $(TIMEOBJ)
endif

#BUILD_SUFFIX := $(ARCH)
OBJDIR := objs-$(ARCH)
OBJDIR_PRECOMPILED := objs_precompiled-$(ARCH)
EXENAME = cardioid

DFLAGS += -DDiff_Weight_Type_Single

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
mpip:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/cardioid-$(ARCH)-mpip \
	BUILD_SUFFIX=$(ARCH)-mpip OBJDIR=objs-$(ARCH)-mpip \
	CFLAGS="$(CFLAGS_MPIP)" CXXFLAGS="$(CXXFLAGS_MPIP)" \
	LDFLAGS="$(LDFLAGS_MPIP)"  EXENAME=cardioid-$(ARCH)-mpip
flops:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/cardioid-$(ARCH)-hpm \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH)-hpm \
	CFLAGS="$(CFLAGS_OPT) -DNTIMING -DHPM" CXXFLAGS="$(CXXFLAGS_OPT) -DNTIMING -DHPM" \
	LDFLAGS="-L/usr/local/tools/mpitrace/lib -lmpihpm_smp $(LDFLAGS_OPT)"  EXENAME=cardioid-$(ARCH)-hpm 
oss:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/cardioid-$(ARCH)-oss \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT)" CXXFLAGS="$(CXXFLAGS_OPT)" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=cardioid-$(ARCH)-oss \
   LD="/usr/global/tools/openspeedshop/oss-dev/bgq/oss202/bgq/bin/osslink -c pcsamp $(CXX)"

$(OBJDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(BINDIR)/$(EXENAME): $(CARDIOIDOBJS) Version.cc
	@if [ -d $(OBJDIR_PRECOMPILED) ]; then  cp $(OBJDIR_PRECOMPILED)/*  $(OBJDIR); fi
	$(LD) Version.cc $(DFLAGS) -o $(BINDIR)/$(EXENAME) $(CARDIOIDOBJS) $(LDFLAGS) \
	-DSVNVERSION="$(SVNVERSION)" \
	-DCXXFLAGS="$(CXXFLAGS)" \
	-DCFLAGS="$(CFLAGS)" \
	-DLDFLAGS="$(LDFLAGS)" \
	-DPATH="$(PWD)" \
	-DTARGET=$(EXENAME) \
	-DARCH=$(ARCH) 


testLoadBalancer:
	@$(MAKE) --no-print-directory makedirs testLoadBalancer-target \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT)" CXXFLAGS="$(CXXFLAGS_OPT)" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=testLoadBalancer-$(ARCH)

testLoadBalancer-target: $(CARDIOIDOBJS) $(OBJDIR)/testLoadBalancer.o
	$(LD) $(DFLAGS) -o $(BINDIR)/$(EXENAME) $(OBJDIR)/testLoadBalancer.o $(CARDIOIDOBJS) $(LDFLAGS)

testGridRouter: makedirs $(CARDIOIDOBJS) $(OBJDIR)/testGridRouter.o
	$(LD) $(DFLAGS) -o testGridRouter $(OBJDIR)/testGridRouter.o $(CARDIOIDOBJS) $(LDFLAGS)

compareSnapshots: 
	@$(MAKE) --no-print-directory makedirs compareSnapshots-target \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT)" CXXFLAGS="$(CXXFLAGS_OPT)" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=compareSnapshots-$(ARCH)

compareSnapshots-target: $(HEARTOBJS) $(OBJDIR)/compareSnapshots.o Version.cc
	$(LD) Version.cc $(DFLAGS) -o $(BINDIR)/$(EXENAME) $(HEARTOBJS) $(OBJDIR)/compareSnapshots.o $(LDFLAGS)

$(DDCMD_FILES):
	./mkLinks_ddcMD.sh $@

ddcMD_dist:
	./mkDist_ddcMD.sh $(DDCMD_FILES)

makedirs:
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR) ; fi
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR) ; fi

#------------------------------------------------------------------------------
# generate dependencies in makefile: use -Y to avoid library header files
# that are likely to be different on other platforms.
 depend :
	makedepend -f .depend_$(BUILD_SUFFIX) -p'$(OBJDIR)/' -o '.o' -Y -- $(CFLAGS) $(CXXFLAGS) -- *.cc *.c *.hh *.h 2> /dev/null

.depend_$(BUILD_SUFFIX): $(CARDIOIDSRC)
	@touch .depend_$(BUILD_SUFFIX)
	@$(MAKE) depend BUILD_SUFFIX=$(BUILD_SUFFIX)

#------------------------------------------------------------------------------
#  Cleanup object files
clean :
	rm -rf $(OBJDIR)/*.o
	rm -f *.lst
	rm -rf $(BINDIR)/$(EXENAME)-$(ARCH)*.*
	rm -f $(BINDIR)/$(EXENAME)-$(ARCH)*
	rm -f .depend_$(ARCH)*
	rm -f $(BINDIR)/testLoadBalancer-$(ARCH)
	rm -f $(BINDIR)/compareSnapshots-$(ARCH)
	rm -f testGridRouter

distclean: clean
	rm -f $(DDCMD_FILES)
	rm -rf ddcMD_files
	rm -rf objs-*
	rm -f .depend*
	rm -rf $(BINDIR)

ifdef BUILD_SUFFIX
  -include .depend_$(BUILD_SUFFIX)
endif
