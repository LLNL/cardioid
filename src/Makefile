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
	stringUtils.cc \
	readCellList.cc \
	PioHeaderData.cc \
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
	BoundingBox.cc \
	initializeSimulate.cc \
	Simulate.cc \
	initializeAnatomy.cc setConductivity.cc assignCellsToTasks.cc \
	getRemoteCells.cc \
	simulationLoop.cc \
	diffusionFactory.cc reactionFactory.cc \
	stimulusFactory.cc sensorFactory.cc drugFactory.cc \
   Drug.cc DrugChannel.cc \
	FibreConductivity.cc JHUConductivity.cc\
	DiffusionUtils.cc \
	FGRDiffusion.cc \
	FGRDiffusionOMP.cc \
	FGRDiffusionThreads.cc \
	FGRDiffusionOverlap.cc \
	FGRDiffusionStrip.cc \
	FGRUtils.cc \
	Reaction.cc \
	ReactionFHN.cc \
	TT04_CellML_Reaction.cc \
	TT06_CellML_Reaction.cc  \
	TT06Dev_Reaction.cc \
	TT06_RRG_Reaction.cc \
	OHaraRudy_Reaction.cc \
	SimpleOHaraRudy.cc \
	Grandi_Reaction.cc \
	SimpleGrandi.cc \
	TT06Tau.cc \
	TT06Func.cc \
	pade.cc \
   TT04_CellML_Endo.cc TT04_CellML_Mid.cc TT04_CellML_Epi.cc \
   TT06_CellML_Endo.cc TT06_CellML_Mid.cc TT06_CellML_Epi.cc \
	TT06_RRG.cc \
   OHaraRudyScan.cc \
   GrandiScan.cc  \
	TestStimulus.cc PointStimulus.cc BoxStimulus.cc \
	PeriodicPulse.cc RandomPulse.cc \
	PointListSensor.cc ActivationTimeSensor.cc \
	ActivationAndRecoverySensor.cc \
	VoronoiCoarsening.cc DataVoronoiCoarsening.cc GradientVoronoiCoarsening.cc \
	CaAverageSensor.cc \
	ECGSensor.cc \
	MinMaxSensor.cc ProcBox.cc MaxDVSensor.cc DVThreshSensor.cc \
   StateVariableSensor.cc

OHaraRudySRC= OHaraRudy.c OHaraRudy_Voltage.c OHaraRudy_Concent.c OHaraRudy_INaCai.c OHaraRudy_INaCass.c  \
              OHaraRudy_INaK.c OHaraRudy_INab.c OHaraRudy_ICab.c OHaraRudy_IKb.c \
              OHaraRudy_IpCa.c OHaraRudy_INaFast.c OHaraRudy_INaL.c OHaraRudy_Ito.c \
              OHaraRudy_IKr.c OHaraRudy_IKs.c OHaraRudy_IK1.c OHaraRudy_ICa.c null_INull.c\
              OHaraRudy_Fluxes.c OHaraRudy_CaMKtrap.c reversalPotentials.c
RTYSC14ASRC = OHaraRudyMod_INaFast.c RTYSC14A_IKr.c  MYBGBKC_INa.c

GrandiSRC=    Grandi.c         Grandi_IClCa.c   Grandi_IKs.c     Grandi_INaL.c    Grandi_Revpots.c \
              Grandi_Concent.c Grandi_IClb.c    Grandi_IKur.c    Grandi_INab.c    Grandi_Voltage.c \
              Grandi_Fluxes.c  Grandi_IK1.c     Grandi_INCX.c    Grandi_INull.c \
              Grandi_ICa.c     Grandi_IKp.c     Grandi_INa.c     Grandi_IpCa.c \
              Grandi_ICab.c    Grandi_IKr.c     Grandi_INaK.c    Grandi_Ito.c

SPISRC=	 spi_impl.c
SIMDSRC=  clooper.c TT06GatesSimd.c TT06GatesSimd2.c TT06Gates.c TT06NonGates.c TT06NonGatesSimd.c workBound.c svd.c $(OHaraRudySRC) $(RTYSC14ASRC) $(GrandiSRC)

CONSTRUCTED_INCLUDES=\
OHaraRudy_Voltage.h OHaraRudy_Ito.h OHaraRudy_IpCa.h OHaraRudy_INab.h OHaraRudy_INaL.h OHaraRudy_INaK.h OHaraRudy_INaCass.h \
OHaraRudy_INaCai.h OHaraRudy_IKs.h OHaraRudy_IKr.h OHaraRudy_IKb.h OHaraRudy_IK1.h OHaraRudy_ICab.h OHaraRudy_ICa.h \
OHaraRudy_Fluxes.h OHaraRudy_Concent.h OHaraRudy_CaMKtrap.h OHaraRudy_INaFast.h OHaraRudyMod_INaFast.h RTYSC14A_IKr.h \
MYBGBKC_INa.h null_INull.h

GRANDI_CONSTRUCTED_INCLUDES=\
Grandi_ICa.h Grandi_IKp.h Grandi_INa.h Grandi_IpCa.h Grandi_ICab.h Grandi_IKr.h Grandi_INaK.h Grandi_Ito.h \
Grandi_IClCa.h Grandi_IKs.h Grandi_INaL.h Grandi_Voltage.h Grandi_Concent.h Grandi_IClb.h Grandi_IKur.h Grandi_INab.h \
Grandi_Fluxes.h Grandi_IK1.h Grandi_INCX.h Grandi_INull.h

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
	CFLAGS="$(CFLAGS_OPT) -DNDEBUG -DNTIMING -DHPM" CXXFLAGS="$(CXXFLAGS_OPT) -DNDEBUG -DNTIMING -DHPM" \
	LDFLAGS="-L/usr/local/tools/mpitrace/lib -lmpihpm_smp $(LDFLAGS_OPT)"  EXENAME=cardioid-$(ARCH)-hpm 
oss:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/cardioid-$(ARCH)-oss \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT)" CXXFLAGS="$(CXXFLAGS_OPT)" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=cardioid-$(ARCH)-oss \
   LD="/usr/global/tools/openspeedshop/oss-dev/bgq/oss202/bgq/bin/osslink -c pcsamp $(CXX)"

GENGETOPT ?= gengetopt
gengetopt:
	$(GENGETOPT) -i singleCellOptions.ggo -F singleCellOptions 

$(OBJDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

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

#OHaraRudy.cc : OHaraRudy.h  $(CONSTRUCTED_INCLUDES) mkCompInclude   mkCompInclude.awk

$(CONSTRUCTED_INCLUDES) :  mkCompInclude mkCompInclude.awk 
	mkCompInclude $@

#Grandi.cc : Grandi.h  $(GRANDI_CONSTRUCTED_INCLUDES) mkCompInclude_Grandi   mkCompInclude_Grandi.awk

$(GRANDI_CONSTRUCTED_INCLUDES) :  mkCompInclude_Grandi mkCompInclude_Grandi.awk 
	mkCompInclude_Grandi $@

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

collectGradientTraces: 
	@$(MAKE) --no-print-directory makedirs collectGradientTraces-target \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT)" CXXFLAGS="$(CXXFLAGS_OPT)" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=collectGradientTraces-$(ARCH)

collectGradientTraces-target: $(HEARTOBJS) $(OBJDIR)/collectGradientTraces.o Version.cc
	$(LD) Version.cc $(DFLAGS) -o $(BINDIR)/$(EXENAME) $(HEARTOBJS) $(OBJDIR)/collectGradientTraces.o $(LDFLAGS)

singleCell: 
	@$(MAKE) --no-print-directory makedirs $@-target \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT)" CXXFLAGS="$(CXXFLAGS_OPT)" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=$@-$(ARCH)

singleCell-target: $(HEARTOBJS) $(OBJDIR)/singleCell.o $(OBJDIR)/singleCellOptions.o singleCellOptions.h Version.cc
	$(LD) Version.cc $(DFLAGS) -o $(BINDIR)/$(EXENAME) $(HEARTOBJS) $(OBJDIR)/singleCell.o $(OBJDIR)/singleCellOptions.o $(LDFLAGS)

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
	./mkCompInclude $(CONSTRUCTED_INCLUDES)
	makedepend -f .depend_$(BUILD_SUFFIX) -p'$(OBJDIR)/' -o '.o' -Y -- $(CFLAGS) $(CXXFLAGS) -- *.cc *.c *.hh *.h 2> /dev/null

	./mkCompInclude_Grandi $(GRANDI_CONSTRUCTED_INCLUDES)
	makedepend -f .depend_$(BUILD_SUFFIX) -p'$(OBJDIR)/' -o '.o' -Y -- $(CFLAGS) $(CXXFLAGS) -- *.cc *.c *.hh *.h 2> /dev/null

.depend_$(BUILD_SUFFIX): $(CARDIOIDSRC)
	@touch .depend_$(BUILD_SUFFIX)
	@$(MAKE) depend BUILD_SUFFIX=$(BUILD_SUFFIX)
	rm -r 

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
	rm -f $(BINDIR)/collectGradientTraces-$(ARCH)
	rm -f testGridRouter

distclean: clean
	rm -f $(DDCMD_FILES)
	rm -f $(CONSTRUCTED_INCLUDES)
	rm -rf ddcMD_files
	rm -rf objs-*
	rm -f .depend*
	rm -rf $(BINDIR)
	rm -rf html
	rm -fr latex

ifdef BUILD_SUFFIX
  -include .depend_$(BUILD_SUFFIX)
endif

