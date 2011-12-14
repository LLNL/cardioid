# Get name of architecture-specific include file from hostname
include Makefile.arch

.PHONY: DEFAULT makedirs clean distclean depend

# If defined, use architecture file set with ARCH variable
ARCH ?= $(ARCHGUESS)
include $(ARCH).mk
#------------------------------------------------------------------------------
#

HEARTSRC=cardioid.cc \
	 Vector.cc \
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
	 Saleheen98Diffusion.cc \
	 ReactionFHN.cc \
	 TT04_bbReaction.cc TT04_CellML_Reaction.cc \
	 TT06_CellML_Reaction.cc \
	 IBM_tenTusscher04.cc IBM_tenTusscher04_endoLUT.cc \
	 IBM_tenTusscher04_epiLUT.cc IBM_tenTusscher04_midLUT.cc \
	 TT04_CellML_Endo.cc TT04_CellML_Mid.cc TT04_CellML_Epi.cc \
	 TT06_CellML_Endo.cc TT06_CellML_Mid.cc TT06_CellML_Epi.cc \
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

OBJDIR = objs-$(ARCH)
INSTALLDIR = ../bin
EXENAME = cardioid-$(ARCH)
CXXFLAGS += -DARCH='"$(ARCH)"'

DEFAULT: $(DDCMD_FILES) makedirs $(EXENAME)

$(OBJDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(EXENAME): $(HEARTOBJS)
	$(LD) $(DFLAGS) -o $(EXENAME) $(HEARTOBJS) $(LDFLAGS)

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

#------------------------------------------------------------------------------
# generate dependencies in makefile: use -Y to avoid library header files
# that are likely to be different on other platforms.
 depend :
	@makedepend -f .depend_$(ARCH) -p'$(OBJDIR)/' -o '.o' -Y -D$(PLT) *.cc *.c *.hh *.h 2> /dev/null

.depend_$(ARCH): $(HEARTSRC)
	@touch .depend_$(ARCH)
	@$(MAKE) depend

#	makedepend -Y -D$(PLT) *.[cCh]
#makedepend -D$(PLT) -I/usr/include/CC *.[cC]
#------------------------------------------------------------------------------
#  Cleanup object files
clean :
	rm -rf $(OBJDIR)/*.o
	rm -f $(EXENAME)
	rm -f .depend_$(ARCH)*
	rm -f testLoadBalancer
	rm -f testGridRouter

distclean: clean
	rm -f $(DDCMD_FILES)
	rm -rf ddcMD_files
	rm -rf $(OBJDIR)
	rm -f .depend*

-include .depend_$(ARCH)
