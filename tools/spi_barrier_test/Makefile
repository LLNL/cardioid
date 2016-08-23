# Get name of architecture-specific include file from hostname
SHELL = /bin/bash

.PHONY: DEFAULT COMMON_TARGETS makedirs clean distclean depend

include Makefile.arch
# If defined, use architecture file set with ARCH variable
ARCH ?= $(ARCHGUESS)
include ../../trunk/src/arch/$(ARCH).mk

SVNVERSION := $(shell svnversion ./ )

#------------------------------------------------------------------------------
#

BINDIR = .

HEARTSRC= test.cc \
	 CommTable.cc 

SPISRC=	 spi_impl.c

DDCMD_FILES = 

DDCMDSRC = $(filter %.c, $(DDCMD_FILES))

HEARTOBJS = $(HEARTSRC:%.cc=$(OBJDIR)/%.o)
HEARTOBJS += $(DDCMDSRC:%.c=$(OBJDIR)/%.o)
HEARTOBJS += $(SPISRC:%.c=$(OBJDIR)/%.o)
ifdef TIMEOBJ
  HEARTOBJS += $(TIMEOBJ)
endif

#BUILD_SUFFIX := $(ARCH)
#OBJDIR := objs-$(BUILD_SUFFIX)
EXENAME = spi_test


DEFAULT: opt

COMMON_TARGETS: $(DDCMD_FILES) makedirs

opt:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/$(EXENAME)-$(ARCH) \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT)" CXXFLAGS="$(CXXFLAGS_OPT)" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=$(EXENAME)-$(ARCH)
debug:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/$(EXENAME)-$(ARCH)-debug \
	BUILD_SUFFIX=$(ARCH)-debug OBJDIR=objs-$(ARCH)-debug \
	CFLAGS="$(CFLAGS_DEBUG)" CXXFLAGS="$(CXXFLAGS_DEBUG)" \
	LDFLAGS="$(LDFLAGS_DEBUG)" EXENAME=$(EXENAME)-$(ARCH)-debug 
profile:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/$(EXENAME)-$(ARCH)-prof \
	BUILD_SUFFIX=$(ARCH)-prof OBJDIR=objs-$(ARCH)-prof \
	CFLAGS="$(CFLAGS_PROF)" CXXFLAGS="$(CXXFLAGS_PROF)" \
	LDFLAGS="$(LDFLAGS_PROF)"  EXENAME=$(EXENAME)-$(ARCH)-prof
mpip:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/$(EXENAME)-$(ARCH)-mpip \
	BUILD_SUFFIX=$(ARCH)-mpip OBJDIR=objs-$(ARCH)-mpip \
	CFLAGS="$(CFLAGS_MPIP)" CXXFLAGS="$(CXXFLAGS_MPIP)" \
	LDFLAGS="$(LDFLAGS_MPIP)"  EXENAME=$(EXENAME)-$(ARCH)-mpip
time:
	@$(MAKE) --no-print-directory COMMON_TARGETS timebase \
	$(BINDIR)/$(EXENAME)-$(ARCH)-time \
	BUILD_SUFFIX=$(ARCH) OBJDIR=objs-$(ARCH) \
	CFLAGS="$(CFLAGS_OPT) -DTIMING" CXXFLAGS="$(CXXFLAGS_OPT) -DTIMING" \
	LDFLAGS="$(LDFLAGS_OPT)"  EXENAME=$(EXENAME)-$(ARCH)-time \
	TIMEOBJ=objs-$(ARCH)/timebase.o

timebase: timebase.c
	/bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gcc -c -O timebase.c -o $(OBJDIR)/timebase.o

$(OBJDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(BINDIR)/$(EXENAME): $(HEARTOBJS)
	$(LD) $(DFLAGS) -o $(BINDIR)/$(EXENAME) $(HEARTOBJS) $(LDFLAGS) \
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

testLoadBalancer-target: $(HEARTOBJS) $(OBJDIR)/testLoadBalancer.o
	$(LD) $(DFLAGS) -o $(BINDIR)/$(EXENAME) $(OBJDIR)/testLoadBalancer.o $(HEARTOBJS) $(LDFLAGS)


testGridRouter: makedirs $(HEARTOBJS) $(OBJDIR)/testGridRouter.o
	$(LD) $(DFLAGS) -o testGridRouter $(OBJDIR)/testGridRouter.o $(HEARTOBJS) $(LDFLAGS)

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

.depend_$(BUILD_SUFFIX): $(HEARTSRC)
	@touch .depend_$(BUILD_SUFFIX)
	@$(MAKE) depend BUILD_SUFFIX=$(BUILD_SUFFIX)

#------------------------------------------------------------------------------
#  Cleanup object files
clean :
	rm -rf objs-$(ARCH)/*.o
	rm -rf $(BINDIR)/$(EXENAME)-$(ARCH)*.*
	rm -f $(BINDIR)/$(EXENAME)-$(ARCH)*
	rm -f .depend_$(ARCH)*
	rm -f $(BINDIR)/testLoadBalancer-$(ARCH)
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
