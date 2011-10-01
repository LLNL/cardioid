# Get name of architecture-specific include file from hostname
include Makefile.arch

# If defined, use architecture file set with ARCH variable
ARCH ?= $(ARCHGUESS)
include $(ARCH).mk
#------------------------------------------------------------------------------
#

HEARTSRC=TT04Model.C SimpleInputParser.C ProcessGrid3D.C DataGrid3D.C \
	 GDLoadBalancer.C CommTable.C

HEARTOBJS = $(HEARTSRC:%.C=$(OBJDIR)/%.o)

OBJDIR = objs-$(ARCH)
INSTALLDIR = ../bin
EXENAME = bigHeart-$(ARCH)
CXXFLAGS += -DARCH='"$(ARCH)"'

DEFAULT: $(EXENAME)

$(OBJDIR)/%.o: %.C
	$(CXX) $(CXXFLAGS) -o $@ -c $^

$(OBJDIR)/bigHeart.o: bigHeart.C
	$(CXX) $(CXXFLAGS) -o $@ -c $^

$(EXENAME): makedirs depend $(HEARTOBJS) $(OBJDIR)/bigHeart.o
	$(LD) $(DFLAGS) -o $(EXENAME) $(OBJDIR)/bigHeart.o $(HEARTOBJS) $(LDFLAGS)

testLoadBalancer: makedirs $(OBJDIR)/GDLoadBalancer.o $(OBJDIR)/testLoadBalancer.o
	$(LD) $(DFLAGS) -o testLoadBalancer $(OBJDIR)/testLoadBalancer.o $(HEARTOBJS) $(LDFLAGS)

makedirs:
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR) ; fi
	@if [ ! -e .depend_$(ARCH) ]; then touch .depend_$(ARCH) ; fi

#------------------------------------------------------------------------------
# generate dependencies in makefile: use -Y to avoid library header files
# that are likely to be different on other platforms.
 depend :
	@makedepend -f .depend_$(ARCH) -p'($(OBJDIR)/' -o '.o)' -Y -D$(PLT) *.[cCh] 2> /dev/null

.depend_$(ARCH): $(HEARTSRC)
	@touch .depend_$(ARCH)
	@$(MAKE) depend

#	makedepend -Y -D$(PLT) *.[cCh]
#makedepend -D$(PLT) -I/usr/include/CC *.[cC]
#------------------------------------------------------------------------------
#  Cleanup object files
 clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(EXENAME)
	rmdir $(OBJDIR)
	rm -f .depend_$(ARCH)
	rm -f testLoadBalancer

-include .depend_$(ARCH)
