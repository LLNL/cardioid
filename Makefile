include Makefile.arch
ARCHES ?= $(ARCHGUESS)
TESTS ?= seq
PYTHON ?= python

test: $(foreach ARCH,$(ARCHES),$(foreach TEST,$(TESTS),test-$(ARCH)-$(TEST)))
build: $(foreach ARCH,$(ARCHES),build-$(ARCH))
clean: $(foreach ARCH,$(ARCHES),$(foreach TEST,$(TESTS),clean-$(ARCH)-$(TEST)))

mkdir-%:
	@mkdir -p $(patsubst mkdir-%,build/%,$@)
.PRECIOUS: build/%/CMakeCache.txt
build/%/CMakeCache.txt: mkdir-%
	(cd $(patsubst mkdir-%,build/%,$<) && cmake -C $(patsubst mkdir-%,../../arch/%.txt,$<) ../..)
build-%: build/%/CMakeCache.txt
	$(MAKE) --no-print-directory -C build/$(patsubst build-%,%,$@)
cleanbuild-%:
	$(MAKE) --no-print-directory -C build/$(patsubst cleanbuild-%,%,$@) clean

define TEST_RULE
.PRECIOUS: test/scratch/$1/$2/results.tap
test/scratch/$1/$2/results.tap : build-$1
	mkdir -p test/scratch/$1/$2
	$$(PYTHON) test/runner.py -b $1 -p $2 >| $$@ 2>&1
test-$1-$2 : test/scratch/$1/$2/results.tap
	@cat $$<
clean-$1-$2: cleanbuild-$1
	rm -f test/scratch/$1/$2/results.tap
endef

ALL_ARCH=$(sort $(patsubst arch/%.txt,%,$(wildcard arch/*.txt)) $(ARCHES))
ALL_PROFILE=$(patsubst test/profiles/%,%,$(wildcard test/profiles/*))

$(foreach arch,$(ALL_ARCH),$(foreach profile,$(ALL_PROFILE),$(eval $(call TEST_RULE,$(arch),$(profile)))))

