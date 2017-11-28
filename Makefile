
ARCHES ?= $(shell $(MAKE) --no-print-directory -C src getArch)
TESTS ?= seq
PYTHON ?= python

test: $(foreach ARCH,$(ARCHES),$(foreach TEST,$(TESTS),test-$(ARCH)-$(TEST)))
build: $(foreach ARCH,$(ARCHES),build-$(ARCH))
clean: $(foreach ARCH,$(ARCHES),$(foreach TEST,$(TESTS),clean-$(ARCH)-$(TEST)))

build-%:
	$(MAKE) -C src ARCH=$(patsubst build-%,%,$@) opt singleCell
cleanbuild-%:
	$(MAKE) -C src ARCH=$(patsubst cleanbuild-%,%,$@) clean

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

ALL_ARCH=$(patsubst src/arch/%.mk,%,$(wildcard src/arch/*.mk))
ALL_PROFILE=$(patsubst test/profiles/%,%,$(wildcard test/profiles/*))

$(foreach arch,$(ALL_ARCH),$(foreach profile,$(ALL_PROFILE),$(eval $(call TEST_RULE,$(arch),$(profile)))))


