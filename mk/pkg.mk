include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

all: build

NETLIB_DIR = http://www.netlib.org
Q_DIR = https://raw.githubusercontent.com/KxSystems/kdb/f4fb81b7dc7448c73d94ae9dff0610fd71c7fe93
GITHUB_DIR = https://github.com/$(1)/$(2)/archive/$(3)

.PRECIOUS: ../download/%.part
../download/%.part:
	mkdir -p $(dir $@)
	$(call fetch,$@,$($*.url))

../download/%: | ../download/%.part
	if $(if $($*.sha256),[ $($*.sha256) = `<$| $(sha256)` ],:); then \
	    mv -- $| $@; else echo; \
	    echo 'checksum mismatch'; \
	    echo 'try deleting $| so it can be downloaded again'; \
	    echo; ! :; fi


extract: .extracted
.extracted:
	$(MAKE) do_extract
	touch $@
do_extract:

patch: .patched
.patched: .extracted
	$(MAKE) do_patch
	touch $@
do_patch:

configure: .configured
.configured: .patched
	$(MAKE) do_configure
	touch $@
do_configure:

build: .built
.built: .configured
	$(MAKE) do_build
	touch $@
do_build:

install: .installed
.installed: .built
	mkdir -p -- ../include ../lib
	$(MAKE) do_install
	touch $@
do_install:

clean: do_clean
	rm -f .extracted .patched .configured .built .installed
do_clean:

.PHONY: $(foreach t,extract patch configure build install clean,$(t) do_$(t))
