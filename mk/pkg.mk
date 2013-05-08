include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

all: build

NETLIB_DIR = http://www.netlib.org
Q_DIR = http://kx.com/q

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

patch: extract .patched
.patched:
	$(MAKE) do_patch
	touch $@
do_patch:

configure: patch .configured
.configured:
	$(MAKE) do_configure
	touch $@
do_configure:

build: configure .built
.built:
	$(MAKE) do_build
	touch $@
do_build:

install: build .installed
.installed:
	mkdir -p -- ../include ../lib
	$(MAKE) do_install
	touch $@
do_install:

clean: do_clean
	rm -f .extracted .patched .configured .built .installed
do_clean:
