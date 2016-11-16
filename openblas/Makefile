include ../mk/pkg.mk

openblas-0.2.19.tgz.url = $(call GITHUB_DIR,xianyi,OpenBLAS,v0.2.19.tar.gz)
openblas-0.2.19.tgz.sha256 = \
    9c40b5e4970f27c5f6911cb0a28aa26b6c83f17418b69f8e5a116bb983ca8557

ifneq ($(BUILD_OPENBLAS),)

do_extract: ../download/openblas-0.2.19.tgz
	tar xzf $<
	ln -sf OpenBLAS-0.2.19 work

# For static architecture we can skip single-precision and complex routines to
# reduce build time. For dynamic architecture all the routines are in a single
# table which would be difficult to untangle. See kernel/setparam_*.c.

do_patch:
	patch -p0 -d work <buffers.patch
	patch -p0 -d work <make.patch
ifneq ($(BUILD_OPENBLAS), dynamic)
	patch -p0 -d work <only_double.patch
	$(SEDI) -e '/^[[:space:]].*\([scz]blat\|[SCZ]BLAT\)/d' \
	    -e '/^[scz]blat/d;/^level/s/[scz]blat[^[:space:]]*//g' \
	    work/test/Makefile
endif

# OpenBLAS has a limit on the number of fixed-size buffers it can allocate
# internally. Both its own threads and q -s threads use up extra buffers and
# OpenBLAS crashes if the limit is exceeded. Therefore, set it to a high number
# with NUM_THREADS. Each buffer is 16-32 MB and they are never freed, but the
# total number allocated should not exceed peak usage plus number of cores. See
# driver/others/memory.c.

do_configure:
	echo "# Settings from qml:"             >>work/Makefile.rule
ifeq ($(BUILD_OPENBLAS), dynamic)
	echo DYNAMIC_ARCH=1                     >>work/Makefile.rule
else
	echo DYNAMIC_ARCH=0                     >>work/Makefile.rule
ifneq ($(BUILD_OPENBLAS), native)
	echo TARGET=$(BUILD_OPENBLAS) | \
	    tr '[:lower:]' '[:upper:]'          >>work/Makefile.rule
endif
endif
	echo CC='$(CC)'                         >>work/Makefile.rule
	echo FC='$(FC)'                         >>work/Makefile.rule
	echo HOSTCC='$(XCC)'                    >>work/Makefile.rule
	echo BINARY=$(if $(KXARCH:%32=),64,32)  >>work/Makefile.rule
	echo USE_THREAD=0                       >>work/Makefile.rule
	echo NUM_THREADS=64                     >>work/Makefile.rule
	echo NO_SHARED=1                        >>work/Makefile.rule
	echo NO_CBLAS=1                         >>work/Makefile.rule
	echo COMMON_OPT='$(FLAGS)'              >>work/Makefile.rule
	echo FCOMMON_OPT='$(FFLAGS)'            >>work/Makefile.rule
	echo CCOMMON_OPT='$(CFLAGS)'            >>work/Makefile.rule
	echo QML_TOOLPREFIX='$(TOOLPREFIX)'     >>work/Makefile.rule
	
# LDFLAGS is used for tests; -static helps when MinGW DLLs are not in PATH
	echo LDFLAGS='$(LDFLAGS) $(if $(WINDOWS),-static)' >>work/Makefile.rule

do_build:
	$(MAKE) -C work

do_install:
	cp work/libopenblas.a ../lib/

do_clean:
	rm -rf OpenBLAS-[0-9]* work

endif # BUILD_OPENBLAS
