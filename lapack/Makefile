include ../mk/pkg.mk

lapack-3.6.0.tgz.url = $(NETLIB_DIR)/lapack/lapack-3.6.0.tgz
lapack-3.6.0.tgz.sha256 = \
    a9a0082c918fe14e377bbd570057616768dca76cbdc713457d8199aaa233ffc3

# BUILD_BLAS sets BUILD_LAPACK too
ifneq ($(BUILD_LAPACK),)

do_extract: ../download/lapack-3.6.0.tgz
	tar xzf $<
	ln -sf lapack-3.6.0 work


do_patch: make.patch
	cat $^ | patch -p1 -d work
	find work/TESTING -name Makefile -exec $(SEDI) \
	    '/^include.*\/make\.inc/{p;s/.*/OPTS=$$(OPTS_TESTING)/;}' {} +


do_configure:
	cp work/make.inc.example work/make.inc
	echo FORTRAN='$(FC)'                           >>work/make.inc
	echo OPTS='$(FLAGS) $(FFLAGS)'                 >>work/make.inc
	echo CFLAGS='$(FLAGS) $(CFLAGS)'               >>work/make.inc
	echo CC='$(CC)'                                >>work/make.inc
	echo ARCH='$(AR)'                              >>work/make.inc
	echo RANLIB='$(RANLIB)'                        >>work/make.inc
	
ifeq ($(BUILD_BLAS),)
	echo BLASLIB='$(LIBS_BLAS)'                    >>work/make.inc
endif
	
# Test routines use huge arrays that overflow the stack with -frecursive.
# Alternative fix is to increase stack size.
	echo OPTS_TESTING='$(FLAGS) $(FFLAGS_NTHREAD)' >>work/make.inc
	
ifneq ($(WINDOWS),)
# MinGW libgfortran.dll is probably not in PATH; this simplifies things
	echo LOADOPTS+=-static                         >>work/make.inc
endif


do_build:
ifneq ($(BUILD_BLAS),)
	$(MAKE) -C work blaslib
	$(MAKE) -C work blas_testing
endif
	$(MAKE) -C work lapacklib
	$(MAKE) -C work lapack_testing


install_libs:
ifneq ($(BUILD_BLAS),)
	cp work/librefblas.a ../lib/
endif
	cp work/liblapack.a ../lib/

endif # BUILD_LAPACK


install_libs:
do_install: install_libs
	cp lapack.h ../include/

do_clean:
	rm -rf lapack-[0-9]* work
