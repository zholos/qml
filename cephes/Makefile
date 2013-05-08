include ../mk/pkg.mk

cprob.tgz.url = $(NETLIB_DIR)/cephes/cprob.tgz
cprob.tgz.sha256 = \
    fae184467171ac96a1d7c6378618dcedeb2c8f04ff1237de8f77c16986c9111c

do_extract: ../download/cprob.tgz
	mkdir -p work
	tar xzf $< -C work

do_patch: cprob.patch
	patch -p1 -d work <$<

# Variables set taking implicit rules into account
do_configure:
	rm -f work/make.inc
	echo CC='$(CC)'                  >>work/make.inc
	echo CPPFLAGS=                   >>work/make.inc
	echo CFLAGS='$(FLAGS) $(CFLAGS)' >>work/make.inc
	echo AR='$(AR)'                  >>work/make.inc
	echo RANLIB='$(RANLIB)'          >>work/make.inc

do_build:
	$(MAKE) -C work -f cprob.mak

do_install:
	cp cprob.h ../include/
	cp work/libprob.a ../lib/

do_clean:
	rm -rf work
