include ../mk/pkg.mk

conmax.f.url = $(NETLIB_DIR)/opt/conmax.f
conmax.f.sha256 = \
    253ed0e323892896cdf066e7cb5c42b41ff49dbf7aa3283875998dd260753738

do_extract: ../download/conmax.f
	mkdir -p work
	cp $< work/

do_patch: conmax.sed
	$(SEDI) -f conmax.sed work/conmax.f

do_build:
	cd work && $(FC) $(FLAGS) $(FFLAGS) -c conmax.f
	cd work && $(AR) r libconmax.a conmax.o
	cd work && $(RANLIB) libconmax.a

do_install:
	cp conmax.h ../include/
	cp work/libconmax.a ../lib/

do_clean:
	rm -rf work
