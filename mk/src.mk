include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

VERSION := 0.5.3-$(shell date +%Y%m%d)

OBJS := const.o alloc.o util.o opt.o \
        libm.o cephes.o lapack.o conmin.o conmax.o nlopt.o
INCLUDES := alloc.h util.h opt.h wrap.h conmin.h conmax.h nlopt.h

CFLAGS += -std=gnu99 -Wall -Wextra -Wno-missing-field-initializers
DEFINES = -DQML_VERSION=$(VERSION) -DKXARCH=$(KXARCH) -DKXVER=$(KXVER)

all: build

build: qml.$(DLLEXT)

%.o: %.c $(INCLUDES)
	$(CC) $(FLAGS) \
	    $(CFLAGS) \
	    $(DEFINES) \
	    -I../include -c -o $@ $<

qml.$(DLLEXT): $(OBJS) qml.symlist qml.mapfile
	$(CC) $(FLAGS) -shared -o $@ $(OBJS) \
	    -L../lib -lprob -lconmax -lnlopt \
	    $(LDFLAGS) \
	    $(if $(BUILD_LAPACK),-llapack,$(LIBS_LAPACK)) \
	    $(if $(BUILD_BLAS),-lrefblas,$(LIBS_BLAS)) \
	    $(call ld_static,$(LIBS_FORTRAN)) \
	    -lm \
	    $(if $(WINDOWS),-lq) \
	    $(call ld_export,qml)

qml.symlist: $(OBJS)
	$(call nm_exports,$(OBJS)) | sed -n 's/^qml_/_&/p' >$@.tmp
	$(if $(BUILD_LAPACK),,echo _xerbla_ >>$@.tmp)
	mv $@.tmp $@

qml.mapfile: qml.symlist
	echo "{ global:"             >$@.tmp
	sed 's/^_/    /;s/$$/;/' $< >>$@.tmp
	echo "  local: *; };"       >>$@.tmp
	mv $@.tmp $@


ifneq ($(QHOME),)
# don't export QHOME because q below should get it from the general environment
install: build
	cp qml.$(DLLEXT) '$(QHOME)'/$(KXARCH)/
	cp qml.q         '$(QHOME)'/
	
uninstall:
	rm -f -- '$(QHOME)'/$(KXARCH)/qml.$(DLLEXT)
	rm -f -- '$(QHOME)'/qml.q
endif

test: build
	q test.q -cd -s 16 </dev/null

clean_src:
clean_misc:
clean: clean_src clean_misc
	rm -f qml.so qml.dll qml.symlist qml.mapfile
	rm -f $(OBJS)


SRC = $(patsubst %.o,%.c,$(OBJS)) $(INCLUDES)

define rule_cp
$(1): $(2)
	cp -- $(2) $(1)
endef

define rule_clean_src
clean_src:
	rm -f -- $$(SRC)
endef

define copy_src
$(foreach src,$(SRC),$(eval $(call rule_cp,$(src),$(addprefix $(1)/,$(src)))))\
$(eval $(rule_clean_src))
endef
