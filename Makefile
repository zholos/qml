include mk/common.mk

PACKAGES := cephes lapack conmax nlopt q
BUILDS := src debug
OTHERS := coverage

$(foreach dir,$(PACKAGES) $(BUILDS) $(OTHERS),\
    $(foreach target,build install uninstall test clean,$(eval \
        $(dir)_$(target): ; make -C $(dir) $(target))))

all: build

build: $(addsuffix _install,$(PACKAGES)) src_build
	@echo build complete.

install: build src_install
	@echo installed.

uninstall: src_uninstall
	@echo removed.

test: build $(addsuffix _test,$(BUILDS))
	@echo all tests passed.

clean: $(addsuffix _clean,$(PACKAGES) $(BUILDS) $(OTHERS))
	rm -rf -- include/ lib/
	$(MAKE) -f mk/test.mk clean
