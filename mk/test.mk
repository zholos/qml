include $(dir $(lastword $(MAKEFILE_LIST)))common.mk


define conftest.c/
int main() { return 0; }
endef

define conftest.f/
       program main
       end
endef

define conftest.c/shared
int x;
int f() { return x; }
int g() { return x; }
endef

define conftest.c/blas
int dgemm_();
int main() {
    dgemm_();
    return 0;
}
endef

define conftest.c/lapack
int dgelsd_();
int main() {
    dgelsd_();
    return 0;
}
endef

define conftest.f/builtin
       subroutine f(a, b, c, i, j, k)
       character a, b, c
       integer   i, j, k
       c = a // b
       k = i ** j
       end
endef

conftest.c/builtin = $(conftest.c/)


export conftest_c = $(conftest.c/$(CONFTEST))
export conftest_f = $(conftest.f/$(CONFTEST))

conftest.c: force
	echo "$$conftest_c" >$@

conftest.f: force
	echo "$$conftest_f" >$@


LINK_FLAGS = $(LDFLAGS) $(LIBS_LAPACK) $(LIBS_BLAS) $(LIBS_FORTRAN) -lm

test/c_version:
	$(CC) --version

test/c_compile: conftest.c
	$(CC) -Wimplicit -Werror -c -o conftest.o $< \
	    $(FLAGS) $(CFLAGS)

test/c_link: conftest.c
	$(CC) -Wimplicit -Werror -o conftest.exe $< \
	    $(FLAGS) $(CFLAGS) $(LINK_FLAGS)


test/xc_version:
	$(XCC) --version

test/xc_compile: conftest.c
	$(XCC) -c -o conftest.o $<

test/xc_link: conftest.c
	$(XCC) -o conftest.exe $<

test/xc_run: test/xc_link
	./conftest.exe


test/f_compile: conftest.f
	$(FC) -Werror -c -o conftest.o $< \
	    $(FLAGS) $(FFLAGS)

test/f_link: conftest.f
	$(FC) -Werror -o conftest.exe $< \
	    $(FLAGS) $(FFLAGS) $(LINK_FLAGS)

test/c_and_f_link: test/c_link test/f_link

test/f_compile_c_link: conftest.f conftest.c
	$(FC) -Werror -c -o conftest.o $< \
	    $(FLAGS) $(FFLAGS)
	$(CC) -Werror -o conftest.exe conftest.o conftest.c \
	    $(FLAGS) $(CFLAGS) $(LINK_FLAGS)


test/shared_link: conftest.c
	echo _f                          >conftest.symlist
	echo "{ global: f; local: *; };" >conftest.mapfile
	$(CC) -Werror -shared -o conftest.$(DLLEXT) $< \
	    $(FLAGS) $(CFLAGS) $(LINK_FLAGS) $(call ld_export,conftest)

test/ld_static: conftest.f conftest.c
	$(FC) -Werror -c -o conftest.o $< \
	    $(FLAGS) $(FFLAGS)
	$(CC) -Werror -shared -o conftest.$(DLLEXT) conftest.o conftest.c \
	    $(FLAGS) $(CFLAGS) $(call ld_static,$(LINK_FLAGS))

test/ld_export: test/shared_link
ifeq ($(WINDOWS),)
	[ f = "$$($(call nm_exports,-D conftest.$(DLLEXT)) | grep '^[fg]$$')" ]
else
# nm doesn't work this way on Windows, so just assume the mapfile worked
endif

test/no_cygwin: test/shared_link
	ldd conftest.$(DLLEXT) | \
	    { ! grep -q '[[:space:]]cygwin[0-9]*[.]dll[[:space:]]'; }


fetch-test/curl = curl --version
fetch-test/wget = wget --version

test/fetch:
	$(fetch-test/$(FETCH))

test/sha256:
	[ ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad = \
	    "$$(printf abc | 2>/dev/null $(sha256))" ]

test/patch:
	patch </dev/null


clean:
	rm -f conftest.*
