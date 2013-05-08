all:

force:

LD      = $(TOOLPREFIX)ld
AR      = $(TOOLPREFIX)ar
RANLIB  = $(TOOLPREFIX)ranlib
DLLTOOL = $(TOOLPREFIX)dlltool
NM      = $(TOOLPREFIX)nm
STRIP   = $(TOOLPREFIX)strip

# GNU and BSD sed differ in argument parsing. This way is portable.
SEDI = sed -i.tmp

# Superfluous leading underscore appears sometimes under Windows
sed_s+ = [[:space:]]\{1,\}
nm_exports = $(NM) $(1) | sed -n \
    's/^[[:alnum:]]*$(sed_s+)T$(sed_s+)_*\([[:alnum:]_]\{1,\}\)$$/\1/p'

ld_export/-exported_symbols_list = $(LD_EXPORT) $(1).symlist
ld_export/-Wl                    = $(LD_EXPORT),$(1).mapfile
ld_export = $(ld_export/$(patsubst -Wl%,-Wl,$(LD_EXPORT)))

ld_static/-Wl,-Bstatic = -Wl,-Bstatic $(1) -Wl,-Bdynamic
ld_static/ = $(1)
ld_static = $(ld_static/$(LD_STATIC))


fetch/curl   = curl -RLo $(1) $(2)
fetch/wget   = wget -O $(1) $(2)
fetch/manual = [ -e '$(1)' ] || \
    { echo; \
      echo 'please fetch $(1) from $(2)'; \
      echo 'or install curl or wget and run ./configure again'; \
      echo; ! :; }

fetch = $(fetch/$(FETCH))


sha256/sha256    = sha256
sha256/sha256sum = sha256sum | cut -b 1-64

sha256 = $(sha256/$(SHA256))


WINDOWS = $(filter w%,$(KXARCH))
DLLEXT = $(if $(WINDOWS),dll,so)

ifndef QML_CONFIGURE
ifeq ($(wildcard $(dir $(lastword $(MAKEFILE_LIST)))../config.mk),)
$(error run ./configure first)
endif
include $(dir $(lastword $(MAKEFILE_LIST)))../config.mk
endif
