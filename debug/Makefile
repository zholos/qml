include ../mk/src.mk

build: .q_files

$(call copy_src,../src)

CFLAGS += -O0 -Wno-parentheses
FLAGS += -ggdb -ftrapv
DEFINES += -DDEBUG

.q_files: debug.patch ../src/qml.q ../src/test.q
# patch is sensitive to line endings, so convert them
	sed 's/[[:space:]]$$//' ../src/qml.q  >qml.q
	sed 's/[[:space:]]$$//' ../src/test.q >test.q
	patch -p0 <$<
	touch $@

valgrind: all
	valgrind --suppressions=q.supp --leak-check=full --show-reachable=yes \
	    --track-origins=yes --trace-children=yes q test.q -cd </dev/null

clean_misc:
	rm -f qml.q qml.q.orig test.q test.q.orig .q_files
	rm -f callgrind.out.*
