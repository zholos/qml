include ../mk/src.mk

build: .q_files

$(call copy_src,../debug)

CFLAGS += -O0 -Wno-parentheses
FLAGS += -ggdb --coverage
DEFINES += -DDEBUG -DCOVERAGE

.q_files: ../debug/.q_files
	cp ../debug/qml.q ../debug/test.q .
	touch $@

coverage: all
	rm -f -- *.gcda *.gcov
	-q test.q -cd
	$(foreach c, \
	    $(patsubst %.o,%.c,$(OBJS)) \
	    $(filter-out wrap.h,$(INCLUDES)),gcov -bc $(c) &&) :

clean_misc:
	rm -f qml.q test.q .q_files
	rm -f -- *.gcno *.gcda *.gcov
	rm -f callgrind.out.*
