#ifndef QML_SRC_OPT_H
#define QML_SRC_OPT_H

#include <k.h>

// options argument can be in one of these formats:
// - symbol atom or vector for flags, e.g. `quiet`full
// - dictionary, e.g. `quiet`iter!1b,42
// - getopt-style list, e.g. `quiet`iter,42

// null-terminated array of optn defines available options
struct optn {
    const char* s;
    H t; // 0 for flag, -KI, or -KF
};

// parsed option values are stored in corresponding array of optv
// must be initialized to default values before calling take_opt()
// note: flag stored in i and can have values other than 0 and 1 (e.g.
// .qml.minx[`quiet`iter!42 10000]); passing it through bool function argument
// doesn't reliably constrain its values as bool is sometimes defined as int
union optv {
    I i; /* first member */
    F f;
};

// returns 0 on failure
int take_opt(K x, const struct optn* n, union optv* v);


#endif // QML_SRC_OPT_H
