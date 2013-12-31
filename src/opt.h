#ifndef QML_SRC_OPT_H
#define QML_SRC_OPT_H

#include <k.h>

struct optn {
    const char* s;
    H t; // 0 for boolean, -KI, or -KF
};

union optv {
    I i; /* first member */
    F f;
};

// Fills in an optv array
// 0 on failure
int take_opt(K x, const struct optn* n, union optv* v);


#endif // QML_SRC_OPT_H
