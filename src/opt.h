#ifndef OPT_H
#define OPT_H

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


#endif // OPT_H
