#ifndef QML_SRC_ALLOC_H
#define QML_SRC_ALLOC_H

#include "util.h"


// On failure sets n=0, err="wsfull" and returns NULL.
// m == n or 0 < m < wi
I* alloc_I(I* n, S* err);
F* alloc_F(I* n, S* err);
I* alloc_II(I* n, I m, S* err);
F* alloc_FF(I* n, I m, S* err);

void free_I(I* p);
void free_F(F* p);


static inline I
add_size(I a, I m, I n) {
    if (!likely(!n || m <= (wi - a) / n))
        return wi;
    return a + m * n;
}


#endif // QML_SRC_ALLOC_H
