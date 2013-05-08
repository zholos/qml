#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "alloc.h"


#ifdef DEBUG

#define alloc_max 32
static int alloc_i = 1, alloc_fail[alloc_max] = { 1 };

K
qml_debug_alloc(K x) {
    (void)x;
    int i;
    for (i = alloc_max - 1; i >= 0; i--)
        if (i >= alloc_i)
            alloc_fail[i] = 1;
        else if (!(alloc_fail[i] ^= 1))
            break;

    alloc_i = 1;
    return kb(i >= 0);
}

void*
alloc(size_t n) {
    if (!n)
        return malloc(n);

    if (!always(alloc_i < alloc_max))
        abort();

    // Only active inside while[.qml.debug_alloc[]] loop, where !alloc_fail[0].
    if (!alloc_fail[0] && alloc_fail[alloc_i++])
        return NULL;

    void* r = malloc(n);
    if (!always(r))
        abort(); // only fail for scheduled reasons
    return r;
}

K
qml_debug_coverage(K x) {
    (void)x;
    free(alloc(0));
    return NULL;
}

#else
#define alloc malloc
#endif



// malloc may or may not return null on n = 0. Although it's likely that there
// is another error pending, don't try to emit wsfull in this case.
#define alloc_TT(T)                                           \
T*                                                            \
alloc_##T##T(I* n, I m, S* err) {                             \
    if (!likely(*n && m))                                     \
        return NULL;                                          \
                                                              \
    assert(*n >= 0 && m > 0);                                 \
    if (!likely(*n < wi && *n <= SIZE_MAX / sizeof(T) / m)) { \
        *n = 0;                                               \
        if (!*err) *err = "limit";                            \
        return NULL;                                          \
    }                                                         \
                                                              \
    T* r = alloc(*n * m * sizeof(T));                         \
    if (!likely(r)) {                                         \
        *n = 0;                                               \
        if (!*err) *err = "wsfull";                           \
    }                                                         \
    return r;                                                 \
}

alloc_TT(I)
alloc_TT(F)


#define alloc_T(T)                  \
T*                                  \
alloc_##T(I* n, S* err) {           \
    return alloc_##T##T(n, 1, err); \
}

alloc_T(I)
alloc_T(F)


void free_I(I* p) { free(p); }
void free_F(F* p) { free(p); }
