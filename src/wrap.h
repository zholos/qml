#ifndef QML_SRC_WRAP_H
#define QML_SRC_WRAP_H

#include "util.h"

// x is scalar or vector, call is an expression of v, result is returned
#define wrap_call(x, call)                 \
    do {                                   \
        if (likely((x)->t == -KF)) {       \
            F v = (x)->f;                  \
            return kf((call));             \
        }                                  \
        wrap_call_vector(x, call);         \
        if (compatible_f((x))) {           \
            F v = convert_f((x));          \
            return kf((call));             \
        }                                  \
        return krr("type");                \
    } while (0)

// vector case separated to avoid dead code when atom case is already handled
#define wrap_call_vector(x, call)          \
    do {                                   \
        K wc__x = convert_F((x));          \
        if (wc__x) {                       \
            K wc__r = ktn(KF, wc__x->n);   \
            repeat (wc__i, wc__x->n) {     \
                F v = kF(wc__x)[wc__i];    \
                kF(wc__r)[wc__i] = (call); \
            }                              \
            q0(wc__x);                     \
            return wc__r;                  \
        }                                  \
    } while (0)


// (atom/vector F)
#define wrap_F(name, func) \
K                          \
qml_##name(K x) {          \
    wrap_call(x, func(v)); \
}

#define wrap_F_(func) wrap_F(func, func)


// (atom/vector F, atom/vector F)
#define wrap_FF(name, func)                 \
K                                           \
qml_##name(K x, K y) {                      \
    if (likely(xt == -KF && y->t == -KF))   \
        return kf(func(xf, y->f));          \
    if (compatible_f(x)) {                  \
        F u = convert_f(x);                 \
        wrap_call(y, func(u, v));           \
    }                                       \
    if (compatible_f(y)) {                  \
        F u = convert_f(y);                 \
        assert(!compatible_f(x));           \
        wrap_call_vector(x, func(v, u));    \
    }                                       \
    check_type(x = convert_F(x),);          \
    check_type(y = convert_F(y), q0(x));    \
    check_length(xn == y->n, q0(y); q0(x)); \
    K r = ktn(KF, xn);                      \
    repeat (i, xn)                          \
        kF(r)[i] = func(xF[i], kF(y)[i]);   \
    q0(y); q0(x);                           \
    return r;                               \
}

#define wrap_FF_(func) wrap_FF(func, func)


// (atom I, atom/vector F)
#define wrap_iF(name, func)       \
K                                 \
qml_##name(K x, K y) {            \
    check_type(compatible_i(x),); \
    I u = convert_i(x);           \
    wrap_call(y, func(u, v));     \
}

// (atom F, atom/vector F)
#define wrap_fF(name, func)       \
K                                 \
qml_##name(K x, K y) {            \
    check_type(compatible_f(x),); \
    F u = convert_f(x);           \
    wrap_call(y, func(u, v));     \
}

// (atom F, atom F, atom/vector F)
#define wrap_iiF(name, func)                         \
K                                                    \
qml_##name(K x, K y, K z) {                          \
    check_type(compatible_i(x) && compatible_i(y),); \
    I t = convert_f(x), u = convert_f(y);            \
    wrap_call(z, func(t, u, v));                     \
}

// (atom F, atom F, atom/vector F)
#define wrap_ffF(name, func)                         \
K                                                    \
qml_##name(K x, K y, K z) {                          \
    check_type(compatible_f(x) && compatible_f(y),); \
    F t = convert_f(x), u = convert_f(y);            \
    wrap_call(z, func(t, u, v));                     \
}


#endif // QML_SRC_WRAP_H
