#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "util.h"


int
compatible_i(K x) {
    return xt >= -KJ && xt <= -KH || xt == -KB;
}


// Must handle default j type in kdb+ 3.0.
// This is mostly for d.f., so 0Wi and 0Wj will give almost the same results.
static inline I
i_from_j(J x) {
    return x == nj ? ni : x > wi ? wi : x < -wi ? -wi : x;
}

static inline I
i_from_h(H x) {
    return x == nh ? ni : x;
}

I
convert_i(K x) {
    switch (xt) {
    case -KJ:
        return i_from_j(xj);
    case -KI:
        return xi;
    case -KH:
        return i_from_h(xh);
    case -KB:
        return xg;
    }
    return ni;
}



int
compatible_f(K x) {
    return xt >= -KF && xt <= -KH || xt == -KB;
}


static inline F
f_from_j(J x) {
    return x == nj ? nf : x == wj ? wf : x == -wj ? -wf : x;
}

static inline F
f_from_i(I x) {
    return x == ni ? nf : x == wi ? wf : x == -wi ? -wf : x;
}

static inline F
f_from_h(H x) {
    return x == nh ? nf : x == wh ? wf : x == -wh ? -wf : x;
}

F
convert_f(K x) {
    switch (xt) {
    case -KF:
        return xf;
    case -KE:
        return xe;
    case -KJ:
        return f_from_j(xj);
    case -KI:
        return f_from_i(xi);
    case -KH:
        return f_from_h(xh);
    case -KB:
        return xg;
    }
    return nf;
}


// returns KF; or NULL
K
convert_F(K x) {
    if (likely(qt(x) == KF))
        // Not certain that x == r1(x), assume in most cases it is
        return r1(x);

    K r;
    switch (qt(x)) {
    case KE:
        r = ktn(KF, qn(x));
        repeat (i, qn(x))
            qF(r, i) = qE(x, i);
        break;
    case KJ:
        r = ktn(KF, qn(x));
        repeat (i, qn(x))
            qF(r, i) = f_from_j(qJ(x, i));
        break;
    case KI:
        r = ktn(KF, qn(x));
        repeat (i, qn(x))
            qF(r, i) = f_from_i(qI(x, i));
        break;
    case KH:
        r = ktn(KF, qn(x));
        repeat (i, qn(x))
            qF(r, i) = f_from_h(qH(x, i));
        break;
    case KB:
        r = ktn(KF, qn(x));
        repeat (i, qn(x))
            qF(r, i) = qB(x, i);
        break;
    case 0:
        r = ktn(KF, qn(x));
        repeat (i, qn(x))
            if (likely(qt(qK(x, i)) == -KF))
                qF(r, i) = qf(qK(x, i));
            else if (compatible_f(qK(x, i)))
                qF(r, i) = convert_f(qK(x, i));
            else
                return q0(r);
        break;
    default:
        return NULL;
    }
    return r;
}


// returns list of KF; or NULL
K
convert_FF(K x) {
    if (qt(x) != 0)
        return NULL;

    repeat (j, qn(x))
        if (!likely(qt(qK(x, j)) == KF))
            goto copy;
    return r1(x);

copy:;
    K r = ktn(0, xn);
    int failed = 0;
    repeat (j, xn)
        failed |= !likely(qK(r, j) = convert_F(qK(x, j)));
    if (failed)
        return q0(r);
    return r;
}


// like "f"$x; or NULL
// typically used for mixed lists of KF and -KF (complex numbers) or
// relatively small multidimensional lists of parameters
K
convert_FFF(K x) {
    if (qt(x) != 0) {
        if (qt(x) == -KF || qt(x) == KF)
            return r1(x);
        else if (compatible_f(x))
            return kf(convert_f(x));
        else
            return convert_F(x);
    }

    repeat (i, qn(x)) {
        K v = qK(x, i);
        if (qt(v) != -KF && qt(v) != KF)
            goto copy;
    }
    // list of KF shouldn't normally appear, so don't try to convert to vector,
    // but do convert () to 0#0.
    if (qn(x))
        return r1(x);

copy:;
    K r = ktn(0, qn(x));
    I vector = 1;
    repeat (i, qn(x)) {
        K v = convert_FFF(qK(x, i));
        if (!v)
            return q0(r);
        assert(qt(v) == -KF || qt(v) == KF || qt(v) == 0);
        vector &= qt(v) == -KF;
        qK(r, i) = v;
    }
    if (vector) {
        K u = ktn(KF, qn(r));
        repeat (i, qn(r))
            qF(u, i) = qf(qK(r, i));
        q0(r);
        return u;
    } else
        return r;
}



int
item_I(I* r, K x, I i) {
    assert(has_n(x) && i < qn(x));
    switch (xt) {
    case KJ:
        *r = i_from_j(qJ(x, i));
        return 1;
    case KI:
        *r = qI(x, i);
        return 1;
    case KH:
        *r = i_from_h(qH(x, i));
        return 1;
    case KB:
        *r = qB(x, i);
        return 1;
    case 0:
        if (compatible_i(qK(x, i))) {
            *r = convert_i(qK(x, i));
            return 1;
        } else
            return 0;
    default:
        return 0;
    }
}


int
item_F(F* r, K x, I i) {
    assert(has_n(x) && i < qn(x));
    switch (xt) {
    case KF:
        *r = qF(x, i);
        return 1;
    case KE:
        *r = qE(x, i);
        return 1;
    case KJ:
        *r = f_from_j(qJ(x, i));
        return 1;
    case KI:
        *r = f_from_i(qI(x, i));
        return 1;
    case KH:
        *r = f_from_h(qH(x, i));
        return 1;
    case KB:
        *r = qB(x, i);
        return 1;
    case 0:
        if (compatible_f(qK(x, i))) {
            *r = convert_f(qK(x, i));
            return 1;
        } else
            return 0;
    default:
        return 0;
    }
}


K
make_F_null(I n) {
    K x = ktn(KF, n);
    repeat (i, n)
        xF[i] = nf;
    return x;
}


K
make_F(const F* a, I n) {
    K x = ktn(KF, n);
    memcpy(xF, a, n * sizeof(F));
    return x;
}


K
new_D()
{
    return xD(ktn(KS, 0), ktn(0,  0));
}

void
append_D(K x, S k, K v)
{
    js(&xK[0], ss(k));
    jk(&xK[1], v);
}
