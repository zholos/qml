#include <string.h>

#include "conmin.h"

#include "alloc.h"
#include "conmax.h"
#include "opt.h"


static I
write_param(K x, F* a)
{
    switch (qt(x)) {
    case -KF:
        if (a)
            *a = qf(x);
        return 1;
    case KF:
        if (a)
            memcpy(a, kF(x), qn(x) * sizeof(F));
        return qn(x);
    case 0:;
        I n = 0;
        repeat (i, qn(x))
            n = add_size(n, write_param(qK(x, i), a ? a + n : a), 1);
        return n;
    default:
        assert(0);
        return 0;
    }
}


F*
take_param(K x, I* n, S* err) {
    *n = write_param(x, NULL);
    F* a = alloc_F(n, err);
    write_param(x, a);
    return a;
}


F*
make_param(K x, F* param, K* r) {
    switch (qt(x)) {
    case -KF:
        *r = kf(*param++);
        break;
    case KF:
        *r = make_F(param, qn(x));
        param += qn(x);
        break;
    case 0:
        *r = ktn(0, qn(x));
        repeat (i, qn(x))
            param = make_param(qK(x, i), param, &qK(*r, i));
        break;
    default:
        assert(0);
    }
    return param;
}



struct k0 no_error_;
struct k0 empty_con_ = { 0 };



static const struct optn solve_opt[] = {
    { "iter",   -KI },
    { "tol",    -KF },
    { "steps",  -KI },
    { "slp",      0 },
    { "rk",       0 },
    { "full",     0 },
    { "quiet",    0 },
    { NULL }
};

K
qml_solvex(K opts, K x, K y)
{
    union optv v[] = { { 1000 }, { .f = -1 }, { -1 },
                       { 0 }, { 0 }, { 0 }, { 0 } };
    if (!take_opt(opts, solve_opt, v))
        return krr("opt");
    return solvemin(NULL, x, y, v[0].i, v[1].f, v[2].i,
                    v[3].i, v[4].i, 0, v[5].i, v[6].i);
}

K
qml_solve(K x, K y)
{
    return qml_solvex(empty_con, x, y);
}

K
qml_minx(K opts, K x, K y)
{
    union optv v[] = { { 1000 }, { .f = -1 }, { -1 },
                       { 0 }, { 0 }, { 0 }, { 0 } };
    if (!take_opt(opts, solve_opt, v))
        return krr("opt");
    return solvemin(x, empty_con, y, v[0].i, v[1].f, v[2].i,
                    v[3].i, v[4].i, 0, v[5].i, v[6].i);
}

K
qml_min(K x, K y)
{
    return qml_minx(empty_con, x, y);
}

static const struct optn conmin_opt[] = {
    { "iter",   -KI },
    { "tol",    -KF },
    { "steps",  -KI },
    { "slp",      0 },
    { "rk",       0 },
    { "lincon",   0 },
    { "full",     0 },
    { "quiet",    0 },
    { NULL }
};

K
qml_conminx(K opts, K x, K y, K z)
{
    union optv v[] = { { 1000 }, { .f = -1 }, { -1 },
                       { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
    if (!take_opt(opts, conmin_opt, v))
        return krr("opt");
    return solvemin(x, y, z, v[0].i, v[1].f, v[2].i,
                    v[3].i, v[4].i, v[5].i, v[6].i, v[7].i);
}

K
qml_conmin(K x, K y, K z)
{
    return qml_conminx(empty_con, x, y, z);
}

static const struct optn rootline_opt[] = {
    { "iter", -KI },
    { "tol",  -KF },
    { "full",   0 },
    { "quiet",  0 },
    { NULL }
};

K
qml_rootx(K opts, K x, K y)
{
    union optv v[] = { { 100 }, { .f = -1 }, { 0 }, { 0 } };
    if (!take_opt(opts, rootline_opt, v))
        return krr("opt");
    return root(x, y, v[0].i, v[1].f, v[2].i, v[3].i);
}

K
qml_root(K x, K y)
{
    return qml_rootx(empty_con, x, y);
}

K
qml_linex(K opts, K x, K y, K z) {
    union optv v[] = { { 100 }, { .f = -1 }, { 0 }, { 0 } };
    if (!take_opt(opts, rootline_opt, v))
        return krr("opt");
    return line(x, y, z, v[0].i, v[1].f, v[2].i, v[3].i);
}

K
qml_line(K x, K y, K z)
{
    return qml_linex(empty_con, x, y, z);
}
