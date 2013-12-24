#include <float.h>
#include <math.h>
#include <string.h>

#include "conmin.h"

#include "alloc.h"
#include "conmax.h"
#include "nlopt.h"
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


F
call_param(struct call_info* info, int sign, K f, F* param)
{
    if (info->error != no_error)
        return 0;

    K a;
    if (info->arg)
        a = knk(1, kf(info->base + info->arg * *param));
    else
        make_param(info->start, param, &a);

    K x = dot(f, a);
    q0(a);
    if (x && compatible_f(x)) {
        F v = convert_f(x);
        q0(x);
        if (isnan(v)) // protect optimization routines from NaNs
            return wf; // this is -wf for ">=" constraints
        return sign * v;
    }

    // function didn't return a float as we'd hoped
    if (!x || qt(x) == -128)
        info->error = x;
    else {
        info->error = krr(callable(x) ? "rank" : "type");
        q0(x);
    }
    return 0;
}


F
eval_param(struct eval_info* info,
           I which, F* param, I n, F* grad, int grad_step, I* contyp_)
{
    K f;
    int sign;
    I contyp;
    assert(callable(info->con) || qt(info->con) == 0);
    if (info->fun) // line, min or conmin
        if (which == 0) { // objective function
            f = info->fun;
            sign = 1;
            contyp = 1;
        } else { // constraints
            f = qt(info->con) ? info->con : qK(info->con, which-1);
            sign = -1;
            contyp = info->contyp;
        }
    else { // root or solve
        f = qt(info->con) ? info->con : qK(info->con, which);
        sign = info->con_sign;
        contyp = 2;
    }

    F v = call_param(&info->call, sign, f, param);
    if (grad) {
        // CONMAX can do this automatically, but we've already set up the call
        if (contyp == -1) // linear function
            // don't need a centered difference for linear constraints
            repeat (i, n) {
                F p = param[i];
                param[i] = p + 1;
                *grad = call_param(&info->call, sign, f, param) - v;
                grad += grad_step;
                param[i] = p;
            }
        else { // nonlinear function
            // same algorithm and step as in CONMAX
            const F h = sqrt(DBL_EPSILON / FLT_RADIX);
            repeat (i, n) {
                F p = param[i], p1 = p + h, p2 = p - h, v1, v2;
                param[i] = p1; v1 = call_param(&info->call, sign, f, param);
                param[i] = p2; v2 = call_param(&info->call, sign, f, param);
                *grad = (v1 - v2) / (p1 - p2);
                grad += grad_step;
                param[i] = p;
            }
        }
    }

    if (info->call.error != no_error)
        contyp = 0; // deactivate all constraints to exit sooner

    if (contyp_)
        *contyp_ = contyp;
    return v;
}


struct k0 no_error_;
struct k0 empty_con_ = { 0 };



static const struct optn solve_opt[] = {
    [0] = { "iter",   -KI },
    [1] = { "tol",    -KF },
    [2] = { "steps",  -KI },
    [3] = { "slp",      0 },
    [4] = { "rk",       0 },
    [5] = { "full",     0 },
    [6] = { "quiet",    0 },
          { NULL }
};

K
qml_solvex(K opts, K x, K y)
{
    union optv v[] = { { -1 }, { .f = -1 }, { -1 },
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

static const struct optn min_opt[] = {
    [0] = { "iter",   -KI },
    [1] = { "tol",    -KF },
    [2] = { "steps",  -KI },
    [3] = { "slp",      0 },
    [4] = { "rk",       0 },
    [5] = { "nm",       0 },
    [6] = { "sbplx",    0 },
    [7] = { "full",     0 },
    [8] = { "quiet",    0 },
          { NULL }
};

K
qml_minx(K opts, K x, K y)
{
    union optv v[] = { { -1 }, { .f = -1 }, { -1 },
                       { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
    if (!take_opt(opts, min_opt, v))
        return krr("opt");
    if (v[5].i || v[6].i) {
        if (v[5].i && v[6].i || v[1].f >= 0 || v[2].i >= 0 || v[3].i || v[4].i)
            return krr("opt");
        return nloptmin(x, empty_con, y, v[0].i, v[1].f,
                        v[6].i, 0, v[7].i, v[8].i);
    } else
        return solvemin(x, empty_con, y, v[0].i, v[1].f, v[2].i,
                        v[3].i, v[4].i, 0, v[7].i, v[8].i);
}

K
qml_min(K x, K y)
{
    return qml_minx(empty_con, x, y);
}

static const struct optn conmin_opt[] = {
    [0] = { "iter",   -KI },
    [1] = { "tol",    -KF },
    [2] = { "steps",  -KI },
    [3] = { "slp",      0 },
    [4] = { "rk",       0 },
    [5] = { "cobyla",   0 },
    [6] = { "lincon",   0 },
    [7] = { "full",     0 },
    [8] = { "quiet",    0 },
          { NULL }
};

K
qml_conminx(K opts, K x, K y, K z)
{
    union optv v[] = { { -1 }, { .f = -1 }, { -1 },
                       { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
    if (!take_opt(opts, conmin_opt, v))
        return krr("opt");
    if (v[5].i) {
        if (v[2].i >= 0 || v[3].i || v[4].i)
            return krr("opt");
        return nloptmin(x, y, z, v[0].i, v[1].f,
                        0, v[6].i, v[7].i, v[8].i);
    } else
        return solvemin(x, y, z, v[0].i, v[1].f, v[2].i,
                        v[3].i, v[4].i, v[6].i, v[7].i, v[8].i);
}

K
qml_conmin(K x, K y, K z)
{
    return qml_conminx(empty_con, x, y, z);
}

static const struct optn rootline_opt[] = {
    [0] = { "iter", -KI },
    [1] = { "tol",  -KF },
    [2] = { "full",   0 },
    [3] = { "quiet",  0 },
          { NULL }
};

K
qml_rootx(K opts, K x, K y)
{
    union optv v[] = { { -1 }, { .f = -1 }, { 0 }, { 0 } };
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
    union optv v[] = { { -1 }, { .f = -1 }, { 0 }, { 0 } };
    if (!take_opt(opts, rootline_opt, v))
        return krr("opt");
    return line(x, y, z, v[0].i, v[1].f, v[2].i, v[3].i);
}

K
qml_line(K x, K y, K z)
{
    return qml_linex(empty_con, x, y, z);
}
