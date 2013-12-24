#include <float.h>

#include <neldermead.h>
#include <cobyla.h>

#include "nlopt.h"

#include "alloc.h"
#include "conmin.h"


static F
objective(unsigned n, const double* param, double* grad, void* info)
{
    // Casting away constness is fine because param was allocated with malloc.
    // eval_param needs to modify param but then reverts it.
    return eval_param(info, 0, (double*)param, n, grad, 1, NULL);
}

static void
constraints(unsigned m, double* result, unsigned n, const double* param,
            double* grad, void* info)
{
    for (unsigned i = 0; i < m; i++) {
        result[i] = eval_param(info, i+1, (double*)param, n, grad, 1, NULL);
        if (grad)
            grad += n;
    }
}


K
nloptmin(K fun, K con, K start_, I maxiter, F tolcon,
         I method, I lincon, I full, I quiet)
{
    //         fun       con               method
    // min     callable  empty list        0=nl, 1=sbplx
    // conmin  callable  callable or list  0=cobyla

    if (!callable(fun) || !callable(con) && qt(con) != 0)
        return krr("type");

    K start = convert_FFF(start_);
    if (!start)
        return krr("type");

    S err = NULL;
    F minf = wf;
    I nparm;
    F* param = take_param(start, &nparm, &err);

    if (nparm > 10000) // approximate limit before malloc overflow in nldrmd.c
        if (!err) err = "limit";

    I ncon = qt(con) ? 1 : qn(con);
    if (ncon > 10000) // with nparm limit, approximate limit for cobyla.c
        if (!err) err = "limit";

    F* xstep = alloc_F(&nparm, &err);
    F* lbound = alloc_F(&nparm, &err);
    F* ubound = alloc_F(&nparm, &err);
    F* xtol_abs = alloc_F(&nparm, &err);
    repeat (i, nparm) {
        xstep[i] = fmax(1, fabs(param[i])); // also see
                                            // nlopt_set_default_initial_step
        lbound[i] = -wf;
        ubound[i] = wf;
        xtol_abs[i] = 0;
    }

    F* fc_tol = alloc_F(&ncon, &err);
    repeat (i, ncon)
        // irrelevant for COBYLA because used only in minf_max check
        fc_tol[i] = 0;

    struct eval_info info;
    info.call.arg = qt(start) == -KF ? 1 : 0;
    info.call.base = 0;
    info.call.start = start;
    info.call.error = no_error;
    info.fun = fun;
    info.con = con;
    info.contyp = lincon ? -1 : -2;
    info.con_sign = 1;

    if (maxiter < 0) {
        maxiter = 10000; // option default
        if (con != empty_con)
            maxiter *= 100; // iterations are much faster in COBYLA
    }

    nlopt_stopping stop = {
        .n = nparm,
        .minf_max = -wf,
        .maxeval = max_i(1, maxiter), // 0 means ignore
        .xtol_abs = xtol_abs
    };

    if (con != empty_con) {
        if (!(tolcon >= 0))
            tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
        stop.xtol_rel = DBL_EPSILON; // COBYLA calculates rho from this
    }

    if (err) {
        info.call.error = krr(err);
        goto skip_call;
    }

    int result;

    // protect optimization routines from NaNs
    repeat (i, nparm)
        if (isnan(param[i])) {
            result = NLOPT_SUCCESS;
            // minf is set to some default value that's not important here
            // sig = "nan" will be set below
            goto skip_call;
        };

    if (con != empty_con) {
        nlopt_constraint fc = {
            .m = ncon,
            .mf = constraints,
            .f_data = &info,
            .tol = fc_tol
        };

        result = cobyla_minimize(nparm, objective, &info,
                                 1, &fc, 0, NULL, lbound, ubound,
                                 param, &minf, &stop, xstep);
    } else
        if (method)
            result = sbplx_minimize(nparm, objective, &info, lbound, ubound,
                                    param, &minf, xstep, &stop);
        else
            result = nldrmd_minimize(nparm, objective, &info, lbound, ubound,
                                     param, &minf, xstep, &stop);

skip_call:
    free_F(fc_tol);
    free_F(xtol_abs);
    free_F(ubound);
    free_F(lbound);
    free_F(xstep);

    if (info.call.error != no_error)
        goto skip_result;

    // any subsequent errors can be silenced
    S sig = NULL;
    switch (result) {
        case NLOPT_FTOL_REACHED:
        case NLOPT_XTOL_REACHED:
        case NLOPT_SUCCESS:
            repeat (i, nparm)
                if (isnan(param[i]))
                    sig = "nan";
            break;
        case NLOPT_MAXEVAL_REACHED:
            sig = "iter";
            break;
        case NLOPT_INVALID_ARGS:
            sig = "qml_assert";
            break;
        case NLOPT_ROUNDOFF_LIMITED:
            sig = "round";
            break;
        case NLOPT_OUT_OF_MEMORY:
            sig = "wsfull";
            break;
        default:
            sig = "fail";
    }

    F con_err = 0;
    K cons;
    if (con != empty_con) {
        if (full)
            cons = ktn(KF, ncon);

        repeat (i, ncon) {
            F e = eval_param(&info, i+1, param, nparm, NULL, 1, NULL);
            con_err = fmax(con_err, e);
            if (full)
                qF(cons, i) = fabs(e) <= tolcon ? 0 : -e;
        }
        if (!sig && con_err > tolcon)
            sig = "feas";
    }

    if (!quiet && sig) {
        if (con != empty_con && full)
            q0(cons);
        info.call.error = krr(sig);
        goto skip_result;
    }

    assert(info.call.error == no_error);

    K x, last;
    if (sig) {
        if (full)
            make_param(start, param, &last);
        repeat (i, nparm)
            param[i] = nf;
    }
    make_param(start, param, &x);

    if (full) {
        //          normal          sig
        // min     `x`f     `iter  `x`last`f         `iter`sig
        // conmin  `x`f`cons`iter  `x`last`f`cons`err`iter`sig

        K d = new_D();
        append_D(d, "x", x);
        if (sig)
            append_D(d, "last", last);
        append_D(d, "f", kf(minf));
        if (con != empty_con) {
            append_D(d, "cons", cons);
            if (sig)
                append_D(d, "err", kf(con_err));
        }
        append_D(d, "iter", ki(stop.nevals));
        if (sig)
            append_D(d, "sig", ks(sig));
        x = d;
    }

skip_result:
    free_F(param);
    q0(start);

    if (info.call.error != no_error)
        return info.call.error;

    return x;
}
