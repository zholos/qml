#include <float.h>
#include <math.h>

#include <conmax.h>

#include "conmax.h"

#include "alloc.h"
#include "conmin.h"


double
d1mach_(int* i) {
    // This should only be called by CONMAX, and only with this argument.
    if (*i != 3)
        abort();

    // DBL_EPSILON is pow(FLT_RADIX, 1-DBL_MANT_DIG),
    // CONMAX needs   pow(FLT_RADIX,  -DBL_MANT_DIG).
    return DBL_EPSILON / FLT_RADIX;
}


int
fnset_(I* nparm, I* numgr, F* pttbl, F* param,
       I* ipt, I* indfn, I* icntyp, F* confun)
{
    I which = *ipt - 1;
    confun += which;
    *confun = eval_param((struct eval_info*)pttbl,
                         which, param, *nparm,
                         *indfn ? confun + *numgr : NULL, *numgr,
                         icntyp + which);
    return 0;
}


K
solvemin(K fun, K con, K start_, I maxiter, F tolcon, I steps,
         int slp, int rk, int lincon, int full, int quiet)
{
    //          fun        con               constraints
    // solve    NULL       callable or list  n type 2
    // min      callable   empty list        1 type 1
    // conmin   callable   callable or list  1 type 1, n type -1

    if (fun && !callable(fun) || !callable(con) && qt(con) != 0)
        return krr("type");

    if (slp && (rk || steps > 0))
        return krr("opt");

    I ifun, numgr = qt(con) ? 1 : qn(con);
    if (fun) { // min or conmin
        numgr = add_size(numgr, 1, 1);
        ifun = 1;
    } else // solve
        ifun = numgr;

    K start = convert_FFF(start_);
    if (!start)
        return krr("type");

    S err = NULL;
    I nparm;
    F* param = take_param(start, &nparm, &err);

    F* vfun = alloc_F(&ifun, &err);
    repeat (i, ifun)
        vfun[i] = 0;

    I lerror = add_size(3, numgr, 1);
    F* error = alloc_F(&lerror, &err);

    I lwrk = add_size(add_size(13, numgr, add_size(11, nparm, 4)), 
                                   nparm, add_size(27, nparm, 2));
    F* work = alloc_F(&lwrk, &err);

    I liwrk = add_size(add_size(3, nparm, 7), numgr, 7);
    I* iwork = alloc_I(&liwrk, &err);

    struct eval_info info;
    info.call.arg = qt(start) == -KF ? 1 : 0;
    info.call.base = 0;
    info.call.start = start;
    info.call.error = no_error;
    info.fun = fun;
    info.con = con;
    info.contyp = lincon ? -1 : -2;
    info.con_sign = 1;

    // not safe to make the call in case of error because some arrays have an
    // assumed minimum size
    if (err) {
        info.call.error = krr(err);
        goto skip_call;
    }

    I ioptn = 200;
    if (maxiter < 0)
        maxiter = 1000; // option default
    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    work[1] = tolcon;
    if (steps > 0) {
        iwork[1] = steps;
        ioptn += 100;
    }
    ioptn += slp ? 1000 : rk ? 2000 : 0;

    I iter;
    conmax_(&ioptn, &nparm, &numgr, &maxiter, vfun, &ifun, (F*)&info,
            iwork, &liwrk, work, &lwrk, &iter, param, error);

skip_call:
    free_I(iwork);
    free_F(work);
    free_F(vfun);

    // either err was set above, or error in user function
    if (info.call.error != no_error)
        goto skip_result;

    S sig = NULL; // different from err, may be returned instead of raised
    if (iter >= maxiter)
        sig = "iter";
    else if (iter < 0) {
        // objective isn't evaluated in this case
        if (fun)
            error[numgr] = eval_param(&info, 0, param, nparm, NULL, 1, NULL);
        sig = "feas";
    } else if (!fun && !(fabs(error[numgr]) <= tolcon))
        sig = "feas";
    else
        repeat (i, nparm)
            if (isnan(param[i]))
                sig = "nan";

    if (!quiet && sig) {
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
        // solve   `x       `iter  `x`last       `err`iter`sig
        // min     `x`f     `iter  `x`last`f         `iter`sig
        // conmin  `x`f`cons`iter  `x`last`f`cons`err`iter`sig

        K d = new_D();
        append_D(d, "x", x);
        if (sig)
            append_D(d, "last", last);

        if (fun) {
            append_D(d, "f", kf(error[numgr]));
            if (con != empty_con) {
                 K cons = ktn(KF, numgr-1);
                 repeat (i, numgr-1) {
                     F e = error[i+1];
                     qF(cons, i) = fabs(e) <= tolcon ? 0 : -e;
                 }
                 append_D(d, "cons", cons);
            }
        }

        if (sig && con != empty_con)
            append_D(d, "err", kf(error[numgr + (!fun ? 0 : lincon ? 1 : 2)]));

        append_D(d, "iter", ki(max_i(0, iter)));
        if (sig)
            append_D(d, "sig", ks(sig));

        x = d;
    }

skip_result:
    free_F(error);
    free_F(param);
    q0(start);

    if (info.call.error != no_error)
        return info.call.error;
    return x;
}



K
root(K fun, K start, I maxiter, F tolcon, int full, int quiet)
{
    if (!callable(fun) || !has_n(start))
        return krr("type");
    if (qn(start) != 2)
        return krr("length");
    F p1, p2;
    if (!item_F(&p1, start, 0) || !item_F(&p2, start, 1))
        return krr("type");

    if (p1 > p2)
        swap_f(&p1, &p2);

    struct eval_info info;
    info.call.arg = 1;
    info.call.base = 0;
    info.call.error = no_error;
    info.fun = NULL; // root/solve flag
    info.con = fun;
    info.contyp = -2;
    info.con_sign = 1;

    F f1 = call_param(&info.call, 1, fun, &p1);
    F f2 = call_param(&info.call, 1, fun, &p2);
    // info.call.error possibly set

    if (maxiter < 0)
        maxiter = 100; // option default
    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    if (f1 < -tolcon && f2 > tolcon) {
        f1 = -f1;
        f2 = -f2;
        info.con_sign = -1;
    }
    int sig_sign = !(f1 > tolcon && f2 < -tolcon);

    I nsrch, ioptn = 0, nparm = 1, numgr = 1, ifun = 1,
      iphse = 0, iwork[17], liwrk = 17, lwrk = 6;
    F dvec = 1, cfun = 0, zwork = 0, work[6], parwrk, err1[4];
    iwork[15] = -2;

    muller_(&maxiter, &nsrch, &ioptn, &nparm, &numgr, &dvec, &cfun, &ifun,
            (F*)&info, &zwork, &tolcon, &iphse, iwork, &liwrk, work, &lwrk,
            &parwrk, err1, &p1, &f1, &p2, &f2);

    if (info.call.error != no_error)
        return info.call.error;

    S sig = NULL;
    if (nsrch >= maxiter)
        sig = "iter";
    else if (!(fabs(f2) <= tolcon))
        sig = "feas";
    else if (isnan(p2))
        sig = "nan";

    if (sig && sig_sign)
        sig = "sign";

    if (!quiet && sig)
        return krr(sig);

    K x = sig ? kf(nf) : kf(p2);
    if (full) {
        //        normal   error
        // root  `x`iter  `x`last`err`iter`sig
        K d = new_D();
        append_D(d, "x", x);
        if (sig) {
            append_D(d, "last", kf(p2));
            append_D(d, "err", kf(fabs(f2)));
        }
        append_D(d, "iter", ki(nsrch));
        if (sig)
            append_D(d, "sig", ks(sig));
        x = d;
    }
    return x;
}


K
line(K fun, K base, K start, I maxiter, F tolcon, int full, int quiet)
{
    if (!callable(fun) || !compatible_f(base) || !compatible_f(start))
        return krr("type");

    struct eval_info info;
    info.call.arg = 1;
    info.call.base = convert_f(base);
    info.call.error = no_error;
    info.fun = fun;
    info.con = fun; // min/line flag

    F projct = convert_f(start) - info.call.base;
    if (projct < 0) {
        projct = -projct;
        info.call.arg = -1;
    }

    // Call sequence based on sample driver from conmax.f

    if (maxiter < 0)
        maxiter = 100; // option default
    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    F prjlim = wf;
    F tol1 = 100 * DBL_EPSILON;
    if (tol1 > tolcon)
        tol1 = tolcon;

    I nadd = maxiter, lims1 = maxiter, ioptn = 0, numgr = 1,
      nparm = 1, ifun = 1, mact = 1, iact = 1, iphse = 0, itypm1 = 0,
      itypm2 = 0, iwork[17], liwrk = 17, lwrk = 42, nsrch;
    F cx[2], cfun = 0, param = 0, error[4], rchdwn = 2, unit = 1, rchin = 2,
      work[42], err1[4], parprj, emin, emin1, parser;
    cx[0] = 1;
    iwork[6] = 1;

    searsl_(&maxiter, &nadd, &lims1, &ioptn, &numgr, &nparm, &prjlim, &tol1, cx,
            &cfun, &ifun, (F*)&info,
            &param, error,
            &rchdwn, &mact, &iact, &iphse, &unit,
            &tolcon, &rchin, &itypm1, &itypm2,
            iwork, &liwrk, work, &lwrk, err1, &parprj, &projct,
            &emin, &emin1, &parser, &nsrch);

    if (info.call.error != no_error)
        return info.call.error;

    S sig = NULL;
    if (nsrch >= maxiter || nsrch >= lims1)
        sig = "iter";
    else if (isnan(projct))
        sig = "nan";

    if (!quiet && sig)
        return krr(sig);

    projct = info.call.base + info.call.arg * projct;
    K x = sig ? kf(nf) : kf(projct);
    if (full) {
        //        normal     error
        // line  `x`f`iter  `x`last`f`iter`sig
        K d = new_D();
        append_D(d, "x", x);
        if (sig)
             append_D(d, "last", kf(projct));
        append_D(d, "f", kf(emin));
        append_D(d, "iter", ki(nsrch));
        if (sig)
            append_D(d, "sig", ks(sig));
        x = d;
    }
    return x;
}
