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


// The first member must be of the same type as pttbl. This allows casting a
// structure pointer to a first-member pointer, passing it through CONMAX, and
// casting it back to get the original pointer in a well-defined manner.
struct fnset_info {
    F base; /* first member */
    K fun, con, start;
    I contyp;
    int neg;
    int con_neg;
    K error;
};

static F
fnset_call(struct fnset_info* info, K f, int neg, F* param)
{
    if (info->error != no_error)
        return 0;

    K a;
    if (!info->start && info->fun) // line
        a = knk(1, kf(info->base + (info->neg ? -*param : *param)));
    else if (!info->start || info->start->t == -KF) // root, or scalar param
        a = knk(1, kf(*param));
    else
        make_param(info->start, param, &a);

    K x = dot(f, a);
    q0(a);
    if (x && compatible_f(x)) {
        F v = convert_f(x);
        q0(x);
        if (isnan(v)) // protect optimization routines from NaNs
            return wf; // this is -wf for ">=" constraints
        return neg ? -v : v;
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

int
fnset_(I* nparm, I* numgr, F* pttbl, F* param,
       I* ipt, I* indfn, I* icntyp, F* confun)
{
    struct fnset_info* info = (struct fnset_info*)pttbl;

    icntyp += *ipt - 1;
    confun += *ipt - 1;

    K f;
    int neg;
    assert(callable(info->con) || qt(info->con) == 0);
    if (info->fun) // line, min or conmin
        if (*ipt == 1) { // objective function
            f = info->fun;
            neg = 0;
            *icntyp = 1;
        } else { // constraints
            f = qt(info->con) ? info->con : qK(info->con, *ipt-2);
            neg = 1;
            *icntyp = info->contyp;
        }
    else { // root or solve
        f = qt(info->con) ? info->con : qK(info->con, *ipt-1);
        neg = info->con_neg;
        *icntyp = 2;
    }

    F v = *confun = fnset_call(info, f, neg, param);
    if (*indfn) {
        I m = *numgr, n = *nparm;
        if (*icntyp == -1) // linear function
            repeat (i, n) {
                F p = param[i];
                param[i] = p + 1;
                *(confun += m) = fnset_call(info, f, neg, param) - v;
                param[i] = p;
            }
        else { // nonlinear function
            F h = sqrt(DBL_EPSILON / FLT_RADIX);
            repeat (i, n) {
                F p = param[i], p1, p2, v1, v2;
                param[i] = p1 = p + h; v1 = fnset_call(info, f, neg, param);
                param[i] = p2 = p - h; v2 = fnset_call(info, f, neg, param);
                *(confun += m) = (v1 - v2) / (p1 - p2);
                param[i] = p;
            }
        }
    }

    if (info->error != no_error) // deactivate all constraints to exit sooner
        *icntyp = 0;

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

    struct fnset_info info;
    info.fun = fun;
    info.con = con;
    info.start = start;
    info.contyp = lincon ? -1 : -2;
    info.con_neg = 0;
    info.error = no_error;

    // not safe to make the call in case of error because some arrays have an
    // assumed minimum size
    if (err) {
        info.error = krr(err);
        goto skip_call;
    }

    I ioptn = 200;
    I itlim = max_i(0, maxiter);
    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    work[1] = tolcon;
    if (steps > 0) {
        iwork[1] = steps;
        ioptn += 100;
    }
    ioptn += slp ? 1000 : rk ? 2000 : 0;

    I iter;
    conmax_(&ioptn, &nparm, &numgr, &itlim, vfun, &ifun, (F*)&info,
            iwork, &liwrk, work, &lwrk, &iter, param, error);

skip_call:
    free_I(iwork);
    free_F(work);
    free_F(vfun);

    // either err was set above, or error in user function
    if (info.error != no_error)
        goto skip_result;

    S sig = NULL; // different from err, may be returned instead of raised
    if (iter >= maxiter)
        sig = "iter";
    else if (iter < 0)
        sig = "feas";
    else if (!fun && !(error[numgr] >= -tolcon && error[numgr] <= tolcon))
        sig = "feas";
    else
        repeat (i, nparm)
            if (isnan(param[i]))
                sig = "nan";

    if (!quiet && sig) {
        info.error = krr(sig);
        goto skip_result;
    }

    assert(info.error == no_error);

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
                     qF(cons, i) = e >= -tolcon && e <= tolcon ? 0 : -e;
                 }
                 append_D(d, "cons", cons);
            }
        }

        if (sig && con != empty_con)
            append_D(d, "err", kf(error[numgr + (!fun ? 0 : lincon ? 1 : 2)]));

        append_D(d, "iter", ki(iter >= 0 ? iter : 0));
        if (sig)
            append_D(d, "sig", ks(sig));

        x = d;
    }

skip_result:
    free_F(error);
    free_F(param);
    q0(start);

    if (info.error != no_error)
        return info.error;
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

    struct fnset_info info;
    info.fun = NULL; // root/solve flag
    info.con = fun;
    info.start = NULL; // root/line flag
    info.contyp = -2;
    info.con_neg = 0;
    info.error = no_error;

    F f1 = fnset_call(&info, fun, 0, &p1);
    F f2 = fnset_call(&info, fun, 0, &p2);
    // info.error possibly set

    if (maxiter < 0)
        maxiter = 0;
    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    if (f1 < -tolcon && f2 > tolcon) {
        f1 = -f1;
        f2 = -f2;
        info.con_neg = 1;
    }
    int sig_sign = !(f1 > tolcon && f2 < -tolcon);

    I limmul = maxiter, nsrch, ioptn = 0, nparm = 1, numgr = 1, ifun = 1,
      iphse = 0, iwork[17], liwrk = 17, lwrk = 6;
    F dvec = 1, cfun = 0, zwork = 0, work[6], parwrk, err1[4];
    iwork[15] = -2;

    muller_(&limmul, &nsrch, &ioptn, &nparm, &numgr, &dvec, &cfun, &ifun,
            (F*)&info, &zwork, &tolcon, &iphse, iwork, &liwrk, work, &lwrk,
            &parwrk, err1, &p1, &f1, &p2, &f2);

    if (info.error != no_error)
        return info.error;

    S sig = NULL;
    if (nsrch >= maxiter)
        sig = "iter";
    else if (!(f2 >= -tolcon && f2 <= tolcon))
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

    struct fnset_info info;
    info.base = convert_f(base);
    info.fun = fun;
    info.con = fun; // min/line flag
    info.start = NULL; // root/line flag
    info.neg = 0;
    info.error = no_error;

    F projct = convert_f(start) - info.base;
    if (projct < 0) {
        projct = -projct;
        info.neg = 1;
    }

    // Call sequence based on sample driver from conmax.f

    if (maxiter < 0)
        maxiter = 0;
    if (!(tolcon >= 0))
        tolcon = sqrt(DBL_EPSILON / FLT_RADIX);
    F prjlim = wf;
    F tol1 = 100 * DBL_EPSILON;
    if (tol1 > tolcon)
        tol1 = tolcon;

    I initlm = maxiter, nadd = initlm, lims1 = initlm, ioptn = 0, numgr = 1,
      nparm = 1, ifun = 1, mact = 1, iact = 1, iphse = 0, itypm1 = 0,
      itypm2 = 0, iwork[17], liwrk = 17, lwrk = 42, nsrch;
    F cx[2], cfun = 0, param = 0, error[4], rchdwn = 2, unit = 1, rchin = 2,
      work[42], err1[4], parprj, emin, emin1, parser;
    cx[0] = 1;
    iwork[6] = 1;

    searsl_(&initlm, &nadd, &lims1, &ioptn, &numgr, &nparm, &prjlim, &tol1, cx,
            &cfun, &ifun, (F*)&info,
            &param, error,
            &rchdwn, &mact, &iact, &iphse, &unit,
            &tolcon, &rchin, &itypm1, &itypm2,
            iwork, &liwrk, work, &lwrk, err1, &parprj, &projct,
            &emin, &emin1, &parser, &nsrch);

    if (info.error != no_error)
        return info.error;

    S sig = NULL;
    if (nsrch >= maxiter || nsrch >= lims1)
        sig = "iter";
    else if (isnan(projct))
        sig = "nan";

    if (!quiet && sig)
        return krr(sig);

    projct = info.base + (info.neg ? -projct : projct);
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
