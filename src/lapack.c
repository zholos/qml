#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <lapack.h>

#include "alloc.h"
#include "opt.h"
#include "util.h"


// Override error-handling routine. This is called by LAPACK before returning an
// info value indicating an error. We check the info value and ignore the call.
int
xerbla_(char* srname, int* info) {
    return 0;
}


static void
check_info(int info, S* err) {
    if (info < 0)
        if (!*err)
            *err = always(info == -7) ? "wsfull" : "qml_assert";
}


static I
take_maxwork(I info, F maxwork) {
    return likely(!info && maxwork >= 0 && maxwork < wi) ? (I)maxwork : wi;
}


// If the malloc fails, uses one item on the stack because LAPACK needs at least
// one item even for invalid requests. The array must be freed in the same
// block.
#define alloc_W(lwork_)                 \
    I lwork = (lwork_);                 \
    F w__min[1],                        \
      *w__full = alloc_F(&lwork, &err), \
      *w = lwork ? w__full : w__min

#define free_W() \
    free_F(w__full)


static F*
take_square_matrix(K x_, I* n, int* triangular, S* err) {
    K x = convert_FF(x_);
    if (!likely(x && (*n = qnw(x)))) {
        *n = 0;
        if (!*err) *err = "type";
    }

    F* a = alloc_FF(n, *n, err);
    repeat (j, *n)
        if (qn(qK(x, j)) != *n) {
            *n = 0;
            if (!*err) *err = "length";
        }

    if (triangular) {
        int upper = 1, lower = 1;
        repeat (j, *n)
            repeat (i, *n)
                if ((a[j + i * *n] = qF(qK(x, j), i)))
                    if (i < j) // non-zero below diagonal, so not upper
                        upper = 0;
                    else if (i > j) // non-zero above diagonal, so not lower
                        lower = 0;
        *triangular = upper ? 1 : lower ? -1 : 0;
    } else
        repeat (j, *n)
            repeat (i, *n)
                a[j + i * *n] = qF(qK(x, j), i);

    if (x) q0(x);
    return a;
}


static int alloc_square_; // flag for passing instead of b_column
#define alloc_square (&alloc_square_)

// ldr = m is allowed, otherwise *ldr must be set on input
static F*
take_matrix(K x_, I* ldr, I* m, I* n, int* column, S* err) {
    if (column && column != alloc_square) {
        K x = convert_F(x_);
        if (x) {
            *m = qnw(x), *n = 1;
            *ldr = max_i(*ldr, *m); // if ldr == m, *ldr is set above
            F* a = alloc_F(ldr, err);
            if (!*ldr)
                *m = 0;
            memcpy(a, xF, *m * sizeof(F));
            q0(x);
            *column = 1;
            return a;
        }
        *column = 0;
    }

    K x = convert_FF(x_);
    if (!likely(x && (*m = qnw(x)) && (*n = qnw(qK(x, 0))))) {
        *m = *n = 0;
        if (!*err) *err = "type";
    }
    repeat (j, *m)
        if (qn(qK(x, j)) != *n) {
            *m = *n = 0;
            if (!*err) *err = "length";
        }
    *ldr = max_i(*ldr, *m); // if ldr == m, *ldr is set above

    I nm = column == alloc_square ? max_i(*m, *n) : *n;
    F* a = alloc_FF(&nm, *ldr, err);
    if (!nm)
        *m = *n = 0;
    repeat (j, *m)
        repeat (i, *n)
            a[j + i * *ldr] = qF(qK(x, j), i);
    if (x) q0(x);
    return a;
}


#define make_transposed (-1)
#define make_upper      (-2)
#define make_lower      (-3)

// if a is null matrix is filled with NaNs
static K
make_matrix(const F* a, I ldr, I m, I n, int column) {
    if (column == make_transposed) {
        if (!a)
            return make_matrix(NULL, 0, n, m, 0);
        K r = ktn(0, n);
        repeat (i, n)
            qK(r, i) = make_F(a + i*ldr, m);
        return r;
    }

    if (column == make_upper || column == make_lower) {
        if (!a)
            return make_matrix(NULL, 0, m, n, 0);
        K r = ktn(0, m);
        repeat (j, m) {
            K x = ktn(KF, n);
            repeat (i, n)
                if (column == make_upper)
                    qF(x, i) = i < j ? 0 : a[j + ldr*i];
                else
                    qF(x, i) = i < j ? a[j + ldr*i] : i == j ? 1 : 0;
            qK(r, j) = x;
        }
        return r;
    }

    if (column) {
        // assert(n == 1);
        if (a)
            return make_F(a, m);
        else
            return make_F_null(m);
    }

    K r = ktn(0, m);
    if (a)
        // Most common case
        // TODO: cache-optimized transpose
        repeat (j, m) {
            K x = ktn(KF, n);
            repeat (i, n)
                qF(x, i) = a[j + i*ldr];
            qK(r, j) = x;
        }
    else {
        K x = make_F_null(n);
        repeat (j, m)
            qK(r, j) = r1(x);
        q0(x);
    }
    return r;
}


static K
make_complex(F a, F b) {
    if (b == 0) // in case isnan(b), return complex pair
        return kf(a);
    else {
        K x = ktn(KF, 2);
        qF(x, 0) = a;
        qF(x, 1) = b;
        return x;
    }
}


// a and b can be null
// returns xt==0 when complex
static K
make_complex_vector(const F* a, const F* b, int n) {
    if (b)
        repeat (i, n)
            if (!(b[i] == 0))
                goto complex;

    if (a) {
        K x = ktn(KF, n);
        repeat (i, n)
            qF(x, i) = a[i];
        return x;
    } else
        return make_F_null(n);

complex:;
    K x = ktn(0, n);
    repeat (i, n)
        qK(x, i) = make_complex(a[i], b[i]);
    return x;
}



// Matrix determinant
K
qml_mdet(K x) {
    int triangular;
    I n, info;
    S err = NULL;

    F* a = take_square_matrix(x, &n, &triangular, &err);
    F r = 1;

    if (!triangular) {
        I* ipiv = alloc_I(&n, &err);

        dgetrf_(&n, &n, a, &n, ipiv, &info);
        check_info(info, &err);
        // info>0 indicates singularity, but we don't treat this specially

        repeat (i, n)
            if (ipiv[i]-1 != i)
                r = -r;
        free_I(ipiv);
    }

    repeat (i, n)
        r *= a[i + i*n];
    free_F(a);

    if (err) return krr(err);
    return kf(r);
}


// Matrix inverse
K
qml_minv(K x) {
    I n, info;
    S err = NULL;

    F* a = take_square_matrix(x, &n, NULL, &err);
    I* ipiv = alloc_I(&n, &err);

    dgetrf_(&n, &n, a, &n, ipiv, &info);
    check_info(info, &err);
    // info > 0 indicates singularity

    if (!info) {
        I lwork_query = -1;
        F maxwork;
        dgetri_(&n, NULL, &n, NULL, &maxwork, &lwork_query, &info);
        check_info(info, &err);

        alloc_W(take_maxwork(info, maxwork));

        dgetri_(&n, a, &n, ipiv, w, &lwork, &info);
        check_info(info, &err);
        // info > 0 indicates singularity

        free_W();
    }

    free_I(ipiv);

    x = make_matrix(info ? NULL : a, n, n, n, 0);
    free_F(a);

    return check_err(x, err);
}


// Matrix multiply
K
qml_mm(K x, K y) {
    int b_column;
    I a_m, a_n, b_m, b_n;
    S err = NULL;

    F* a = take_matrix(x, &a_m, &a_m, &a_n, NULL, &err);
    F* b = take_matrix(y, &b_m, &b_m, &b_n, &b_column, &err);
    if (a_n != b_m)
        if (!err) err = "length";

    F* r = alloc_FF(&b_n, a_m, &err);

    if (!err) {
        int i_1 = 1;
        double f_0 = 0, f_1 = 1;
        if (b_n == 1)
            dgemv_("N", &a_m, &a_n, &f_1, a, &a_m, b, &i_1, &f_0, r, &i_1);
        else if (a_m == 1)
            dgemv_("T", &b_m, &b_n, &f_1, b, &b_m, a, &i_1, &f_0, r, &i_1);
        else
            dgemm_("N", "N", &a_m, &b_n, &a_n,
                &f_1, a, &a_m, b, &b_m, &f_0, r, &a_m);
    }

    free_F(b);
    free_F(a);
    x = make_matrix(err ? NULL : r, a_m, a_m, b_n, b_column);
    free_F(r);

    return check_err(x, err);
}


// Matrix substitution solve
K
qml_ms(K x, K y) {
    int a_triangular, b_column;
    I a_n, b_m, b_n;
    S err = NULL;

    F* a = take_square_matrix(x, &a_n, &a_triangular, &err);
    if (!a_triangular)
        if (!err) err = "domain";

    F* b = take_matrix(y, &b_m, &b_m, &b_n, &b_column, &err);
    if (a_n != b_m)
        if (!err) err = "length";

    I info;
    for (info = a_n; info; info--)
        if (!(a[(info-1) + (info-1)*a_n] != 0)) // consider NaN
            break; // info > 0 indicates singularity

    if (!err && !info) {
        int i_1 = 1;
        double f_1 = 1;
        if (b_n == 1)
            dtrsv_(a_triangular > 0 ? "U" : "L", "N", "N",
                   &a_n, a, &a_n, b, &i_1); 
        else
            dtrsm_("L", a_triangular > 0 ? "U" : "L", "N", "N",
                   &b_m, &b_n, &f_1, a, &a_n, b, &b_m);
    }

    free_F(a);
    x = make_matrix(info ? NULL : b, b_m, b_m, b_n, b_column);
    free_F(b);

    return check_err(x, err);
}


// Eigenvalues and eigenvectors
K
qml_mevu(K x) {
    I n, info;
    S err = NULL;

    F* a = take_square_matrix(x, &n, NULL, &err);

    I lwork_query = -1;
    F maxwork;
    dgeev_("N", "V", &n, NULL, &n, NULL, NULL,
           NULL, &n, NULL, &n, &maxwork, &lwork_query, &info);
    check_info(info, &err);

    alloc_W(take_maxwork(info, maxwork));
    F* lr = alloc_F (&n,    &err);
    F* li = alloc_F (&n,    &err);
    F* ev = alloc_FF(&n, n, &err);

    info = -1;
    dgeev_("N", "V", &n, a, &n, lr, li, NULL, &n, ev, &n, w, &lwork, &info);
    check_info(info, &err);
    free_W();
    free_F(a);

    x = ktn(0, 2);
    qK(x, 0) = make_complex_vector(info ? NULL : lr, info ? NULL : li, n);

    if (qt(qK(x, 0))) // no complex elements (also when info)
        qK(x, 1) = make_matrix(info ? NULL : ev, n, n, n, make_transposed);
    else {
        assert(!info);
        qK(x, 1) = ktn(0, n);
        repeat (j, n) {
            if (li[j] != 0 && j+1 < n) {
                K q1 = qK(qK(x, 1), j)   = ktn(0, n);
                K q2 = qK(qK(x, 1), j+1) = ktn(0, n);
                repeat (i, n) {
                    F a = ev[i + j*n], b = ev[i + (j+1)*n];
                    qK(q1, i) = make_complex(a,  b);
                    qK(q2, i) = make_complex(a, -b);
                }
                j++;
            } else
                qK(qK(x, 1), j) = make_F(ev + j*n, n);
        }
    }

    free_F(ev);
    free_F(li);
    free_F(lr);

    return check_err(x, err);
}


// Cholesky decomposition
K
qml_mchol(K x) {
    I n, info;
    S err = NULL;

    F* a = take_square_matrix(x, &n, NULL, &err);

    dpotrf_("U", &n, a, &n, &info);
    check_info(info, &err);

    x = make_matrix(info ? NULL : a, n, n, n, make_upper);
    free_F(a);

    return check_err(x, err);
}


// QR factorization
static K
mqr(K x, const int pivot) {
    I n, m, info;
    S err = NULL;

    F* a = take_matrix(x, &m, &m, &n, alloc_square, &err);
    I min = min_i(m, n);

    I lwork_query = -1;
    F maxwork;
    if (pivot)
        dgeqp3_(&m, &n, NULL, &m, NULL, NULL, &maxwork, &lwork_query, &info);
    else
        dgeqrf_(&m, &n, NULL, &m, NULL,       &maxwork, &lwork_query, &info);
    check_info(info, &err);
    I lwork_ge = take_maxwork(info, maxwork);

    dorgqr_(&m, &m, &min, NULL, &m, NULL, &maxwork, &lwork_query, &info);
    check_info(info, &err);
    I lwork_or = take_maxwork(info, maxwork);

    alloc_W(max_i(lwork_ge, lwork_or));
    F* tau = alloc_F(&min, &err);
    if (!min)
        n = 0;

    K p;
    if (pivot) {
        I* ipiv = alloc_I(&n, &err);
        if (!n)
            min = 0;

        repeat (i, n)
            ipiv[i] = 0;
        dgeqp3_(&m, &n, a, &m, ipiv, tau, w, &lwork, &info);

        p = ktn(KL, n);
        repeat (i, n)
            qL(p, i) = ipiv[i] - 1;
        free_I(ipiv);
    } else
        if (lwork)
            dgeqrf_(&m, &n, a, &m,   tau, w, &lwork, &info);
    check_info(info, &err);

    K r = make_matrix(a, m, m, n, make_upper);

    dorgqr_(&m, &m, &min, a, &m, tau, w, &lwork, &info);
    check_info(info, &err);
    free_W();
    free_F(tau);

    x = ktn(0, 2);
    qK(x, 0) = make_matrix(a, m, m, m, 0);
    qK(x, 1) = r;
    if (pivot)
        jk(&x, p);
    free_F(a);

    return check_err(x, err);
}

K
qml_mqr(K x) {
    return mqr(x, 0);
}

K
qml_mqrp(K x) {
    return mqr(x, 1);
}


// LUP factorization
K
qml_mlup(K x) {
    I m, n, info;
    S err = NULL;

    F* a = take_matrix(x, &m, &m, &n, NULL, &err);

    I min = min_i(m, n);
    I* ipiv = alloc_I(&min, &err);
    if (!min)
        n = 0;

    dgetrf_(&m, &n, a, &m, ipiv, &info);
    check_info(info, &err);

    x = ktn(0, 3);
    qK(x, 0) = make_matrix(a, m, m, min, make_lower);
    qK(x, 1) = make_matrix(a, m, min, n, make_upper);
    L* q = kL(qK(x, 2) = ktn(KL, m));
    repeat (i, m)
        q[i] = i;
    repeat (i, min)
        swap_l(q + i, q + ipiv[i]-1);

    free_I(ipiv);
    free_F(a);

    return check_err(x, err);
}


// Singular value decomposition
K
qml_msvd(K x) {
    I m, n, info;
    S err = NULL;

    F* a = take_matrix(x, &m, &m, &n, NULL, &err);
    I min = min_i(m, n);

    I lwork_query = -1;
    F maxwork;
    dgesdd_("A", &m, &n, NULL, &m, NULL, NULL, &m, NULL, &n,
            &maxwork, &lwork_query, NULL, &info);
    check_info(info, &err);

    alloc_W(take_maxwork(info, maxwork));
    F* s  = alloc_F (&min,    &err);
    F* u  = alloc_FF(&m, m,   &err);
    F* vt = alloc_FF(&n, n,   &err);
    I* iw = alloc_II(&min, 8, &err);

    if (!min)
        m = 0;

    dgesdd_("A", &m, &n, a, &m, s, u, &m, vt, &n, w, &lwork, iw, &info);
    check_info(info, &err);
    free_W();
    free_I(iw);
    free_F(a);

    x = ktn(0, 3);
    qK(x, 0) = make_matrix(info ? NULL : u, m, m, m, 0);
    qK(x, 1) = ktn(0, m);
    repeat (j, m) {
        K q = qK(qK(x, 1), j) = ktn(KF, n);
        repeat (i, n)
            qF(q, i) = 0;
        if (j < n)
            qF(q, j) = info ? nf : s[j];
    }
    qK(x, 2) = make_matrix(info ? NULL : vt, n, n, n, make_transposed);
    free_F(vt);
    free_F(u);
    free_F(s);

    return check_err(x, err);
}


// Polynomial root finding
K
qml_poly(K x_) {
    I n, info;
    S err = NULL;

    K x = convert_F(x_);
    if (!x)
        return krr("type");
    if (!(n = qnw(x))) {
        q0(x);
        return krr("length");
    }
    if (n >= wi)
        if (!err) err = "limit";
    n--;

    F lc = qF(x, 0);
    if (lc == 0)
        if (!err) err = "roots";

    F* lr = alloc_F(&n, &err);
    F* li = alloc_F(&n, &err);
    if (!n)
        goto done;

    I lwork_query = -1;
    F maxwork;
    dgeev_("N", "N", &n, NULL, &n, NULL, NULL, NULL, &n, NULL, &n,
           &maxwork, &lwork_query, &info);
    check_info(info, &err);

    alloc_W(take_maxwork(info, maxwork));
    F* cm = alloc_FF(&n, n, &err);

    /* make companion matrix */
    repeat (i, n-1)
        repeat (j, n)
            cm[j + n*i] = j == i + 1;
    repeat (j, n)
        cm[j + n*(n-1)] = qF(x, n-j) / -lc;

    dgeev_("N", "N", &n, cm, &n, lr, li, NULL, &n, NULL, &n, w, &lwork, &info);
    check_info(info, &err);

    free_W();
    free_F(cm);

    if (info) {
        if (!err) err = "roots"; // "Failures are rare." - dhseqr.f
        n = 0;
    }

done:
    q0(x);

    x = make_complex_vector(lr, li, n);
    free_F(li);
    free_F(lr);

    return check_err(x, err);
}


// Linear equations
static K
mls(K x, K y, int equi) {
    int b_column;
    I a_n, b_m, b_n, info;
    S err = NULL;

    F* a = take_square_matrix(x, &a_n, NULL, &err);
    F* b = take_matrix(y, &b_m, &b_m, &b_n, &b_column, &err);
    if (a_n != b_m)
        if (!err) err = "length";

    I* ipiv = alloc_I(&a_n, &err);

    if (equi) {
        F* bi = b;
        b     = alloc_FF(&a_n, b_n, &err);

        F* af = alloc_FF(&a_n, a_n, &err);
        F* r  = alloc_F (&a_n,      &err);
        F* c  = alloc_F (&a_n,      &err);
        F* w  = alloc_FF(&a_n, 4,   &err);
        F* ferr = alloc_F(&b_n, &err);
        F* berr = alloc_F(&b_n, &err);

        I* iw = alloc_I(&a_n, &err);

        F rcond;
        char equed;

        dgesvx_("E", "N", &a_n, &b_n, a, &a_n, af, &a_n, ipiv, &equed, r, c,
                bi, &b_m, b, &b_m, &rcond, ferr, berr, w, iw, &info);
        free_I(iw);
        free_F(berr);
        free_F(ferr);
        free_F(w);
        free_F(c);
        free_F(r);
        free_F(af);
        free_F(bi);
    } else
        dgesv_(&a_n, &b_n, a, &a_n, ipiv, b, &b_m, &info);

    check_info(info, &err);
    free(ipiv);
    free(a);

    x = make_matrix(info ? NULL : b, b_m, b_m, b_n, b_column);
    free(b);
    return check_err(x, err);
}

static const struct optn ls_opt[] = {
    [0] = { "equi", 0 },
          { NULL }
};

K
qml_mlsx(K opts, K x, K y) {
    union optv v[] = { { 0 } };
    if (!take_opt(opts, ls_opt, v))
        return krr("opt");
    return mls(x, y, v[0].i);
}

K
qml_mls(K x, K y) {
    return mls(x, y, 0);
}


// Linear least squares
static K
mlsq(K x, K y, int svd) {
    int b_column;
    I a_m, a_n, b_m, b_n, info;
    S err = NULL;

    F* a = take_matrix(x, &a_m, &a_m, &a_n, NULL, &err);
    I ldb = a_n;
    F* b = take_matrix(y, &ldb, &b_m, &b_n, &b_column, &err);
    if (a_m != b_m)
        if (!err) err = "length";

    I lwork_query = -1;
    I liwork = wi; // bug0038 in LAPACK <3.2.2 does not initialize this
    I rank;
    F maxwork;
    F rcond = -1;
    if (svd)
        dgelsd_(&a_m, &a_n, &b_n, NULL, &a_m, NULL, &ldb,
                NULL, &rcond, &rank, &maxwork, &lwork_query, &liwork, &info);
    else
        dgels_("N", &a_m, &a_n, &b_n, NULL, &a_m, NULL, &ldb,
               &maxwork, &lwork_query, &info);
    check_info(info, &err);

    alloc_W(take_maxwork(info, maxwork));
    if (svd && !info) {
        I min = min_i(a_m, a_n);
        F* s = alloc_F(&min, &err);
        if (!min)
            a_n = 0;
        I* iw = alloc_I(&liwork, &err);

        // iw access pattern is non-trivial, so don't call if error
        if (!err)
            dgelsd_(&a_m, &a_n, &b_n, a, &a_m, b, &ldb,
                    s, &rcond, &rank, w, &lwork, iw, &info);
        free(iw);
        free(s);
    } else
        dgels_("N", &a_m, &a_n, &b_n, a, &a_m, b, &ldb,
               w, &lwork, &info);

    check_info(info, &err);
    free_W();
    free(a);

    x = make_matrix(info ? NULL : b, ldb, a_n, b_n, b_column);
    free(b);
    return check_err(x, err);
}

static const struct optn lsq_opt[] = {
    [0] = { "svd", 0 },
          { NULL }
};

K
qml_mlsqx(K opts, K x, K y) {
    union optv v[] = { { 0 } };
    if (!take_opt(opts, lsq_opt, v))
        return krr("opt");
    return mlsq(x, y, v[0].i);
}

K
qml_mlsq(K x, K y) {
    return mlsq(x, y, 0);
}
