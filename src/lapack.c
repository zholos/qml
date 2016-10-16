#include <stddef.h>
#include <stdlib.h>

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


// column argument values
// note: some functions do arithmetic on these, e.g. !!column to distinguish
// make_column from make_general
enum {
    make_upper = -1, // upper triangular
    make_lower = -2, // lower triangular with 1 on diagonal
    make_general = 0,
    make_column = 1
};

// check if x is a q matrix, return number of rows (m) and columns (n)
// returned x must be passed to copy_matrix or released with if (x) q0(x);
// if check fails, err is set and m=0 - so copy_matrix() may still be called
// if column is set, accepts vector as m*1 matrix
static K
check_matrix(K x_, I* m, I* n, int* column, S* err) {
    if (column) {
        K x = convert_F(x_);
        if (x) {
            *n = 1;
            if (!likely(*m = qnw(x)))
                if (!*err) *err = "type";
            *column = make_column;
            return x;
        }
        *column = make_general;
    }

    K x = convert_FF(x_);
    if (!likely(x && (*m = qnw(x)) && (*n = qnw(qK(x, 0))))) {
        *m = *n = 0;
        if (!*err) *err = "type";
    }
    repeati (j, *m)
        if (qn(qK(x, j)) != *n) {
            *m = *n = 0;
            if (!*err) *err = "length";
        }
    return x;
}

// copy data from output of check_matrix into pre-allocated array
// flip=0: convert row-major x to column-major a, transposing data in memory
// flip=1: don't transpose data in memory, which is faster but flips matrix
// lda, m and n refer to output array shape regardless of flip
// column=1: requires (flip?m:n)=1, so for flip=1 swap m and n from check_matrix
static void
copy_matrix(K x, F* a, I lda, I m, I n, int column, int flip)
{
    assert(a || !m || !n);
    if (column && x) {
        if (flip)
            swap_i(&m, &n);
        assert(n <= 1);
        copy_F(a, 0, qrF(x, m), 0, m);
    } else
        repeati (i, n)
            if (!flip)
                repeati (j, m)
                    a[j + i*lda] = qF(qK(x, j), i);
            else
                copy_F(a, i*lda, qrF(qK(x, i), m), 0, m);
    if (x) q0(x);
}


// if triangular is set, returns 1 if matrix is upper, -1 if lower, else 0
// upper and lower refer to output matrix shape regardless of flip
// triangular without flip isn't used so isn't implemented
static F*
take_square_matrix(K x, I* n, int* triangular, int flip, S* err) {
    assert(flip || !triangular);
    I m;
    x = check_matrix(x, &m, n, NULL, err);
    if (!likely(m == *n)) {
        *n = 0;
        if (!*err) *err = "length";
    }
    // don't use m past here

    F* a = alloc_FF(n, *n, err);

    if (triangular) {
        int upper = 1, lower = 1;
        repeati (i, *n)
            repeati (j, *n)
                if ((a[j + i * *n] = qF(qK(x, i), j)))
                    if (i < j) // non-zero below diagonal, so not upper
                        upper = 0;
                    else if (i > j) // non-zero above diagonal, so not lower
                        lower = 0;
        *triangular = upper ? 1 : lower ? -1 : 0;
        if (x) q0(x);
    } else
        copy_matrix(x, a, *n, *n, *n, 0, flip);

    return a;
}


// m and n refer to output array shape regardless of flip
// if column is set and true, (flip?m:n)=1
// on error a = NULL and m = n = 0; the latter protects from functions here and
// in LAPACK (e.g. dgeqrf) accessing a based on one of m and n != 0
static F*
take_matrix(K x, I* m, I* n, int* column, int flip, S* err) {
    x = check_matrix(x, m, n, column, err);
    F* a = alloc_FF(m, *n, err); // works for vector too (*n==1)
    if (!*m)
        *n = 0;
    if (flip)
        swap_i(m, n);
    copy_matrix(x, a, *m, *m, *n, column && *column, flip);
    // note: if this function ever rounds up lda, keep lda=1 for m=1 so that
    // a column or row matrix can be flipped by swapping dimensions
    return a;
}


// if a is null matrix is filled with NaNs
// column=make_upper/make_lower refers to output matrix layout (after flip)
static K
make_matrix(const F* a, I lda, I m, I n, int column, int flip) {
    if (flip)
        swap_i(&m, &n);
    // m and n now refer to x shape

    if (!likely(a))
        if (column > 0) // may also be make_upper or make_lower here
            // n == 1, unless result is being discarded ('length in qml_mm)
            return make_F_null(m);
        else {
            K r = ktn(0, m);
            K x = make_F_null(n);
            repeati (j, m)
                qK(r, j) = r1(x);
            q0(x);
            return r;
        }

    if (!column) { // most common case
        K r = ktn(0, m);
        repeati (j, m)
            if (!flip) {
                // TODO: cache-optimized transpose
                K x = ktn(KF, n);
                repeati (i, n)
                    qF(x, i) = a[j + i*lda];
                qK(r, j) = x;
            } else
                qK(r, j) = make_F(a, j*lda, n);
        return r;
    }

    if (column == make_column) {
        assert(n <= 1);
        return make_F(a, 0, m);
    }

    assert(column == make_upper || column == make_lower);

    K r = ktn(0, m);
    repeati (j, m) {
        K x = ktn(KF, n);
        if (column == make_upper) {
            repeati (i, min_i(j, n))
                qF(x, i) = 0;
            if (!flip)
                repeati_ (i, j, n)
                    qF(x, i) = a[j + lda*i];
            else
                copy_F(qrF(x, n), j, a, j + j*lda, max_i(0, n-j));
        } else {
            if (!flip)
                repeati (i, min_i(j, n))
                    qF(x, i) = a[j + lda*i];
            else
                copy_F(qrF(x, n), 0, a, j*lda, min_i(j, n));
            if (j < n)
                qF(x, j) = 1;
            repeati_ (i, j+1, n)
                qF(x, i) = 0;
        }
        qK(r, j) = x;
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
make_complex_vector(const F* a, const F* b, I n) {
    if (b)
        repeati (i, n)
            if (!(b[i] == 0))
                goto complex;

    if (a) {
        K x = ktn(KF, n);
        repeati (i, n)
            qF(x, i) = a[i];
        return x;
    } else
        return make_F_null(n);

complex:;
    K x = ktn(0, n);
    repeati (i, n)
        qK(x, i) = make_complex(a[i], b[i]);
    return x;
}



// Matrix determinant
K
qml_mdet(K x) {
    int triangular;
    I n, info;
    S err = NULL;

    // |A|=|A'|, so flipping input doesn't affect output
    F* a = take_square_matrix(x, &n, &triangular, 1, &err);
    F r = 1;

    if (!triangular) {
        I* ipiv = alloc_I(&n, &err);

        dgetrf_(&n, &n, a, &n, ipiv, &info);
        check_info(info, &err);
        // info>0 indicates singularity, but we don't treat this specially

        repeati (i, n)
            if (ipiv[i]-1 != i)
                r = -r;
        free_I(ipiv);
    }

    repeati (i, n)
        r *= a[i + i*n];
    free_F(a);

    return check_err(kf(r), err);
}


// Matrix inverse
K
qml_minv(K x) {
    I n, info;
    S err = NULL;

    // inv(A)'=inv(A'), so flip both input and output
    F* a = take_square_matrix(x, &n, NULL, 1, &err);
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

    x = make_matrix(info ? NULL : a, n, n, n, make_general, 1);
    free_F(a);

    return check_err(x, err);
}


// Matrix multiply
static K
mm(K x, K y, int lflip, int rflip) {
    int b_column = 0;
    I a_m, a_n, b_m, b_n;
    S err = NULL;

    // (AB)'=B'A', so flip inputs, multiply in reverse order, and flip output
    // rflip can't accept a column because result wouldn't be a column
    F* a = take_matrix(x, &a_m, &a_n, NULL, 1, &err);
    F* b = take_matrix(y, &b_m, &b_n, rflip ? NULL : &b_column, 1, &err);
    // option flips are applied to this: r(m * n) = b(b_m * b_n) a(a_m * a_n)
    I m = rflip ? b_n : b_m;
    I n = lflip ? a_m : a_n;
    I k = lflip ? a_n : a_m;
    if (k != (rflip ? b_m :b_n))
        if (!err) err = "length";

    F* r = alloc_FF(&m, n, &err);

    // BLAS doesn't check for zero size like LAPACK, so skip call on error
    if (!err) {
        int i_1 = 1;
        double f_0 = 0, f_1 = 1;
        if (m == 1)
            dgemv_(lflip?"N":"T", &a_m, &a_n, &f_1, a, &a_m, b, &i_1,
                   &f_0, r, &i_1);
        else if (n == 1)
            dgemv_(rflip?"T":"N", &b_m, &b_n, &f_1, b, &b_m, a, &i_1,
                   &f_0, r, &i_1);
        else
            dgemm_(rflip?"T":"N", lflip?"T":"N", &m, &n, &k,
                   &f_1, b, &b_m, a, &a_m, &f_0, r, &m);
    }

    free_F(b);
    free_F(a);
    x = make_matrix(err ? NULL : r, m, m, n, b_column, 1);
    free_F(r);

    return check_err(x, err);
}

static const struct optn mm_opt[] = {
    [0] = { "lflip", 0 },
    [1] = { "rflip", 0 },
          { NULL }
};

K
qml_mmx(K opts, K x, K y) {
    union optv v[] = { { 0 }, { 0 } };
    if (!take_opt(opts, mm_opt, v))
        return krr("opt");
    return mm(x, y, v[0].i, v[1].i);
}

K
qml_mm(K x, K y) {
    return mm(x, y, 0, 0);
}


// Matrix substitution solve
K
qml_ms(K x, K y) {
    int a_triangular, b_column;
    I a_n, b_m, b_n;
    S err = NULL;

    // AX=B <=> X'A'=B', so flip inputs and outputs and specify this equation
    // type (A on right side) in dtrsm
    F* a = take_square_matrix(x, &a_n, &a_triangular, 1, &err);
    if (!a_triangular)
        if (!err) err = "domain";

    F* b = take_matrix(y, &b_m, &b_n, &b_column, 1, &err);
    if (a_n != b_n)
        if (!err) err = "length";

    I info;
    for (info = a_n; info; info--)
        if (!(a[(info-1) + (info-1)*a_n] != 0)) // consider NaN
            break; // info > 0 indicates singularity

    if (!err && !info) {
        int i_1 = 1;
        double f_1 = 1;
        if (b_m == 1)
            dtrsv_(a_triangular > 0 ? "U" : "L", "T", "N",
                   &a_n, a, &a_n, b, &i_1); 
        else
            dtrsm_("R", a_triangular > 0 ? "U" : "L", "N", "N",
                   &b_m, &b_n, &f_1, a, &a_n, b, &b_m);
    }

    free_F(a);
    x = make_matrix(info ? NULL : b, b_m, b_m, b_n, b_column, 1);
    free_F(b);

    return check_err(x, err);
}


// Eigenvalues and eigenvectors
K
qml_mevu(K x) {
    I n, info;
    S err = NULL;

    // A v = lambda v <=> v* A' = lambda* v*
    // Therefore, flip input, query left eigenvectors instead of right (which
    // returns v, not v*), and apply conjugate transpose to eigenvalues. Then
    // swap complex conjugate pairs of lambda and corresponding pairs of v to
    // restore original order (positive imaginary part first). Combined, this is
    // a no-op for lambda, leaving a swap for pairs of v.
    F* a = take_square_matrix(x, &n, NULL, 1, &err);

    I lwork_query = -1;
    F maxwork;
    dgeev_("V", "N", &n, NULL, &n, NULL, NULL,
           NULL, &n, NULL, &n, &maxwork, &lwork_query, &info);
    check_info(info, &err);

    alloc_W(take_maxwork(info, maxwork));
    F* lr = alloc_F (&n,    &err);
    F* li = alloc_F (&n,    &err);
    F* ev = alloc_FF(&n, n, &err);

    info = -1;
    dgeev_("V", "N", &n, a, &n, lr, li, ev, &n, NULL, &n, w, &lwork, &info);
    check_info(info, &err);
    free_W();
    free_F(a);

    x = ktn(0, 2);
    qK(x, 0) = make_complex_vector(info ? NULL : lr, info ? NULL : li, n);

    if (qt(qK(x, 0))) // no complex elements (also when info)
        qK(x, 1) = make_matrix(info ? NULL : ev, n, n, n, make_general, 1);
    else {
        assert(!info);
        qK(x, 1) = ktn(0, n);
        repeati (j, n) {
            if (li[j] != 0 && j+1 < n) {
                K q1 = qK(qK(x, 1), j)   = ktn(0, n);
                K q2 = qK(qK(x, 1), j+1) = ktn(0, n);
                repeati (i, n) {
                    F a = ev[i + j*n], b = ev[i + (j+1)*n];
                    qK(q1, i) = make_complex(a, -b); // swap pair like lambda
                    qK(q2, i) = make_complex(a,  b);
                }
                j++;
            } else
                qK(qK(x, 1), j) = make_F(ev, j*n, n);
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

    // Operation requires A=A', so flipping input doesn't affect result.
    // For output, find lower factor and flip it to get expected upper factor.
    F* a = take_square_matrix(x, &n, NULL, 1, &err);

    // with "L" uses lower triangular part of a, which is upper triangular part
    // of input
    dpotrf_("L", &n, a, &n, &info);
    check_info(info, &err);

    x = make_matrix(info ? NULL : a, n, n, n, make_upper, 1);
    free_F(a);

    return check_err(x, err);
}


// QR factorization
static K
mqr(K x, const int pivot) {
    I n, m, info;
    S err = NULL;

    // allocate m-by-max(m,n) matrix to fit m-by-m Q later
    x = check_matrix(x, &m, &n, NULL, &err);
    I nm = max_i(m, n);
    F* a = alloc_FF(&nm, m, &err);
    if (!nm)
        m = n = 0;
    copy_matrix(x, a, m, m, n, 0, 0);

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
        m = n = 0;

    K p;
    if (pivot) {
        I* ipiv = alloc_I(&n, &err);
        if (!n)
            min = m = 0;

        repeati (i, n)
            ipiv[i] = 0;
        dgeqp3_(&m, &n, a, &m, ipiv, tau, w, &lwork, &info);

        p = ktn(KL, n);
        repeati (i, n)
            qL(p, i) = ipiv[i] - 1;
        free_I(ipiv);
    } else
        // This could be replaced by dgelqf (LQ factorization) to avoid flips,
        // but then there would be less common code with pivoting version.
        if (lwork)
            dgeqrf_(&m, &n, a, &m,   tau, w, &lwork, &info);
    check_info(info, &err);

    K r = make_matrix(a, m, m, n, make_upper, 0);

    dorgqr_(&m, &m, &min, a, &m, tau, w, &lwork, &info);
    check_info(info, &err);
    free_W();
    free_F(tau);

    x = ktn(0, 2);
    qK(x, 0) = make_matrix(a, m, m, m, make_general, 0);
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

    F* a = take_matrix(x, &m, &n, NULL, 0, &err);

    I min = min_i(m, n);
    I* ipiv = alloc_I(&min, &err);
    if (!min)
        m = n = 0;

    dgetrf_(&m, &n, a, &m, ipiv, &info);
    check_info(info, &err);

    x = ktn(0, 3);
    qK(x, 0) = make_matrix(a, m, m, min, make_lower, 0);
    qK(x, 1) = make_matrix(a, m, min, n, make_upper, 0);
    L* q = kL(qK(x, 2) = ktn(KL, m));
    repeati (i, m)
        q[i] = i;
    repeati (i, min)
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

    // (U Sigma V)'=V'Sigma'U', so flip input and reverse and flip outputs
    F* a = take_matrix(x, &m, &n, NULL, 1, &err);
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
        m = n = 0;

    dgesdd_("A", &m, &n, a, &m, s, u, &m, vt, &n, w, &lwork, iw, &info);
    check_info(info, &err);
    free_W();
    free_I(iw);
    free_F(a);

    x = ktn(0, 3);
    qK(x, 0) = make_matrix(info ? NULL : vt, n, n, n, make_general, 1);
    qK(x, 1) = ktn(0, n);
    repeati (i, n) {
        K q = qK(qK(x, 1), i) = ktn(KF, m);
        repeati (j, m)
            qF(q, j) = 0;
        if (i < m)
            qF(q, i) = info ? nf : s[i];
    }
    qK(x, 2) = make_matrix(info ? NULL : u, m, m, m, make_general, 0);
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
    repeati (i, n-1)
        repeati (j, n)
            cm[j + n*i] = j == i + 1;
    repeati (j, n)
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
mls(K x, K y, int equi, int flip) {
    int b_column;
    I a_n, b_m, b_n, info;
    S err = NULL;

    // with equi always flip A because dgesvx can accept transposed data
    F* a = take_square_matrix(x, &a_n, NULL, flip || equi, &err);
    F* b = take_matrix(y, &b_m, &b_n, &b_column, flip, &err);
    if (flip && b_column)
        swap_i(&b_m, &b_n); // turn row back into column and don't flip result
    if (a_n != b_m)
        if (!err) err = "length";

    I* ipiv = alloc_I(&a_n, &err);

    if (equi) {
        F* bi = b;
        b     = alloc_FF(&b_m, b_n, &err);

        F* af = alloc_FF(&a_n, a_n, &err);
        F* r  = alloc_F (&a_n,      &err);
        F* c  = alloc_F (&a_n,      &err);
        F* w  = alloc_FF(&a_n, 4,   &err);
        F* ferr = alloc_F(&b_n, &err);
        F* berr = alloc_F(&b_n, &err);

        I* iw = alloc_I(&a_n, &err);

        F rcond;
        char equed;

        dgesvx_("E", flip?"N":"T", &a_n, &b_n, a, &a_n, af, &a_n, ipiv, &equed,
                r, c, bi, &b_m, b, &b_m, &rcond, ferr, berr, w, iw, &info);
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

    x = make_matrix(info ? NULL : b, b_m, b_m, b_n, b_column,
                    flip && !b_column);
    free(b);
    return check_err(x, err);
}

static const struct optn ls_opt[] = {
    [0] = { "equi", 0 },
    [1] = { "flip", 0 },
          { NULL }
};

K
qml_mlsx(K opts, K x, K y) {
    union optv v[] = { { 0 }, { 0 } };
    if (!take_opt(opts, ls_opt, v))
        return krr("opt");
    return mls(x, y, v[0].i, v[1].i);
}

K
qml_mls(K x, K y) {
    return mls(x, y, 0, 0);
}


// Linear least squares
static K
mlsq(K x, K y, int svd, int flip) {
    int b_column;
    I a_m, a_n, b_m, b_n, info;
    S err = NULL;

    // without svd always flip A because dgels can accept transposed data
    F* a = take_matrix(x, &a_m, &a_n, NULL, flip || !svd, &err);
    I eqns = !svd && !flip ? a_n : a_m;
    I vars = !svd && !flip ? a_m : a_n;

    // b has a_m rows in input but a_n rows in output, so allocate larger array
    y = check_matrix(y, &b_m, &b_n, &b_column, &err);
    if (flip && !b_column)
        swap_i(&b_m, &b_n); // flip right matrix and corresponding result
        // note: compared to mls() the swap goes in the other branch because
        // there take_matrix() is called and does another swap
    I ldb = max_i(vars, b_m);
    F* b = alloc_FF(&ldb, b_n, &err);
    if (!ldb)
        b_n = b_m = 0;
    copy_matrix(y, b, ldb, b_m, b_n, b_column, flip && !b_column);

    if (eqns != b_m) {
        // avoid accessing uninitialized part of b
        ldb = eqns = vars = a_m = a_n = b_m = 0;
        if (!err) err = "length";
    }

    I lwork_query = -1;
    I liwork = wi; // bug0038 in LAPACK <3.2.2 does not initialize this
    I rank;
    F maxwork;
    F rcond = -1;
    if (svd)
        dgelsd_(&a_m, &a_n, &b_n, NULL, &a_m, NULL, &ldb,
                NULL, &rcond, &rank, &maxwork, &lwork_query, &liwork, &info);
    else
        dgels_(flip?"N":"T", &a_m, &a_n, &b_n, NULL, &a_m, NULL, &ldb,
               &maxwork, &lwork_query, &info);
    check_info(info, &err);

    alloc_W(take_maxwork(info, maxwork));
    if (svd && !info) {
        I min = min_i(a_m, a_n);
        F* s = alloc_F(&min, &err);
        if (!min)
            a_m = a_n = 0;
        I* iw = alloc_I(&liwork, &err);

        // can't convey failed alloc of iw (liwork=0), so simply skip the call
        if (liwork)
            dgelsd_(&a_m, &a_n, &b_n, a, &a_m, b, &ldb,
                    s, &rcond, &rank, w, &lwork, iw, &info);
        free(iw);
        free(s);
    } else
        dgels_(flip?"N":"T", &a_m, &a_n, &b_n, a, &a_m, b, &ldb,
               w, &lwork, &info);

    check_info(info, &err);
    free_W();
    free(a);

    x = make_matrix(info ? NULL : b, ldb, vars, b_n, b_column,
                    flip && !b_column);
    free(b);
    return check_err(x, err);
}

static const struct optn lsq_opt[] = {
    [0] = { "svd",  0 },
    [1] = { "flip", 0 },
          { NULL }
};

K
qml_mlsqx(K opts, K x, K y) {
    union optv v[] = { { 0 }, { 0 } };
    if (!take_opt(opts, lsq_opt, v))
        return krr("opt");
    return mlsq(x, y, v[0].i, v[1].i);
}

K
qml_mlsq(K x, K y) {
    return mlsq(x, y, 0, 0);
}


// undocumented no-op function to benchmark matrix format conversion
static K
mnoop(K x, int square, int triangular_, int mark, int upper, int lower,
      int flip)
{
    int triangular = 0;
    int column = upper ? make_upper : lower ? make_lower : 0;
    I m, n;
    S err = NULL;
    F* a;
    if (square) {
        a = take_square_matrix(
            x, &n, triangular_ ? &triangular : NULL, flip, &err);
        m = n;
    } else
        a = take_matrix(x, &m, &n, column ? NULL : &column, flip, &err);

    if (mark)
        repeati (i, n)
            repeati (j, m) {
                F* v = &a[j+i*m];
                *v = *v*10000 + ((j+1)*100 + (i+1))*(*v<0?-1:1);
            }

    x = make_matrix(a, m, m, n, column, flip);
    free_F(a);
    if (triangular_) {
        K d = new_D();
        append_D(d, "x", x);
        append_D(d, "triangular", ks(triangular > 0 ? "upper" :
                                     triangular < 0 ? "lower" : ""));
        x = d;
    }
    return check_err(x, err);
}

static const struct optn noop_opt[] = {
    [0] = { "square", 0 },
    [1] = { "triangular", 0 },
    [2] = { "mark", 0 },
    [3] = { "upper", 0 },
    [4] = { "lower", 0 },
    [5] = { "flip", 0 },
          { NULL }
};

K
qml_mnoopx(K opts, K x) {
    union optv v[] = { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
    if (!take_opt(opts, noop_opt, v))
        return krr("opt");
    if (v[1].i && (!v[0].i || !v[5].i) || v[3].i && v[4].i)
        return krr("opt");
    return mnoop(x, v[0].i, v[1].i, v[2].i, v[3].i, v[4].i, v[5].i);
}

K
qml_mnoop(K x) {
    return mnoop(x, 0, 0, 0, 0, 0, 0);
}
