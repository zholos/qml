#!/usr/bin/env python
from __future__ import division

import sympy as sp
import sympy.mpmath as mp
from sympy import Matrix, S, I

from qform import *


def hilbert_matrix(n, m):
    return Matrix(n, m, lambda i, j: S(1)/(1+i+j))

def random_matrix(n, m, scale=1):
    return Matrix(n, m, lambda i, j:
        sp.Rational(sp.tan((1+i*n+j)*scale).evalf(17)).limit_denominator(10))

subjects = [
    Matrix([[42]]),
    Matrix([[1, 2], [-3, 4]]),
    Matrix([[5, -6]]),
    Matrix([[7], [8]]),
    Matrix([[1, 2, 3], [4, -5, 6], [7, 8, 9]]),
    Matrix([[1, 2, 3], [4, -5, 6]]),
    Matrix([[1, 2], [4, -5], [7, 8]]),
    hilbert_matrix(5, 5),
    hilbert_matrix(3, 5),
    hilbert_matrix(5, 3),
    random_matrix(7, 7),
    random_matrix(7, 4, sp.sqrt(2)),
    random_matrix(4, 7, sp.sqrt(3)),
    Matrix([[1, 3], [2, 6]]),
    Matrix([[1, 3, 4], [2, 6, 8], [-3, -9, -12]]),
    Matrix([[1, 3, 4], [2, 6, 8], [-3, -9, 12]]),
    Matrix([[1,0,0],[1,1,0],[1,0,1],[0,1,1],[0,0,1]]) * hilbert_matrix(3, 5),
    random_matrix(7, 3, sp.sqrt(2)) * sp.eye(3).row_join(sp.ones(3, 1)),
    Matrix([[1,0],[0,1],[1,1],[1,-1]]) * random_matrix(2, 7, sp.sqrt(3))
]

large_subjects = [
    "{til[x 1]+/:1+til x 0} 100 100",
    "{til[x 1]+/:1+til x 0} 20 100",
    "{til[x 1]+/:1+til x 0} 100 20",
    "{y#tan x*1+til prd y}[1] 100 100",
    "{y#tan x*1+til prd y}[sqrt 2] 20 100",
    "{y#tan x*1+til prd y}[sqrt 3] 100 20"
]

eigenvalue_subjects = [
    Matrix([[42]]),
    Matrix([[1, 2], [-3, 4]]),
    Matrix([[1, 2], [3, 4]]),
    Matrix([[5, 6], [0, 5]]),
    Matrix([[5, 0], [0, 5]]),
    Matrix([[1, 3], [2, 6]]),
    Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
    Matrix([[1, 2, 3], [0, 5, 6], [7, 8, 9]]),
    Matrix([[1, 2, 3], [0, 5, -6], [-7, 8, -9]]),
    Matrix([[3, 0, 2], [0, 5, 0], [1, 0, 4]]),
    Matrix([[1, 2, 0, 0], [3, 4, 0, 0], [0, 0, 1, 2], [0, 0, 3, 4]]),
    hilbert_matrix(5, 5),
    random_matrix(7, 7)
]

cholesky_subjects = [
    Matrix([[42]]),
    Matrix([[1, 2], [2, 5]]),
    Matrix([[1, 1, 3], [1, 5, 13], [3, 13, 98]]),
    Matrix([[3, -2, 3], [-2, 7, -1], [3, -1, 5]]),
    hilbert_matrix(5, 5),
    Matrix([[1, 2], [2, -5]]),
    Matrix([[1, -2], [-2, 4]]),
    Matrix([[1, 2], [-3, 4]])
]

def N(A):
    for r in A.atoms(sp.Pow):
        if r.exp == S(1)/2 and r.base.is_Rational and r.base > 1000000000:
            return mp.matrix(A.evalf(mp.mp.dps))
    return A

def test_mev():
    output("""\
    reim:{$[0>type x;1 0*x;2=count x;x;'`]};
    mc:{((x[0]*y 0)-x[1]*y 1;(x[0]*y 1)+x[1]*y 0)};
    mmc:{((.qml.mm[x 0]y 0)-.qml.mm[x 1]y 1;(.qml.mm[x 0]y 1)+.qml.mm[x 1]y 0)};
    mev_:{[b;x]
        if[2<>count wv:.qml.mev x;'`length];
        if[not all over prec>=abs
            mmc[flip vc;flip(flip')(reim'')flip x]-
            flip(w:reim'[wv 0])mc'vc:(flip')(reim'')(v:wv 1);'`check];
        / Normalize sign; LAPACK already normalized to real
        v*:1-2*0>{x a?max a:abs x}each vc[;0];
        (?'[prec>=abs w[;1];w[;0];w];?'[b;v;0n])};""")

    for A in eigenvalue_subjects:
        if A.rows <= 3:
            V = []
            for w, n, r in A.eigenvects():
                w = sp.simplify(sp.expand_complex(w))
                if len(r) == 1:
                    r = r[0]
                    r = sp.simplify(sp.expand_complex(r))
                    r = r.normalized() / sp.sign(max(r, key=abs))
                    r = sp.simplify(sp.expand_complex(r))
                else:
                    r = None
                V.extend([(w, r)]*n)
            V.sort(key=lambda (x, _): (-abs(x), -sp.im(x)))
        else:
            Am = mp.matrix(A)
            # extra precision for complex pairs to be equal in sort
            with mp.extradps(mp.mp.dps):
                W, R = mp.eig(Am)
            V = []
            for w, r in zip(W, (R.column(i) for i in range(R.cols))):
                w = mp.chop(w)
                with mp.extradps(mp.mp.dps):
                    _, S, _ = mp.svd(Am - w*mp.eye(A.rows))
                if sum(x == 0 for x in mp.chop(S)) == 1:
                    # nullity 1, so normalized eigenvector is unique
                    r /= mp.norm(r) * mp.sign(max(r, key=abs))
                    r = mp.chop(r)
                else:
                    r = None
                V.append((w, r))
            V.sort(key=lambda (x, _): (-abs(x), -x.imag))
        W, R = zip(*V)
        test("mev_[%sb" % "".join("0" if r is None else "1" for r in R), A,
             (W, [r if r is None else list(r) for r in R]), complex_pair=True)

def test_mchol():
    for A in cholesky_subjects:
        if A.is_hermitian and all(x > 0 for x in A.berkowitz_minors()):
            R = A.cholesky().T
        else:
            R = A * sp.nan
        test("mchol", A, R)

def qrp(A):
    Q = sp.zeros(A.rows, 0)
    P = []
    J = list(range(A.cols))
    while J:
        U = A - Q*Q.T*A
        j = max(J, key=lambda j: (U.col(j).norm(), -j))
        J.remove(j)
        P.append(j)
        a = U.col(j)
        if a.norm():
            Q = Q.row_join(a.normalized())
    R = Q.T * A.extract(range(A.rows), P)
    assert Q.T*Q == sp.eye(Q.cols) and R.is_upper
    assert A.extract(range(A.rows), P) == Q * R
    return Q, R, P

def test_mqr(pivot=False):
    if not pivot:
        output("""\
    mzero:{all all each prec>=0w^abs x};
    mortho:{
        mzero[.qml.diag[count[m]#1.]-m:.qml.mm[x] flip x] and
        mzero .qml.diag[count[n]#1.]-n:.qml.mm[flip x] x};
    mqr_:{[k;x]
        $[2<>count qr:.qml.mqr x;::;
          not (.qml.mdim[q:qr 0]~2#d 0) and .qml.mdim[r:qr 1]~d:.qml.mdim x;::;
          not mortho[q] and mzero (d[1]&til d 0)#'r;::;
          not mzero x-.qml.mm[q] r;::;
          null k;1b;(s*/:k#/:q;(s:1-2*0>.qml.mdiag k#r)*k#k#/:r)]};""")
    else:
        output("""\
    mqrp_:{[k;x]
        $[3<>count qrp:.qml.mqrp x;::;
          not (.qml.mdim[q:qrp 0]~2#d 0) and (.qml.mdim[r:qrp 1]~d)
            and asc[p:qrp 2]~til last d:.qml.mdim x;::;
          not mortho[q] and mzero (d[1]&til d 0)#'r;::;
          not mzero (x@\:p)-.qml.mm[q] r;::;
          null k;1b;(s*/:k#/:q;(s:1-2*0>.qml.mdiag k#r)*k#k#/:r;k#p)]};""")

    for A in subjects:
        rank = A.rank()
        if not pivot:
            if A.rows == A.cols == rank:
                Q, R = A.QRdecomposition()
                test("mqr_[%s" % rank, A, (N(Q), N(R)))
            else:
                test("mqr_[0N", A, qstr("1b"))
        else:
            Q, R, P = qrp(A)
            R = R.extract(range(R.rows), range(rank))
            P = P[:rank]
            test("mqrp_[%s" % rank, A, (N(Q), N(R), P))

    reps(250)
    for Aq in large_subjects:
        test("%s[0N" % ("mqr_" if not pivot else "mqrp_"), qstr(Aq), qstr("1b"))
    reps(10000)

def test_mlup():
    output("""\
    mlup_b:{
        $[3<>count lup:.qml.mlup x;::;
          not (.qml.mdim[l:lup 0]~mins d) and (.qml.mdim[u:lup 1]~(min d;d 1))
            and asc[p:lup 2]~til first d:.qml.mdim x;::;
          not mzero[(i+1)_'l] and mzero[til[min d]#'u]
            and all[prec>=abs 1-l@'i:til d 0] and mzero x[p]-.qml.mm[l] u;::;
          1]};""")

    for A in subjects:
        test("mlup_b", A, qstr("1b"))

    reps(250)
    for Aq in large_subjects:
        test("mlup_b", qstr(Aq), qstr("1b"))
    reps(10000)

def test_msvd():
    output("""\
    msvd_:{[b;x]
        $[3<>count usv:.qml.msvd x;::;
          not (.qml.mdim[u:usv 0]~2#d 0) and (.qml.mdim[v:usv 2]~2#d 1) and
            .qml.mdim[s:usv 1]~d:.qml.mdim x;::;
          not mortho[u] and mortho[v] and
            all[0<=f:.qml.mdiag s] and mzero s _'til d 0;::;
          not mzero x-.qml.mm[u] .qml.mm[s] flip v;::;
          b;1b;(p*/:m#/:u;m#m#/:s;
              (p:(f>.qml.eps*f[0]*max d)*1-2*0>m#u 0)*/:(m:min d)#/:v)]};""")

    for A in subjects:
        U, S, V = map(mp.matrix, mp.svd(mp.matrix(A)))
        m = min(A.rows, A.cols)
        p = mp.diag([mp.sign(s * u) for s, u in zip(mp.chop(S), U[0, :])])
        U *= p
        S = mp.diag(S)
        V = V.T * p
        test("msvd_[0b", A, (U, S, V))

    reps(250)
    for Aq in large_subjects:
        test("msvd_[1b", qstr(Aq), qstr("1b"))
    reps(10000)


poly_subjects = [
    [1],
    [S(1)/2, 21],
    [1, 3, 2],
    [1, 2, 3],
    [5, 0, 0, 0, 0, 1],
    [10, 0, 0, 0, 75, 420],
    (lambda n: [sp.binomial(n, i) for i in range(n+1)])(2),
    list(random_matrix(25, 1)),
    [1, -1, 1+I],
    [sp.sqrt(-2), sp.sqrt(2)-sp.sqrt(-2),2+I,-3*I],
    list(random_matrix(25, 2)*Matrix([1, I]))
]

def test_poly():
    output("""\
    poly_:{
        r:{$[prec>=abs(x:reim x)1;x 0;x]} each .qml.poly x;
        r iasc {x:reim x;(abs x 1;neg x 1;x 0)} each r};""")

    for A in poly_subjects:
        A = sp.Poly(A, sp.Symbol("x"))
        if A.degree() <= 3 and all(a.is_real for a in A.all_coeffs()):
            R = []
            for r in sp.roots(A, multiple=True):
                r = sp.simplify(sp.expand_complex(r))
                R.append(r)
            R.sort(key=lambda x: (abs(sp.im(x)), -sp.im(x), sp.re(x)))
        else:
            R = mp.polyroots([mp.mpc(*(a.evalf(mp.mp.dps)).as_real_imag())
                              for a in A.all_coeffs()])
            R.sort(key=lambda x: (abs(x.imag), -x.imag, x.real))
        assert len(R) == A.degree()
        test("poly_", A.all_coeffs(), R, complex_pair=True)


def tests():
    prec("1e-9")
    test_mev()
    test_mchol()
    test_mqr(False)
    test_mqr(True)
    test_mlup()
    test_msvd()
    prec("1e-8")
    test_poly()


if __name__ == "__main__":
    tests()
