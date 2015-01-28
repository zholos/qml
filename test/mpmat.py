#!/usr/bin/env python
from __future__ import division

import sympy as sp
import sympy.mpmath as mp
from sympy import Matrix, S

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


def tests():
    test_msvd()


if __name__ == "__main__":
    tests()
