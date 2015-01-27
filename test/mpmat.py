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
    Matrix([[1, 3], [2, 6]])
]

def test_msvd():
    output("""\
    msvd_:{
        $[3<>count usv:.qml.msvd x;::;
          not (.qml.mdim[u:usv 0]~2#d 0) and (.qml.mdim[v:usv 2]~2#d 1) and
            .qml.mdim[s:usv 1]~d:.qml.mdim x;::;
          not mortho[u] and mortho[v] and
            all[0<=.qml.mdiag s] and mzero s _'til d 0;::;
          not mzero x-.qml.mm[u] .qml.mm[s] flip v;::;
          (p*/:m#/:u;m#m#/:s;(p:1-2*0>m#u 0)*/:(m:min d)#/:v)]};""")

    for A in subjects:
        U, S, V = map(mp.matrix, mp.svd(mp.matrix(A)))
        m = min(A.rows, A.cols)
        p = mp.diag([mp.sign(u) for u in U[0, :]])
        U *= p
        S = mp.diag(S)
        V = V.T * p
        test("msvd_", A, (U, S, V))


def tests():
    test_msvd()


if __name__ == "__main__":
    tests()
