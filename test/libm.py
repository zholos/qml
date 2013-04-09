#!/usr/bin/env python
from __future__ import division

import sympy as sp
from sympy import S, FiniteSet

from qform import *


posarg = FiniteSet(0, S(1)/4, S(1)/3, S(1)/2, S(3)/4, 1, 2, 3)

def test_pow():
    def emit_root(n, rational):
        arg = posarg + FiniteSet(8, 9)
        if n % 2:
            arg += FiniteSet(*(-x for x in arg))
        func = {2: "sqrt", 3: "cbrt"}[n]
        if not rational:
            # Later testcases can use .qml.sqrt to represent irrational numbers,
            # but while testing .qml.sqrt itself just check the inverse.
            output("    %s_i:{%s:.qml.%s x};" % (func, "*".join(["x"]*n), func))
        for x in sorted(arg):
            y = sp.real_root(x, n)
            if bool(sp.ask(sp.Q.rational(y))) == bool(rational):
                if rational:
                    test(func, x, y)
                else:
                    test("%s_i" % func, x, x)

    for rational in (True, False):
        for n in (2, 3):
            emit_root(n, rational)

def tests():
    test_pow()


if __name__ == "__main__":
    tests()
