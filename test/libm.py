#!/usr/bin/env python
from __future__ import division

import sympy as sp
from sympy import S, pi, FiniteSet, Interval

from qform import *


posarg = FiniteSet(0, S(1)/4, S(1)/3, S(1)/2, S(3)/4, 1, 2, 3)
exparg = posarg + FiniteSet(*(-x for x in posarg))

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

    for x in sorted(exparg):
        test("exp", x, sp.exp(x), no_pow=True)
    for x in sorted(exparg):
        test("pow", pi, x, pi**x, no_pow=True)

def test_trig():
    arg = FiniteSet(0, *(pi/x for x in (12, 6, 4, 3, 2, 1)))
    arg += FiniteSet(*sum(((pi-x, pi+x, 2*pi-x) for x in arg), ()))
    arg += FiniteSet(*(-x for x in arg))

    def emit(name, func):
        for x in sorted(arg):
            if func(x) != sp.zoo:
                test(name, x, sp.simplify(sp.sqrtdenest(func(x))))
        test(name, S.NaN, S.NaN)

    def aemit(aname, func, domain):
        for x, y in sorted((sp.simplify(sp.sqrtdenest(func(x))), x)
                           for x in arg & domain):
            test(aname, x, y)
        test(aname, S.NaN, S.NaN)

    emit("sin", sp.sin)
    emit("cos", sp.cos)
    emit("tan", sp.tan)
    aemit("asin", sp.sin, Interval(-pi/2, pi/2))
    aemit("acos", sp.cos, Interval(0, pi))
    aemit("atan", sp.tan, Interval(-pi/2, pi/2, True, True))

    for x in sorted(arg & Interval(-pi, pi, True, True)):
        a = sp.simplify(sp.sqrtdenest(sp.sin(x)))
        b = sp.simplify(sp.sqrtdenest(sp.cos(x)))
        f = sp.sqrt(3) if \
            sp.sqrt(3) == abs(sp.simplify(a/b).as_numer_denom()[0]) else 1
        if b:
            test("atan2", sp.simplify(a/abs(b))*f, sp.sign(b)*f, x)
        if a and a != b:
            test("atan2", sp.sign(a)*f, sp.simplify(b/abs(a))*f, x)
    test("atan2", S.NaN, 1, S.NaN)
    test("atan2", 1, S.NaN, S.NaN)

def test_trigh():
    def emit(name, func):
        for x in sorted(exparg):
            test(name, x, sp.ratsimp(sp.radsimp(func(x).rewrite(sp.exp))),
                 no_trigh=True)

    emit("sinh", sp.sinh)
    emit("cosh", sp.cosh)
    emit("tanh", sp.tanh)

    def aemit(aname, func, even):
        for x in sorted(exparg):
            if not even or x >= 0:
                test(aname, func(x), x)

    aemit("asinh", sp.sinh, False)
    aemit("acosh", sp.cosh, True)
    aemit("atanh", sp.tanh, False)

def test_extra():
    m1arg = FiniteSet(*(x / max(exparg) for x in exparg))
    for x in sorted(m1arg):
        test("expm1", x, sp.exp(x)-1)

    def emit_log(name, base, floor=False):
        for x in sorted(exparg):
            test(name, base**x, sp.floor(x) if floor else x)

    emit_log("log", sp.E)
    emit_log("log10", 10)
    emit_log("logb", 2, True)

    for x in sorted(m1arg):
        test("log1p", qstr(".qml.expm1 %s" % qform(x)), x)

    for x in sorted(posarg):
        y = sum(tuple(sorted(posarg))[-2:]) - x
        test("hypot", x, y, sp.sqrt(x**2 + y**2))

    roundarg = FiniteSet(*sum(((x-42, x+42) for x in posarg if x < 1), ()))
    def emit_round(name, func):
        for x in sorted(roundarg):
            test(name, x, func(x))

    emit_round("floor", sp.floor)
    emit_round("ceil", sp.ceiling)

    for x in sorted(exparg):
        test("fabs", x, abs(x))

    for x in sorted(exparg):
        y = S(1)/3
        test("fmod", x, y, x % y - y if x < 0 and x % y else x % y)


def tests():
    test_pow()
    test_trig()
    test_trigh()
    test_extra()


if __name__ == "__main__":
    tests()
