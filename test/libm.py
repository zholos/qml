#!/usr/bin/env python2
from __future__ import division

import sympy as sp
import sympy.stats as st
import mpmath as mp
from sympy import S, pi, oo, FiniteSet, ProductSet, Interval
mpf = mp.mpf

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


def N(x):
    if x.is_Rational and abs(x.p * x.q) < 2**63-1:
        return x
    else:
        return mpf(x.evalf(mp.mp.dps))

def test_hyper():
    for x in sorted(exparg):
        test("erf", x, N(sp.erf(x)))
    for x in sorted(exparg):
        test("erfc", x, N(sp.erfc(x)))

    gamarg = FiniteSet(*(x+S(1)/12 for x in exparg))
    betarg = ProductSet(gamarg, gamarg)
    for x in sorted(gamarg):
        test("lgamma", x, N(sp.log(abs(sp.gamma(x)))))
    for x in sorted(gamarg):
        test("gamma", x, N(sp.gamma(x)))
    for x, y in sorted(betarg, key=lambda (x, y): (y, x)):
        test("beta", x, y, N(sp.beta(x, y)))

    pgamarg = FiniteSet(S(1)/12, S(1)/3, S(3)/2, 5)
    pgamargp = ProductSet(gamarg & Interval(0, oo, True), pgamarg)
    for a, x in sorted(pgamargp):
        test("pgamma", a, x, N(sp.lowergamma(a, x)))
    for a, x in sorted(pgamargp):
        test("pgammac", a, x, N(sp.uppergamma(a, x)))
    for a, x in sorted(pgamargp):
        test("pgammar", a, x, N(sp.lowergamma(a, x)/sp.gamma(a)))
    for a, x in sorted(pgamargp):
        test("pgammarc", a, x, N(sp.uppergamma(a, x)/sp.gamma(a)))
    for a, x in sorted(pgamargp):
        test("ipgammarc", a, N(sp.uppergamma(a, x)/sp.gamma(a)), x)

    pbetargp = [(a, b, x) for a, b, x in ProductSet(betarg, pgamarg)
                if a > 0 and b > 0 and x < 1]
    pbetargp.sort(key=lambda (a, b, x): (b, a, x))
    for a, b, x in pbetargp:
        test("pbeta", a, b, x, mp.betainc(mpf(a), mpf(b), x2=mpf(x)))
    for a, b, x in pbetargp:
        test("pbetar", a, b, x, mp.betainc(mpf(a), mpf(b), x2=mpf(x),
                                           regularized=True))
    for a, b, x in pbetargp:
        test("ipbetar", a, b, mp.betainc(mpf(a), mpf(b), x2=mpf(x),
                                         regularized=True), x)

    for x in sorted(posarg):
        test("j0", x, N(sp.besselj(0, x)))
    for x in sorted(posarg):
        test("j1", x, N(sp.besselj(1, x)))
    for x in sorted(posarg-FiniteSet(0)):
        test("y0", x, N(sp.bessely(0, x)))
    for x in sorted(posarg-FiniteSet(0)):
        test("y1", x, N(sp.bessely(1, x)))


def test_prob():
    def emit(name, iname, cdf, args, no_small=False):
        V = []
        for arg in sorted(args):
            y = cdf(*arg)
            if isinstance(y, mpf):
                e = sp.nsimplify(y, rational=True)
                if e.is_Rational and e.q <= 1000 and \
                        mp.almosteq(mp.mpf(e), y, 1e-25):
                    y = e
            else:
                y = N(y)
            V.append(arg + (y,))
        for v in V:
            if name:
                test(name, *v)
        for v in V:
            if iname and (not no_small or 1/1000 <= v[-1] <= 999/1000):
                test(iname, *(v[:-2] + v[:-3:-1]))

    x = sp.Symbol("x")
    emit("ncdf", "nicdf",
         sp.Lambda(x, st.cdf(st.Normal("X", 0, 1))(x)), zip(exparg))
    # using cdf() for anything more complex is too slow

    df = FiniteSet(1, S(3)/2, 2, S(5)/2, 5, 25)
    emit("c2cdf", "c2icdf",
         lambda k, x: sp.lowergamma(k/2, x/2)/sp.gamma(k/2),
         ProductSet(df, posarg), no_small=True)

    dfint = df & sp.fancysets.Naturals()
    def cdf(k, x):
        k, x = map(mpf, (k, x))
        return .5 + .5*mp.sign(x)*mp.betainc(k/2, .5, x1=1/(1+x**2/k),
                                             regularized=True)
    emit("stcdf", "sticdf", cdf, ProductSet(dfint, exparg))

    def cdf(d1, d2, x):
        d1, d2, x = map(mpf, (d1, d2, x))
        return mp.betainc(d1/2, d2/2, x2=x/(x+d2/d1), regularized=True)

    emit("fcdf", "ficdf", cdf, ProductSet(dfint, dfint, posarg))

    kth = ProductSet(sp.ImageSet(lambda x: x/5, df),
                     posarg - FiniteSet(0))
    emit("gcdf", "gicdf",
         lambda k, th, x: sp.lowergamma(k, x/th)/sp.gamma(k),
         ProductSet(kth, posarg), no_small=True)

    karg = FiniteSet(0, 1, 2, 5, 10, 15, 40)
    knparg = [(k, n, p) for k, n, p
              in ProductSet(karg, karg, posarg & Interval(0, 1, True, True))
              if k <= n and n > 0]
    def cdf(k, n, p):
        return st.P(st.Binomial("X", n, p) <= k)
    emit("bncdf", "bnicdf", cdf, knparg, no_small=True)

    def cdf(k, lamda):
        return sp.uppergamma(k+1, lamda)/sp.gamma(k+1)
    emit("pscdf", "psicdf", cdf,
         ProductSet(karg, posarg + karg - FiniteSet(0)), no_small=True)

    x, i = sp.symbols("x i")
    def smcdf(n, e):
        return 1-sp.Sum(sp.binomial(n, i)*e*(e+i/n)**(i-1)*(1-e-i/n)**(n-i),
                        (i, 0, sp.floor(n*(1-e)))).doit()
    kcdf = sp.Lambda(x,
        sp.sqrt(2*pi)/x*sp.Sum(sp.exp(-pi**2/8*(2*i-1)**2/x**2), (i, 1, oo)))
    smarg = ProductSet(karg - FiniteSet(0), posarg & Interval(0, 1, True, True))
    karg = FiniteSet(S(1)/100, S(1)/10) + (posarg & Interval(S(1)/4, oo, True))

    for n, e in sorted(smarg):
        test("smcdf", n, e, N(smcdf(n, e)))
    prec("1e-10")
    for x in sorted(karg):
        test("kcdf", x, N(kcdf(x)))
    prec("1e-9")
    for n, e in sorted(smarg):
        p = smcdf(n, e)
        if p < S(9)/10:
            test("smicdf", n, N(p), e)
    prec("1e-6")
    for x in sorted(karg):
        p = kcdf(x)
        if N(p) > S(10)**-8:
            test("kicdf", N(p), x)


def tests():
    test_pow()
    test_trig()
    test_trigh()
    test_extra()
    test_hyper()
    test_prob()


if __name__ == "__main__":
    tests()
