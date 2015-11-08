#!/usr/bin/env python
from __future__ import division

from fractions import Fraction
from decimal import Decimal

from qform import *


def qforms(x):
    if isinstance(x, str):
        return x
    else:
        return qform(x)


def test_root_opt():
    output("""\
    rootx_opt:{
        f:{(x-3)*x-5};
        $[x=0;.qml.root[f;4 7];
          x=1;(::;0<)@'.qml.rootx[`full;f;4,7.]`x`iter;
          x=2;(null;::;::;::;`sign=)@'
                  .qml.rootx[`full`quiet`iter,200,`tol,0;f;6 7.]
                  `x`last`err`iter`sig;
          x=3;y~@[(),.qml.root[;-1 1]@;z;`$];
          x=4;y~@[(),.qml.root[::;]@;z;`$];
          x=5;y~@[(),.qml.rootx[;::;-1 1]@;z;`$];
          '`]};
    test["rootx_opt[0;::;::]";"5"];
    test["rootx_opt[1;::;::]";"5 1"];
    test["rootx_opt[2;::;::]";"1 7 8 200 1"];
    test["rootx_opt[3;`type;`]";"1"];
    test["rootx_opt[3;`type;()]";"1"];
    test["rootx_opt[3;`type;{0},{0}]";"1"];
    test["rootx_opt[3;`foo;{'`foo}]";"1"];
    test["rootx_opt[3;`type;{()}]";"1"];
    test["rootx_opt[3;`rank;{y}]";"1"];
    test["rootx_opt[4;`type;0]";"1"];
    test["rootx_opt[4;`length;()]";"1"];
    test["rootx_opt[4;`length;1#0]";"1"];
    test["rootx_opt[4;`length;3#0]";"1"];
    test["rootx_opt[4;`type;`,0]";"1"];
    test["rootx_opt[4;`type;0,`]";"1"];
    test["rootx_opt[4;`sign;1 1]";"1"];
    test["rootx_opt[5;`opt;`foo]";"1"];""")

def test_root():
    output("""\
    root_n:{[norm;f;x0]
        norm .qml.rootx[`quiet;f;x0]};
    solve_n:{[norm;f;x0]
        norm .qml.solvex[`quiet;f;x0]};""")

    def emit(func, alt_x0, x, norm = "::", root = True, solve = True):
        alt_funcs = [func, "'[neg;]" + func, "'[{x*x};]" + func]
        if root:
            for func in alt_funcs[:2]:
                for i, j in ((0, 2), (3, 0), (2, 1), (1, 3)):
                    output("    test[\"root_n[%s;%s;(%s;%s)]\";\"%s\"];" %
                        (norm, func, qforms(alt_x0[i]), qforms(alt_x0[j]),
                         qforms(x)))
        if solve:
            for func in alt_funcs:
                for x0 in alt_x0:
                    output("    test[\"solve_n[%s;%s;%s]\";\"%s\"];" %
                        (norm, func, qforms(x0), qforms(x)))

    emit("(42-)",
         [0, 1, 100, 200],
         42)
    emit("{(x-100)*x-3}",
         [-2, 0, 5, 10],
         3)
    emit("{(x-2)*5+(4*x)+x*x}",
         [-100, 0, 5, 100],
         2,
         solve = False) # 1
    emit("{.5+.qml.sin x}",
         [-5, 0, 5, 10],
         "pi%3",
         norm = "{abs(1.5*pi)-x mod 2*pi}")
    emit("{.qml.sin[x]+(x%4)+1%7}",
         [-10, -5, 0, 5],
         Decimal("-0.11448565782234166611"),
         solve = False)
    emit("{(x*.qml.exp x)-2}",
         [-Fraction(1, 2), 0, 1, Fraction(3, 2)],
         Decimal("0.85260550201372549135"))
    emit("(abs@.qml.e-)",
         [-2, -1, 3, 4],
         "e",
         root = False)
    emit("{-1e-1+.qml.atan[x]+pi%2}",
         [0, -1, -25, -100],
         Decimal("-9.9666444232592378598"))
    emit("{$[0<=x-:1;.qml.sqrt[x]-1e-1;1%x]}",
         [-10, 0, 2, 10],
         Fraction(101, 100),
         solve = False)
    emit("{.qml.log 1-abs 42-x}",
         [Fraction(209, 5), Fraction(293, 7), Fraction(337, 8), Fraction(169, 4)],
         42,
         root = False)


def test_line_opt():
    output("""\
    linex_opt:{
        f:{(x-3)*x-5};
        $[x=0;.qml.line[f;0;1];
          x=1;(::;::;0<)@'.qml.linex[enlist`full;f;0.;1]`x`f`iter;
          x=2;{(null y 0;x[y 1]-y 2;0<y 3;`iter=y 4)}[f]
                  .qml.linex[`full`quiet`iter`tol!1 1 3 0;f;0.;1.]
                  `x`last`f`iter`sig;
          x=3;y~@[(),.qml.linex[;f;0;1]@;z;`$];
          x=4;y~@[(),.qml.line[;0;1]@;z;`$];
          x=5;y~@[(),.qml.line[abs;;].;z;`$];
          x=6;y~@[(),.qml.linex[;abs;0;1]@;z;`$];
          '`]};
    test["linex_opt[0;::;::]";"4"];
    test["linex_opt[1;::;::]";"4 -1 1"];
    test["linex_opt[2;::;::]";"1 0 1 1"];
    test["linex_opt[3;`iter;`iter,1]";"1"];
    test["linex_opt[4;`type;`]";"1"];
    test["linex_opt[4;`type;()]";"1"];
    test["linex_opt[4;`type;{0},{0}]";"1"];
    test["linex_opt[4;`foo;{'`foo}]";"1"];
    test["linex_opt[4;`type;{()}]";"1"];
    test["linex_opt[4;`rank;{y}]";"1"];
    test["linex_opt[5;`type;`,0]";"1"];
    test["linex_opt[5;`type;0,`]";"1"];
    test["linex_opt[5;`type;(`;enlist 0)]";"1"];
    test["linex_opt[5;`type;(enlist 0;`)]";"1"];
    test["linex_opt[6;`opt;`foo]";"1"];""")


def test_line():
    output("""\
    line_n:{[norm;f;base;x0]
        norm .qml.linex[`quiet;f;base;x0]};
    minx_n:{[norm;opt;f;x0]
        norm .qml.minx[`quiet,opt;f;x0]};""")

    def emit(func, alt_x0, x, norm = "::"):
        for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3),
                     (3, 2), (3, 1), (3, 0), (2, 1), (2, 0)):
            output("    test[\"line_n[%s;%s;%s;%s]\";\"%s\"];" %
                (norm, func, qforms(alt_x0[i]), qforms(alt_x0[j]), qforms(x)))
        for x0 in alt_x0:
            for opts in ["()", "`slp", "`nm", "`sbplx"]:
                output("    test[\"minx_n[%s;%s;%s;%s]\";\"%s\"];" %
                    (norm, opts, func, qforms(x0), qforms(x)))

    emit("{(x-2)*x+1}",
         [-1, 0, 1, 2],
         Fraction(1, 2))
    emit("{(x*(x+1)*x+3)-(x<-2)*8+x*x*x}",
         [-2, -1, 0, 3],
         "(1%3)*-4+.qml.sqrt 7")
    emit("{1+x*-4+(x*.qml.sqrt[2]-4)+x*x*x}",
         [-3, -1, 2, 5],
         ".qml.sqrt 2")
    emit(".qml.cos",
         [-4, 0, 1, 7],
         "pi",
         norm = "mod[;2*pi]")
    emit("{.qml.sin[x]+x*x%5}",
         [-5, -3, 3, 5],
         Decimal("-1.1105105035811121290"))
    emit("{x*.qml.exp[x]-2}",
         [-1, 0, 1, 2],
         Decimal("0.37482252818362338162"))
    emit("abs@.qml.e-",
         [-2, -1, 3, 4],
         "e")
    emit("{-1%1+1e3*.qml.pow[x-1e4;2]}",
         [0, 1, 19999, 20000],
         10000)
    emit("{-1%x-x<1}",
         [1, "1+100*.qml.eps", 2, 3],
         1)
    emit("{-1e300|.qml.log abs 42-x*x}",
         [-10, 0, 50, 100],
         ".qml.sqrt 42",
         norm = "abs")


def test_solve_opt():
    output("""\
    solvex_opt:{
        f:{20-x-2*y},{-10-(3*x)+4*y};
        $[x=0;.qml.solve[f;0 0];
          x=1;(::;0<)@'
                  .qml.solvex[`full`rk`steps`tol`iter!1,1b,10,.1,100;f;0.,0]
                  `x`iter;
          x=2;(all null@;::;prec>abs@;0<;`feas=)@'
                  .qml.solvex[`full`quiet`slp`tol,0;f;0.,0]
                  `x`last`err`iter`sig;
          x=3;y~@[(),.qml.solve[{1};]@;z;`$];
          x=4;y~@[(),.qml.solve[;0]@;z;`$];
          x=5;y~@[(),.qml.solvex[;{0};0]@;z;`$];
          '`]};
    test["solvex_opt[0;::;::]";"6 -7"];
    test["solvex_opt[1;::;::]";"(6 -7;1)"];
    test["solvex_opt[2;::;::]";"(1;6 -7;1;1;1)"];
    test["solvex_opt[3;`feas;enlist 0#0]";"1"];
    test["solvex_opt[4;`type;`]";"1"];
    test["solvex_opt[4;`type;`,0]";"1"];
    test["solvex_opt[4;`type;{0},`]";"1"];
    test["solvex_opt[4;`type;({0};`,0)]";"1"];
    test["solvex_opt[4;`type;({0};{0},`)]";"1"];
    test["solvex_opt[4;`type;({0};`,{0})]";"1"];
    test["solvex_opt[4;`foo;{'`foo}]";"1"];
    test["solvex_opt[4;`type;{()}]";"1"];
    test["solvex_opt[4;`rank;{y}]";"1"];
    test["solvex_opt[5;`opt;`foo]";"1"];
    test["solvex_opt[5;`opt;`slp`rk]";"1"];
    test["solvex_opt[5;`opt;`slp`steps,10]";"1"];""")


def test_solve():
    def emit(funcs, dfuncs, alt_x0, x, comment = ""):
        if isinstance(funcs, str):
            assert not dfuncs
            alt_funcs = [
                funcs,
                "'[neg;]each reverse" + funcs,
                "raze .[;(::;0);'[neg;]]reverse each 2 cut" + funcs
            ]
        else:
            if len(funcs) == 0:
                alt_funcs = ["()"]
            elif len(funcs) == 1:
                alt_funcs = [
                    funcs[0],
                    "'[neg;]%s" % funcs[0]
                ]
            else:
                alt_funcs = [
                    ",".join(funcs),
                    "'[neg;]each %s" % ",".join(funcs[::-1]),
                    "(%s)" % ";".join(
                        ("%s" if i % 2 else "'[neg;]%s") %
                        funcs[min(i ^ 1, len(funcs)-1)]
                        for i in range(len(funcs)))
                ]
            if dfuncs:
                assert len(funcs) == len(dfuncs) != 0
                if len(funcs) == 1:
                    alt_funcs += [
                        "enlist%s,%s" % (funcs[0], dfuncs[0]),
                        "enlist('[neg;]each %s,%s)" % (funcs[0], dfuncs[0])
                    ]
                else:
                    alt_funcs += [
                        "(%s)" % ";".join(
                            map(",".join, zip(funcs, dfuncs))),
                        "('[neg;]'')(%s)" % ";".join(
                            map(",".join, reversed(zip(funcs, dfuncs)))),
                        "(%s)" % ";".join(
                            ("%s" if i % 2 else "'[neg;]each %s") %
                            "%s,%s" % (funcs[min(i ^ 1, len(funcs)-1)],
                                       dfuncs[min(i ^ 1, len(funcs)-1)])
                            for i in range(len(funcs)))
                    ]

        if comment:
            comment = " / " + comment

        for funcs in alt_funcs:
            for x0 in alt_x0:
                output("    test[\"solve_n[::;%s;%s]\";\"%s\"];%s" %
                    (funcs, qforms(x0), qforms(x), comment))

    emit(["{0}", "{0}"], None, ["enlist 0#0."], [[]])
    emit([], None, [[6, 7]], [6, 7])
    emit(["(42-)"], ["{-1}"], [0], 42)
    emit(["(1-.qml.log@)", "{(x-.qml.e)*x-1}"],
         None,
         [-1],
         None,
         comment = "was invalid read")
    emit(["(1-.qml.log@)", "{(x-.qml.e)*x-1}"],
         ["(-1%)", "{-1-.qml.e-2*x}"],
         [1, 2, 3],
         "e")
    emit(["{25-(x*x)+y*y}", "{64-x*7+3*y}"],
         ["{-2*(x;y)}", "{(-7-3*y;-3*x)}"],
         [[0, 0], [5, 10], [7, -15]],
         [4, 3])
    emit(["{25-(x*x)+y*y}", "{y-1+x%.qml.pow[16;1%x]}"],
         ["{-2*(x;y)}", "{y;((-1-.qml.log[16]%x)%.qml.pow[16;1%x];1)}"],
         [[1, 0], [2, -1], [10, 10]],
         [4, 3])
    emit(["{y-.qml.sin x}", "{y-.qml.sin[x]+1+x%7}"],
         ["{y;(neg .qml.cos x;1)}", "{y;(neg .qml.cos[x]+1%7;1)}"],
         [[0, 0], [30, 0], [-20, -10]],
         "(-7;.qml.sin -7)")
    emit(["{3+y-y*x*x}", "{5-(x*x)+y*y}", "{y-.qml.exp x-2}"],
         ["{(-2*x*y;1-x*x)}", "{-2*(x;y)}", "{y;(neg .qml.exp x-2;1)}"],
         [[1, 0], [4, -2], [10, 2]],
         [2, 1])
    emit(["{1-z-y}", "{2-x-(4*y)-3*z}", "{3-(4*z)-y+4*x}"],
         ["{z;0 1 -1}", "{z;-1 4 -3}", "{z;4 1 -4}"],
         [[0, 0, 0], [5, -3, 1], [-1000, 2000, 4000]],
         [4, 5, 6])
    emit(["{-6-x+y+z}", "{50-(x*x)+(y*y)+z*z}", "{60-x*y*z}", "{32-y*z-x}"],
         ["{z;-1 -1 -1}", "{-2*(x;y;z)}", "{neg(y*z;x*z;x*y)}", "{(y;x-z;neg y)}"],
         [[1, -3, -3], [5, 0, -5], [12, -8, 4]],
         [3, -4, -5]);
    emit("{1-last x},{sum .qml.pow[x;til 6]*}each -1 1 -2 2 -5",
         None,
         [[[0] * 6], [[-100] * 6], [[0, 100, -200, -200, 100, 0]]],
         [[20, 4, -25, -5, 5, 1]])
    emit(["{y;.qml.log[x]-(.qml.asin[.qml.sin[pi*x]]%pi)+x*(1+2*.qml.log 1.5)%3}", "{sum x*x:21 -17-(1 -2.;3 4.)mmu y}"],
         None,
         [[Fraction(7, 4), [1, 2]], [Fraction(5, 4), [2, -4]], [Fraction(1, 10), [2, 0]]],
         [Fraction(3, 2), [5, -8]])


def test_min_opt():
    output("""\
    minx_opt:{[x;opt;y;z]
        f:{(x*x)+(2*y*y)-x*y+1};
        m:$[opt~();.qml.min;.qml.minx opt];
        $[x=0;.qml.min[f;0 0];
          x=1;(::;::;0<)@'.qml.minx[opt,`full;f;0 0]`x`f`iter;
          x=2;(all null@;::;::;::;`iter=)@'
                  .qml.minx[opt,`full`quiet`iter,0;f;1 -1]
                  `x`last`f`iter`sig;
          x=3;(0#0.)~last .qml.min[{y;0};(0;())];
          x=4;y~@[(),m[abs;]@;z;`$];
          x=5;y~@[(),m[;0]@;z;`$];
          x=6;y~@[(),m .;z;`$];
          x=7;y~@[(),.qml.minx[;{0};0]@;z;`$];
          '`]};
    test["minx_opt[0;();::;::]";"4 1%7"];
    test["minx_opt[1;();::;::]";    "(4 1%7;-2%7;1)"];
    test["minx_opt[1;`nm;::;::]";   "(4 1%7;-2%7;1)"];
    test["minx_opt[1;`sbplx;::;::]";"(4 1%7;-2%7;1)"];
    test["minx_opt[2;();::;::]";    "(1;1 -1;3;0;1)"];
    test["minx_opt[2;`nm;::;::]";   "(1;1 -1;3;1;1)"];
    test["minx_opt[2;`sbplx;::;::]";"(1;1 -1;3;1;1)"];
    test["minx_opt[3;();::;::]";"1"];""")

    def emit(err, arg):
        output('    test["minx_opt[%s;%s;%s;%s]";"1"];' %
               (which, opt, err, arg))

    which = 4
    for opt in ["()", "`nm", "`sbplx"]:
        emit("`type", "`")
        emit("`nan", "0n")

    for opt in ["`nm", "`sbplx"]:
        emit("`limit", "10001#0")
        emit("`limit", "100 cut 10001#0")

    which = 5
    for opt in ["()", "`nm", "`sbplx"]:
        emit("`type", "`")
        emit("`type", "()")
        emit("`type", "`,0")
        emit("`type", "{},`")
        emit("`type", "`,{}")
        emit("`foo", "{'`foo}")
        emit("`type", "{()}")
        emit("`rank", "{y}")

    which = 6
    opt = "()"
    emit("`rank", "({0},{y};0)")
    emit("`type", "({0},{`};0)")
    emit("`type", "({0},{0};1#0)")
    emit("`type", "({0},{1#0};0)")
    emit("`type", "({0},{2#0};1#0)")
    emit("`type", "({y;0},{y;3#0};2#0)")
    emit("`type", "({y;0},{y;enlist 2#0};2#0)")
    emit("`type", "({0},{2#0};enlist 2#0)")

    which = 7
    opt = "::"
    emit("`opt", "`foo")
    emit("`opt", "`slp`rk")
    emit("`opt", "`slp`steps,10")
    emit("`opt", "`nm`sbplx")
    emit("`opt", "`nm`slp")
    emit("`opt", "`nm`rk")
    emit("`opt", "`nm`tol,1e-6")
    emit("`opt", "`nm`steps,10")


def test_min():
    def emit(func, dfunc, alt_x0, x, more = False):
        alt_funcs = [(False, func)]
        if dfunc:
            alt_funcs += [(True, "%s,%s" % (func, dfunc))]
        for opts in ["", "`nm", "`sbplx"]:
            if more:
                opts += "`iter,100000"
            for grad, func in alt_funcs:
                for x0 in alt_x0:
                    if grad and ("`nm" in opts or "`sbplx" in opts):
                        continue
                    output("    test[\"minx_n[::;%s;%s;%s]\";\"%s\"];" %
                        (opts or "()", func, qforms(x0), qforms(x)))

    emit("{-1}", None, ["enlist 0#0."], [[]])
    emit("{(a*a:x-1)+10*b*b:y-1+x*x}",
         "{((2*x-1)-40*x*y-1+x*x;20*y-1+x*x)}",
         [[0, 0], [5, 5], [-10, 10]],
         [1, 2])
    emit("{(x*x*x)+(y*y*y)+(2*x*y)+.qml.pow[5;neg x]+.qml.pow[4;neg y]}",
         "{((3*x*x)+(2*y)-.qml.log[5]*.qml.pow[5;neg x];(2*x)+(3*y*y)-.qml.log[4]*.qml.pow[4;neg y])}",
         [[0, 0], [-10, 0], [10, -20]],
         [Decimal("0.34343411653951616845"), Decimal("0.28609347201908362344")])
    emit("{(x*x*x)+(y*y*y)+(2*x*y)+.qml.pow[5;neg x]+.qml.pow[4;neg y]+z*z+y}",
         "{((3*x*x)+(2*y)-.qml.log[5]*.qml.pow[5;neg x];(2*x)+(3*y*y)+z-.qml.log[4]*.qml.pow[4;neg y];y+2*z)}",
         [[0, 0, 0], [1, 1, -1], [10, -20, 30]],
         [Decimal("0.29002456462855436461"), Decimal("0.37840366475888345284"), Decimal("-0.18920183237944172642")])
    emit("{[x;y;z;u]sum x*x:(x-3;(9*y)-x*x;z-x-y;u-z*z)}",
         "{[x;y;z;u]sum 2*(1 0 0 0;(-2*x;9;0;0);-1 1 1 0;(0;0;-2*z;1))*(x-3;(9*y)-x*x;z-x-y;u-z*z)}",
         [[0, 0, 0, 0], [-1, -2, -3, -4], [500, -1000, 1000, -500]],
         [3, 1, 2, 4],
         more = True)
    emit("{{sum{x*x}sum[z*x]-y}[(x;x*x);.qml.sin[10*x]+2*.qml.atan 10*x:til[11]%10]}[]",
         "{{enlist x mmu 2*sum[x*z]-y}[(x;x*x);.qml.sin[10*x]+2*.qml.atan 10*x:til[11]%10]}[]",
         [[[0, 0]], [[10, -5]], [[-100, 50]]],
         [[Decimal("9.4011843228733705351"), Decimal("-6.7445676817128794990")]],
         more = True)



def test_conmin_opt():
    output("""\
    conminx_opt:{[x;opt;y;z]
        f:{(x*x)+(y*y)+(2*x*y)-2*x};
        c:{y-2+x},{y-1};
        m:$[opt~();.qml.conmin;.qml.conminx opt];
        $[x=0;.qml.conmin[f;c;0 0];
          x=1;{(::;::;::;0<)@'x`x`f`cons`iter}
                  .qml.conminx[`full`quiet`slp`lincon`tol,.1;f;c;0 0];
          x=2;{(all null@;::;::;::;::;::;`iter=)@'x`x`last`f`cons`err`iter`sig}
                  .qml.conminx[`full`quiet`rk`steps`iter!1 1 1 3 0;f;c;-.5 -1];
          x=3;{(z[`f]-x . l;0<>z`f;z[`err]+min[y .\:l:z`last];`feas=z`sig)}[f;c1]
                  .qml.conminx[opt,`full`quiet;f;c1:c,{neg y};0 0];
          x=4;y~@[(),m[;{x};0]@;z;`$];
          x=5;y~@[(),m[{x};;0n]@;z;`$];
          x=6;y~@[(),m[{x};;0]@;z;`$];
          x=7;y~@[(),m .;z;`$];
          x=8;y~@[(),m[{0};{-1},{1};]@;z;`$];
          x=9;y~@[(),.qml.conminx[;{0};();0]@;z;`$];
          '`]};
    test["conminx_opt[0;();::;::]";"-3 5%4"];
    test["conminx_opt[1;();::;::]";"(-3 5%4;7%4;0 1%4;1)"];
    test["conminx_opt[2;();::;::]";"(1;-.5 -1;3.25;-2.5 -2;2.5;0;1)"];
    test["conminx_opt[3;();::;::]";     "0 1 0 1"];
    test["conminx_opt[3;`cobyla;::;::]";"0 1 0 1"];""")

    def emit(err, arg):
        output('    test["conminx_opt[%s;%s;%s;%s]";"1"];' %
               (which, opt, err, arg))

    which = 4
    for opt in ["()", "`cobyla"]:
        emit("`type", "`")
        emit("`type", "()")
        emit("`type", "`,0")
        emit("`type", "{},`")
        emit("`type", "`,{}")
        emit("`foo", "{'`foo}")
        emit("`type", "{()}")
        emit("`rank", "{y}")

    which = 5
    for opt in ["`cobyla"]:
        emit("`nan", "{x}")
        emit("`nan", "()")

    which = 6
    for opt in ["()", "`cobyla"]:
        emit("`type", "`")
        emit("`type", "{0},`,{0}")
        emit("`type", "({0};`,0)")
        emit("`type", "({0};{0},`)")
        emit("`type", "({0};`,{0})")
        emit("`foo", "{'`foo}")
        emit("`foo", "{0},{'`foo},{0}")
        emit("`type", "{()}")
        emit("`type", "{0},{()},{0}")
        emit("`rank", "{y}")
        emit("`rank", "{0},{y},{0}")

    for opt in ["`cobyla"]:
        emit("`limit", "10001#{0}")
        emit("`limit", "10001#enlist{0},{0}")

    which = 7
    for opt in ["()"]:
        emit("`rank", "({0};enlist{0},{y};0)")
        emit("`type", "({0};enlist{0},{`};0)")
        emit("`type", "({0};enlist{0},{0};1#0)")
        emit("`type", "({0};enlist{0},{1#0};0)")
        emit("`type", "({0};enlist{0},{2#0};1#0)")
        emit("`type", "({y;0};enlist{y;0},{y;3#0};2#0)")
        emit("`type", "({y;0};enlist{y;0},{y;enlist 2#0};2#0)")
        emit("`type", "({0};enlist{0},{2#0};enlist 2#0)")

    which = 8
    for opt in ["()", "`cobyla", "`cobyla`full"]:
        emit("`feas", "enlist 0#0")

    for opt in ["`cobyla"]:
        emit("`limit", "10001#0")
        emit("`limit", "100 cut 10001#0")

    which = 9
    opt = "::"
    emit("`opt", "`foo")
    emit("`opt", "`slp`rk")
    emit("`opt", "`slp`steps,10")
    emit("`opt", "`cobyla`slp")
    emit("`opt", "`cobyla`rk")
    emit("`opt", "`cobyla`steps,10")


def test_conmin():
    def emit(func, dfunc, cons, dcons, alt_x0, x, linear = False, more = False):
        opts = "`quiet"
        alt_opts = [opts]
        if linear:
            alt_opts += [opts + "`lincon",
                         opts + "`slp"]
        alt_opts += [opts + "`cobyla"]

        alt_funcs = [(False, func)]
        if dfunc:
            alt_funcs += [(True, "%s,%s" % (func, dfunc))]

        if isinstance(cons, str):
            assert not dcons
            alt_cons = [(False, cons),
                        (False, "reverse" + cons)]
        elif cons:
            alt_cons = [(False, ",".join(cons))]
            if len(cons) > 1:
                alt_cons += [(False, ",".join(cons[::-1]))]
            if dcons:
                assert len(cons) == len(dcons)
                if len(cons) == 1:
                    alt_cons += [
                        (True, "enlist%s,%s" % (cons[0], dcons[0]))]
                else:
                    alt_cons += [
                        (True, "(%s)" % ";".join(
                            map(",".join, zip(cons, dcons)))),
                        (True, "(%s)" % ";".join(
                            "%s,%s" % (cons[i], dcons[i]) if i % 2 else cons[i]
                            for i in reversed(range(len(cons)))))
                    ]
        else:
            alt_cons = [(False, "()")]

        for cgrad, cons in alt_cons:
            for fgrad, func in alt_funcs:
                for x0 in alt_x0:
                    for opts in alt_opts:
                        if (fgrad or cgrad) and ("`nm" in opts or
                                                 "`cobyla" in opts):
                            continue
                        if more:
                            opts += "`iter,100000"
                            if "`cobyla" in opts:
                                opts += "00,`tol,1e-6"
                        output(
                            "    test[\".qml.conminx[%s;%s;%s;%s]\";\"%s\"];" %
                            (opts, func, cons, qforms(x0), qforms(x)))

    emit("{-1}", None, ["{1}"], None, ["enlist 0#0."], [[]])
    emit("{(a*a:x-1)+10*b*b:y-1+x*x}",
         "{((2*x-1)-40*x*y-1+x*x;20*y-1+x*x)}",
         [],
         [],
         [[0, 0], [5, 5], [-10, 10]],
         [1, 2])
    emit("{(x*x*x*x)+y*(3*x*1+x*6)+y*y*y-8}",
         "{((4*x*x*x)+y*3*1+x*12;(3*x*1+x*6)+4*y*y*y-6)}",
         ["{x-y+1}"],
         ["{y;1 -1}"],
         [[0, 0], [15, 5], [-10, -5]],
         [Decimal("6.1546405692889590154"), Decimal("-4.1526183884860082232")],
         linear = True)
    emit("{neg(5*x)+(3*y)+7*z}",
         "{z;-5 -3 -7}",
         ["{z;10-(2*x)+4*y}", "{15-(3*z)-y}", "{9-x+y+z}", "{z;x}", "{z;y}", "{z}"],
         ["{z;-2 -4 0}", "{z;0 1 -3}", "{z;-1 -1 -1}", "{z;1 0 0}", "{z;0 1 0}", "{z;0 0 1}"],
         [[0, 0, 0], [-8, -10, -3], [200, 400, 100]],
         [4, 0, 5],
         linear = True)
    emit("{[x;y;z;u]sum x*x:(x-3;(9*y)-x*x;z-x-y;u-z*z)}",
         "{[x;y;z;u]sum 2*(1 0 0 0;(-2*x;9;0;0);-1 1 1 0;(0;0;-2*z;1))*(x-3;(9*y)-x*x;z-x-y;u-z*z)}",
         ["{[x;y;z;u]sum -1 -1 1 1%v+1e-99%v:(x;y;z;u)}", "{[x;y;z;u]min(x;y;z;u)}"],
         ["{[x;y;z;u] -1 -1 1 1*(1e-99-v)%w*w:1e-99+v*:v:(x;y;z;u)}", "{[x;y;z;u]v=min v:(x;y;z;u)}"],
         [[1, 1, 1, 1], [2, -3, 4, -5], [1000, 1000, 500, 500]],
         [Decimal("3.1355504511518109809"), Decimal("1.1046418527109901434"), Decimal("1.4278718594741174466"), Decimal("1.9089393921649469239")],
         more = True)
    emit("{(x mmu (.01 .02 -.003;.02 .2 -.004;-.003 -.004 .001) mmu x)-.7*.1 .25 .05 mmu x:0.+x}",
         "{enlist(2*(.01 .02 -.003;.02 .2 -.004;-.003 -.004 .001)mmu 0.+x)-.7*.1 .25 .05}",
         ["{x 0}", "{x 1}", "{x 2}", "{sum(1%3)-x}", "{sum x-1%3}"],
         ["{enlist 1 0 0}", "{enlist 0 1 0}", "{enlist 0 0 1}", "{enlist -1 -1 -1}", "{enlist 1 1 1}"],
         [[[0, 0, 0]], [[0, Fraction(1, 2), Fraction(1, 2)]], [[Fraction(371, 3), Fraction(-217, 48), Fraction(4445, 12)]]],
         [[Fraction(3, 4), Fraction(1, 4), 0]],
         linear = True)
    emit("{neg x+2*y}",
         "{y;-1 -2}",
         ["{25-(x*x)+y*y}", "{1+(3*x)-x*y*(3+1%.qml.sqrt 11)%.qml.sqrt 14}", "{(67+2*.qml.sqrt 154)-(x*x)+(2*x*y)+4*y*y}", "{x-.5*y}"],
         ["{-2*(x;y)}", "{3 0-(y;x)*(3+1%.qml.sqrt 11)%.qml.sqrt 14}", "{(67+2*.qml.sqrt 154)-2*(x+y;x+4*y)}", "{y;1 -.5}"],
         [[0, 0], [-2, -8], [-7, 15]],
         ".qml.sqrt 11 14")
    emit("{y;0}",
         "{y;(0;0 0)}",
         "{x,(''[neg;x])}{y;.qml.log[x]-(.qml.asin[.qml.sin[pi*x]]%pi)+x*(1+2*.qml.log 1.5)%3},{sum{x*x}21 -17-(1 -2.;3 4.)mmu y}",
         None,
         [[Fraction(7, 4), [1, 3]], [Fraction(5, 4), [2, -4]], [Fraction(1, 10), [2, 0]]],
         [Fraction(3, 2), [5, -8]],
         more = True)


def tests():
    reps(25)
    prec("1e-5")
    test_root_opt()
    test_root()
    prec("1e-6")
    test_line_opt()
    test_line()
    prec("1e-4")
    test_solve_opt()
    test_solve()
    test_min_opt()
    test_min()
    prec("1e-3")
    test_conmin_opt()
    test_conmin()


if __name__ == "__main__":
    tests()
