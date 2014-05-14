#!/usr/bin/env python
from __future__ import division

from fractions import Fraction
from decimal import Decimal

from qform import qform


def qforms(x):
    if isinstance(x, str):
        return x
    else:
        return qform(x)


def test_root():
    output("""\
    rootx_opt:{
        f:{(x-3)*x-5};
        $[x=0;.qml.root[f;4 7];
          x=1;(::;0<)@'.qml.rootx[`full;f;4,7.]`x`iter;
          x=2;(null;::;::;::;`sign=)@'
                  .qml.rootx[`full`quiet`iter,200,`tol,0;f;6 7.]
                  `x`last`err`iter`sig;
          x=3;`type~@[.qml.root[;0 0];0;`$];
          x=4;`type`length`type`type`sign~@[.qml.root[::;];;`$] each
                  (0;enlist 0;`,0;0,`;1 1);
          '`]};
    test["rootx_opt 0";"5"];
    test["rootx_opt 1";"5 1"];
    test["rootx_opt 2";"1 7 8 200 1"];
    test["rootx_opt 3";"1"];
    test["rootx_opt 4";"1"];
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


def test_line():
    output("""\
    linex_opt:{
        f:{(x-3)*x-5};
        $[x=0;.qml.line[f;0;1];
          x=1;(::;::;0<)@'.qml.linex[enlist`full;f;0.;1]`x`f`iter;
          x=2;{(null y 0;x[y 1]-y 2;0<y 3;`iter=y 4)}[f]
                  .qml.linex[`full`quiet`iter`tol!1 1 3 0;f;0.;1.]
                  `x`last`f`iter`sig;
          x=3;`type`type`type~.[.qml.line;;`$] each @[(abs;0;0);;:;`]'[til 3];
          '`]};
    test["linex_opt 0";"4"];
    test["linex_opt 1";"4 -1 1"];
    test["linex_opt 2";"1 0 1 1"];
    test["linex_opt 3";"1"];
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


def test_solve():
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
          '`]};
    test["solvex_opt 0";"6 -7"];
    test["solvex_opt 1";"(6 -7;1)"];
    test["solvex_opt 2";"(1;6 -7;1;1;1)"];""")

    def emit(funcs, alt_x0, x, comment = ""):
        if isinstance(funcs, str):
            alt_funcs = [funcs,
                         "'[neg;]each reverse" + funcs,
                         "raze .[;(::;0);'[neg;]]reverse each 2 cut" + funcs]
        elif len(funcs) > 1:
            alt_funcs = [",".join(funcs),
                         "'[neg;]each " + ",".join(funcs[::-1]),
                         "(%s)" % ";".join([("" if i % 2 else "'[neg;]") +
                                            funcs[min(i ^ 1, len(funcs) - 1)]
                                            for i in range(len(funcs))])]

        if comment:
            comment = " / " + comment

        for funcs in alt_funcs:
            for x0 in alt_x0:
                output("    test[\"solve_n[::;%s;%s]\";\"%s\"];%s" %
                    (funcs, qforms(x0), qforms(x), comment))

    emit(["(1-.qml.log@)", "{(x-.qml.e)*x-1}"],
         [-1],
         None,
         comment = "was invalid read")
    emit(["(1-.qml.log@)", "{(x-.qml.e)*x-1}"],
         [1, 2, 3],
         "e")
    emit(["{25-(x*x)+y*y}", "{64-x*7+3*y}"],
         [[0, 0], [5, 10], [7, -15]],
         [4, 3])
    emit(["{25-(x*x)+y*y}", "{y-1+x%.qml.pow[16;1%x]}"],
         [[1, 0], [2, -1], [10, 10]],
         [4, 3])
    emit(["{y-.qml.sin x}", "{y-.qml.sin[x]+1+x%7}"],
         [[0, 0], [30, 0], [-20, -10]],
         "(-7;.qml.sin -7)")
    emit(["{3+y-y*x*x}", "{5-(x*x)+y*y}", "{y-.qml.exp x-2}"],
         [[1, 0], [4, -2], [10, 2]],
         [2, 1])
    emit(["{1-z-y}", "{2-x-(4*y)-3*z}", "{3-(4*z)-y+4*x}"],
         [[0, 0, 0], [5, -3, 1], [-1000, 2000, 4000]],
         [4, 5, 6])
    emit(["{-6-x+y+z}", "{50-(x*x)+(y*y)+z*z}", "{60-x*y*z}", "{32-y*z-x}"],
         [[1, -3, -3], [5, 0, -5], [12, -8, 4]],
         [3, -4, -5]);
    emit("{1-last x},{sum .qml.pow[x;til 6]*}each -1 1 -2 2 -5",
         [[[0] * 6], [[-100] * 6], [[0, 100, -200, -200, 100, 0]]],
         [[20, 4, -25, -5, 5, 1]])
    emit(["{y;.qml.log[x]-(.qml.asin[.qml.sin[pi*x]]%pi)+x*(1+2*.qml.log 1.5)%3}", "{sum x*x:21 -17-(1 -2.;3 4.)mmu y}"],
         [[Fraction(7, 4), [1, 2]], [Fraction(5, 4), [2, -4]], [Fraction(1, 10), [2, 0]]],
         [Fraction(3, 2), [5, -8]])


def test_min():
    output("""\
    minx_opt:{
        f:{(x*x)+(2*y*y)-x*y+1};
        $[x=0;.qml.min[f;0 0];
          x=1;(::;::;0<)@'.qml.minx[y,`full;f;0 0]`x`f`iter;
          x=2;(all null@;::;::;::;`iter=)@'
                  .qml.minx[y,`full`quiet`iter,0;f;1 -1]
                  `x`last`f`iter`sig;
          x=3;`type`nan~@[.qml.minx[y;abs;];;`$] each `,0n;
          x=4;(0#0.)~last .qml.min[{y;0};(0;())];
          '`]};
    test["minx_opt[0;()]";"4 1%7"];
    test["minx_opt[1;()]";    "(4 1%7;-2%7;1)"];
    test["minx_opt[1;`nm]";   "(4 1%7;-2%7;1)"];
    test["minx_opt[1;`sbplx]";"(4 1%7;-2%7;1)"];
    test["minx_opt[2;()]";    "(1;1 -1;3;0;1)"];
    test["minx_opt[2;`nm]";   "(1;1 -1;3;1;1)"];
    test["minx_opt[2;`sbplx]";"(1;1 -1;3;1;1)"];
    test["minx_opt[3;()]";"1"];
    test["minx_opt[3;`nm`full]";"1"];
    test["minx_opt[4;()]";"1"];""")

    def emit(func, alt_x0, x, more = False):
        for opts in ["", "`nm", "`sbplx"]:
            if more:
                opts += "`iter,100000"
            for x0 in alt_x0:
                output("    test[\"minx_n[::;%s;%s;%s]\";\"%s\"];" %
                    (opts or "()", func, qforms(x0), qforms(x)))

    emit("{(a*a:x-1)+10*b*b:y-1+x*x}",
         [[0, 0], [5, 5], [-10, 10]],
         [1, 2])
    emit("{(x*x*x)+(y*y*y)+(2*x*y)+.qml.pow[5;neg x]+.qml.pow[4;neg y]}",
         [[0, 0], [-10, 0], [10, -20]],
         [Decimal("0.34343411653951616845"), Decimal("0.28609347201908362344")])
    emit("{(x*x*x)+(y*y*y)+(2*x*y)+.qml.pow[5;neg x]+.qml.pow[4;neg y]+z*z+y}",
         [[0, 0, 0], [1, 1, -1], [10, -20, 30]],
         [Decimal("0.29002456462855436461"), Decimal("0.37840366475888345284"), Decimal("-0.18920183237944172642")])
    emit("{[x;y;z;u]sum x*x:(x-3;(9*y)-x*x;z-x-y;u-z*z)}",
         [[0, 0, 0, 0], [-1, -2, -3, -4], [500, -1000, 1000, -500]],
         [3, 1, 2, 4],
         more = True)
    emit("{{sum{x*x}sum[z*x]-y}[(x;x*x);.qml.sin[10*x]+2*.qml.atan 10*x:til[11]%10]}[]",
         [[[0, 0]], [[10, -5]], [[-100, 50]]],
         [[Decimal("9.4011843228733705351"), Decimal("-6.7445676817128794990")]],
         more = True)


def test_conmin():
    output("""\
    conminx_opt:{
        f:{(x*x)+(y*y)+(2*x*y)-2*x};
        c:{y-2+x},{y-1};
        $[x=0;.qml.conmin[f;c;0 0];
          x=1;(::;::;::;0<)@'
                  .qml.conminx[`full`quiet`slp`lincon`tol,.1;f;c;0 0]
                  `x`f`cons`iter;
          x=2;(all null@;::;::;::;::;::;`iter=)@'
                  .qml.conminx[`full`quiet`rk`steps`iter!1 1 1 3 0;f;c;-.5 -1]
                  `x`last`f`cons`err`iter`sig;
          x=3;{(z[`f]-x . l;0<>z`f;z[`err]+min[y .\:l:z`last];`feas=z`sig)}[f;c1]
                  .qml.conminx[y,`full`quiet;f;c1:c,{neg y};0 0];
          x=4;`nan`nan~@[.qml.conminx[y,`full;{x};;0n];;`$] each ({x};());
          x=5;`type`foo`type`rank~@[.qml.conmin[;{x};0];;`$] each
                  0,{'`foo},{()},{y};
          '`]};
    test["conminx_opt[0;()]";"-3 5%4"];
    test["conminx_opt[1;()]";"(-3 5%4;7%4;0 1%4;1)"];
    test["conminx_opt[2;()]";"(1;-.5 -1;3.25;-2.5 -2;2.5;0;1)"];
    test["conminx_opt[3;()]";     "0 1 0 1"];
    test["conminx_opt[3;`cobyla]";"0 1 0 1"];
    test["conminx_opt[4;`cobyla]";"1"];
    test["conminx_opt[5;()]";"1"];""")

    def emit(func, cons, alt_x0, x, linear = False, more = False):
        opts = "`quiet"
        alt_opts = [opts]
        if linear:
            alt_opts += [opts + "`lincon",
                         opts + "`slp"]
        alt_opts += [opts + "`cobyla"]

        if isinstance(cons, str):
            alt_cons = [cons,
                        "reverse" + cons]
        elif cons:
            alt_cons = [",".join(cons)]
            if len(cons) > 1:
                alt_cons.append(",".join(cons[::-1]))
        else:
            alt_cons = ["()"]

        for cons in alt_cons:
            for x0 in alt_x0:
                for opts in alt_opts:
                    if more:
                        opts += "`iter,100000"
                        if "`cobyla" in opts:
                            opts += "00,`tol,1e-6"
                    output("    test[\".qml.conminx[%s;%s;%s;%s]\";\"%s\"];" %
                        (opts, func, cons, qforms(x0), qforms(x)))

    emit("{(a*a:x-1)+10*b*b:y-1+x*x}",
         [],
         [[0, 0], [5, 5], [-10, 10]],
         [1, 2])
    emit("{(x*x*x*x)+y*(3*x*1+x*6)+y*y*y-8}",
         ["{x-y+1}"],
         [[0, 0], [15, 5], [-10, -5]],
         [Decimal("6.1546405692889590154"), Decimal("-4.1526183884860082232")],
         linear = True)
    emit("{neg(5*x)+(3*y)+7*z}",
         ["{z;10-(2*x)+4*y}", "{15-(3*z)-y}", "{9-x+y+z}", "{z;x}", "{z;y}", "{z}"],
         [[0, 0, 0], [-8, -10, -3], [200, 400, 100]],
         [4, 0, 5],
         linear = True)
    emit("{[x;y;z;u]sum{x*x}(x-3;(9*y)-x*x;z-x-y;u-z*z)}",
         ["{[x;y;z;u]sum -1 -1 1 1%v+1e-99%v:(x;y;z;u)}", "{[x;y;z;u]min(x;y;z;u)}"],
         [[1, 1, 1, 1], [2, -3, 4, -5], [1000, 1000, 500, 500]],
         [Decimal("3.1355504511518109809"), Decimal("1.1046418527109901434"), Decimal("1.4278718594741174466"), Decimal("1.9089393921649469239")],
         more = True)
    emit("{(x mmu (.01 .02 -.003;.02 .2 -.004;-.003 -.004 .001) mmu x)-.7*.1 .25 .05 mmu x:0.+x}",
         ["{x 0}", "{x 1}", "{x 2}", "{1-sum x}", "{sum[x]-1}"],
         [[[0, 0, 0]], [[0, Fraction(1, 2), Fraction(1, 2)]], [[Fraction(371, 3), Fraction(-217, 48), Fraction(4445, 12)]]],
         [[Fraction(3, 4), Fraction(1, 4), 0]],
         linear = True)
    emit("{neg x+2*y}",
         ["{25-(x*x)+y*y}", "{1+(3*x)-x*y*(3+1%.qml.sqrt 11)%.qml.sqrt 14}", "{(67+2*.qml.sqrt 154)-(x*x)+(2*x*y)+4*y*y},{x-.5*y}"],
         [[0, 0], [-2, -8], [-7, 15]],
         ".qml.sqrt 11 14")
    emit("{y;0}",
         "{x,(''[neg;x])}{y;.qml.log[x]-(.qml.asin[.qml.sin[pi*x]]%pi)+x*(1+2*.qml.log 1.5)%3},{sum{x*x}21 -17-(1 -2.;3 4.)mmu y}",
         [[Fraction(7, 4), [1, 3]], [Fraction(5, 4), [2, -4]], [Fraction(1, 10), [2, 0]]],
         [Fraction(3, 2), [5, -8]],
         more = True)


def test_prec(prec):
    output("    prec:%s;" % prec)

def test_reps(reps):
    output("    reps:%d;" % reps)


def tests():
    test_reps(25)
    test_prec("1e-6")
    test_root()
    test_prec("1e-7")
    test_line()
    test_prec("1e-4")
    test_solve()
    test_min()
    test_conmin()


if __name__ == "__main__":
    def output(s):
        print s
    tests()
