from __future__ import division

import sys
from fractions import Fraction
import decimal
from decimal import Decimal
from types import NoneType
import sympy as sp
import sympy.mpmath as mp
from sympy import S

decimal.getcontext().prec = 50
mp.mp.dps = 50

__all__ = ["qform", "qstr", "output", "prec", "reps", "test"]


def qform(self, no_pow=False, no_trigh=False, complex_pair=False):
    # form types:
    #        "i", "j", "f" - partial scalar atoms
    #   "S", "s" - scalar expression
    #   "V", "v" - list expression
    #    ^.   ^---- doesn't need brackets
    #     `-------- needs brackets as left argument

    # gathers:
    #   *             -> "s"/"S"
    #   ["i"/"j"/"f"] -> "v"
    #   [*]           -> "v"/"V"

    def encode_vector(forms, types):
        if any([t not in ("i", "j", "f") for t in types]):
            raise Exception()
        q = " ".join(forms)
        if "f" in types:
            if "." not in q and "0n" not in forms:
                q += "."
        elif "j" in types:
            q += "j"
        return q

    def encode_item(f, t, left=False):
        if t in ("i", "j", "f"):
            return encode_vector([f], [t])
        elif left and t in ("S", "V"):
            return "(%s)" % f
        else:
            return f

    def rational(n):
        if isinstance(n, (int, long)):
            return n, 1
        if isinstance(n, Fraction):
            return n.numerator, n.denominator
        if isinstance(n, sp.Rational):
            return n.p, n.q

    def number_form(n):
        if rational(n):
            p, q = rational(n)
            if q == 1:
                if abs(n) < 2**31-1:
                    return str(n), "i"
                if abs(n) < 2**63-1:
                    return str(n), "j"
            else:
                forms, types = zip(*map(number_form, (p, q)))
                if all([t in ("i", "j") for t in types]):
                    return "%".join(map(encode_item, forms, types)), "S"
            return number_form(p / Decimal(q))

        if isinstance(n, Decimal):
            if abs(n) < 10 ** Decimal(-20):
                return "0", "f"
            q = abs(n).log10().quantize(1, rounding=decimal.ROUND_FLOOR)
            q = n.quantize(10 ** (q - 20)).as_tuple()
            l = len(q.digits) + q.exponent
            q = "-" * q.sign + \
                ''.join(map(str, q.digits[:max(0, l)])) + \
                "." + "0" * max(0, -l) + \
                ''.join(map(str, q.digits[max(0, l):]))
            return q.rstrip("0").rstrip("."), "f"

        if isinstance(n, (mp.mpf, mp.mpc)) and not n.imag:
            return number_form(Decimal(mp.nstr(n.real, mp.mp.dps)))

        if complex_pair and isinstance(n, mp.mpc):
            return list_form([n.real, n.imag])

        if n is None:
            return "0n", "f"

        raise Exception()

    def expr_form(e, left=False):
        def op(name, a, b):
            return "%s%s%s" % (
                encode_item(*expr_form(a, left=True), left=True), name,
                encode_item(*expr_form(b))), "S"
        def func(name, x):
            if left:
                return "%s[%s]" % (name, encode_item(*expr_form(x))), "s"
            else:
                return "%s %s" % (name, encode_item(*expr_form(x))), "S"

        if complex_pair:
            a, b = e.as_real_imag()
            if b:
                return "(%s;%s)" % (encode_item(*expr_form(a)),
                                    encode_item(*expr_form(b))), "s"
        a, b = e.as_coeff_mul()
        if not b:
            return number_form(a)
        if abs(a.p) != 1 and abs(sp.Mul(*b).as_numer_denom()[0]) != 1:
            return op("*", a, sp.Mul(*b))
        if a == -1 and b:
            return func("neg", sp.Mul(*b))
        a, b = e.as_coeff_add()
        if a != 0 and b:
            if b[0].as_coeff_mul()[0] < 0:
                return op("-", a, -sp.Add(*b))
            else:
                return op("+", a, sp.Add(*b))
        a, b = e.as_numer_denom()
        if b != 1:
            for r in a.as_coeff_mul()[1]:
                r, p = r.as_base_exp()
                if p.is_Rational and p.p == 1 and p.q > 1:
                    r = r.as_coeff_mul()[0]
                    if r.q == 1 and b.as_coeff_mul()[0] % r.p == 0:
                        a, b = a / r.p**p, b / r.p**p
            if b.as_coeff_mul()[0] != 1:
                if a.as_coeff_mul()[0] == -1 or all(
                        x.as_coeff_mul()[0]<0 for x in a.as_coeff_add()[1]):
                    a, b = -a, -b
            return op("%", a, b)
        if e.func == sp.Add:
            a, b = e.as_two_terms()
            if b.as_coeff_mul()[0] < 0:
                return op("-", a, -b)
            else:
                return op("+", a, b)
        if e.func == sp.Mul:
            a, b = e.as_two_terms()
            return op("*", a, b)
        b, p = e.as_base_exp()
        if p != 1:
            if p.is_Rational:
                if not no_pow and p not in (S(1)/2, S(1)/3, 2, 3):
                    if b == sp.E:
                        return func(".qml.exp", p)
                    else:
                        return ".qml.pow[%s;%s]" % (
                            encode_item(*expr_form(b)),
                            encode_item(*expr_form(p))), "s"
                for x, name in ((3, "cbrt"), (2, "sqrt")):
                    if p.q % x == 0:
                        return func(".qml.%s" % name, e**x)
                if p.q != 1:
                    raise Exception
                assert 1 <= p.p <= 5
                f, t = expr_form(b)
                if t == "S":
                    s = "{"+"*".join(["x"]*p.p)+"} "+encode_item(f, t)
                else:
                    s = "*".join([encode_item(f, t, left=True)]*p.p)
                return s, "S"
            raise Exception
        if e == sp.E:
            return "e", "s"
        if e == sp.pi:
            return "pi", "s"
        if e == S.NaN:
            return "0n", "f"
        if not no_trigh and e.func in (sp.sinh, sp.cosh, sp.tanh):
            return func(".qml.%s" % e.func, *e.args)
        raise Exception

    def common_denominator(f):
        def mul(a, b):
            if isinstance(a, (list, tuple)):
                return [mul(x, b) for x in a]
            else:
                return a * b

        if isinstance(f, (list, tuple)):
            ns, ds, rs = zip(*map(common_denominator, f))
            if None not in ds and len(set(rs)) == 1:
                cd = sp.lcm(map(S, ds))
                return [mul(n, cd // d) for n, d in zip(ns, ds)], cd, rs[0]
        if rational(f):
            return rational(f) + (1,)
        if isinstance(f, sp.Basic):
            a, b = f.as_coeff_mul()
            if a.is_Rational and len(b) == 1 and b[0].is_Pow and \
                    b[0].base.is_Rational and b[0].exp == S(1)/2:
                return a.p, a.q, b[0].base
        return None, None, None

    def list_form(self, no_flip=False):
        if len(self) == 0:
            return "()", "v"

        if len(self) == 1 or len(set(value_form(e) for e in self)) == 1:
            f, t = value_form(self[0])
            if t in ("v", "V"):
                q = "enlist " + f
            else:
                q = encode_item(f, t)
            if len(self) != 1 or t not in ("v", "V"):
                q = str(len(self)) + "#" + q
            return q, "V"

        forms, types = zip(*map(value_form, self))
        if all([t in ("i", "j", "f") for t in types]):
            return encode_vector(forms, types), "v"

        # special case for Hilbert matrix
        if all(rational(e) and rational(e)[0] == 1 for e in self):
            hq = "1%"+encode_item(*list_form([rational(e)[1] for e in self]))
            return hq, "V"

        q, t = "(%s)" % ";".join(map(encode_item, forms, types)), "v"

        n, d, r = common_denominator(self)
        if d > 1 or r > 1:
            e, op = max(((d/sp.sqrt(r), "%"), (sp.sqrt(r)/d, "*")))
            cdq = op.join((encode_item(*list_form(n), left=True),
                           encode_item(*expr_form(e))))
            if len(cdq) < len(q):
                q, t = cdq, "V"

        if not no_flip and all(isinstance(x, (list, tuple))
                               for x in self) and len(set(map(len, self))) == 1:
            fq = "flip "+encode_item(*list_form(list(zip(*self)), no_flip=True))
            if len(fq) < len(q):
                q, t = fq, "V"

        return q, t

    def value_form(self):
        if rational(self) or \
                isinstance(self, (Decimal, mp.mpf, mp.mpc, NoneType)):
            return number_form(self)
        if isinstance(self, (sp.Basic)):
            return expr_form(self)
        if isinstance(self, (list, tuple)):
            return list_form(self)
        if isinstance(self, (sp.MatrixBase, mp.matrix)):
            return value_form(self.tolist())
        if hasattr(self, 'qform'):
            return self.qform(), "S"
        raise Exception()

    return encode_item(*value_form(self))


class qstr:
    def __init__(self, qform):
        self._qform = qform

    def qform(self):
        return self._qform


output_file = None

def output(s):
    (output_file or sys.stdout).write((s + "\n").decode("ascii"))

def prec(prec):
    output("    prec:%s;" % prec)

def reps(reps):
    output("    reps:%d;" % reps)

def test(name, *args, **kwargs):
    comment = kwargs.pop("comment", "")
    if comment:
        comment = " / %s" % comment

    args = [qform(x, **kwargs) for x in args]
    result = args.pop(-1)

    if "[" in name:
        call = ";%s]" % ";".join(args)
    elif len(args) == 1:
        call = " %s" % args[0]
    else:
        call = "[%s]" % ";".join(args)

    if "_" not in name:
        name = ".qml.%s" % name

    output('    test["%s%s";"%s"];%s' % (name, call, result, comment))
