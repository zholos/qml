#!/usr/bin/env python
from __future__ import division

from fractions import Fraction
import math
import operator

from qform import qform


class Matrix:
    def __init__(self, rows):
        self.n = len(rows)
        self.m = len(rows[0])
        self.rows = rows

    def row(self, i):
        return self.rows[i]

    def column(self, j):
        return [self[i][j] for i in range(self.n)]

    def _assert_square(self):
        if self.m != self.n:
            raise Exception()

    def diagonal(self):
        return [self[i][i] for i in range(min(self.n, self.m))]

    def __getitem__(self, i):
        return self.row(i)

    def __cmp__(self, right):
        if isinstance(right, Matrix):
            return cmp(self.rows, right.rows)
        return 1

    def __mul__(self, right):
        if self.m != right.n:
            raise Exception()
        def cell(i, j):
            return sum([self[i][k] * right[k][j] for k in range(self.m)])
        return Matrix([[cell(i, j) for j in range(right.m)]
                                   for i in range(self.n)])

    def transpose(self):
        return Matrix(map(list, zip(*self.rows)))

    def take_lower(self):
        return Matrix([[self[i][j] if j <= i else 0 for j in range(self.m)]
                                                    for i in range(self.n)])

    def take_upper(self):
        return self.transpose().take_lower().transpose()

    def subst_solve(self, right):
        self._assert_square()
        if self.m != right.n:
            raise Exception()
        X = [None] * self.n
        for iteration in range(self.n):
            for i in range(self.n):
                if X[i] is None:
                    s = right[i][:]
                    for j in range(i) + range(i+1, self.n):
                        if self[i][j] != 0:
                            if X[j] is None:
                                break
                            for k in range(right.m):
                                s[k] -= self[i][j] * X[j][k]
                    else:
                        if self[i][i] == 0:
                            raise Exception()
                        X[i] = [s[k] / Fraction(self[i][i])
                                for k in range(right.m)]
                        break
        return Matrix(X)

    def _assert_row_echelon(self):
        zeroes = -1
        for i in range(self.n):
            for j in range(self.m):
                if self[i][j] != 0:
                    if j <= zeroes:
                        raise Exception()
                    zeroes = j
                    break
            else:
                zeroes = self.m

    def _assert_reduced_row_echelon(self):
        self._assert_row_echelon()
        for i0 in range(self.n):
            for j0 in range(self.m):
                if self[i0][j0] != 0:
                    if self[i0][j0] != 1:
                        raise Exception()
                    for i in range(i0) + range(i0+1, self.n):
                        if self[i][j0] != 0:
                            raise Exception()
                    break

    def _gauss(self):
        factor, rows = [1], [self[i][:] for i in range(self.n)]

        def switch(i0, i1):
            if i0 != i1:
                factor[0], rows[i0], rows[i1] = -factor[0], rows[i1], rows[i0]

        def multiply(i0, a):
            factor[0] /= Fraction(a)
            for j in range(self.m):
                rows[i0][j] *= a

        def subtract(i0, i1, a):
            for j in range(self.m):
                rows[i0][j] -= rows[i1][j] * a

        i0 = 0
        for j0 in range(self.m):
            if i0 < self.n:
                c = [abs(rows[i][j0]) for i in range(i0, self.n)]
                switch(i0, i0 + c.index(max(c)))
                if rows[i0][j0] != 0:
                    for i1 in range(i0+1, self.n):
                        subtract(i1, i0, rows[i1][j0] / Fraction(rows[i0][j0]))
                    i0 += 1

        multiply(0, factor[0])
        if factor[0] != 1:
            raise Excetpion()

        matrix = Matrix(rows)
        matrix._assert_row_echelon()
        return matrix

    def _jordan(self):
        self._assert_row_echelon()
        rows = [self[i][:] for i in range(self.n)]

        def divide(i0, a):
            for j in range(self.m):
                rows[i0][j] /= Fraction(a)

        def subtract(i0, i1, a):
            for j in range(self.m):
                rows[i0][j] -= rows[i1][j] * a

        for i0 in range(self.n):
            for j0 in range(self.m):
                if rows[i0][j0] != 0:
                    divide(i0, rows[i0][j0])
                    for i in range(i0):
                        subtract(i, i0, rows[i][j0])
                    break

        matrix = Matrix(rows)
        matrix._assert_reduced_row_echelon()
        return matrix

    def det(self):
        self._assert_square()
        return reduce(operator.mul, self._gauss().diagonal())

    def rank(self):
        g = self._gauss()
        return sum([max([1 if g[i][j] != 0 else 0 for j in range(self.m)])
                                                  for i in range(self.n)])

    def inverse(self):
        self._assert_square()
        if self.det() == 0:
            raise Exception()
        identity = Matrix.identity_matrix(self.n)
        augmented = Matrix([self[i] + identity[i] for i in range(self.n)])
        augmented = augmented._gauss()._jordan()
        inverse = Matrix([augmented[i][self.m:] for i in range(self.n)])
        if self * inverse != identity:
            raise Exception()
        return inverse

    def _rank_factorize(self):
        rre = self._gauss()._jordan()
        c, r = [], []
        for i in range(self.n):
            for j in range(self.m):
                if rre[i][j] != 0:
                    c.append(self.column(j))
                    r.append(rre[i])
                    break
        c, r = Matrix(c).transpose(), Matrix(r)

        if c * r != self:
            raise Exception()
        return c, r

    def pseudo_inverse(self):
        if self.rank() == min(self.n, self.m):
            transpose = self.transpose()
            if self.n < self.m:
                inverse = transpose * (self * transpose).inverse()
            else:
                inverse = (transpose * self).inverse() * transpose
        else:
            c, r = self._rank_factorize()
            inverse = r.pseudo_inverse() * c.pseudo_inverse()

        if self * inverse * self != self or inverse * self * inverse != inverse:
            raise Exception()
        return inverse

    @staticmethod
    def identity_matrix(n):
        return Matrix.diagonal_matrix([1] * n)

    @staticmethod
    def diagonal_matrix(d):
        return Matrix([[0]*i + [d[i]] + [0]*(len(d)-1-i)
                       for i in range(len(d))])

    @staticmethod
    def hilbert_matrix(n, m):
        matrix = Matrix([[Fraction(1, 1+i+j) for j in range(m)]
                                             for i in range(n)])
        matrix.preserve = True
        return matrix

    @staticmethod
    def random_matrix(n, m, scale = 1):
        def cell(i, j):
            return Fraction.from_float(math.tan((1+i*n+j)*scale)). \
                            limit_denominator(10)
        matrix = Matrix([[cell(i, j) for j in range(m)] for i in range(n)])
        matrix.preserve = True
        return matrix

    @staticmethod
    def null_matrix(n, m):
        return Matrix([[None] * m for i in range(n)])

    def qform(self):
        return qform(self.rows, preserve = getattr(self, 'preserve', False))


class Random: # meaning arbitrary, not unpredictable
    def __init__(self):
        self.i = 0

    def __call__(self, x):
        if isinstance(x, (list, tuple)):
            return x[self(len(x))]
        self.i += 1
        return self.i % x

    def pop(self, x):
        return x.pop(self(len(x)))

    def permute(self, x):
        r, x = [], x[:]
        while x:
            r.append(self.pop(x))
        return r


def test_prec(prec):
    output("    prec:%s;" % prec)

def test(name, *args, **kwargs):
    args, result = map(qform, args[:-1]), qform(args[-1])

    if "[" in name:
        call = ";%s]" % ";".join(args)
    elif len(args) == 1:
        call = " %s" % args[0]
    else:
        call = "[%s]" % ";".join(args)

    comment = kwargs.get("comment", "")
    if comment:
        comment = " / %s" % comment

    output('    test[".qml.%s%s";"%s"];%s' % (name, call, result, comment))


subjects = [
    #Matrix([[0]]),
    Matrix([[42]]),
    Matrix([[1, 2], [-3, 4]]),
    Matrix([[5, -6]]),
    Matrix([[7], [8]]),
    Matrix([[1, 2, 3], [4, -5, 6], [7, 8, 9]]),
    Matrix([[1, 2, 3], [4, -5, 6]]),
    Matrix([[1, 2], [4, -5], [7, 8]]),
    Matrix.hilbert_matrix(5, 5),
    Matrix.hilbert_matrix(3, 5),
    Matrix.hilbert_matrix(5, 3),
    Matrix.random_matrix(7, 7),
    Matrix.random_matrix(7, 4, math.sqrt(2)),
    Matrix.random_matrix(4, 7, math.sqrt(3)),
    Matrix([[1, 3], [2, 6]]),
    Matrix.random_matrix(3, 1, math.sqrt(5)),
    Matrix.random_matrix(1, 3, math.sqrt(7)),
    Matrix.random_matrix(4, 1, math.sqrt(11)),
    Matrix.random_matrix(1, 4, math.sqrt(13))
]

def test_diag():
    for d in ([1], [-2, 3], [-4, 5, 6]):
        test("diag", d, Matrix.diagonal_matrix(d))

def test_mdiag():
    for A in subjects:
        test("mdiag", A, A.diagonal())

def test_mdet():
    lower = True
    random = Random()
    for A in subjects:
        if A.m == A.n:
            det = A.det()
            test("mdet", A, det)
            if det:
                L = A.take_lower() if lower else A.take_upper()
                if L != A:
                    det = L.det()
                    test("mdet", L, det)
                    if det:
                        i = random(A.n)
                        L[i][i] = 0
                        if L.det() != 0:
                            raise Exception()
                        test("mdet", L, 0)
                lower = not lower

def test_mrank():
    done = []
    random = Random()
    for A in subjects:
        Z = Matrix([[0] * A.m for i in range(A.n)])
        if Z not in done:
            done.append(Z)
            test("mrank", Z, 0)

        for rank in range(1, A.n+1):
            rows, free = range(A.n), []
            for i in range(rank):
                free.append(random.pop(rows))
            B = []
            for i in range(A.n):
                if i in free:
                    B.append(A[i])
                else:
                    B.append(map(sum, zip(A[random(free)], A[random(free)])))
            B = Matrix(B)
            if B.rank() == rank:
                if B not in done:
                    done.append(B)
                    test("mrank", B, rank)

def test_minv():
    for A in subjects:
        if A.m == A.n:
            if A.det() != 0:
                B = A.inverse()
            else:
                B = Matrix.null_matrix(A.n, A.m)
            test("minv", A, B)

def test_mpinv():
    for A in subjects:
        test("mpinv", A, A.pseudo_inverse())

def test_mm():
    def emit(A, B, C):
        test("mm", A, B, C)
        if B.m == 1:
            test("mm", A, B.column(0), C.column(0))

    for A in subjects:
        for B in subjects:
            if A.m == B.n:
                C = A * B
                emit(A, B, C)

def test_ms():
    def emit(A, B, X):
        test("ms", A, B, X)
        if B.m == 1:
            test("ms", A, B.column(0), X.column(0))

    single_done = 0
    zero_done = 0
    random = Random()
    for A in subjects:
        for B in subjects:
            if A.m == A.n == B.n:
                if A.n == 1:
                    if single_done == 3:
                        continue
                    single_done += 1

                L = [A.take_lower(), A.take_upper()]
                if L[0] == L[1]:
                    L.pop()
                i = random(A.n)
                for Z in L:
                    emit(Z, B, Z.subst_solve(B))
                    if zero_done == 10:
                        continue
                    zero_done += 1
                    Z[i][i] = 0
                    emit(Z, B, Matrix.null_matrix(B.n, B.m))

def test_mls(equi):
    routine = "mlsx[`equi" if equi else "mls"
    def emit(A, B, X):
        test(routine, A, B, X)
        if B.m == 1:
            test(routine, A, B.column(0), X.column(0))

    if not equi:
        test(routine, [[0]], [1], [None], comment="ATLAS dgetf2() return code")

    single_done = 0
    random = Random()

    for A in subjects:
        for B in subjects:
            if A.m == A.n == B.n:
                if A.n == 1:
                    if single_done == 3:
                        continue
                    single_done += 1

                if A.det() != 0:
                    X = A.inverse() * B
                else:
                    X = Matrix.null_matrix(B.n, B.m)

                emit(A, B, X)

                if equi and A.det() != 0 and A.m >= 3 and A.m <= 5:
                    L, R = [Matrix.diagonal_matrix(random.permute(
                            [1000 ** i for i in range(A.m)])) for i in range(2)]

                    emit(L * A,     L * B,               X)
                    emit(    A * R,     B, R.inverse() * X)
                    emit(L * A * R, L * B, R.inverse() * X)

def test_mlsq(svd):
    routine = "mlsqx[`svd" if svd else "mlsq"
    def emit(A, B, X):
        test(routine, A, B, X)
        if B.m == 1:
            test(routine, A, B.column(0), X.column(0))

    single_done = 0
    for A in subjects:
        for B in subjects:
            if A.n == B.n:
                if A.n == 1:
                    if single_done == 3:
                        continue
                    single_done += 1

                X = A.pseudo_inverse() * B

                if not svd and A.rank() != min(A.n, A.m):
                    # Replace matrix by more-obviously-rank-deficient one
                    # as routine is poor at detecting them.
                    A = Matrix.diagonal_matrix([1] * (A.n - 1) + [0]) * A
                    X = Matrix.null_matrix(B.n, B.m)

                emit(A, B, X)

def test_mkron():
    output("""\
    test[".qml.mkron[(1 2;-3 4);-1 2]";"(-1 -2;2 4;3 -4;-6 8)"];
    test[".qml.mkron[(0 1 0;-2.5 0 3);(1 2;-3 4)]";"(0 0 1 2 0 0;0 0 -3 4 0 0;-2.5 -5 0 0 3 6;7.5 -10 0 0 -9 12)"];""")

def tests():
    test_prec("1e-9")
    test_diag()
    test_mdiag()
    test_mdet()
    test_mrank()
    test_minv()
    test_mpinv()
    test_mm()
    test_prec("1e-8")
    test_ms()
    test_mls(False)
    test_mls(True)
    test_prec("1e-7")
    test_mlsq(False)
    test_mlsq(True)
    test_mkron()


if __name__ == "__main__":
    def output(s):
        print s
    tests()
