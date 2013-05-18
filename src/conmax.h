#ifndef QML_SRC_CONMAX_H
#define QML_SRC_CONMAX_H

#include "util.h"

K solvemin(K fun, K con, K start_, I maxiter, F tolcon, I steps,
           int slp, int rk, int lincon, int full, int quiet);

K root(K fun, K start, I maxiter, F tolcon, int full, int quiet);

K line(K fun, K base, K start, I maxiter, F tolcon, int full, int quiet);


#endif // QML_SRC_CONMAX_H
