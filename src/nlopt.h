#ifndef QML_SRC_NLOPT_H
#define QML_SRC_NLOPT_H

#include "util.h"

K nloptmin(K fun, K con, K start_, I maxiter, F tolcon,
           I method, I lincon, I full, I quiet);

#endif // QML_SRC_NLOPT_H
