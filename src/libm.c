#include <math.h>

#include "wrap.h"


wrap_F_ (cos)
wrap_F_ (sin)
wrap_F_ (tan)
wrap_F_ (acos)
wrap_F_ (asin)
wrap_F_ (atan)
wrap_FF_(atan2)
wrap_F_ (cosh)
wrap_F_ (sinh)
wrap_F_ (tanh)
wrap_F_ (acosh)
wrap_F_ (asinh)
wrap_F_ (atanh)
wrap_F_ (exp)
wrap_F_ (log)
wrap_F_ (log10)
wrap_F_ (logb)
wrap_F_ (expm1)
wrap_F_ (log1p)
wrap_FF_(pow)
wrap_F_ (floor)
wrap_F_ (ceil)
wrap_F_ (fabs)
wrap_FF_(fmod)
wrap_F_ (erf)
wrap_F_ (erfc)
wrap_F_ (lgamma)
wrap_F  (gamma, tgamma)
wrap_F_ (j0)
wrap_F_ (j1)
wrap_F_ (y0)
wrap_F_ (y1)
wrap_F_ (sqrt)
wrap_F_ (cbrt)
wrap_FF_(hypot)


static int
sgamma(F x) {
    if (x > 0)
        return 1;
    return fmod(x, 2) > -1 ? -1 : 1;
}

static F
beta(F x, F y) {
    return sgamma(x) * sgamma(y) * sgamma(x + y) *
       exp(lgamma(x) + lgamma(y) - lgamma(x + y));
}

wrap_FF_(beta)
