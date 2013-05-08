#include <cprob.h>

#include "wrap.h"


// chdtri() returns the inverse of the complementary CDF
static F
c2icdf(F x, F y) {
    return chdtri(x, 1 - y);
}

// fdtri() returns the inverse of the complementary CDF
static F
ficdf(I x, I y, F z) {
    return fdtri(x, y, 1 - z);
}

// gdtr() interprets second parameter differently
static F
gcdf(F x, F y, F z) {
    return igam(x, z / y);
}

// no gdtri() function
static F
gicdf(F x, F y, F z) {
    return y * igami(x, 1 - z);
}

// smirnov() returns complementary CDF
static F
smcdf(I x, F y) {
    return 1 - smirnov(x, y);
}

// smirnovi() returns inverse of the complementary CDF
static F
smicdf(I x, F y) {
    return smirnovi(x, 1 - y);
}

// komogorov() returns complementary CDF
static F
kcdf(F x) {
    return 1 - kolmogorov(x);
}

// komogoi() returns inverse of the complementary CDF
static F
kicdf(F x) {
    check(x>=1e-8, nf,); // doesn't converge well for small values
    return kolmogi(1 - x);
}


wrap_fF (pgammar,   igam)
wrap_fF (pgammarc,  igamc)
wrap_fF (ipgammarc, igami)
wrap_ffF(pbetar,    incbet)
wrap_ffF(ipbetar,   incbi)
wrap_F  (ncdf,      ndtr)
wrap_F  (nicdf,     ndtri)
wrap_iF (stcdf,     stdtr)
wrap_iF (sticdf,    stdtri)
wrap_iiF(fcdf,      fdtr)
wrap_iiF(ficdf,     ficdf)
wrap_fF (c2cdf,     chdtr)
wrap_fF (c2icdf,    c2icdf)
wrap_ffF(gcdf,      gcdf)
wrap_ffF(gicdf,     gicdf)
wrap_iiF(bncdf,     bdtr)
wrap_iiF(bnicdf,    bdtri)
wrap_iF (pscdf,     pdtr)
wrap_iF (psicdf,    pdtri)
wrap_iF (smcdf,     smcdf)
wrap_iF (smicdf,    smicdf)
wrap_F_ (kcdf)
wrap_F_ (kicdf)
