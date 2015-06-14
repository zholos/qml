\d .qml

dll:`qml^dll^:`; / optional override
{x[0] set @[dll 2:;1_x;
    {'x}"qml: loading ",string[x 1],"(): ",]} each flip
    {(r;`$"qml_",/:string r:raze x;where count each x)}

    1 2 3 4!(
        `sin`cos`tan`asin`acos`atan`sinh`cosh`tanh`asinh`acosh`atanh,`exp`expm1,
            `log`log10`logb`log1p`sqrt`cbrt`floor`ceil`fabs`erf`erfc,
            `lgamma`gamma`j0`j1`y0`y1`ncdf`nicdf`kcdf`kicdf,
            `mdet`minv`mevu`mchol`mqr`mqrp`mlup`msvd`poly`const;
        `atan2`pow`hypot`fmod`beta`pgammar`pgammarc`ipgammarc`c2cdf`c2icdf,
            `stcdf`sticdf`pscdf`psicdf`smcdf`smicdf`mm`ms`mls`mlsq,
            `solve`min`root;
        `pbetar`ipbetar`fcdf`ficdf`gcdf`gicdf`bncdf`bnicdf`mlsx`mlsqx,
            `solvex`minx`rootx`conmin`line;
        `conminx`linex);

`version`pi`e`eps set'const each til 4;
pgamma:{gamma[x]*pgammar[x;y]};
pgammac:{gamma[x]*pgammarc[x;y]};
pbeta:{beta[x;y]*pbetar[x;y;z]};
diag:{@[count[x]#abs[type x]$0;;:;]'[til count x;x]};
mdim:{(count x;count x 0)};
mdiag:{(n#x)@'til n:min mdim x};
mrank:{sum not (d<eps*d[0]*max mdim x)|0=d:mdiag msvd[x]1};
mpinv:{mm[x 2] flip mm[x 0]
    ?'[(s=0)|s<eps*s[0;0]*max mdim s;s*0;reciprocal s:(x:msvd x)1]};
mev:{x@\:idesc sum each {x*x} first x:mevu x};
mkron:{raze(raze'')(flip')x*\:\:y};
