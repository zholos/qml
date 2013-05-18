#ifndef QML_PKG_CPROB_H
#define QML_PKG_CPROB_H

double igam(double, double);
double igamc(double, double);
double igami(double, double);
double incbet(double, double, double);
double incbi(double, double, double);
double ndtr(double);
double ndtri(double);
double stdtr(int, double);
double stdtri(int, double);
double fdtr(int, int, double);
double fdtrc(int, int, double);
double fdtri(int, int, double);
double chdtrc(double, double);
double chdtr(double, double);
double chdtri(double, double);
double bdtr(int, int, double);
double bdtri(int, int, double);
double pdtr(int, double);
double pdtri(int, double);
double smirnov(int, double);
double kolmogorov(double);
double smirnovi(int, double);
double kolmogi(double);

#endif
