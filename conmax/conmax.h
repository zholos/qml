#ifndef QML_PKG_CONMAX_H
#define QML_PKG_CONMAX_H

int conmax_(int* ioptn, int* nparm, int* numgr, int* itlim,
            double* fun, int* ifun, double* pttbl,
            int* iwork, int* liwrk, double* work, int* lwrk,
            int* iter, double* param, double* error);
int muller_(int* limmul, int* nsrch, int* ioptn, int* nparm, int* numgr,
            double* dvec, double* fun, int* ifun, double* pttbl,
            double* zwork, double* tolcon, int* iphse,
            int* iwork, int* liwrk, double* work, int* lwrk, double* parwrk,
            double* err1, double* p1, double* f1, double* procor, double* emin);
int searsl_(int* initlm, int* nadd, int* lims1, int* ioptn,
            int* numgr, int* nparm, double* prjlim, double* tol1,
            double* x, double* fun, int* ifun, double* pttbl,
            double* param, double* error,
            double* rchdwn, int* mact, int* iact, int* iphse, double* unit,
            double* tolcon, double* rchin, int* itypm1, int* itypm2,
            int* iwork, int* liwrk, double* work, int* lwrk,
            double* err1, double* parprj, double* projct,
            double* emin, double* emin1, double* parser, int* nsrch);

#endif
