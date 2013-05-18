#ifndef QML_PKG_LAPACK_H
#define QML_PKG_LAPACK_H

// Fortran BLAS prototypes
int dgemv_(char* trans, int* m, int* n,
           double* alpha, double* a, int* lda,
           double* x, int* incx, double* beta, double* y, int* incy);
int dgemm_(char* transa, char* transb, int* m, int* n, int* k,
           double* alpha, double* a, int* lda,
           double* b, int* ldb, double* beta, double* c, int* ldc);
int dtrsv_(char* uplo, char* trans, char* diag, int* n,
           double* a, int* lda, double* x, int* incx);
int dtrsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n,
           double* alpha, double* a, int* lda, double* b, int* ldb);

// Fortran LAPACK prototypes
int dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
int dgetri_(int* n, double* a, int* lda, int* ipiv,
            double* work, int* lwork, int* info);
int dgeqrf_(int* m, int* n, double* a, int* lda,
            double* tau, double* work, int* lwork, int* info);
int dgeqp3_(int* m, int* n, double* a, int* lda, int* jpvt,
            double* tau, double* work, int* lwork, int* info);
int dorgqr_(int* m, int* n, int* k, double* a, int* lda,
            double* tau, double* work, int* lwork, int* info);
int dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
int dgesdd_(char* jobz, int* m, int* n, double* a, int* lda,
            double* s, double* u, int* ldu, double* vt, int* ldvt,
            double* work, int* lwork, int* iwork, int* info);
int dgesv_ (int* n, int* nrhs, double* a, int* lda, int* ipiv,
            double* b, int* ldb, int* info);
int dgesvx_(char* fact, char* trans, int* n, int* nrhs,
            double* a, int* lda, double* af, int* ldaf, int* ipiv,
            char* equed, double* r, double* c,
            double* b, int* ldb, double* x, int* ldx,
            double* rcond, double* ferr, double* berr,
            double* work, int* iwork, int* info);
int dgels_ (char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
            double* b, int* ldb, double* work, int* lwork, int* info);
int dgelsd_(int* m, int* n, int* nrhs,
            double* a, int* lda, double* b, int* ldb,
            double* s, double* rcond, int* rank,
            double* work, int* lwork, int *iwork, int* info);
int dgeev_ (char* jobvl, char* jobvr, int* n, double* a, int* lda,
            double* wr, double* wi_,
            double* vl, int* ldvl, double* vr, int* ldvr,
            double* work, int* lwork, int* info);
int zgeev_ (char* jobvl, char* jobvr, int* n, double* a, int* lda,
            double* w,
            double* vl, int* ldvl, double* vr, int* ldvr,
            double* work, int* lwork, double* rwork, int* info);

#endif
