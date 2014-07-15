/ Different ways to solve a system of linear equations

\l qml.q

/ generate system of equations
n:1000;
while[0=.qml.mdet A:n cut .qml.nicdf .005+(n*n)?.99];
b:.qml.nicdf .005+n?.99;


/ solve by matrix inversion
/   Ax=b  =>  x=inv(A)b
1"inversion: ";
\t x0:.qml.mm[.qml.minv A]b;


/ solve by QR factorization
/   Ax=b  <=>  QRx=b  =>  Rx=Q'b (R is triangular, solved by back substitution)
1"QR:        ";
\t x1:.qml.ms[QR 1].qml.mm[flip (QR:.qml.mqr A)0]b;

/ solve by QR factorization with column pivoting
/   APP'x=b  <=>  QRP'x=b  <=>  RP'x=Q'b  =>  Ry=Q'b, x=Py
1"QRP:       ";
\t x2:.qml.ms[QRP 1;.qml.mm[flip QRP 0]b]iasc(QRP:.qml.mqrp A)2;

/ solve by LUP factorization
/   Ax=b  <=>  PAx=Pb  <=>  LUx=Pb  <=>  Ly=Pb, Ux=y (L and U are triangular)
1"LUP:       ";
\t x3:.qml.ms[LUP 1].qml.ms[LUP 0]b(LUP:.qml.mlup A)2;

/ solve by Cholesky factorization
/   Ax=b  =>  A'Ax=A'b  <=>  R'Rx=A'b  <=>  R'y=A'b, Rx=y (R is triangular)
1"Cholesky:  ";
\t x4:.qml.ms[R].qml.ms[flip R:.qml.mchol .qml.mm[A_;A]].qml.mm[A_:flip A]b;

/ solve by singular value decomposition
/   Ax=b  <=>  USV'x=b  =>  x=V inv(S)U'b (S is diagonal)
1"SVD:       ";
\t x5:.qml.mm[SVD 2].qml.mm[flip SVD 0;b]%.qml.mdiag(SVD:.qml.msvd A)1;


/ solve with appropriate function
1"mls:       ";
\t x6:.qml.mls[A;b];

/ solve with appropriate function, equilibrating
1"mlsx`equi: ";
\t x7:.qml.mlsx[`equi][A;b];


/ solve as least-squares problem
1"mlsq:      ";
\t x8:.qml.mlsq[A;b];

/ solve as least-squares problem using SVD algorithm
1"mlsq`svd:  ";
\t x9:.qml.mlsqx[`svd][A;b];


/ check results
if[1e-7<{max -1+(b%x)|x%b}.qml.mm[A]x0;'`incorrect];
if[1e-5<max{max -1+(x0%x)|x%x0}each(x1;x2;x3;x4;x5;x6;x7;x8;x9);'`different];
