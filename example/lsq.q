/ Different ways to solve a linear least squares problem

\l qml.q

/ generate overdetermined system of equations
n:500;  / variables
m:2000; / equations
A:n cut .qml.nicdf .005+(m*n)?.99;
b:.qml.nicdf .005+m?.99;


/ OLS estimator
1"OLS:           ";
\t x0:.qml.mm[.qml.minv .qml.mm[flip A]A].qml.mm[A_:flip A]b;


/ pseudoinverse
1"pseudoinverse: ";
\t x1:.qml.mm[.qml.mpinv A]b;

/ pseudoinverse by QR factorization
1"QR:            ";
\t x2:.qml.mm[.qml.ms[R].qml.ms[flip R:count[A 0]#.qml.mqr[A]1]flip A]b;


/ solve with appropriate function
1"mlsq:          ";
\t x3:.qml.mlsq[A;b];

/ solve with appropriate function using SVD algorithm
1"mlsq`svd:      ";
\t x4:.qml.mlsqx[`svd][A;b];


/ check results
if[1e-7<{max -1+(y%x)|x%y}[.qml.mm[A_].qml.mm[A]x0;.qml.mm[A_:flip A]b];
   '`incorrect];
if[1e-5<max{max -1+(x0%x)|x%x0}each(x1;x2;x3;x4);'`different];
