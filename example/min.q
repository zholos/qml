/ Different optimization methods

\l qml.q

/ generate system of equations (as in ls.q but much smaller)
n:15;
while[0=.qml.mdet A:n cut .qml.nicdf .005+(n*n)?.99];
b:.qml.nicdf .005+n?.99;

/ solve with appropriate function for comparison
1"mls:        ";
\t x0:.qml.mls[A;b];

/ solve with generic function
1"solve:      ";
\t x1:first .qml.solve[b{x-y$}'A;enlist n#0];

/ solve with optimization as least-squares problem
opts:`iter,1000000;
1"min:        ";
\t x2:first .qml.minx[opts][{x$x}b-A mmu;enlist n#0];
1"minx`nm:    ";
\t x3:first .qml.minx[`nm,opts][{x$x}b-A mmu;enlist n#0];
1"minx`sbplx: ";
\t x4:first .qml.minx[`sbplx,opts][{x$x}b-A mmu;enlist n#0];

/ check results
if[1e-7<{max -1+(b%x)|x%b}A mmu x0;'`incorrect];
if[1e-3<max{max -1+(x0%x)|x%x0}each(x1;x2;x3;x4);'`different];

-1"";

/ generate overdetermined system of equations (as in lsq.q but much smaller)
n:50;  / variables
m:200; / equations
A:n cut .qml.nicdf .005+(m*n)?.99;
b:.qml.nicdf .005+m?.99;

/ solve with appropriate function for comparison
1"mlsq:       ";
\t x0:.qml.mlsq[A;b];

/ solve with optimization
1"min:        ";
\t x1:first .qml.minx[opts][{x$x}b-A mmu;enlist n#0];
1"minx`nm:    ";
\t x2:first .qml.minx[`nm,opts][{x$x}b-A mmu;enlist n#0];
1"minx`sbplx: ";
\t x3:first .qml.minx[`sbplx,opts][{x$x}b-A mmu;enlist n#0];

/ check results
if[1e-7<{max -1+(y%x)|x%y}[flip[A]mmu A mmu x0;flip[A]mmu b];
   '`incorrect];
if[1e-3<max{max -1+(x0%x)|x%x0}each(x1;x2;x3);'`different];
