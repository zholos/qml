/ Portfolio optimization

\l qml.q

/ generate expected returns (r) and covariance matrix (S)
n:40;
bates:{[a;b;m;n]avg a+n?'m#b-a};
rho:{1^x^flip x} (sums[0,til n-1]_bates[-.2;.9;5](n*n-1)div 2)@\:til n;
sigma:bates[0;12.;3] n;
S:rho*sigma*/:sigma;
r:bates[-1;1.1;7] n;

/ find a portfolio on the efficient frontier (parametrized by q)
q:1.;
f:{(x mmu S mmu x)-q*r mmu x};
cons:{@[;x]}'[til n],{sum[x]-1},{1-sum x}; / 0<=x, sum[x]=1

1"conmax:        ";
\t x0:first .qml.conminx[`lincon;f;cons;enlist n#0];

/ analytic gradients
df:{enlist (2*S mmu x)-q*r};
dcons:{:[;enlist x]}each .qml.diag[n#1.],n#/:1 -1.;

1"conmax':       ";
\t x1:first .qml.conminx[`lincon;(f;df);flip(cons;dcons);enlist n#0];

/ substitute equality constraint
f1:{(x mmu S mmu x)-q*r mmu x,:1-sum x};
cons1:{@[;x]}'[til n-1],{1-sum x};

df1:{enlist {(-1_x)-last x} (2*S mmu x,1-sum x)-q*r};
dcons1:{:[;enlist x]}each .qml.diag[(n-1)#1.],enlist (n-1)#-1.

1"conmax(subs):  ";
\t x2:{x,1-sum x} first .qml.conminx[`lincon;f1;cons1;enlist (n-1)#0];
1"conmax(subs)': ";
\t x3:{x,1-sum x} first .qml.conminx[`lincon;f1,df1;flip(cons1;dcons1);enlist (n-1)#0];

1"cobyla:        ";
\t x4:first .qml.conminx[`lincon`cobyla;f;cons;enlist n#0];
1"cobyla(subs):  ";
\t x5:{x,1-sum x} first .qml.conminx[`lincon`cobyla;f1;cons1;enlist (n-1)#0];

/ check results
if[1e-5<max{max abs x-x0}each(x1;x2;x3;x4;x5);'`different];
