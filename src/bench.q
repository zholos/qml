if[any"-cd"~/:.z.x;
    .qml.dll:` sv hsym[`$system"cd"],`qml]; / full path stops Windows searching
raw:any"-raw"~/:.z.x;
speed:$[any"-fast"~/:.z.x;0;any"-slow"~/:.z.x;2;1];
patterns:.z.x 1+where"-like"~/:.z.x;

\l qml.q
-1"qml ",string .qml.version;

/ \t is also susceptible to wall clock jumps, so use .z.p for higher precision.
/ If call is fast, iterate it multiple times between measurements - for at least
/ one second. Then take mean of multiple samples.

iters:0D.1 0D00:00:01 0D00:00:02 speed;
sample:3 5 10 speed;

time:{
    / Evaluate each argument, leaving only the top-level call for timing.
    if[10h=type x;
        if[not(0h=type x)and -11h=type first x:parse x;'`type];
        x:eval each x];
    / Multiply iterations by 10 until it takes about 50 ms, then compute
    / approximate iterations for 1 second.
    {if[(iters div 20)>t:x y;:.z.s[x;10*y]];
        if[iters>t;t:x y*:1+"j"$iters%t];
        y,(t%0D.001)%y}[{.z.p-do[z;x . y].z.p}[x 0;1_x];1j]};

bench:{
    if[not $[count patterns;any x like/:patterns;1b];:(::)];
    1 x,$[raw;"";(count[x]_47#" ")],"\t";
    if[@[{get x;0b};x;{-1"'",x;1b}];:(::)]; / prime it and skip on error
    t:flip time each sample#enlist x;
    i:.qml.sticdf[sample-1;.975]*dev[t 1]%sqrt sample; / 95% confidence interval
    -1 $[raw;
        string[avg t 1],"\t",string i;
        .Q.f[6;avg t 1]," +- ",.Q.f[6;i]," ms (",string["j"$avg t 0]," iters)"];
    };

n:3 10 100 1000;
A:n!{x cut sin sqrt 1+til x*x} each n;
B:n!{x cut cos sqrt 1+til x*x} each n;
a:first each A;
b:{x[;0]}each B;
AS:{.qml.mm[x] flip x} each A;
AU:{.qml.mlup[x]1} each A;
ml:(4 3;3 4;6 10;11 5;60 100;110 50;600 1000;1100 500);
mr:(3 2;4 2;10 4;5 9;100 40;50 90;1000 400;500 900);
A_,:{reverse'[key x]!flip'[value x]}A_:ml!{x#sin sqrt 1+til prd x} each ml;
B_,:{reverse'[key x]!flip'[value x]}B_:mr!{x#sin sqrt 1+til prd x} each mr;
b_:{first'[key x]!value[x][;;0]}B_;

bench1:{bench each (-3!'y)sv\:"{}"vs x};
bench2:{bench each (raze raze flip ("{}"vs x;)@)each (-3!'y)(;;"")'(-3!'z)};

bench1[".qml.mdet A {}";n];
bench1[".qml.minv A {}";n];
bench1[".qml.dot[a {};b {}]";n];
bench1[".qml.mm[A {};B {}]";n];
bench1[".qml.mm[A {};b {}]";n];
bench1[".qml.mm[1#A {};B {}]";n];
bench2[".qml.mm[A_ {};B_ {}]";ml;mr];
bench1[".qml.mmx[`lflip;A {};B {}]";n];
bench1[".qml.mmx[`lflip;A {};b {}]";n];
bench1[".qml.mmx[`rflip;A {};B {}]";n];
bench1[".qml.mmx[`lflip`rflip;A {};B {}]";n];
bench1[".qml.ms[AU {};B {}]";n];
bench1[".qml.ms[AU {};b {}]";n];
bench1[".qml.mevu A {}";n];
bench1[".qml.mchol AS {}";n];
bench1[".qml.mqr A {}";n];
bench1[".qml.mqr A_ {}";ml];
bench1[".qml.mqrp A {}";n];
bench1[".qml.mqrp A_ {}";ml];
bench1[".qml.mlup A {}";n];
bench1[".qml.mlup A_ {}";ml];
bench1[".qml.msvd A {}";n];
bench1[".qml.msvd A_ {}";ml];
bench1[".qml.mls[A {};B {}]";n];
bench1[".qml.mls[A {};b {}]";n];
bench1[".qml.mlsx[`equi;A {};B {}]";n];
bench1[".qml.mlsx[`equi;A {};b {}]";n];
bench1[".qml.mlsx[`flip;A {};B {}]";n];
bench1[".qml.mlsx[`flip;A {};b {}]";n];
bench1[".qml.mlsx[`equi`flip;A {};B {}]";n];
bench1[".qml.mlsx[`equi`flip;A {};b {}]";n];
bench1[".qml.mlsq[A {};B {}]";n];
bench1[".qml.mlsq[A {};b {}]";n];
bench2[".qml.mlsq[A_ {};B_ {}]";reverse'[ml];mr];
bench2[".qml.mlsq[A_ {};b_ {}]";reverse'[ml];mr[;0]];
bench1[".qml.mlsqx[`svd;A {};B {}]";n];
bench1[".qml.mlsqx[`svd;A {};b {}]";n];
bench2[".qml.mlsqx[`svd;A_ {};B_ {}]";reverse'[ml];mr];
bench2[".qml.mlsqx[`svd;A_ {};b_ {}]";reverse'[ml];mr[;0]];
bench2[".qml.mlsqx[`flip;A_ {};B_ {}]";ml;reverse'[mr]];
bench2[".qml.mlsqx[`flip;A_ {};b_ {}]";ml;mr[;0]];
bench1[".qml.mnoopx[`square;A {}]";n];
bench1[".qml.mnoopx[`square`flip;A {}]";n];
bench1[".qml.mnoopx[`square`flip`triangular;A {}]";n];
bench1[".qml.mnoop A_ {}";ml];
bench1[".qml.mnoopx[`flip;A_ {}]";ml];
bench1[".qml.mnoopx[`square`upper;A {}]";n];
bench1[".qml.mnoopx[`square`lower;A {}]";n];
bench1[".qml.mnoopx[`square`upper`flip;A {}]";n];
bench1[".qml.mnoopx[`square`lower`flip;A {}]";n];
bench1[".qml.mnoopx[`upper;A_ {}]";ml];
bench1[".qml.mnoopx[`lower;A_ {}]";ml];
bench1[".qml.mnoopx[`upper`flip;A_ {}]";ml];
bench1[".qml.mnoopx[`lower`flip;A_ {}]";ml];

\\
