option show_stats 1;
option gentimes 1;
option knitro_options 'alg=1 outlev=3 opttol=1.0e-8 feastol=1.0e-8 maxit=250';

param nmkt ; # number of markets
param nprod ; # number of brands per market
param ngroup ; # number of category nests

# Define sets/ indices
set M := 1..nmkt ;  # index set of market
set J ordered := 1..nprod ;  # index set of brand (products), excluding outside good
set J1 = J union {nprod+1};
set GROUP ordered= 1..ngroup;

## Define input data:
param Q {m in M, j in J1} ;
param Q0 {m in M} = Q[m,nprod+1];
param AV {m in M, j in J};
param Market {m in M} = sum{j in J1} Q[m,j];

#param D0 {j in J, l in L};
param category {j in J, g in GROUP}, integer;
param group {j in J}, integer;
param groupsales {m in M, g in GROUP}  = sum{j in J: Q[m,j] > 0} Q[m,j] * category[j,g];
param groupavail {m in M, g in GROUP}  = max{j in J: group[j] = g} AV[m,j];
param Mtotal = sum {m in M} Q0[m];

#### Define variables
var LAMBDA {g in GROUP} >=.01, <=0.99;
var dj {j in J};
var xi {m in M};

set AVAILABLE := {m in M, j in J : AV[m,j] > 0};

# var utils {j in J, m in M} = dj[j] 
# Define Equations / Variables
var LIV{m in M, g in GROUP : groupavail[m,g] > 0} = log(sum{j in J : AV[m,j] > 0} category[j,g] *exp(dj[j]))
var IVS{m in M} = 1.0 + sum{g in GROUP: groupavail[m,g] > 0} exp(LIV[m,g] * LAMBDAONE);
var inside{m in M, j in J: Q[m,j] > 0} =  Q[m,j] * (dj[j]- sum{g in GROUP: groupavail[m,g] > 0} LIV[m,g] * category[j,g]);
var groupshare {m in M, g in GROUP : groupavail[m,g] > 0}  = LIV[m,g] * LAMBDAONE - log(IVS[m]);
var outside{ m in M} = -log(IVS[m]);

var IV{m in M, g in GROUP} = sum{j in J: (m,j) in AVAILABLE} exp(utils[j,m])*category[j,g];
var IVL{m in M, g in GROUP}= IV[m,g]^LAMBDAONE;
var denom{m in M} = (xi[m] + sum{g in GROUP} IVL[m,g]);
var pjt {(m,j) in AVAILABLE} = (exp(utils[j,m]) / IV[m,group[j]]) * IVL[m,group[j]] / denom[m];
var pin{m in M} = sum{j in J: (m,j) in AVAILABLE} pjt[m,j];

maximize MSL:  (sum{m in M, g in GROUP: groupavail[m,g] > 0} groupsales[m,g] * groupshare[m,g]) + (sum { m in M} outside[m] * Q0[m]) + (sum{m in M, j in J: Q[m,j] > 0} inside[m,j]); 
subject to addup: sum{m in M} xi[m] = 0;
