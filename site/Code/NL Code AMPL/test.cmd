reset;
model VendingNL.mod;

data vending_data_full.dat
solve;
display sum{m in M} xi[m];

table Utils2 OUT "utils_full.tab" :  [J], dj;
write table Utils2;
table Lambda2 OUT "lambda_full.tab" : [G], LAMBDA;
write table Lambda2;
table xi2 OUT "xi_full.tab" : [M], xi;
write table xi2;
table share2 OUT "share_full.tab": [ M, J], {(m,j) in AVAILABLE} pjt[m,j];
write table share2;

option solver gjh;
option gjh_options 'sparse';
write "hessian_full", solve;

