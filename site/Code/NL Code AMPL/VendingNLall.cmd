model VendingNL.mod;
data vending_data_ignore.dat
solve;
display sum{m in M} xi[m];

table Utils3 OUT "utils_ig.tab" :  [J], dj;
write table Utils3;
table Lambda3 OUT "lambda_ig.tab" : [G], LAMBDA;
write table Lambda3;
table xi3 OUT "xi_ig.tab" : [M], xi;
write table xi3;
table share3 OUT "share_ig.tab": [ M, J], {(m,j) in AVAILABLE} pjt[m,j];
write table share3;


option solver gjh;
option gjh_options 'sparse';
write "bhessianignore";
solve;


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
write "bhessianfull";
solve;

