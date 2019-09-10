model VendingNL-EM.mod;

data vending_em5.dat
solve;
display {g in GROUP} LAMBDA[g];

table Utils2 OUT "utils_em5.tab" :  [J], dj;
write table Utils2;
table Lambda2 OUT "lambda_em5.tab" : [G], LAMBDA;
write table Lambda2;
table xi OUT "xi_em5.tab" : [NFE], xi2;
write table xi;

#option solver gjh;
#option gjh_options 'sparse';
#write "hessian_em", solve;

