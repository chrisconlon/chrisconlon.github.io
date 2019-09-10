include cat_ignore;
param COEF= 39601;

for {i in 1..COEF} {
   print {j in 1..COEF} H[i,j] > hessian_ignore.txt;
}


