clear
load dsarray2-mach.mat;
for k=1:length(dsarray)
        N(k) = size(dsarray(k).avail,2);
end
[J K] = size(dsarray(2).Xmat);

% NESTED LOGIT ADJUSTMENT
K = max(dsarray(2).group);

NN = [0 cumsum([N])];
a1=[1:NN(2)]';
a2=repvec([NN(2)+1:NN(3)],2);
a3=repvec([NN(3)+1:NN(4)],4);
a4=repvec([NN(4)+1:NN(5)],8);
mktid=[a1;a2;a3;a4];


% load the parameter starting values and expand for the extra EM observations

dj_f=dlmread('utils_1.tab','\t',2,1);
lambda_f=dlmread('lambda_1.tab','\t',2,1);
xi_f=dlmread('xi_1.tab','\t',2,1);
xi=accumarray([mktid],xi_f,[max(mktid) 1],@max);

t0 = [dj_f; lambda_f; xi];

ds3=estepn(dsarray,t0);

% machine fe correction
ds3.mktid=ds3.mach_id(ds3.mktid);

save esteppedN-2.mat ds3 t0;
dstoampl(ds3,'vending_em2.dat');
