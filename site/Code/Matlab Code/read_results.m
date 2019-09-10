% Read in dataset structure and parameters

load('~/demand_final/paper-results/results-em.mat','t0');
t_em=t0;
load('~/demand_final/paper-results/results-ignore.mat','t0');
t_ig=t0;
load('~/demand_final/paper-results/results-full.mat','t0');
t_full = t0;

dj_f=dlmread('~/demand_final/NL_single/utils_full1.tab','\t',2,1);
xi_f=dlmread('~/demand_final/NL_single/xi_full1.tab','\t',2,1);
lambda_f = ones(5,1)*[0.582637 ];
tn_f1 = [dj_f; lambda_f; xi_f;];

dj_ig=dlmread('~/demand_final/NL_single/utils_ig1A.tab','\t',2,1);
xi_ig=dlmread('~/demand_final/NL_single/xi_ig1A.tab','\t',2,1);
lambda_ig = ones(5,1)*[0.697498];
tn_ig1 = [dj_ig; lambda_ig; xi_ig;];

dj_em=dlmread('~/demand_final/NL_single/utils_em4.tab','\t',2,1);
xi_em=dlmread('~/demand_final/NL_single/xi_em4.tab','\t',2,1);
lambda_em = ones(5,1)*[0.704003];
tn_em1 = [dj_em; lambda_em; xi_em;];


dj_f=dlmread('~/demand_final/NL_category/utils_full.tab','\t',2,1);
lambda_f=dlmread('~/demand_final/NL_category/lambda_full.tab','\t',2,1);
xi_f=dlmread('~/demand_final/NL_category/xi_full.tab','\t',2,1);
tn_f = [dj_f; lambda_f; xi_f;];

dj_ig=dlmread('~/demand_final/NL_category/utils_ig.tab','\t',2,1);
lambda_ig=dlmread('~/demand_final/NL_category/lambda_ig.tab','\t',2,1);
xi_ig=dlmread('~/demand_final/NL_category/xi_ig.tab','\t',2,1);
tn_ig = [dj_ig; lambda_ig; xi_ig;];

dj_em	 =dlmread('~/demand_final/NL_category/utils_em5.tab','\t',2,1);
lambda_em=dlmread('~/demand_final/NL_category/lambda_em5.tab','\t',2,1);
xi_em	 =dlmread('~/demand_final/NL_category/xi_em5.tab','\t',2,1);
tn_em = [dj_em; lambda_em; xi_em;];

xi_em_rc = t_em(47:end);
xi_ig_rc = t_ig(47:end);
xi_f_rc = t_full(47:end);

%save '~/demand_final/paper-results/results-all.mat' t_full t_ig t_em tn_f1 tn_ig1 tn_em1 tn_f tn_ig tn_em;

save '~/demand_final/paper-results/xi-all.mat' xi_em_rc xi_ig_rc xi_f_rc xi_f xi_em xi_ig;
