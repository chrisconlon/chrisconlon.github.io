%{
clear
load ~/demand_final/paper-results/results-full
load ~/demand_final/paper-results/results-all
sen_f=nl_stderr(ds3,tn_f);
clear ds3 t0
save results-all

clear
load ~/demand_final/paper-results/results-ignore
load ~/demand_final/paper-results/results-all
sen_ig=nl_stderr(ds3,tn_ig);
clear ds3 t0
save results-all

clear
load ~/demand_final/paper-results/results-emN
load ~/demand_final/paper-results/results-all
sen_em=nl_stderr(ds3,tn_em);
clear ds3 t0
save results-all
%}

clear
load ~/demand_final/paper-results/results-full
load ~/demand_final/paper-results/results-all
sen_f1=nl_stderr(ds3,tn_f1);
clear ds3 t0
save results-all

clear
load ~/demand_final/paper-results/results-ignore
load ~/demand_final/paper-results/results-all
sen_ig1=nl_stderr(ds3,tn_ig1);
clear ds3 t0
save results-all

clear
load ~/demand_final/paper-results/results-emN
load ~/demand_final/paper-results/results-all
sen_em1=nl_stderr(ds3,tn_em1);
clear ds3 t0
save results-all




