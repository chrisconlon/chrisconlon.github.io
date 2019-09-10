clear
load dsarray2-mach

% load the parameter starting values and expand for the extra EM observations
load ~/demand_final/mach_fe/outfull.mat

[J K] = size(dsarray(1).Xmat);
theta1=t1(1:J+K);

xi = t1(J+K+1:end);
a=[dsarray(1).mkt_id; dsarray(2).mkt_id; dsarray(3).mkt_id; dsarray(4).mkt_id];
xi2=xi(a);

ds3=estep2(dsarray,[theta1;xi2]);
ds3.mktid = ds3.mktid(ds3.mach_id);
t0=t1;

save estepped2.mat ds3 t0;

