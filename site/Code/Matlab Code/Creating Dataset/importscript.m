clear
c = importdata('formatlab.csv');
raw.availid2=c.data(:,5);
raw.availid = c.data(:,6);
raw.availid2(isnan(raw.availid2)) = raw.availid(isnan(raw.availid2));
raw.pid=statagroup(c.data(:,1));
raw.qvec = c.data(:,3);
raw.Mvec = c.data(:,4);
raw.availidfull=c.data(:,7);
raw.id = c.data(:,end);
J=max(raw.pid);
clear c;

a = importdata('availsets.csv');
availold = sparse(a.data(:,2:end));
[ai,aj,as]=find(availold);
raw.avail=sparse((a.data(ai,1)),aj,as);
[N Ja] = size(raw.avail);
if(Ja ~= J),
    disp(['Mismatch between sales product number ' num2str(J) ' and availability ' num2str(Ja)])
end
clear a;

b = importdata('product-data.csv');
raw.class = b.data(:,1);
if(max(statagroup(b.data(:,2))) ~=J)
    disp(['Mismatch between number of products in product data ' num2str(max(statagroup(b.data(:,2)))) ' and sales data ' 
num2str(J)]) 
end
X2 = mynorm2(b.data(:,end-6:end));
[coeff,score,latent]=princomp(X2(:,[1:3 5:7]));
cumsum(latent.^2./(sum(latent.^2)))
raw.X2=score(:,1:3);
[raw.id raw.avs] = statagroup(raw.id);

save rawnew.mat raw;
tic
clear
load rawnew
importvending3(raw,'full');
importvending3(raw,'ignore');
%importvending(raw,'em');
t=toc
