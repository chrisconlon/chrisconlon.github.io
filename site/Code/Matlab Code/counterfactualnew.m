clear
% set Constants
M= 3360; whichxi= 5007;

% Read in product data
names=strtrim(cellstr(importdata('prodnames.txt')));
names=names(2:end);

PDATA = csvread('table1.csv',1,1);
pid = PDATA(:,8);
category = PDATA(:,1)+1;
exclude =  [9 24 28 32 36 37 38 92 97 998];
stockouts = [25  65  27  31  17  15  5  12  95  66];
av=~ismember(pid,exclude);
av1=~ismember(pid,exclude) & ~ismember(pid,stockouts);
prods=find(av);

%Sort order for products
[y,order]= sortrows([PDATA(av,1) -PDATA(av,9)]);
names2=names(av);
% Order of Included 35 Products
names2(order)

% excluded Products
names(~av)

% Simulated Stockout Products
names(ismember(pid,stockouts))


% load data and results
load ~/demand_final/paper-results/results-full
load ~/demand_final/paper-results/results-all

%load ~/demand_final/mach_fe/ds-full-mach
%load ~/demand_final/mach_fe/results-all

J=44; L=2; G=5;

% all xi_t
xi_ignore 	= t_ig(J+L+1:end);
nmkt=length(xi_ignore);

xi_full 	= t_full(J+L+nmkt+1:end);
xi_em 		= t_em(J+L+nmkt+1:end);

xi_nfull 	= tn_f(J+G+nmkt+1:end);
xi_nignore 	= tn_ig(J+G+1:end);
xi_nem 		= tn_em(J+G+nmkt+1:end);

[mean(xi_full) mean(xi_ignore) mean(xi_em)]
[mean(t_full(J+L+1:end)) mean(t_ig(J+L+1:end)) mean(t_em(J+L+1:end)) ]

%whichxi=findmedian(xi_full);

xi_full 	= quantile(xi_full,1/2);
xi_ignore 	= quantile(xi_ignore,1/2);
xi_em 		= quantile(xi_em,1/2);

xi_nfull 	= quantile(xi_nfull,1/2);
xi_nignore 	= quantile(xi_nignore,1/2);
xi_nem 		= quantile(xi_nem,1/2);

t_full 	= t_full(1:(J+L))
t_ig 	= t_ig(1:(J+L))
t_em 	= t_em(1:(J+L))

tn_f 	= tn_f(1:(J+G))
tn_ig 	= tn_ig(1:(J+G))
tn_em 	= tn_em(1:(J+G))

% Best substitutes
ds3.avail = sparse(repmat(double(av),[1 length(prods)]));
for j=1:length(prods),
	ds3.avail(prods(j),j)=0;
end
ds3.avail = [ds3.avail av av1];
n = size(ds3.avail,2);
ds3 = rmfield(ds3, 'mktid');
%ds3.mktid = ones(n,1);
ds3.q0t = 1*ones(n,1);
ds3.nmiss =1*ones(n,1);
ds3.qmat = sparse([double(ds3.avail); ones(1,n)]);

theta=[t_em(1:J+L); xi_em*ones(n,1)];
thetaN=[tn_em(1:J+G); xi_nem*ones(n,1)];

P = choiceprob(ds3,theta);
Pn=choiceprobN(ds3,thetaN);

P0 = repmat(P(1:J,end-1),[1 length(prods)]);

Pstar=nantozero(P(1:J,1:end-2)-P0);

[y,i]=max(Pstar);

% Best Substitutes 
[names(prods) names(i)]


format bank
% Weekly Sales
ds3.avail = sparse([av  av1]);
n = size(ds3.avail,2)
ds3.qmat = sparse([double(ds3.avail); ones(1,n)]);

P_em 	= choiceprob(ds3,[t_em; ones(n,1)*xi_em]);
P_full 	= choiceprob(ds3,[t_full; ones(n,1)*xi_full]);
P_ignore= choiceprob(ds3,[t_ig; ones(n,1)*xi_ignore]);

P_nignore= choiceprobN(ds3,[tn_ig; ones(n,1)*xi_nignore]);
P_nfull= choiceprobN(ds3,[tn_f; ones(n,1)*xi_nfull]);
P_nem= choiceprobN(ds3,[tn_em; ones(n,1)*xi_nem]);

P_before = [P_full(:,1) P_ignore(:,1) P_em(:,1) P_nfull(:,1) P_nignore(:,1) P_nem(:,1)];
P_after = [P_full(:,2) P_ignore(:,2) P_em(:,2)  P_nfull(:,2) P_nignore(:,2) P_nem(:,2)];


% Table 6 : Weekly Sales
weekly = P_before(av,:)*M;
weekly = [weekly(:,1:3) weekly(:,3)-weekly(:,1) weekly(:,4:6) weekly(:,6)-weekly(:,4)];
weekly(order,:)

% Table 7: Weekly Stockout Cost
group = category(av);
sales=(P_after(av,:)-P_before(av,:))*M;
sales(order,:)

% Table 8: Economic Impact
% Lost to Stockouts
[i,j,s] = find(sales .* (sales < 0));
part1=[accumarray([group(i) j], s); sum(sales.* (sales < 0))]

% Gained by Substitutes
[i2,j2,s2] = find(sales .* (sales > 0));
part2=[accumarray([group(i2) j2], s2); sum(sales.* (sales > 0))]

%Change in Overall Sales
[i3,j3,s3] = find(sales);
part3=[accumarray([group(i3) j3], s3); sum(sales)]

% Pct inside
inside=part2*(-100)./part1

% Pct Change in Profit
table1=importdata('table1.csv');
margin=table1.data(:,4)-table1.data(:,5);
profits=repmat(margin(av),[ 1 6]).*(P_after(av,:) - P_before(av,:))*M;
[i4,j4,s4] = find(profits);
part5=[accumarray([group(i4) j4], s4); sum(profits)]

% Regression Results
format short
s=regstats(weekly(:,4),table1.data(av,3),'linear',{'rsquare','beta','tstat'})
% beta / se / t-stat
r1=[s.tstat.beta s.tstat.se]'
s.rsquare'

s=regstats(weekly(:,8),table1.data(av,3),'linear',{'rsquare','beta','tstat'})
% beta / se / t-stat
r2=[s.tstat.beta s.tstat.se]';
s.rsquare'
[r1(:) r2(:)]



% Best substitutes by rank figure   
% percentagechnge
result=(Pstar(av,:)./P0(av,:))*100
resultsorted=sort(result,'descend');
resultsorted=resultsorted(1:34,:);
barh(median(resultsorted,2))
xlabel('Sales Change in Percent (for Substitute)','FontSize',18)
ylabel('Rank of Substitute','FontSize',18)
