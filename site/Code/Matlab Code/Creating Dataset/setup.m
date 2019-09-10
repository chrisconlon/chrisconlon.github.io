clear
data = csvread('salesdata.csv',1,0);
id = data(:,1);
Q=sparse(data(:,2:end));
data2= csvread('dexavail.csv',1,0);
id2 = data2(:,1);
av = sparse(data2(:,2:end));
datafull= csvread('fullavail.csv',1,0);
id3 = datafull(:,1);
avfull = sparse(datafull(:,3:end)>0);

if isequal(id,id2)
    disp('IDs are sorted')
end

if isequal(id,id3)
    disp('Fulll avail IDs are sorted')
end

if isequal(sum(Q(:)), sum(Q(av > 0)))
    disp('No sales when unavailable')    
end

bigq =full(sum(Q,2));
data3= csvread('dexreaddata.csv',1,0);
rate=bigq./max(data3(:,3),1);
mach_id = statagroup(data3(:,2));

dataX = csvread('characteristics.csv',1,0);
%X = dataX(:,end-3:end);
X = zscore(dataX(:,[2 6 ]));

group = csvread('class.csv',1,2)+1;

% minimum period of 90 minutes
M=max(data3(:,3),1.5)*20;
disp(['Number of too small markets: ' num2str(sum(bigq > M))])

% outside good share
q0 = M-bigq;

% adjust some borderline cases
idx1=find(q0 < 0 & q0 > -10) ;
M(idx1) = M(idx1) - q0(idx1);
disp(['Number of too small markets (adjusted): ' num2str(sum(bigq > M))])
q0 = M-bigq;
[q0(q0 < 0) M(q0<0)]                                          

idx=find(q0>=0);
av=av(idx,:);
avfull=avfull(idx,:);
Q=Q(idx,:);
q0=q0(idx);
M = M(idx);

s0 = q0./M;
disp(['Mean outside good share: ' num2str(mean(s0(q0>=0)))])
disp(['quantile outside good share: ' num2str(quantile(s0(q0>=0),[.1 .25 .5 .75 .9]))])

[v,w] = nwspgr('KPN',size(X,2),8);
N=length(av);
ns = length(w);
[J,K] = size(X);
Z=permute(repmat(X,[1 1 ns]) .* permute(repmat(v,[1 1 J]),[3 2 1]),[2 1 3]);

ds3.avail = av';
ds3.availf = avfull';
ds3.Xmat = X;
ds3.group = group;
ds3.qmat = [Q q0]';
ds3.v = v;
ds3.w = w;
ds3.q0t = q0;
ds3.ns =ns;
ds3.Z = Z;
ds3.mach_id = mach_id;

t0= ones(N+J+K,1);

% base case
save em-new0.mat ds3 t0;

clear
load em-new0;
M=sum(ds3.qmat);
bigQ=sum(ds3.qmat(1:end-1,:));
mis = find(max(ds3.avail )==2);
obs = find(max(ds3.avail )==1);

[i2,j2,focalsales]=find((ds3.avail == 2) .* ds3.qmat(1:end-1,:));
M2=M(j2);
nmiss=histc(j2,1:length(M));
histc(nmiss,[0:max(nmiss)])
full([[0:max(nmiss)]' histc(nmiss,[0:max(nmiss)])  accumarray(nmiss+1,M') accumarray(nmiss+1,bigQ') accumarray(nmiss+1,accumarray(j2,focalsales,size(nmiss)))])
cutoff = 43581;

% find missing data (unobserved stockouts)
[i2,j2,focalsales]=find((ds3.avail == 2) .* ds3.qmat(1:end-1,:));
focalsales = focalsales(j2 < cutoff);
i2 = i2(j2 < cutoff);
j2 = j2(j2 < cutoff);
M2=M(j2);
nmiss=histc(j2,1:length(M));
missing_N=histc(nmiss,[0:max(nmiss)])
NN=missing_N(1:4)'*(2.^[0:3])';
full([[0:max(nmiss)]' histc(nmiss,[0:max(nmiss)])  accumarray(nmiss+1,M') accumarray(nmiss+1,bigQ') accumarray(nmiss+1,accumarray(j2,focalsales,size(nmiss)))])

% split data by # of stockouts
dsarray=splitds(ds3,nmiss);

% restrict to three stockouts or less
dsarray=dsarray(1:4);

% cut out some very large markets with many stockouts from the 3 SO case
exclusions=[21    22    23    28    37    40    49    70    96   113   138   139   149   160   166 168   170 179   186];
keep = ~sparse(exclusions,ones(size(exclusions)),ones(size(exclusions)),size(dsarray(4).avail,2),1);
dsarray(4).avail = dsarray(4).avail(:,keep);
dsarray(4).availf = dsarray(4).availf(:,keep);
dsarray(4).qmat = dsarray(4).qmat(:,keep);
dsarray(4).mach_id = dsarray(4).mach_id(keep);

% correctly identify which observations are from which stockout scenario
dsarray=adjust_mask(dsarray);
%dsarray=collapse_markets(dsarray);

dsarray=dsarray(2:end);
save dsarray2-allmissing.mat dsarray;

%ds3=flattends(dsarray);
%t0=ones(sum(size(ds3.Xmat)) + length(ds3.avail),1);
save ds-em-none.mat ds3 t0;


%ds3=dsarray(1);
%ds3.mktid = ds3.mach_id;
%NFE = max(ds3.mktid);
%t0=ones(sum(size(ds3.Xmat)) + NFE,1);
%save ds-ignore-none.mat ds3 t0;

ds3=flattends(dsarray);
%ds3.mktid = ds3.mach_id;
NFE = length(ds3.mach_id);
t0=ones(sum(size(ds3.Xmat)) + NFE,1);
ds3.avail=ds3.availf > 0;
save ds-full-allmissing.mat ds3 t0;



