function [ds]=estep(dsarray,t0)

% generic stuff for all cases
for k=1:length(dsarray)
        N(k) = size(dsarray(k).avail,2);
end
[J K] = size(dsarray(2).Xmat);
NN = [0 cumsum([N])];

% extract parameters;
theta1=t0(1:J+K);
xi = t0(J+K+1:end);
if(length(xi) ~=NN(end))
	error('Wrong number of fixed effects!')
end

% one stockout
myds=dsarray(2);
xi2 = xi(NN(2)+1:NN(3));
thet=[theta1; repmat(xi2,[2 1]) ];
[J T] = size(myds.avail);

% find stockouts and tag them;
[focal,j,s] = find(myds.avail==2);
AV0 = myds.avail > 0;
AV1 = myds.avail ==1;
myds.avail = [AV0 AV1] >0;

QQ = myds.qmat;
myds.qmat=[myds.qmat myds.qmat];

P=choiceprob(myds,thet);
myds.qmat=QQ;
Parray = matsplit(P',T);

% after stockout
P0 = squeeze(Parray(:,:,1))';
P1 = squeeze(Parray(:,:,2))';
M=full(sum(myds.qmat,1))';

ystar=zeros(J+1,2*T);
A=[];
tic
for t=1:T,
    qvec = full(QQ(:,t));
    [y1,y2] = onestockout(M(t),focal(t),P0(:,t),P1(:,t),qvec);
    ystar(:,(2*t-1): 2*t)= [y1 y2]  ;
    diff(:,t) = QQ(:,t) - y1 -y2;
end
toc
dsarray(2).avail=reshape(([AV0;AV1] >0),size(AV0,1),size(AV0,2)*2);
dsarray(2).qmat = ystar;
d1 =  sum(dsarray(2).qmat(:)) - sum(QQ(:))

% two stockout case
myds = dsarray(3);
N= size(myds.avail,2);
[J T] = size(myds.avail);
xi3=xi(NN(3)+1:NN(4));
thet = [theta1; repmat(xi3,[4 1])];

M = full(sum(myds.qmat))';

clear focal
for t=1:T,
    qvec = myds.qmat(:,t);
    focal(t,:) = find(myds.avail(:,t)==2)';
    fsales(t,:) = full(myds.qmat(focal(t,:),t));
end

% Construct four availability regimes
AV0 = myds.avail > 0;
AV1 = myds.avail ==1;
focal1 = sparse([1:T],focal(:,1),ones(T,1),T, J);
focal2 = sparse([1:T],focal(:,2),ones(T,1),T, J);
AVA = AV0 .* ~focal1';
AVB = AV0 .* ~focal2';

[AV0; AVA; AVB; AV1];
%Compute Probabilities
myds=dsarray(3); 
myds.avail = [AV0 AVA AVB AV1] >0;
QQ = myds.qmat;
myds.qmat = repmat(myds.qmat,[1 4]);
P = choiceprob(myds,thet);
myds.qmat = QQ;
Parray = matsplit(P',N);


qhat = zeros(J+1,4,T);
for t=1:T,
%disp(['Iteration ' num2str(t)] )
    tic
    qvec = myds.qmat(:,t);
    Pslice = squeeze(Parray(t,:,:));
    fslice = focal(t,:);
    qhat(:,:,t)=twoSO(qvec,Pslice,fslice);
    timing(t) = toc;
    diff(t) = sum(sum(qhat(:,:,t))) - M(t);
end
d2=full(max(abs(diff(1:T))))
dsarray(3).qmat=reshape(qhat,[size(qhat,1), size(qhat,2)*size(qhat,3)]);
dsarray(3).avail=reshape([AV0; AVA; AVB; AV1]>0,size(AV0,1),size(AV0,2)*4);

% three stockout case
myds = dsarray(4);
N= size(myds.avail,2);
[J T] = size(myds.avail);
xi4 = xi(NN(4)+1:NN(5));
thet= [theta1; repmat(xi4,[8 1])];
M = full(sum(myds.qmat))';

clear focal
for t=1:T,
    qvec = myds.qmat(:,t);
    focal(t,:) = find(myds.avail(:,t)==2)';
end

% Construct eight availability regimes
AV0 = myds.avail > 0;
AV1 = myds.avail ==1;
focal1 = sparse([1:T],focal(:,1),ones(T,1),T, J);
focal2 = sparse([1:T],focal(:,2),ones(T,1),T, J);
focal3 = sparse([1:T],focal(:,3),ones(T,1),T, J);
AVA = AV0 .* ~focal1';
AVB = AV0 .* ~focal2';
AVC = AV0 .* ~focal3';
AVAB = AV0 .* ~focal1' .*~focal2';
AVAC = AV0 .* ~focal1' .*~focal3';
AVBC = AV0 .* ~focal2' .*~focal3';

%Compute Probabilities
myds=dsarray(4);
myds.avail=[AV0 AVA AVB AVC AVAB AVAC AVBC AV1]>0;
QQ = myds.qmat;
myds.qmat =repmat(QQ,[1 8]);
P = choiceprob(myds,thet);
myds.qmat =QQ;
Parray = matsplit(P',N);


qhat = zeros(J+1,8,T);
for t=1:T,
disp(['Iteration ' num2str(t)] )
    tic
    qvec = myds.qmat(:,t);
    Pslice = squeeze(Parray(t,:,:));
    fslice = focal(t,:);
    qhat(:,:,t)=threeSO(qvec,Pslice,fslice);
    timing3(t) = toc;
    diff3(t) = sum(sum(qhat(:,:,t))) - M(t);
end
d3=max(abs(diff3(:)))
dsarray(4).qmat=reshape(qhat,[size(qhat,1), size(qhat,2)*size(qhat,3)]);
a3part=[ones(7,1) diag(-1*ones(7,1))   ];
dsarray(4).avail=reshape(([AV0;AVA;AVB;AVC; AVAB; AVAC; AVBC; AV1] >0),size(AV0,1),size(AV0,2)*8) ;
ds=flattends(dsarray);

a1=[1:NN(2)]'; 
a2=repvec([NN(2)+1:NN(3)],2);
a3=repvec([NN(3)+1:NN(4)],4);
a4=repvec([NN(4)+1:NN(5)],8);
ds.mktid=[a1;a2;a3;a4];

end

