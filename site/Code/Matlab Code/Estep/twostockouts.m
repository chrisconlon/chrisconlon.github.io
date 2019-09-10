function [alpha,beta,p] = twostockouts(qvec,Phat,fslice)

fsales=full(qvec(fslice));
Q0=full(sum(qvec)-sum(fsales));
P0 = Phat(:,1);
PA = Phat(:,2);

% Presume A before B: 
% [q_0^0, q_0^a, q_b^0]


% Presume A before B
% [q_0^0, q_0^a, q_b^0]
[X1,X2,X3] = ndgrid([0:Q0],[0:Q0],[0:fsales(2)-1]);
ind=(X1+X2<Q0);
A=[X1(ind) X2(ind) X3(ind)];
beta = A(:,3)/fsales(2);
alpha=A(:,[1:2])./Q0;
alpha(:,3) = 1-sum(alpha,2);

[n k]=size(A);

% construct second part
p_ba=PA(fslice(2));
p2=nbinpdf(A(:,2),fsales(2)-A(:,end),p_ba);

% construct first part
p_0 = P0([fslice(1) fslice(2)]);
p_0 = [p_0; 1-sum(p_0)];
q_0 = [repmat(fsales(1),[n 1]) A(:,end) A(:,1)] ;
p1  = negmnpdf(q_0,p_0,1);

% multiply probabilities together
p= p1.*p2;




%{
% ndgridM imposes constraint on M1 + M2<=M
A1=ndgridM(Q0,2);
[x1,x2]=ndgrid([1:length(A1)],[0:fsales(2)-1]);
A=[A1(x1(:),:) x2(:)];

alpha=A(:,[1:2])./Q0;
alpha(:,3) = 1-sum(alpha,2);
beta = A(:,3)/fsales(2);

[n k]=size(A);

% construct second part
p_ba=PA(fslice(2));
p2=nbinpdf(A(:,2),fsales(2)-A(:,end),p_ba);

if(sum(isnan(p2))>0)
	error('Problem with negative binomial density');
end

% construct first part
p_0 = P0([fslice(1) fslice(2)]);
p_0 = [p_0; 1-sum(p_0)];
q_0 = [repmat(fsales(1),[n 1]) A(:,end) A(:,1)] ;
p1  = negmnpdf(q_0,p_0,1);

% multiply probabilities together
p= p1.*p2;
%}
end
