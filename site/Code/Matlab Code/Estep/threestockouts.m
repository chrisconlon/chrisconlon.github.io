function [alpha,beta,p] = threestockouts(qvec,P,fslice)

fsales=full(qvec(fslice));
Q0=full(sum(qvec)-sum(fsales));

P0 = P(:,1);
PA = P(:,2);
PAB = P(:,3);

% Presume A before B
% [q_0^0, q_0^a, q_0^ab, q_b^0, q_c^0, q_c^A]
A1=ndgridM(Q0,3,0);
A2=ndgridM(fsales(3),2,1);
[x1,x2,x3]=ndgrid([1:size(A1,1)],[0:fsales(2)-1],[1:size(A2,1)]);
A=[A1(x1(:),:) x2(:) A2(x3(:),:)];
[n k]=size(A);
clear x1 x2 x3 A1 A2

% construct third part
p_cab=PAB(fslice(3));
p3=nbinpdf(A(:,3),fsales(3)-A(:,5) -A(:,6),p_cab);

% construct second part
p_a = PA([fslice(2) fslice(3)]);
p_a = [p_a; 1-sum(p_a)];
q_a = [fsales(2) - A(:,4), A(:,6) A(:,2)];
p2  = negmnpdf(q_a,p_a,1);
clear p_a q_a

% construct first part
p_0 = P0([fslice(1) fslice(2) fslice(3)]);
p_0 = [p_0; 1-sum(p_0)];
q_0 = [repmat(fsales(1),[n 1]) A(:,4) A(:,5) A(:,1)];
p1  =  negmnpdf(q_0,p_0,1);
clear p_0 q_0

% multiply probabilities together
p= p1.*p2.*p3;

% for non-stockout products
alpha=A(:,[1:3])./Q0;
alpha(:,4) = 1-sum(alpha,2);

% for stocked out products
beta(:,1)=A(:,4)./fsales(2);
beta(:,2)=A(:,5)./fsales(3);
beta(:,3)=A(:,6)./fsales(3);

end
