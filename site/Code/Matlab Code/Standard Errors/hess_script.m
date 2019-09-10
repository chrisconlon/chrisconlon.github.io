clear
X=dlmread('hessian-ampl.txt','\t',1,0);
H=sparse(X(:,1),X(:,2),-X(:,3));
a=abs(diag(H))==0;
%a(46) = 1;
H1=H(~a,:);
H2=H1(:,~a);

%[L,p]=chol(H);
num2str(blockinverse(H2,45))
