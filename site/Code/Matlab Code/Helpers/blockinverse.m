% this returns only the sqrt of the diagonal matrix of the top left block of the matrix
function [y]=blockinverse(X,m)
	A=X(1:m,1:m);
	B=X(1:m,m+1:end);
	C=X(m+1:end,1:m);
	D=X(m+1:end,m+1:end);
	n=length(D);

	%Dinv = D\speye(n);
	Dinv = spdiags(1./diag(D),0,n,n);
	Ainv=(A-B*Dinv*C)\speye(m);
	y=full(sqrt(diag(Ainv)));
end

