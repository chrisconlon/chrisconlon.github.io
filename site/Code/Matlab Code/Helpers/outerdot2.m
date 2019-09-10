% This takes the outer product of each vector in the N-D array X and then
% the inner product of that with Y.

%
% Suppose that X is N x M x S and Y is S x 1
% The resulting product Z will be N x M x M
% This works by blocking on 

function B=outerdot2(X,Y)
[N M K] = size(X);
B=zeros(M,M,N);
for i=1:N,
     XX = squeeze(X(i,:,:));
     B(:,:,i) = dot2(outer(XX,XX,1,1),Y,3,1);
end
B=permute(B,[3 1 2]);
end
