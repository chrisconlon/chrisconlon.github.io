function y = matsplit(X,n)
    [N M] = size(X);
    if(mod(N,n)==0),
      K=N./n;
      y=zeros(n, M, K);
       for i=1:K,
           y(:,:,i)=X((i-1)*n+1:i*n,:);
       end
    end
end

