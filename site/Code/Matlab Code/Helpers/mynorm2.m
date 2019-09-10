function  A =mynorm2 (X)
	X=bsxfun(@minus,X,mean(X));
	A=bsxfun(@times,X,1./max(X));
end

