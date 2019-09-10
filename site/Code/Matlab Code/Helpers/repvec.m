function [y] = repvec(x,k)
	a = repmat(x,[k 1]);
	y=a(:);
end

