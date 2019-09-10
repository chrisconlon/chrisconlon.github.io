%returns position of median element of vector (sort of)
%
function [ind]=findmedian(x)
	n=length(x);
	med = floor(1/3*n);
	[y,i] =sort(x);
	ind=i(med);	
end

