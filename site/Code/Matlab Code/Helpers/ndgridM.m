% Let M be a scalar from [0:M] such that we construct [0:M]^k gridpoints.
% where the sum of the gridpoints is less than M.
%
% This only works for k=2 and k=3

function [A] = ndgridM(M,k,offset)
	if(nargin < 3)
		offset=0;
	end
	if(k ~= 2 & k~= 3),
		error('K must be 2 or 3 only')
	end
	if(k==2 )
		[x1,x2]=ndgrid([0:M],[0:M]);
		ind = find(x1+x2 <= M-offset);
		A=[x1(ind) x2(ind)];
	end
	if(k==3)
		[x1,x2,x3]=ndgrid([0:M],[0:M],[0:M]);
		ind = find(x1+x2+x3<= M-offset);
		A=[x1(ind) x2(ind) x3(ind)];
	end
end


