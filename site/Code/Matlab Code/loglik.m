function [ll]=loglik(ds3,t0)
	% define sizes
	[J N] = size(ds3.avail);
	[J2 K] = size(ds3.Xmat);
	[NS L] =size(ds3.v);

	%extract parameters
	delta = t0(1:J);
	sigma = t0(J+1:J+L);
	xi    = t0(J+L+1:end);
	xi2   = xi(ds3.mktid);

	muij = squeeze(dot2(ds3.Z,sigma));
	utils = exp(repmat(delta,[1 NS]) + muij);
	denom = 1./(multiprod(ds3.avail',utils) + repmat(exp(xi2),[1 NS]));

	for t=1:N,
		avslice = ds3.avail(:,t);
		dslice = denom(t,:);
            	p(:,t)=((avslice*dslice).*utils)*ds3.w;
	end
	lp=poslog([p; 1-sum(p)]);
	ll=lp(:)'*ds3.qmat(:);

end

