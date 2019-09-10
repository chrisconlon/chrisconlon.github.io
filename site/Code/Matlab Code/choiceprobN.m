function [pjt,ll] = choiceprobN(ds3,t)
	[J N ] = size(ds3.avail);
	G = max(ds3.group);
	
	dj = t(1:J);
	lam = t(J+1:J+G);
	xi = t(J+G+1:end);
	if(~isfield(ds3,'mktid'))
		mktid=[1:N];	
	else
		mktid=ds3.mktid;
	end
	xi = xi(mktid);
	if(length(xi) ~= N),
		disp('Wrong dimensional parameter vector!');
	end
	
	% construct inclusive values
	edelta=exp(dj);
	[pid, id] = find(ds3.avail);
	IV=accumarray([id ds3.group(pid)],edelta(pid).*exp(xi(id)));

	% construct group shares
	IVL = IV.^repmat(lam',[size(IV,1),1]);
 	denom=1.0+sum(IVL,2);
    	nterm = (IVL ./ IV)./repmat(denom,[1, G]);

	% fix cases where whole category is unavailable
	nterm(isnan(nterm)) = 1;

	% construct within group shares
	indterm = (exp(xi)*exp(dj)') .* ds3.avail';
	
	pjt = indterm.*nterm(:,ds3.group);
	pjt=full([pjt 1-sum(pjt,2)]');

	ll=poslog(pjt(:))'*ds3.qmat(:);
	
end

