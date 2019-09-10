function [se] =nl_stderr(ds3,t)

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
    	Pg = IVL./repmat(denom,[1, G]);

	% construct within group shares
	Pkg = nantozero((exp(xi)*exp(dj)')./IV(:,ds3.group)) .* ds3.avail';
	Pk  = Pkg .* Pg(:,ds3.group);
	
	[mkt,pid,sales] = find(ds3.qmat(1:J,:)');
	subs=sub2ind(size(ds3.qmat(1:J,:)'),mkt,pid);
	catsubs=sub2ind(size(IV),mkt,ds3.group(pid));
	q0t = full(ds3.qmat(J+1,:))';
	
	% SCORE D_J part
	part1	=repmat(lam(ds3.group)',[N 1]).*(Pkg-Pk)-Pkg;
	score_dj=sparse(1:length(mkt),pid,ones(length(mkt),1)) + part1(mkt,:);
	score0_dj = -repmat(lam(ds3.group)',[N 1]).*Pk;

	% xi_t part
	xipart=repmat(lam',[N 1])-repmat(Pg*lam,[1 G]);
	score_xi= xipart(catsubs);
	score0_xi = -(Pg*lam);

	%lambda part
	lam_part=Pg.*IV;
	score_lam=sparse(1:length(mkt),ds3.group(pid),IV(catsubs)) - lam_part(mkt,:);
	score_lam1=IV(catsubs) - sum(lam_part(mkt,:),2);
	score0_lam = -lam_part;
	score0_lam1 = -sum(lam_part,2);

	score = full([score_dj score_lam1; score0_dj score0_lam1]);
	score_xi2 = [score_xi; score0_xi];

	[T Nparam]= size(score);
	A=dot2([sales;q0t],outerprodall(score),1,3);

	% Accumulate row-wise from long format for B
	FEvec = [mkt; [1:N]'];
	SparseMatrix = bsxfun(@eq,sparse(1:N).',FEvec.');
	B=SparseMatrix*(score.*repmat(score_xi2.*[sales;q0t],[1 Nparam]));
	D = q0t.*(score0_xi.^2) + accumarray(mkt,sales.*(score_xi.^2));

	se=sqrt(diag((A-B'*sparse(diag(1./D))*B )\speye(size(A))));

end
