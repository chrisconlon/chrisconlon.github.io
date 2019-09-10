function [qhat]=sloweintegral(alphas,choiceprobs,pweights,qvec)
	[N K] = size(alphas);
	[J K2]= size(choiceprobs);
	N2 = length(pweights);

	if(N ~= N2 | K ~= K2)
		error('Input must be of correct length')
	end

	eprobs = zeros(J,K);

	% only care about probs for nonzero sales
	ind = find(qvec);
	for i=1:length(ind),
		j=ind(i);
                pslice=choiceprobs(j,:);
		b=alphas*pslice';
%		eprobs(j,:)=pweights'*nantozero((alphas.*repmat(pslice,[N 1])) ./b(:,ones(K,1)));
		for k= 1:K,
			eprobs(j,k)=(pweights'* nantozero(alphas(:,k)./b)) *choiceprobs(j,k);
		end
	end
	
	% multiply by q_jt
	qhat=full(repmat(qvec,[1 K]).* nantozero(eprobs));
end

