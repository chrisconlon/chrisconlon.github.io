function [dsarray]=adjust_mask(dsarray)

for(k=1:length(dsarray))
        N(k) = length(dsarray(k).avail) * 2^(k-1);
end

NN=[0 cumsum(N)];
mask=zeros(sum(N),length(dsarray));

for(k=1:length(dsarray)),
	mask([NN(k)+1:NN(k+1)],k)=1;
	dsarray(k).mask = mask(:,k) >0;
end


end

