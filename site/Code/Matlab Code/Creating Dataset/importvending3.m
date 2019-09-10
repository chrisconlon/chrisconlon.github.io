function importvending3(raw,cut)
pid2 = raw.pid; 
qvec = raw.qvec; 
Mvec = raw.Mvec;
id=raw.id;

[J,K] =size(raw.X2);

av1=raw.avail(sub2ind(size(raw.avail),raw.availid,raw.pid));
av2=raw.avail(sub2ind(size(raw.avail),raw.availid2,raw.pid));
stockout=av2>av1;

totstockouts = accumarray(raw.id,raw.qvec.*stockout);;

switch lower(cut)
case{'em'}
    % end of period availability
    av = raw.availid;
    savename = 'dataset-em.mat';
    disp(['Creating em  data as ' savename])
    % adjust stocked out products to be before stockout event
case{'full'}
    % use the full availability choice set (from last restock) always
    av = raw.availidfull;
    savename = 'dataset-full.mat';
    disp(['Creating assume full availability data as ' savename])
case{'ignore'}
    % discard periods with missing data (ignore)
    index=(raw.availid == raw.availid2);
    av=raw.availid(index); pid2=raw.pid(index); qvec = raw.qvec(index);  Mvec=Mvec(index);
    id=statagroup(id(index));
    savename = 'dataset-ignore.mat';
    disp(['Creating ignore periods with missing data as ' savename])
end        
disp(['Missing availids (should be zero): ' num2str(sum(isnan(av)))])

N=max(id);
minM=accumarray(id,ceil(Mvec),[],@min);
qmat = sparse(accumarray([id,pid2],qvec,[N J]));
bigq = accumarray(id,qvec,[N 1]); 

avs=accumarray(id,av,[],@max);
disp(['Markets where M is too small (should be blank): ' num2str(sum(bigq > minM))])
disp(['Total number of consumers (including outside good) ' num2str(sum(minM))])

ds2.q0t = minM-bigq;
ds2.qmat = [qmat ds2.q0t];
ds2.Xmat = raw.X2;
ds2.avail = raw.avail(avs,:);
ds2.cat = statagroup(raw.class);

save(savename,'ds2');
end
