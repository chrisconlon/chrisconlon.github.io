function [ds] =flattends(dsarray)

nds = length(dsarray);
	avail=[];
	availf=[];
	qmat=[];
	nmiss=[];
	mach_id =[];
	for(k=1:nds),
		[J N] = size(dsarray(k).avail);
		avail=[avail dsarray(k).avail];
		availf=[availf dsarray(k).availf];
		qmat=[qmat dsarray(k).qmat];
		nmiss=[nmiss; dsarray(k).nmiss*ones(N,1)];
		mach_id=[mach_id; dsarray(k).mach_id];
	end
	ds.avail=double(avail >0);
	ds.availf=double(availf >0);
	ds.Xmat=dsarray(1).Xmat;
	ds.qmat=qmat;
	ds.v=dsarray(1).v;
	ds.w=dsarray(1).w;
 	ds.q0t=qmat(end,:)';
	ds.ns=dsarray(1).ns;
	ds.Z=dsarray(1).Z;
	ds.nmiss=nmiss;
	ds.group = dsarray(1).group;
	ds.mach_id = mach_id;
end

