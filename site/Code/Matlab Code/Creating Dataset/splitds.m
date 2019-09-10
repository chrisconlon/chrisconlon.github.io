function [dsarray] = splitds(ds,index)
    ind = unique(index)
    
    for i=1:length(ind),
        dsarray(i) = dsslice(ds,index==ind(i));
    end
    for i=1:length(ind),
        dsarray(i).nmiss = ind(i);
    end
    return
    
function [dss] = dsslice(ds,mymask)
    dss.avail  = ds.avail(:,mymask);
    dss.availf = ds.availf(:,mymask);
    dss.Xmat = ds.Xmat;
    dss.group = ds.group;
    dss.qmat = ds.qmat(:,mymask);
    dss.v = ds.v;
    dss.w = ds.w;
    dss.q0t = ds.q0t(mymask);
    dss.ns = ds.ns;
    dss.Z = ds.Z;
    dss.mask = mymask;
    dss.mach_id = ds.mach_id(mymask);
    return;
end
end
