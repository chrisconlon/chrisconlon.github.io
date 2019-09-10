% collapses markets to the avail-mach_id level (0 stockout case only)
function [dsarray]=collapse_markets(dsarray)
	[C,IA,IC ] = unique([dsarray(1).avail'],'rows');
	dsarray(1).avail = dsarray(1).availf(:,IA);
	dsarray(1).availf = dsarray(1).availf(:,IA);
	[pid mkt sales] = find(dsarray(1).qmat);
	dsarray(1).qmat=accumarray([pid IC(mkt)],sales);
	dsarray(1).q0t = accumarray(IC,dsarray(1).q0t);
	dsarray(1).mach_id = dsarray(1).mach_id(IA);

	for(k = 1:length(dsarray)),
		dsarray(k).mkt_id = dsarray(k).mach_id;
	end
end

