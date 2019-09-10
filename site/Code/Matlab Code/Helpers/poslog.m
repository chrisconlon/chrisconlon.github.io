function poslog(x)
	ind = find(x);
	x(ind) = log(x(ind));
	y = x;
end

