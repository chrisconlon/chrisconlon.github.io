function y=nantozero(x)
	y=x;
	y(isnan(x))= 0;
end

