function matlabtoampl(X,vname,fname)
if nargin < 3,
    fname = strcat(vname,'.txt');
end
[n k] = size(X);
disp([vname num2str(k)])
Z=[[1:n]' X];
a= [1:k];

fid=fopen(fname,'w');
fprintf(fid,'param %s ', vname);
if(k > 1),
	fprintf(fid,': ');
	fprintf(fid,' %d ', a);
end
	fprintf(fid,' := ');
for i = 1:n,
    fprintf(fid,'\n');
    fprintf(fid,'%f\t', Z(i,:));
end
fprintf(fid,' ; \n \n');
fclose(fid);
end


