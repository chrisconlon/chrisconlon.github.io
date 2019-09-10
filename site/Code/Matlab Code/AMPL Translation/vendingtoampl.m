function vendingtoampl(ds2,degree)

if nargin < 2,
	degree=5;
end
[nmkt, J] = size(ds2.avail);
[J,K] = size(ds2.Xmat);
[n w] = nwspgr('KPN',K,degree);
[ns k] = size(n);

ds2.v=n;
ds2.w=w;
ds2.ns = ns;

Q = ds2.qmat(:,1:J+1);
matlabtoampl(full(Q),'Q');
matlabtoampl(ds2.cat,'cat');
matlabtoampl(full(ds2.avail),'AV');
matlabtoampl(n,'v');
matlabtoampl(w,'w');
matlabtoampl(ds2.Xmat,'X2');


fid=fopen('header.txt','w');
fprintf(fid,'param ns := %d; \n', ns);
fprintf(fid,'param nmkt := %d; \n', nmkt);
fprintf(fid,'param nprod := %d; \n', J);
fprintf(fid,'param nk := %d; \n \n', K);
fclose(fid);

t0 = [-5*ones(44,1); ones(nmkt+K,1)];
save ds-temp.mat ds2 t0;

! cat header.txt Q.txt AV.txt X2.txt v.txt w.txt > vending.dat
! rm Q.txt cat.txt AV.txt v.txt X2.txt w.txt header.txt
end
