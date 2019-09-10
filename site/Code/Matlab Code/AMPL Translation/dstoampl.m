function dstoampl(ds3,myfname)

	[J N ] = size(ds3.avail);
	NFE = max(ds3.mktid);

	group = ds3.group;
	category = full(sparse([1:44]',ds3.group,ones(44,1)));

	matlabtoampl(full(category),'category');
	matlabtoampl(full(group),'group');
	matlabtoampl(full(ds3.avail'),'AV');
	matlabtoampl(full(ds3.qmat'),'Q');
	matlabtoampl(full(ds3.mktid),'mktid');


	fid = fopen('params.txt','wt');
	fprintf(fid,'param nmkt := %g;\n', N);
	fprintf(fid,'param nprod := %g;\n', J);
	fprintf(fid,'param NFE := %g;\n', NFE);
	fprintf(fid,'param ngroup := %g;\n\n\n',max(group));
	fclose(fid);

	!cat params.txt Q.txt AV.txt category.txt group.txt mktid.txt > vending_data.dat
	copyfile('vending_data.dat',myfname);

	%Deleting these files
	!rm -f params.txt Q.txt AV.txt category.txt group.txt mktid.txt;
end
