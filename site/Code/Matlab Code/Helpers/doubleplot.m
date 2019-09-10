function doublelot(data,n)
data1=data(1:n);
data2=data(n+1:end);
[f,xi]=ksdensity(data1,'npoints',length(unique(data1)),'width',.3);
[f2,xi2]=ksdensity(data2,xi,'npoints',length(unique(data1)),'width',.3);
plot(xi,[f; f2])
legend('No stockouts','Stockouts')
end