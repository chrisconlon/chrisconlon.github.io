function [f,vals] = statagroup(x)
[n k] = size(x);
[vals,i,j] = unique(x,'rows');
vals=1:length(vals);
vals=vals';
f=vals(j);
end

