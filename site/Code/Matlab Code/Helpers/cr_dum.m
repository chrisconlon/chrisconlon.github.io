function f=cr_dum(myvar)
    [B,I,J] = unique(statagroup(myvar));
    f=full(sparse(1:length(J),J,ones(length(J),1)));
end

