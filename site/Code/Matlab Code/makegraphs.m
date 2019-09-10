load results-all

J=44; L=2; G=5;
xi_em=-(t_em(J+L+1:end)-mean(t_em(J+L+1:end)));
xi_full=-t_full(J+L+1:end);
xi_ig=-t_ig(J+L+1:end);

xin_em=tn_em(J+G+1:end);
xin_full=tn_f(J+G+1:end);
xin_ig=tn_ig(J+G+1:end);

nmkt_ignore=length(xi_ig);

%random coeffcients
[mean(xi_full) mean(xi_ig) mean(xi_em)]

% Mean without stockouts / with stockouts
std(xi_em)
[mean(xi_em(1:nmkt_ignore)) mean(xi_em(nmkt_ignore+1:end))]

% Mean without stockouts / with stockouts
std(xi_full)
[mean(xi_full(1:nmkt_ignore)) mean(xi_full(nmkt_ignore+1:end))]


%Nested Logits
[mean(xin_full) mean(xin_ig) mean(xin_em)]

% Mean without stockouts / with stockouts
std(xin_em)
-[mean(xin_em(1:nmkt_ignore)) mean(xin_em(nmkt_ignore+1:end))]

% Mean without stockouts / with stockouts
std(xin_full)
-[mean(xin_full(1:nmkt_ignore)) mean(xin_full(nmkt_ignore+1:end))]

%doubleplot(xi_em,nmkt_ignore)
%xlim([-10 10]);
%title('Fixed Effect distribution: RC-EM Specification')

%doubleplot(xi_full,nmkt_ignore)
%xlim([-10 10]);
%title('Fixed Effect distribution: RC-Full Specification')

%doubleplot(xin_em,nmkt_ignore)
%xlim([-10 10]);
%title('Fixed Effect distribution: NL-EM Specification')

doubleplot(xin_full,nmkt_ignore);
xlim([-10 10]);
%title('Fixed Effect distribution: NL-Full Specification')
