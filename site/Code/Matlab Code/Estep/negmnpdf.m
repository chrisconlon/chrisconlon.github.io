% mnpdf Multinomial probability density function (pdf).
% Y = mnpdf(X,PROB,FOCAL) returns the pdf for the multinomial distribution with
% probabilities PROB, evaluated at each row of X. X and PROB are M-by-K
% matrices or 1-by-K vectors, where K is the number of multinomial bins or
% categories.  Each row of PROB must sum to one, and the sample sizes for
% each observation (rows of X) are given by the row sums SUM(X,2).  Y is
% a M-by-1 vector, and mnpdf computes each row of Y using the corresponding
% rows of the inputs, or replicates them if needed. And FOCAL indicates the
% position of the focal product.
% 
% Example:
%  Generate a random vector with sample size 20 and probabilities P and
%  compute the multinomial pdf of X with probabilities P
%  P=[0.3,0.7];
%  X=mnrnd(20,P);
%  Y=mnpdf(X,P,1);

function Y=negmnpdf(X,prob,focal)
    if(length(X) > 1e9),
        [A,I,J] = unique(X,'rows');
        Y=mnpdf(A,prob).*(A(:,focal)./sum(A,2));  
        Y=Y(J,:);
    else
        Y=mnpdf(X,prob).*(X(:,focal)./sum(X,2));  
    end
    
end

