function [y1, y2] = onestockout(M,focal, P0, P1 ,qvec)
    % sparsity means we only compute for products with nonzero sales
    ind = find(qvec);
    f2 = find(ind==focal);
    
    qvec2 = qvec(ind);
    pvec=nantozero(P1(ind)./P0(ind));
    J = length(pvec);

    % get the sales and p for focal product
    omega = qvec2(f2);
    pfocal = P0(focal);
    
    % Setup the grid and evaluate the pdf
    x = [0:M-omega]';
    pr=nbinpdf(x,omega,pfocal);
    pr = pr./sum(pr);
    
    part1=repmat(x,[1 J]);
    part2=(M-omega - x)* pvec';
    integrand=nantozero(part1./(part1+part2));
    
    % correct the focal product
    integrand(1,f2) = 1;
    
    % Do the integration and correct for <= M
    y1=(pr'*integrand);
    y2 =1-y1;
    
    y1=sparse(ind,ones(J,1),y1,length(P0),1);
    y2=sparse(ind,ones(J,1),y2,length(P0),1);
    y1 = y1.*qvec;
    y2 = y2.*qvec;
    
    %fix the stocked out good --redundant
    y2(focal) =0;
    y1(focal) = qvec(focal);
end
