function [qhat2] = twoSO(qvec,Pslice,fslice)

    Phat=Pslice(:,[1 2 4]);
    Phat2=Pslice(:,[1 3 4]);
    fvends = full(qvec(fslice));

    % construct the grid points A > B and B > A
    [alpha1,beta1,p1]=twostockouts(qvec,Phat,fslice);
    [alpha2,beta2,p2]=twostockouts(qvec,Phat2,fslice([2 1]));

    n1=length(p1);
    n2=length(p2);
    ngrid = n1+n2;

    % Normalize to condition on X <= M
    p=[p1;p2];
    p=p./sum(p);

   
    % put the grid points together
    alphas=[[alpha1(:,1); alpha2(:,1)] [alpha1(:,2); zeros(length(p2),1)] [alpha2(:,2); zeros(length(p1),1)] [alpha1(:,3); alpha2(:,3)]];
    pfocal1 = zeros(1,4);
    pfocal2 = zeros(1,4);

    % adjustments for focal products
    pfocal1=p'*[[ones(length(p1),1) zeros(length(p1),3)]; [beta2 zeros(length(p2),1) 1-beta2 zeros(length(p2),1)]] * full(qvec(fslice(1)));
    pfocal2=p'*[[beta1 1-beta1 zeros(length(p1),2)]; [ones(length(p2),1) zeros(length(p2),3)]] * full(qvec(fslice(2)));

    %do the integral
    qhat2 = sloweintegral(alphas,Pslice,p,qvec);
    qhat2(fslice(1),:) = pfocal1;
    qhat2(fslice(2),:) = pfocal2;

end

