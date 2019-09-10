function [qhat3]=threeSO(qvec,Pslice,fslice)

    fsales = full(qvec(fslice));

    % extract chocie probabilities
    PhatABC=Pslice(:,[1 2 5 8]);
    PhatACB=Pslice(:,[1 2 6 8]);
    PhatBAC=Pslice(:,[1 3 5 8]);
    PhatBCA=Pslice(:,[1 3 7 8]);
    PhatCAB=Pslice(:,[1 4 6 8]);
    PhatCBA=Pslice(:,[1 4 7 8]);

    % generate grid points
    [alphaABC,betaABC,pABC] = threestockouts(qvec,PhatABC,fslice([1 2 3]));
    [alphaACB,betaACB,pACB] = threestockouts(qvec,PhatACB,fslice([1 3 2]));
    [alphaBAC,betaBAC,pBAC] = threestockouts(qvec,PhatBAC,fslice([2 1 3]));
    [alphaBCA,betaBCA,pBCA] = threestockouts(qvec,PhatBCA,fslice([2 3 1]));
    [alphaCAB,betaCAB,pCAB] = threestockouts(qvec,PhatCAB,fslice([3 1 2]));
    [alphaCBA,betaCBA,pCBA] = threestockouts(qvec,PhatCBA,fslice([3 2 1]));

    n1=length(pABC);
    n2=length(pACB);
    n3=length(pBAC);
    n4=length(pBCA);
    n5=length(pCAB);
    n6=length(pCBA);

    p=[pABC; pACB; pBAC; pBCA; pCAB; pCBA];
    p=p./sum(p);
    ngrid =length(p);

    clear pABC pACB pBAC pBCA pCAB pCBA;

    % adjust focal products one at a time
    pfocal1 = zeros(1,8);
    pfocal2 = zeros(1,8);
    pfocal3 = zeros(1,8);

    pfocal1(1) =[ones(n1+n2,1); betaBAC(:,1); betaBCA(:,2); betaCAB(:,1); betaCBA(:,2)]'*p;
    pfocal1(3) =[zeros(n1+n2,1); 1-betaBAC(:,1); betaBCA(:,3); zeros(n5+n6,1)]'*p;
    pfocal1(4) =[zeros(n1+n2+n3+n4,1); 1-betaCAB(:,1); betaCBA(:,3)]'*p;
    pfocal1(7) =[zeros(n1+n2+n3,1); 1-sum(betaBCA(:,[2 3]),2); zeros(n5,1); 1-sum(betaCBA(:,[2 3]),2) ]'*p;

    pfocal2(1) =[betaABC(:,1); betaACB(:,2); ones(n3+n4,1); betaCAB(:,2); betaCBA(:,1)]'*p;
    pfocal2(2) =[1-betaABC(:,1); betaACB(:,3); zeros(n3+n4+n5+n6,1)]'*p;
    pfocal2(4) =[zeros(n1+n2+n3+n4,1); betaCAB(:,3); 1-betaCBA(:,1)]'*p;
    pfocal2(6) =[zeros(n1,1); 1-sum(betaACB(:,[2 3]),2); zeros(n3+n4,1); 1-sum(betaCAB(:,[2 3]),2); zeros(n6,1)]'*p;
    
    pfocal3(1) = [betaABC(:,2); betaACB(:,1); betaBAC(:,2); betaBCA(:,1); ones(n5+n6,1)]'*p;
    pfocal3(2) = [betaABC(:,3); 1-betaACB(:,1); zeros(n3+n4+n5+n6,1)]'*p;
    pfocal3(3) = [zeros(n1+n2,1); betaBAC(:,3); 1-betaBCA(:,1); zeros(n5+n6,1)]'*p;
    pfocal3(5) = [1-sum(betaABC(:,[2 3]),2); zeros(n2,1); 1-sum(betaBAC(:,[2 3]),2); zeros(n4+n5+n6,1)]'*p;

    clear betaABC betaACB betaBAC betaBCA betaCAB betaCBA;
 
    % match grid points to availability regimes
    alphas=zeros(ngrid,8);
    alphas(:,1)=[alphaABC(:,1); alphaACB(:,1); alphaBAC(:,1); alphaBCA(:,1); alphaCAB(:,1); alphaCBA(:,1)];
    alphas(:,2)=[alphaABC(:,2); alphaACB(:,2); zeros(n3+n4+n5+n6,1)];
    alphas(:,3)=[zeros(n1+n2,1); alphaBAC(:,2); alphaBCA(:,2); zeros(n5+n6,1)];
    alphas(:,4)=[zeros(n1+n2+n3+n4,1); alphaCAB(:,2); alphaCBA(:,2)];
    alphas(:,5)=[alphaABC(:,3); zeros(n2,1); alphaBAC(:,3); zeros(n4+n5+n6,1)];
    alphas(:,6)=[zeros(n1,1); alphaACB(:,3); zeros(n3+n4,1); alphaCAB(:,3); zeros(n6,1)];
    alphas(:,7)=[zeros(n1+n2+n3,1); alphaBCA(:,3); zeros(n5,1); alphaCBA(:,3)];
    alphas(:,8)=[alphaABC(:,4); alphaACB(:,4); alphaBAC(:,4); alphaBCA(:,4); alphaCAB(:,4); alphaCBA(:,4)];
    
    clear alphaABC alphaACB alphaBAC alphaBCA alphaCAB alphaCBA pABC pACB pBAC pBCA pCAB pCBA;
    
    qhat3 = sloweintegral(alphas,Pslice,p,qvec);
    
    qhat3(fslice(1),:) = pfocal1 * fsales(1);
    qhat3(fslice(2),:) = pfocal2 * fsales(2);
    qhat3(fslice(3),:) = pfocal3 * fsales(3);

end

