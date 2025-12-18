function [errorparams logmlhres dellogmlh bic]=mlhTPABIsimfitEintParerrorC(initparams,effb,resparamfix,onet3rkk,stateid,fitfu,nintm,drvintv,drvpnt)
%
% [errorparams logmlhres dellogmlh bic] = mlhTPABIsimfitEintParerrorC(initparams,effb,resparamfix,onet3rkk,stateid,fitfu,nintm,drvintv,drvpnt)
%
% Calculate errors (standard deviation) of the parameter determined by mlhTPABIsimfitEintParC.m.  
% Outputs are parameter errors, log likelihood value, log likelihood value
% difference from an instantaneous transition model, and Bayesian
% Information Criterion (BIC) calculated using the optimized parameters.

    function logmlh=mlhTPfitsub(params)
        nparI=(length(params)-3)/3;  % number of parallel pathways
        fretEF=params(1);
        fretEU=params(2);
        ktobrt=params(3);
        frnb0=params(4);        % for a photon count rate of 100 ms-1
        trplottime=params(4+(1:nparI));
        Eint=params(4+(nparI+1:2*nparI));  % Eint
        pI=params(4+(2*nparI+1:3*nparI-1));  % relative fractions of pathways (proportional to the TP times). Note that these are not the same as the relative populations of Is.
        pI=[pI; 1-sum(pI)];

        resparamone=resparamfix;resparamone(3)=resparamone(3)/1000;
        resparamone(1:2)=[fretEF fretEU];

        logmlhfu(1:2)=0;
        stateidvit=1;
        logmlhfu(stateidvit)=mlhrateeigtrpathrevcalvitABIparC([resparamone; effb; ktobrt; frnb0],[trplottime trplottime],burstbint3rall1,cumindex1,indexone1,cntrate1,stateidvit,nintm,Eint,pI);
        stateidvit=2;
        logmlhfu(stateidvit)=mlhrateeigtrpathrevcalvitABIparC([resparamone; effb; ktobrt; frnb0],[trplottime trplottime],burstbint3rall2,cumindex2,indexone2,cntrate2,stateidvit,nintm,Eint,pI);

        if fitfu==1 || fitfu==2, logmlh=-logmlhfu(fitfu); end
        if fitfu==3, logmlh=-sum(logmlhfu); end
    end

stidnum=[0 0];
for k=1:length(onet3rkk)
    stateidone=stateid{k}(1,3);
    stidnum(stateidone)=stidnum(stateidone)+1;
end
burstbint3rall1=[];burstbint3rall2=[];
burstbint3rone1=[];burstbint3rone2=[];
burstphotons1=[];burstphotons2=[];
stidcount=[1 1];
for k=1:length(onet3rkk)
    stateidone=stateid{k}(1,3);
    onet3r=onet3rkk{k};
    if stateidone == 1
        burstbint3rone1=[burstbint3rone1; onet3r];
        burstphotons1=[burstphotons1; length(onet3r(:,1))];
        if mod(stidcount(stateidone),3000) == 0 || stidcount(stateidone) == stidnum(stateidone)
            burstbint3rall1=[burstbint3rall1; burstbint3rone1];
            burstbint3rone1=[];
        end
        stidcount(stateidone)=stidcount(stateidone)+1;
    else
        burstbint3rone2=[burstbint3rone2; onet3r];
        burstphotons2=[burstphotons2; length(onet3r(:,1))];
        if mod(stidcount(stateidone),3000) == 0 || stidcount(stateidone) == stidnum(stateidone)
            burstbint3rall2=[burstbint3rall2; burstbint3rone2];
            burstbint3rone2=[];
        end
        stidcount(stateidone)=stidcount(stateidone)+1;
    end
end

indexone1=1:length(burstphotons1);
cumindex1=[0; cumsum(burstphotons1)];
indexone2=1:length(burstphotons2);
cumindex2=[0; cumsum(burstphotons2)];

cntrate1=zeros(length(indexone1),1);
for k=1:length(indexone1)
    oneburst=burstbint3rall1(cumindex1(indexone1(k))+1:cumindex1(indexone1(k)+1),:);
    photoninterval=diff(oneburst(:,end-2))*1e-4; % 1 ms timeunit
    cntrate1(k)=1/mean(photoninterval);
end
cntrate2=zeros(length(indexone2),1);
for k=1:length(indexone2)
    oneburst=burstbint3rall2(cumindex2(indexone2(k))+1:cumindex2(indexone2(k)+1),:);
    photoninterval=diff(oneburst(:,end-2))*1e-4; % 1 ms timeunit
    cntrate2(k)=1/mean(photoninterval);
end

drvintvs=drvintv*(-drvpnt:drvpnt);
initparamsA=initparams;
if initparamsA(4) > 0.5, initparamsA(4)=1-initparamsA(4); end     % pb

logmlhdiag=zeros(length(initparams),length(drvintvs));
for kkk=1:length(initparams);
    clear logmlhdiag
    for nnn=1:length(drvintvs);
        initparamone=initparams;
        initparamone(kkk)=initparams(kkk)+initparamsA(kkk)*drvintvs(nnn);
        logmlhdiag(nnn)=mlhTPfitsub(initparamone);
    end
    diagcoeff=polyfit(drvintvs,logmlhdiag,2);
    Hessianmat(kkk,kkk)=diagcoeff(1)/initparamsA(kkk)^2;
end

for kkk=1:length(initparams);
    for mmm=kkk+1:length(initparams);
        clear logmlhoffdiag
        for nnn=1:length(drvintvs);
            initparamone=initparams;
            initparamone(kkk)=initparams(kkk)+initparamsA(kkk)*drvintvs(nnn);
            initparamone(mmm)=initparams(mmm)+initparamsA(mmm)*drvintvs(nnn);
            logmlhoffdiag(nnn)=mlhTPfitsub(initparamone);
        end
        offdiagcoeff=polyfit(drvintvs,logmlhoffdiag,2);
        Hessianmat(kkk,mmm)=(offdiagcoeff(1)-Hessianmat(kkk,kkk)*initparamsA(kkk)^2-Hessianmat(mmm,mmm)*initparamsA(mmm)^2)/2/initparamsA(kkk)/initparamsA(mmm);
        Hessianmat(mmm,kkk)=Hessianmat(kkk,mmm);
    end
end
covmat=inv(2*Hessianmat);
errorparams=sqrt(diag(covmat));

logmlhres=-mlhTPfitsub(initparams);
dellogmlh=logmlhres+mlhTPfitsub([initparams(1:length(initparams)-1); 1e-6]);
numphotons=sum(cumindex1(indexone1+1)-cumindex1(indexone1))+sum(cumindex2(indexone2+1)-cumindex2(indexone2));
bic=-2*logmlhres + length(initparams)*log(numphotons);
end