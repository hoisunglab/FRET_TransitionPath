function [errorparams logmlhres dellogmlh bic]=mlhTPABIfit1TPindierror(initparams,effb,resparamfix,oneburst,stateid,fitparamid,trcc,drvintv,drvpnt)
%
% [errorparams logmlhres dellogmlh bic] = mlhTPABIfit1TPindierror(initparams,effb,resparamfix,oneburst,stateid,fitparamid,trcc,drvintv,drvpnt)
%
% Calculate errors (standard deviation) of the parameter determined by mlhTPABIfit1TPindi.m.  
% Outputs are parameter errors, log likelihood value, log likelihood value
% difference from an instantaneous transition model, and Bayesian
% Information Criterion (BIC) calculated using the optimized parameters.

    function logmlh=mlhTPfitsub(params)
        fretEF=params(1);
        fretEU=params(2);
        ktobrt=params(3);
        frnb0=params(4);        % for a photon count rate of 100 ms-1
        trplottime=params(5);
        Eint=params(6);

        resparamone=resparamfix(end-1:end);resparamone(1)=resparamone(1)/1000;
        resparamone=[fretEF; fretEU; resparamone];

        if trcc == 3, effb=fretEU; end      % For DA2 segments

%        logmlhfu=mlhrateeigtrpathrevcalvitABIC([resparamone; effb; ktobrt; frnb0],[trplottime trplottime],oneburst,cumindex,indexone,cntrate,stateid,nintm,Eint);
        logmlhfu=mlhrateeigtrpathrevcalvitintEABI([resparamone; trplottime; effb; ktobrt; frnb0],oneburst,cumindex,indexone,stateid,1,Eint);

        logmlh=-logmlhfu;
    end

burstphotons=size(oneburst,1);
indexone=1:length(burstphotons);
cumindex=[0; cumsum(burstphotons)];

drvintvs=1+drvintv*(-drvpnt:drvpnt);

errorparams=zeros(size(initparams));
logmlhdiag=zeros(length(fitparamid),length(drvintvs));
for kkk=1:length(fitparamid);
    clear logmlhdiag
    for nnn=1:length(drvintvs);
        initparamone=initparams;
        initparamone(fitparamid(kkk))=initparams(fitparamid(kkk))*drvintvs(nnn);
        logmlhdiag(nnn)=mlhTPfitsub(initparamone);
    end
    diagcoeff=polyfit(drvintvs,logmlhdiag,2);
    Hessianmat(kkk,kkk)=diagcoeff(1)/initparams(fitparamid(kkk))^2;
end

for kkk=1:length(fitparamid);
    for mmm=kkk+1:length(fitparamid);
        clear logmlhoffdiag
        for nnn=1:length(drvintvs);
            initparamone=initparams;
            initparamone(fitparamid(kkk))=initparams(fitparamid(kkk))*drvintvs(nnn);initparamone(fitparamid(mmm))=initparams(fitparamid(mmm))*drvintvs(nnn);
            logmlhoffdiag(nnn)=mlhTPfitsub(initparamone);
        end
        offdiagcoeff=polyfit(drvintvs,logmlhoffdiag,2);
        Hessianmat(kkk,mmm)=(offdiagcoeff(1)-Hessianmat(kkk,kkk)*initparams(fitparamid(kkk))^2-Hessianmat(mmm,mmm)*initparams(fitparamid(mmm))^2)/2/initparams(fitparamid(kkk))/initparams(fitparamid(mmm));
        Hessianmat(mmm,kkk)=Hessianmat(kkk,mmm);
    end
end
covmat=inv(2*Hessianmat);
errorparams(fitparamid)=sqrt(diag(covmat));

logmlhres=-mlhTPfitsub(initparams);
dellogmlh=logmlhres+mlhTPfitsub([initparams(1:length(initparams)-2); 1e-6; initparams(end)]);
numphotons=sum(cumindex(indexone+1)-cumindex(indexone));
bic=-2*logmlhres + length(fitparamid)*log(numphotons);
end