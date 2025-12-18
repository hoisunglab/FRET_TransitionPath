function [errorparams logmlhres dellogmlh bic]=mlhTP3cABIfit1TPindierror(initparams,fixparam,resparamfix,oneburst,stateid,fitparamid,drvintv,drvpnt)
%
% [errorparams logmlhres dellogmlh bic] = mlhTP3cABIfit1TPindierror(initparams,fixparam,resparamfix,oneburst,stateid,fitparamid,drvintv,drvpnt)
%
% Calculate errors (standard deviation) of the parameter determined by mlhTP3cABIfit1TPindi.m.  
% Outputs are parameter errors, log likelihood value, log likelihood value
% difference from an instantaneous transition model, and Bayesian
% Information Criterion (BIC) calculated using the optimized parameters.

    function logmlh=mlhTPfitsub(params)
            darkparams=params(5:8);
            ktobrt=darkparams(1:2);  % rate from dark to bright state of A1 and A2
            frnb0=darkparams(3:4);    % A1 and A2 bright population at photon count rate of 100 ms-1
            eff1s1=params(1:2);             % 3 color E1app
            eff2s1=params(3:4);    % 3 color E2app
        
        resparamone=resparamfix(end-1:end);resparamone(1)=resparamone(1)/1000;
        effb1=fixparam{1}(1);  % E1 of A1,A2 dark state (= A1/(A2 + A1 + D))
        effb2=fixparam{1}(2);  % E2 of A1,A2 dark state (= A2/(A2 + A1 + D)) (depends on the A2-protein concentration)
        effs2=fixparam{2};  % Pre-determined E_DA1 (= (A1+A2)/(D+A1+A2)) will be used unless modified below using fitting parameters
        effs3=fixparam{3};  % Pre-determined E_DA2 (= A2/(D+A1+A2)) will be used unless modified below using fitting parameters
        effb12=fixparam{4};  % E12 when no A2: A1 leak into A2 channel
        tpE2c=[];
        if length(fixparam) == 5, tpE2c=fixparam{5}; end

                trplottime=params(9);
                tpE=params(10:11)';  % Eint
                tpE=[tpE tpE2c];
        
                resparamtwo=[eff1s1; eff2s1; effs2(1); eff1s1(2)+eff2s1(2); effs3; resparamone];
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                logmlhfu=mlhrateeigtrpathrevcalvit3cintEABI([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,oneburst,cumindex,indexone,stateid,1,tpE);
        
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
dellogmlh=logmlhres+mlhTPfitsub([initparams(1:length(initparams)-3); 1e-6; initparams(end-1:end)]);
numphotons=sum(cumindex(indexone+1)-cumindex(indexone));
bic=-2*logmlhres + length(fitparamid)*log(numphotons);
end