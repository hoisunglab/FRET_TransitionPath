function [errorparams logmlhres bic]=mlhTP3cABIsimfitEintParerrorC(fitparams,fixparam,resparamfix,onet3rkk,stateid,chans,fitfu,nintm,fitidin,drvintv,drvpnt)
%
% [errorparams logmlhres bic] = mlhTP3cABIsimfitEintParerrorC(fitparams,fixparam,resparamfix,onet3rkk,stateid,chans,fitfu,nintm,fitidin,drvintv,drvpnt)
%
% Calculate errors (standard deviation) of the parameter determined by mlhTP3cABIsimfitEintParC.m.  
% Outputs are parameter errors, log likelihood, and Bayesian Information Criterion (BIC) values.

    function logmlh=mlhTPfitsub(params)
        if any(fitid == [1 4:6])
            switch fitid
                case 1
                    darkparams=params(5:8);
                case 4
                    darkparams=params(6:9);
                case 5
                    darkparams=params(7:10);
                case 6
                    darkparams=params(8:11);
            end
            ktobrt=darkparams(1:2);  % rate from dark to bright state of A1 and A2
            frnb0=darkparams(3:4);    % A1 and A2 bright population at photon count rate of 100 ms-1
            eff1s1=params(1:2);             % 3 color E1app
            eff2s1=params(3:4);    % 3 color E2app
        end
        
        resparamone=resparamfix(end-1:end);resparamone(1)=resparamone(1)/1000;
        effb1=fixparam{1}(1);  % E1 of A1,A2 dark state (= A1/(A2 + A1 + D))
        effb2=fixparam{1}(2);  % E2 of A1,A2 dark state (= A2/(A2 + A1 + D)) (depends on the A2-protein concentration)
        effs2=fixparam{2};  % Pre-determined E_DA1 (= (A1+A2)/(D+A1+A2)) will be used unless modified below using fitting parameters
        effs3=fixparam{3};  % Pre-determined E_DA2 (= A2/(D+A1+A2)) will be used unless modified below using fitting parameters
        effb12=fixparam{4};  % E12 when no A2: A1 leak into A2 channel
        tpE2c=[];
        if length(fixparam) == 5, tpE2c=fixparam{5}; end

        logmlhfu(1:2)=0;
        switch fitid
            case 1
                trplottime=params(8+(1:nparI));
                tpE=params(8+(nparI+1:3*nparI));  % Eint
                tpE=[tpE(1:nparI) tpE(nparI+1:2*nparI) tpE2c];
                pItemp=params(8+(3*nparI+1:4*nparI-1));  % relative fractions of pathways (proportional to the TP times).
                if ~isempty(pItemp)
                    pfactor=cumprod([1; 1-pItemp(1:end-1)]);
                    pI=pItemp.*pfactor; % Convert to fraction pI
                    pI=[pI; 1-sum(pI)];
                else
                    pI=1;   % single pathway
                end
        
                resparamtwo=[eff1s1; eff2s1; effs2(1); eff1s1(2)+eff2s1(2); effs3; resparamone];
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall1{1},cumindex1{1},indexone1{1},cntrate1{1},stateidvit,nintm,tpE,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall2{1},cumindex2{1},indexone2{1},cntrate2{1},stateidvit,nintm,tpE,pI);
            case 2
                trplottime=params(4+(1:nparI));
                Eint=params(4+(nparI+1:2*nparI));  % Eint
                pItemp=params(4+(2*nparI+1:3*nparI-1));  % relative fractions of pathways (proportional to the TP times).
                if ~isempty(pItemp)
                    pfactor=cumprod([1; 1-pItemp(1:end-1)]);
                    pI=pItemp.*pfactor; % Convert to fraction pI
                    pI=[pI; 1-sum(pI)];
                else
                    pI=1;   % single pathway
                end
        
                resparamtwo=[params(1:2); resparamone];
                effb=effb1+effb2;
                ktobrt2=params(3);
                frnb02=params(4);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{2},cumindex1{2},indexone1{2},cntrate1{2},stateidvit,nintm,Eint,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{2},cumindex2{2},indexone2{2},cntrate2{2},stateidvit,nintm,Eint,pI);
            case 3
                trplottime=params(4+(1:nparI));
                Eint=params(4+(nparI+1:2*nparI));  % Eint
                pItemp=params(4+(2*nparI+1:3*nparI-1));  % relative fractions of pathways (proportional to the TP times).
                if ~isempty(pItemp)
                    pfactor=cumprod([1; 1-pItemp(1:end-1)]);
                    pI=pItemp.*pfactor; % Convert to fraction pI
                    pI=[pI; 1-sum(pI)];
                else
                    pI=1;   % single pathway
                end
        
                resparamtwo=[params(1:2); resparamone];
                effb2=resparamtwo(2);
                effb=effb2;
                ktobrt2=params(3);
                frnb02=params(4);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{3},cumindex1{3},indexone1{3},cntrate1{3},stateidvit,nintm,Eint,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{3},cumindex2{3},indexone2{3},cntrate2{3},stateidvit,nintm,Eint,pI);
            case 4
                trplottime=params(9+(1:nparI));
                tpE=params(9+(nparI+1:3*nparI));  % Eint
                Eint=params(9+(3*nparI+1:4*nparI));  % Eint
                tpE=[tpE(1:nparI) tpE(nparI+1:2*nparI) Eint];
                if size(tpE2c,2) == 2, tpE=[tpE tpE2c(:,2)]; end 
                pItemp=params(9+(4*nparI+1:5*nparI-1));  % relative fractions of pathways (proportional to the TP times).
                if ~isempty(pItemp)
                    pfactor=cumprod([1; 1-pItemp(1:end-1)]);
                    pI=pItemp.*pfactor; % Convert to fraction pI
                    pI=[pI; 1-sum(pI)];
                else
                    pI=1;   % single pathway
                end
        
                effs2=[params(5);  eff1s1(2)+eff2s1(2)];
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall1{1},cumindex1{1},indexone1{1},cntrate1{1},stateidvit,nintm,tpE,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall2{1},cumindex2{1},indexone2{1},cntrate2{1},stateidvit,nintm,tpE,pI);

                resparamtwo=[effs2; resparamone];
                effb=effb1+effb2;
                ktobrt2=ktobrt(1);
                frnb02=frnb0(1);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{2},cumindex1{2},indexone1{2},cntrate1{2},stateidvit,nintm,Eint,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{2},cumindex2{2},indexone2{2},cntrate2{2},stateidvit,nintm,Eint,pI);
            case 5
                trplottime=params(10+(1:nparI));
                tpE=params(10+(nparI+1:3*nparI));  % Eint
                Eint=params(10+(3*nparI+1:4*nparI));  % Eint DA2
                tpE=[tpE(1:nparI) tpE(nparI+1:2*nparI) zeros(size(Eint)) Eint];
                if size(tpE2c,2) > 0, tpE(:,3)=tpE2c(:,1); end 
                pItemp=params(10+(4*nparI+1:5*nparI-1));  % relative fractions of pathways (proportional to the TP times).
                if ~isempty(pItemp)
                    pfactor=cumprod([1; 1-pItemp(1:end-1)]);
                    pI=pItemp.*pfactor; % Convert to fraction pI
                    pI=[pI; 1-sum(pI)];
                else
                    pI=1;   % single pathway
                end
        
                effs3=params(5:6); % efficiencies of segment 3 (DA2) from mlhfitparam
                effb2=effs3(2);
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall1{1},cumindex1{1},indexone1{1},cntrate1{1},stateidvit,nintm,tpE,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall2{1},cumindex2{1},indexone2{1},cntrate2{1},stateidvit,nintm,tpE,pI);

                resparamtwo=[effs3; resparamone];
                effb=effb2;
                ktobrt2=ktobrt(2);
                frnb02=frnb0(2);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{3},cumindex1{3},indexone1{3},cntrate1{3},stateidvit,nintm,Eint,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{3},cumindex2{3},indexone2{3},cntrate2{3},stateidvit,nintm,Eint,pI);
            case 6
                trplottime=params(11+(1:nparI));
                tpE=params(11+(nparI+1:3*nparI));  % Eint
                Eint1=params(11+(3*nparI+1:4*nparI));  % Eint DA1
                Eint2=params(11+(4*nparI+1:5*nparI));  % Eint DA2
                tpE=[tpE(1:nparI) tpE(nparI+1:2*nparI) Eint1 Eint2];
                pItemp=params(11+(5*nparI+1:6*nparI-1));  % relative fractions of pathways (proportional to the TP times).
                if ~isempty(pItemp)
                    pfactor=cumprod([1; 1-pItemp(1:end-1)]);
                    pI=pItemp.*pfactor; % Convert to fraction pI
                    pI=[pI; 1-sum(pI)];
                else
                    pI=1;   % single pathway
                end
        
                effs2=[params(5);  eff1s1(2)+eff2s1(2)];
                effs3=params(6:7); % efficiencies of segment 3 (DA2) from mlhfitparam
                effb2=effs3(2);
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall1{1},cumindex1{1},indexone1{1},cntrate1{1},stateidvit,nintm,tpE,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall2{1},cumindex2{1},indexone2{1},cntrate2{1},stateidvit,nintm,tpE,pI);

                resparamtwo=[effs2; resparamone];
                effb=effb1+effb2;
                ktobrt2=ktobrt(1);
                frnb02=frnb0(1);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{2},cumindex1{2},indexone1{2},cntrate1{2},stateidvit,nintm,Eint1,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{2},cumindex2{2},indexone2{2},cntrate2{2},stateidvit,nintm,Eint1,pI);

                resparamtwo=[effs3; resparamone];
                effb=effb2;
                ktobrt2=ktobrt(2);
                frnb02=frnb0(2);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{3},cumindex1{3},indexone1{3},cntrate1{3},stateidvit,nintm,Eint2,pI);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIparC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{3},cumindex2{3},indexone2{3},cntrate2{3},stateidvit,nintm,Eint2,pI);
        end
        
            if fitfu==1 || fitfu==2, logmlh=-logmlhfu(fitfu); end
            if fitfu==3, logmlh=-sum(logmlhfu); end
       
    end

donchan=chans(3);
ac1chan=chans(2);
ac2chan=chans(1);

trccmax=length(onet3rkk);

burstbint3rall1={[],[],[]};burstbint3rall2={[],[],[]};
burstphotons1={[],[],[]};burstphotons2={[],[],[]};
for trcc=1:trccmax
    stidnum=[0 0];
    for k=1:length(onet3rkk{trcc})
        stateidone=stateid{trcc}{k}(1,3);
        stidnum(stateidone)=stidnum(stateidone)+1;
    end

    stidcount=[1 1];
    burstbint3rone1=[];burstbint3rone2=[];
    for k=1:length(onet3rkk{trcc});
        stateidone=stateid{trcc}{k}(1,3);
        onet3r=onet3rkk{trcc}{k};
    
        if trcc == 1
            donid=find(onet3r(:,end) == donchan);
            ac1id=find(onet3r(:,end) == ac1chan);
            onet3r(:,end)=1;
            onet3r(ac1id,end)=2;
            onet3r(donid,end)=3;  % Convert colors for 3-color analysis
        elseif trcc == 2
            donid=find(onet3r(:,end) == donchan);
            onet3r(:,end)=1;
            onet3r(donid,end)=2;
        else
            ac2id=find(onet3r(:,end) == ac2chan);
            onet3r(:,end)=2;
            onet3r(ac2id,end)=1;
        end                        
        
        if stateidone == 1,
            burstbint3rone1=[burstbint3rone1; onet3r];
            burstphotons1{trcc}=[burstphotons1{trcc}; length(onet3r(:,1))];
            if mod(stidcount(stateidone),3000) == 0 || stidcount(stateidone) == stidnum(stateidone)
                burstbint3rall1{trcc}=[burstbint3rall1{trcc}; burstbint3rone1];
                burstbint3rone1=[];
            end
            stidcount(stateidone)=stidcount(stateidone)+1;
        else
            burstbint3rone2=[burstbint3rone2; onet3r];
            burstphotons2{trcc}=[burstphotons2{trcc}; length(onet3r(:,1))];
            if mod(stidcount(stateidone),3000) == 0 || stidcount(stateidone) == stidnum(stateidone)
                burstbint3rall2{trcc}=[burstbint3rall2{trcc}; burstbint3rone2];
                burstbint3rone2=[];
            end
            stidcount(stateidone)=stidcount(stateidone)+1;
        end
    end
end

for trcc=1:trccmax
    indexone1{trcc}=1:length(burstphotons1{trcc});
    cumindex1{trcc}=[0; cumsum(burstphotons1{trcc})];
    indexone2{trcc}=1:length(burstphotons2{trcc});
    cumindex2{trcc}=[0; cumsum(burstphotons2{trcc})];

    cntrate1{trcc}=zeros(length(indexone1{trcc}),1);
    for k=1:length(indexone1{trcc})
        oneburst=burstbint3rall1{trcc}(cumindex1{trcc}(indexone1{trcc}(k))+1:cumindex1{trcc}(indexone1{trcc}(k)+1),:);
        photoninterval=diff(oneburst(:,end-2))*1e-4; % 1 ms timeunit
        cntrate1{trcc}(k)=1/mean(photoninterval);
    end
    cntrate2{trcc}=zeros(length(indexone2{trcc}),1);
    for k=1:length(indexone2{trcc})
        oneburst=burstbint3rall2{trcc}(cumindex2{trcc}(indexone2{trcc}(k))+1:cumindex2{trcc}(indexone2{trcc}(k)+1),:);
        photoninterval=diff(oneburst(:,end-2))*1e-4; % 1 ms timeunit
        cntrate2{trcc}(k)=1/mean(photoninterval);
    end
end

drvintvs=drvintv*(-drvpnt:drvpnt);

nparI=(length(fitparams)-11)/6;  % number of parallel pathways
errorparams=zeros(size(fitparams));
logmlhres=zeros(6,1);
bic=logmlhres;
for fitid=fitidin    % 1: 3 color, 2: 2 color, 3: all
    if ~(fitid == 2 && nparI > 1)
        clear Hessianmat
        fitparamid=find(abs(fitparams(:,fitid)) > 0);
        initparams=fitparams(fitparamid,fitid);
    
        initparamsA=initparams;
        if fitid == 2 || fitid == 3
            if initparamsA(end) > 0.5, initparamsA(end)=1-initparamsA(end); end     % pb, A1 or A2
        else
            if initparamsA(end-1) > 0.5, initparamsA(end-1)=1-initparamsA(end-1); end     % pb, A1
            if initparamsA(end) > 0.5, initparamsA(end)=1-initparamsA(end); end     % pb, A2
        end
    
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
        errorparamsone=sqrt(diag(covmat));
        
        errorparams(fitparamid,fitid)=errorparamsone;
        logmlhres(fitid)=-mlhTPfitsub(initparams);
    
        switch fitid
            case 1
                trcc=1;numphotons=sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
            case 2
                trcc=2;numphotons=sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
            case 3
                trcc=3;numphotons=sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
            case 4
                trcc=1;numphotons=sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
                trcc=2;numphotons=numphotons+sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
            case 5
                trcc=1;numphotons=sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
                trcc=3;numphotons=numphotons+sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
            case 6
                trcc=1;numphotons=sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
                trcc=2;numphotons=numphotons+sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
                trcc=3;numphotons=numphotons+sum(cumindex1{trcc}(indexone1{trcc}+1)-cumindex1{trcc}(indexone1{trcc}))+sum(cumindex2{trcc}(indexone2{trcc}+1)-cumindex2{trcc}(indexone2{trcc}));
        end
        bic(fitid)=-2*logmlhres(fitid) + length(initparams)*log(numphotons);
    end
end     % fitid
end