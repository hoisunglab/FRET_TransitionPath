function respconv=mlhTPABIsimfitEintNintmEC(initparams,LUbounds,effb,resparamfix,onet3rkk,stateid,fitfu,nintm)
%
% resparam = mlhTPABIsimfitEintNintmEC(initparams,LUbounds,effb,resparamfix,onet3rkk,stateid,fitfu,nintm)
%
% Optimize transition path parameters including the transition path time
% (t_TP) and the FRET efficiency of the transition path (E_TP) of an n-step TP model (n = 1, 2, or 3)
% by maximizing likelihood function calculated by mlhrateeigtrpathrevcalvitABIC.m
%
% inputdata (onet3rkk) is a cell array of 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% There are 6 optimization parameters:
% fretEF: Apparent FRET efficiency of State 1 (e.g., folded)
% fretEU: Apparent FRET efficiency of State 2 (unfolded)
% ktobrt: Rate constant of a transition from acceptor dark to bright state
% frnb0: Acceptor bright state population at a reference photon count rate (100 ms-1)
% trplottime: transition path time, t_TP
% Eint: FRET efficiency of the TP, E_TP, which is the same for all TP steps
%
% "effb" is the FRET efficiency of the acceptor dark state.
%
% "resparamfix" has 4 parameters State 1 and State 2 FRET efficiencies,
% which are not used (optimization parameters), relaxation rate, which is
% not so important unless too high because it is scaled by a factor of
% 1000, and fraction of state 1, which does not affect the optimization.
%
% "stateid" is a cell array from which the beginning state is extracted.
%
% "fitfu" selects the likelihood output.
% 1: transitions from 1 to 2
% 2: transitions from 2 to 1
% 3: all transitions
%
% "nintm" is the number of steps of the transition path
%
% Output is a set of optimized parameters.

    function logmlh=mlhTPfitsub(params)
        pconv=diff(LUbounds,1,2)./(1+params.^2)+LUbounds(:,1);
        fretEF=pconv(1);
        fretEU=pconv(2);
        ktobrt=pconv(3);
        frnb0=pconv(4);        % for a photon count rate of 100 ms-1
        trplottime=pconv(5);
        Eint=pconv(6);

        resparamone=resparamfix;resparamone(3)=resparamone(3)/1000;
        resparamone(1:2)=[fretEF fretEU];

        logmlhfu(1:2)=0;
        stateidvit=1;
        logmlhfu(stateidvit)=mlhrateeigtrpathrevcalvitABIC([resparamone; effb; ktobrt; frnb0],[trplottime trplottime],burstbint3rall1,cumindex1,indexone1,cntrate1,stateidvit,nintm,Eint);
        stateidvit=2;
        logmlhfu(stateidvit)=mlhrateeigtrpathrevcalvitABIC([resparamone; effb; ktobrt; frnb0],[trplottime trplottime],burstbint3rall2,cumindex2,indexone2,cntrate2,stateidvit,nintm,Eint);

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

options=optimset('MaxFunEval', 100000, 'MaxIter', 100000);

invinitparams=sqrt(diff(LUbounds,1,2)./(initparams-LUbounds(:,1))-1);
resparams=fminsearch(@mlhTPfitsub,invinitparams,options);
respconv=diff(LUbounds,1,2)./(1+resparams.^2)+LUbounds(:,1);

end