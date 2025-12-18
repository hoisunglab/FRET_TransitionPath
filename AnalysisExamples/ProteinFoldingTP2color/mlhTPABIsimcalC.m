function logmlhres=mlhTPABIsimcalC(params,trplottime,effb,resparamfix,onet3rkk,stateid,nintm,isincE)
%
% logmlhres = mlhTPABIsimcalC(params,trplottime,effb,resparamfix,onet3rkk,stateid,nintm,isincE)
%
% Calculate a likelihood value for transition path data with given
% parameters and an array of transition path times.
%
% inputdata (onet3rkk) is a cell array of 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% Input parameters (params) contain 4 optimization parameters determined by mlhTPABIsimfit0C.
%
% "trplottime" is an array of transition path times for which likelihood values are calculated.
%
% Likelihood value is calculated for the same type of transitions.
%
% "effb" is the FRET efficiency of the acceptor dark state.
%
% "resparamfix" has 4 parameters State 1 and State 2 FRET efficiencies,
% which are replaced by optimization parameters, relaxation rate, which is
% not so important unless too high because it is scaled by a factor of
% 1000, and fraction of state 1, which does not affect the optimization.
%
% "stateid" is a cell array from which the beginning state is extracted.
%
% "nintm" is the number of steps of the transition path
%
% "inincE" determines the FRET efficiency distribution in a transition path (E_TP).
% [] or 1: E_TPs interpolates eff1 and eff2 with the number of steps.
% 2: E_TP = (eff1 + eff2)/2 for all TP steps
% x (0 < x < 1): E_TP = x for all TP steps
%
% Output is a set of likelihood values with different TP times.

        fretEF=params(1);
        fretEU=params(2);

        ktobrt=params(3);
        frnb0=params(4);        % for a photon count rate of 100 ms-1
        
        resparamone=resparamfix;resparamone(3)=resparamone(3)/1000;
        resparamone(1:2)=[fretEF fretEU];

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

        logmlhfu(1:2)={0};
        stateidvit=1;
        logmlhfu{stateidvit}=mlhrateeigtrpathrevcalvitABIC([resparamone; effb; ktobrt; frnb0],trplottime,burstbint3rall1,cumindex1,indexone1,cntrate1,stateidvit,nintm,isincE);
        stateidvit=2;
        logmlhfu{stateidvit}=mlhrateeigtrpathrevcalvitABIC([resparamone; effb; ktobrt; frnb0],trplottime,burstbint3rall2,cumindex2,indexone2,cntrate2,stateidvit,nintm,isincE);

        logmlhres=logmlhfu;
end