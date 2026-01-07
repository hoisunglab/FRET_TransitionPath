function respconv=mlhTP3cABIsimfit0C(initparams,LUbounds,fixparam,resparamfix,onet3rkk,stateid,chans,fitfu,nintm,isincE,fitidin)
%
% resparam = mlhTP3cABIsimfit0C(initparams,LUbounds,fixparam,resparamfix,onet3rkk,stateid,chans,fitfu,nintm,isincE,fitidin)
%
% Optimize transition path parameters with instantaneous transition model (t_TP = 0)
% by maximizing likelihood function calculated by
% mlhrateeigtrpathrevcalvit3cABIparC.m (3-color) and mlhrateeigtrpathrevcalvitABIC.m (2-color).
% This analysis is for a 3-color binding experiment. IDP labeled with donor (D) 
% and acceptor 1 (A1) and binding partner is labeled with acceptor 2 (A2).
%
% inputdata (onet3rkk) is a cell array of 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% Optimization parameters vary depending on the data types analyzed (see "fitidin" below):
% eff1s1: Acceptor 1 (A1) fractions of State 1 (bound) and State 2 (unbound)
% eff2s1: Acceptor 2 (A2) fractions of State 1 and State 2
% effs2: Apparent FRET efficiency E1 for the data with donor and A1 (DA1)
% effs3: Apparent FRET efficiency E2 for the data with donor and A2 (DA2)
% ktobrt: Rate constants of transitions from A1 and A2 dark to bright states
% frnb0: A1 and A2 bright state populations at a reference photon count rate (100 ms-1)
%
% "fixparam" contains donor leak into A1 and A2 channel, A1 leak into A2
% channel, and predetermined DA1 and DA2 FRET efficiencies for 3-color-only analysis.
%
% "resparamfix" has 4 parameters State 1 and State 2 FRET efficiencies,
% which are not used (optimization parameters), relaxation rate, which is
% not so important unless too high because it is scaled by a factor of
% 1000, and fraction of state 1, which does not affect the optimization.
%
% "stateid" is a cell array from which the beginning state is extracted.
%
% "chans" defines the photon detection channels stored in the 5th column of
% the photon trajectory data, onet3rkk.
% chans=[donor A1 A2]
%
% "fitfu" selects the likelihood output.
% 1: transitions from 1 to 2
% 2: transitions from 2 to 1
% 3: all transitions
%
% "nintm" is the number of steps of the transition path
%
% "inincE" determines the FRET efficiency distribution in a transition path (E_TP).
% [] or 1: E_TPs interpolates eff1 and eff2 with the number of steps.
% 2: E_TP = (eff1 + eff2)/2 for all TP steps
% x (0 < x < 1): E_TP = x for all TP steps
%
% "fitidin" defines the analysis type:
% 1: 3-color data, 2: DA1, 3: DA2, 4: 3-color + DA1, 5: 3-color + DA2, 6: 3-color + DA1 + DA2
%
% Output is a set of optimized parameters.

    function logmlh=mlhTPfitsub(params)
        pconv=diff(LUbounds(fitparamid{fitid},:),1,2)./(1+params.^2)+LUbounds(fitparamid{fitid},1);
        if any(fitid == [1 4:6])
            darkparams=pconv(end-3:end);
            pconv(end-3:end)=[];
            ktobrt=darkparams(1:2);  % rate from dark to bright state of A1 and A2
            frnb0=darkparams(3:4);    % A1 and A2 bright population at photon count rate of 100 ms-1
            eff1s1=pconv(1:2);             % 3 color E1app
            eff2s1=pconv(3:4);    % 3 color E2app
        end
        
        resparamone=resparamfix(end-1:end);resparamone(1)=resparamone(1)/1000;
        effb1=fixparam{1}(1);  % E1 of A1,A2 dark state (= A1/(A2 + A1 + D))
        effb2=fixparam{1}(2);  % E2 of A1,A2 dark state (= A2/(A2 + A1 + D)) (depends on the A2-protein concentration)
        effs2=fixparam{2};  % Pre-determined E_DA1 (= (A1+A2)/(D+A1+A2)) will be used unless modified below using fitting parameters
        effs3=fixparam{3};  % Pre-determined E_DA2 (= A2/(D+A1+A2)) will be used unless modified below using fitting parameters
        effb12=fixparam{4};  % E12 when no A2: A1 leak into A2 channel

        trplottime=1e-6;

        logmlhfu(1:2)=0;
        switch fitid
            case 1
                resparamtwo=[eff1s1; eff2s1; effs2(1); eff1s1(2)+eff2s1(2); effs3; resparamone];
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall1{1},cumindex1{1},indexone1{1},cntrate1{1},stateidvit,nintm,[],1);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall2{1},cumindex2{1},indexone2{1},cntrate2{1},stateidvit,nintm,[],1);
            case 2
                resparamtwo=[pconv(1:2); resparamone];
                effb=effb1+effb2;
                ktobrt2=pconv(3);
                frnb02=pconv(4);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{2},cumindex1{2},indexone1{2},cntrate1{2},stateidvit,nintm,isincE);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{2},cumindex2{2},indexone2{2},cntrate2{2},stateidvit,nintm,isincE);
            case 3 % DA2 transitions
                resparamtwo=[pconv(1:2); resparamone];
                effb2=resparamtwo(2);
                effb=effb2;
                ktobrt2=pconv(3);
                frnb02=pconv(4);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{3},cumindex1{3},indexone1{3},cntrate1{3},stateidvit,nintm,isincE);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{3},cumindex2{3},indexone2{3},cntrate2{3},stateidvit,nintm,isincE);
            case 4
                effs2=[pconv(5);  eff1s1(2)+eff2s1(2)];
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall1{1},cumindex1{1},indexone1{1},cntrate1{1},stateidvit,nintm,[],1);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall2{1},cumindex2{1},indexone2{1},cntrate2{1},stateidvit,nintm,[],1);

                resparamtwo=[effs2; resparamone];
                effb=effb1+effb2;
                ktobrt2=ktobrt(1);
                frnb02=frnb0(1);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{2},cumindex1{2},indexone1{2},cntrate1{2},stateidvit,nintm,isincE);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{2},cumindex2{2},indexone2{2},cntrate2{2},stateidvit,nintm,isincE);
            case 5 % 3c + DA2 transitions
                % 3c
                effs3=pconv(5:6); % efficiencies of segment 3 (DA2) from mlhfitparam
                effb2=effs3(2);
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall1{1},cumindex1{1},indexone1{1},cntrate1{1},stateidvit,nintm,[],1);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall2{1},cumindex2{1},indexone2{1},cntrate2{1},stateidvit,nintm,[],1);
                % DA2
                resparamtwo=[effs3; resparamone];
                effb=effb2;
                ktobrt2=ktobrt(2);
                frnb02=frnb0(2);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{3},cumindex1{3},indexone1{3},cntrate1{3},stateidvit,nintm,isincE);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{3},cumindex2{3},indexone2{3},cntrate2{3},stateidvit,nintm,isincE);
            case 6 % 3c + DA1 + DA2 transitions
                % 3c
                effs2=[pconv(5);  eff1s1(2)+eff2s1(2)];
                effs3=pconv(6:7); % efficiencies of segment 3 (DA2) from mlhfitparam
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                effb2=effs3(2);
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall1{1},cumindex1{1},indexone1{1},cntrate1{1},stateidvit,nintm,[],1);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,burstbint3rall2{1},cumindex2{1},indexone2{1},cntrate2{1},stateidvit,nintm,[],1);
                % DA1
                resparamtwo=[effs2; resparamone];
                effb=effb1+effb2;
                ktobrt2=ktobrt(1);
                frnb02=frnb0(1);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{2},cumindex1{2},indexone1{2},cntrate1{2},stateidvit,nintm,isincE);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{2},cumindex2{2},indexone2{2},cntrate2{2},stateidvit,nintm,isincE);
                % DA2
                resparamtwo=[effs3; resparamone];
                effb=effb2;
                ktobrt2=ktobrt(2);
                frnb02=frnb0(2);
                stateidvit=1;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall1{3},cumindex1{3},indexone1{3},cntrate1{3},stateidvit,nintm,isincE);
                stateidvit=2;
                logmlhfu(stateidvit)=logmlhfu(stateidvit)+mlhrateeigtrpathrevcalvitABIC([resparamtwo; effb; ktobrt2; frnb02],[trplottime trplottime],burstbint3rall2{3},cumindex2{3},indexone2{3},cntrate2{3},stateidvit,nintm,isincE);
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
options=optimset('MaxFunEval', 100000, 'MaxIter', 100000);

invinitparams=sqrt(diff(LUbounds,1,2)./(initparams-LUbounds(:,1))-1);
respconv=zeros(length(initparams),6);
fitparamid={[1:4 9:12],[5 6 9 11], [7 8 10 12],[1:5 9:12], [1:4 7:12], [1:5 7:12]};
for fitid=fitidin    % 1: 3 color, 2: DA1, 3: DA2, 4: 3c+DA1, 5: 3c+DA2, 6: all
    fitid
    resparams=fminsearch(@mlhTPfitsub,invinitparams(fitparamid{fitid}),options);
    respconv(fitparamid{fitid},fitid)=diff(LUbounds(fitparamid{fitid},:),1,2)./(1+resparams.^2)+LUbounds(fitparamid{fitid},1);
end

end