function logmlhres=mlhTP3cABIsimintEcalindiC(fitparams,intE3cTPT,fixparam,resparamfix,onet3rkk,stateid,chans,nintm,fitid)
%
% logmlhres = mlhTP3cABIsimintEcalindiC(fitparams,intE3cTPT,fixparam,resparamfix,onet3rkk,stateid,chans,nintm,fitid)
%
% Calculate likelihood values for individual transition path data with given
% parameters and arrays of transition path times (trplottime) and fractions of acceptor 1 and acceptor 2 (intE3c).
%
% inputdata (onet3rkk) is a cell array of 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% Input parameters (fitparams) were determined by mlhTP3cABIsimfit0C.
%
% "trplottime" is an array of transition path times for which likelihood values are calculated.
%
% Likelihood value is calculated for the same type of transitions.
%
% "fixparam" contains donor leak into A1 and A2 channel, A1 leak into A2
% channel, and predetermined DA1 and DA2 FRET efficiencies for 3-color-only analysis.
%
% "resparamfix" has 4 parameters State 1 and State 2 FRET efficiencies,
% which are replaced by optimization parameters, relaxation rate, which is
% not so important unless too high because it is scaled by a factor of
% 1000, and fraction of state 1, which does not affect the optimization.
%
% "stateid" is a cell array from which the beginning state is extracted.
%
% "chans" defines the photon detection channels stored in the 5th column of
% the photon trajectory data, onet3rkk.
% chans=[donor A1 A2]
%
% "nintm" is the number of steps of the transition path
%
% "intE3c" is a cell array of acceptor 1 fractions and acceptor 2 fractions.
%
% "fitidin" defines the analysis type: [1 4 5 6] can be used
% 1: 3-color data, 2: DA1, 3: DA2, 4: 3-color + DA1, 5: 3-color + DA2, 6: 3-color + DA1 + DA2
%
% Output is a set of likelihood values with different TP times and acceptor fractions of individual transitions.

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

logmlhfu={zeros(length(indexone1{1}),size(intE3cTPT,1)),zeros(length(indexone2{1}),size(intE3cTPT,1))};
        fitparamid=find(abs(fitparams(:,fitid)) > 0);
        params=fitparams(fitparamid,fitid);
        if any(fitid == [1 4:6])
            darkparams=params(end-3:end);
            params(end-3:end)=[];
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

        switch fitid
            case 1
                resparamtwo=[eff1s1; eff2s1; effs2(1); eff1s1(2)+eff2s1(2); effs3; resparamone];
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
            case {2, 3}
                error('This function does not calculte two-color (DA1 or DA2) likelihoods.');
            case 4
                effs2=[params(5);  eff1s1(2)+eff2s1(2)];
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
            case 5
                effs3=params(5:6); % efficiencies of segment 3 (DA2) from mlhfitparam
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
            case 6
                effs2=[params(5);  eff1s1(2)+eff2s1(2)];
                effs3=params(6:7); % efficiencies of segment 3 (DA2) from mlhfitparam
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                resparamtwo=[eff1s1; eff2s1; effs2; effs3; resparamone];
        end
        
        trcc=1;
        stateidvit=1;
        for ii=1:length(indexone1{trcc})
            for jj=1:size(intE3cTPT,1)
                tpE=intE3cTPT(jj,1:2);
                trplottime=intE3cTPT(jj,3);
                oneburst=burstbint3rall1{trcc}(cumindex1{trcc}(indexone1{trcc}(ii))+1:cumindex1{trcc}(indexone1{trcc}(ii)+1),:);
                indexone=1;
                cumindex=[0; size(oneburst,1)];
                
                logmlhfu{stateidvit}(ii,jj)=logmlhfu{stateidvit}(ii,jj)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,oneburst,cumindex,indexone,cntrate1{1}(ii),stateidvit,nintm,tpE,1);
            end
        end
        stateidvit=2;
        for ii=1:length(indexone2{trcc})
            for jj=1:size(intE3cTPT,1)
                tpE=intE3cTPT(jj,1:2);
                trplottime=intE3cTPT(jj,3);
                oneburst=burstbint3rall2{trcc}(cumindex2{trcc}(indexone2{trcc}(ii))+1:cumindex2{trcc}(indexone2{trcc}(ii)+1),:);
                indexone=1;
                cumindex=[0; size(oneburst,1)];

                logmlhfu{stateidvit}(ii,jj)=logmlhfu{stateidvit}(ii,jj)+mlhrateeigtrpathrevcalvit3cABIparC([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,oneburst,cumindex,indexone,cntrate2{1}(ii),stateidvit,nintm,tpE,1);
            end
        end
        
%        if fitfu==1 || fitfu==2, logmlh=-logmlhfu(fitfu); end
%        if fitfu==3, logmlh=-sum(logmlhfu); end
        logmlhres=logmlhfu;
end