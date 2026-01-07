function logmlhval=mlhrateeigtrpathrevcalvit3cABIparC(params,TPtime,frburstdata,cumindex,indexone,cntrate,stateid,numstep,tpE,pI)
%
% logmlhval = mlhrateeigtrpathrevcalvit3cABIparC(params,TPtime,frburstdata,cumindex,indexone,cntrate,stateid,numstep,tpE,pI)
%
% Calculate a likelihood value for transition path data with given
% parameters and transition path time.
%
% inputdata (frburstdata) is a 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% Input parameters (params) contain 17 parameters:
%
% eff1: Acceptor 1 fractions of State 1 and State 2 (unbound)
% eff2: Acceptor 2 fractions of State 1 and State 2 (unbound)
% effs2: Apparent FRET efficiency E1 for the data with donor and A1 (DA1)
% effs3: Apparent FRET efficiency E2 for the data with donor and A2 (DA2)
% ratesum: relaxation rate of the two-state kinetic model (between states 1 and 2)
% fr1: fraction of State 1 in the two-state model
% effb1, effb2: Donor leak into A1 and A2 channels
% effb12: A1 leak into A2 channel
% ktobrt: Rate constants of transitions from A1 and A2 dark to bright states
% frnb0: A1 and A2 bright state populations at a reference photon count rate (100 ms-1)
%
% "TPtime" is the transition path time for which the likelihood is calculated.
% So far developed for the calculation of only one TPtime input.
% The number of TPtime input is two but mlhTPgen3cABIcal_MT calculates only for the first value.
%
% Likelihood value is calculated for the same type of transitions.
% "stateid" is the beginning state in the transition. 1: transition 1 --> 2, 2: transition 2 --> 1
%
% "numstep" is the number of steps in the transition path, but so far only
% the one-step model has been developed.
%
% "tpE" determines the acceptor fractions and 2-color FRET efficiencies of the transition paths
% []: 3-color acceptor fractions of the TP interpolate those of State 1 and 2 with the number of steps.
%     2-color FRET efficiencies are the midway values of those of State 1 and 2
% ~[]: 3-color acceptor fractions are defined by tpE
%      2-color FRET efficiencies are also given in tpE or calculated using
%      3-color acceptor fractions.
%
% "pI" is the fractions of the TPs.
%
% Output is a likelihood value.

        % just in case
        offbndid=find(params(1:8) <= 0);params(offbndid)=0.01;
        offbndid=find(params(1:8) >= 1);params(offbndid)=0.99;

        eff1F=params(1);
        eff1U=params(2);
        eff2F=params(3);
        eff2U=params(4);
        if isempty(tpE)
            effsteps1=interp1q([0 1]',[eff1F eff1U]',(1:numstep)'/(numstep+1))';
            effsteps2=interp1q([0 1]',[eff2F eff2U]',(1:numstep)'/(numstep+1))';
        else
%             effsteps1=repmat(tpE(:,1),[1 numstep]);
%             effsteps2=repmat(tpE(:,2),[1 numstep]);
            effsteps1=tpE(:,1)';
            effsteps2=tpE(:,2)';
        end
        nIntm=length(pI);

        effs2F=params(5);
        effs2U=params(6);
        effs3F=params(7);
        effs3U=params(8);
        if isempty(tpE)
            effstepss2=repmat((effs2F+effs2U)/2,[nIntm 1]);
            effstepss3=repmat((effs3F+effs3U)/2,[nIntm 1]);
        else
            effstepss2=effs2U+(effs2F-effs2U)*(sum(tpE,2)-eff1U-eff2U)/(eff1F+eff2F-eff1U-eff2U);
            effstepss3=effs3U+(effs3F-effs3U)*(sum(tpE,2)-eff1U-eff2U)/(eff1F+eff2F-eff1U-eff2U);
            if size(tpE,2) > 2 && tpE(1,3) > 0, effstepss2=tpE(:,3); end
            if size(tpE,2) == 4, effstepss3=tpE(:,4); end
        end
        % just in case
        offbndid=find(effstepss2 <= 0);effstepss2(offbndid)=0.01;
        offbndid=find(effstepss2 >= 1);effstepss2(offbndid)=0.99;
        offbndid=find(effstepss3 <= 0);effstepss3(offbndid)=0.01;
        offbndid=find(effstepss3 >= 1);effstepss3(offbndid)=0.99;
        
        fr1=params(10);
        ratesum=params(9);

%         effb1=params(11);
%         effb2=params(12);
%         effb12=params(13);
%         ktobrt=params(14:15);
%         frnb0=params(16:17);    % for a photon count rate of 100 ms-1

        trsteprate=1./TPtime/2;       
%         if numstep==2, trsteprate=4/3./TPtime; end
%         if numstep==3, trsteprate=5/2./TPtime; end        

        effvec=[eff1F; effsteps1'; eff1U; eff2F; effsteps2'; eff2U; effs2F; effstepss2; effs2U; effs3F; effstepss3; effs3U];
        
        if stateid == 1
            pini=[1; zeros(1+nIntm,1)];
            pfin=[zeros(1+nIntm,1); 1];
        else
            pini=[zeros(1+nIntm,1); 1];
            pfin=[1; zeros(1+nIntm,1)];
        end
            
        ratemat0=zeros(2+nIntm);
        ratemat0(1:end-1,1)=[-1; pI]*2*ratesum*(1-fr1);
        ratemat0(2:end,end)=[pI; -1]*2*ratesum*fr1;
        for ii=1:nIntm
            ratemat0([1 ii+1 end],ii+1)=[1; -2; 1]*trsteprate(ii,1);
        end
        peq=[-1/ratemat0(1,1); pI./trsteprate(:,1); -1/ratemat0(end,end)];
        pinput=[peq./sum(peq) pini pfin];
            
        inputparam=repmat([effvec; params(11:end)],1,2);
        LUboundsdummy=[min(inputparam,[],2)*0.8 max(inputparam,[],2)*1.2]; % dummy LUbounds
        logmlhval=-mlhTPgen3cABIcal_MT(inputparam,LUboundsdummy,frburstdata,cumindex,indexone,cntrate,ratemat0,pinput);