function logmlhval=mlhrateeigtrpathrevcalvit3cintEABI(params,TPtime,frburstdata,cumindex,indexone,stateid,numstep,tpE)
%
% logmlhval = mlhrateeigtrpathrevcalvit3cintEABI(params,TPtime,frburstdata,cumindex,indexone,stateid,numstep,tpE)
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
%
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
        nIntm=1;

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

        effb1=params(11);
        effb2=params(12);
        eff12ub=params(13);
        ktobrt=params(14:15);
        frnb0=params(16:17);    % for a photon count rate of 100 ms-1

        trsteprate=1./TPtime/2;       
%         if numstep==2, trsteprate=4/3./TPtime; end
%         if numstep==3, trsteprate=5/2./TPtime; end        

        if stateid == 1
            pini0=[1; zeros(1+nIntm,1)];
            pfin0=[zeros(1+nIntm,1); 1]';
        else
            pini0=[zeros(1+nIntm,1); 1];
            pfin0=[1; zeros(1+nIntm,1)]';
        end
            
        k=1;
            oneburst=frburstdata(cumindex(indexone(k))+1:cumindex(indexone(k)+1),:);
            photoncolor=oneburst(:,end);
            photoninterval=diff(oneburst(:,end-2))*1e-4;        % in ms
            probonesub=diag(ones(4*(2+nIntm),1));logamp=0;

            cntrate=1/mean(photoninterval);
            ktodark=ktobrt.*(1-frnb0)./frnb0*cntrate/100;
            ratesumb=ktobrt+ktodark;
            frnb=ktobrt./ratesumb;

            peqb1=[frnb(1); 1-frnb(1)];
            peqb2=[frnb(2); 1-frnb(2)];
            peqtemp=[pini0*peqb2(1); pini0*peqb2(2)];
            pini=[peqtemp*peqb1(1); peqtemp*peqb1(2)];
            pfin=[pfin0 pfin0 pfin0 pfin0];
                
            nS=2+nIntm;
            ratemat00=zeros(nS);
            ratemat00(1:2,1)=[-1; 1]*2*ratesum*(1-fr1);
            ratemat00(end-1:end,end)=[1; -1]*2*ratesum*fr1;
            for ii=2:length(ratemat00(:,1))-1;
                ratemat00(ii-1:ii+1,ii)=[1; -2; 1]*trsteprate;
            end
            ratemat0=zeros(4*nS);
            % A2 blinking
            ratemat0(1:nS,1:nS)=ratemat00;
            for ii=1:nS
                ratemat0([ii ii+nS],[ii ii+nS])=ratemat0([ii ii+nS],[ii ii+nS])+ratesumb(2)*[-peqb2(2) peqb2(1); peqb2(2) -peqb2(1)];
            end
            % A1 blinking
            ratemat0(2*nS+1:end,2*nS+1:end)=ratemat0(1:2*nS,1:2*nS);
            for jj=1:2*nS
                ratemat0([jj jj+2*nS],[jj jj+2*nS])=ratemat0([jj jj+2*nS],[jj jj+2*nS])+ratesumb(1)*[-peqb1(2) peqb1(1); peqb1(2) -peqb1(1)];
            end
            
            [eigmat ratemateig]=eig(ratemat0);
            inveigmat=inv(eigmat);

        effvec=[eff1F; effsteps1'; eff1U; eff2F; effsteps2'; eff2U; effs2F; effstepss2; effs2U; effs3F; effstepss3; effs3U];
            
                Emat(:,:,1)=inveigmat*diag([effvec(nS+1:2*nS); effvec(2*nS+1:3*nS)*eff12ub; effvec(3*nS+1:4*nS); effb2*ones(nS,1)])*eigmat;    % 2&3-color A2
                Emat(:,:,2)=inveigmat*diag([effvec(1:nS); effvec(2*nS+1:3*nS)*(1-eff12ub); effb1*(1-effvec(3*nS+1:4*nS)); effb1*ones(nS,1)])*eigmat;    % 2&3-color A1
                Emat(:,:,3)=inveigmat*diag([1-effvec(1:nS)-effvec(nS+1:2*nS); 1-effvec(2*nS+1:3*nS); (1-effvec(3*nS+1:4*nS))*(1-effb1); (1-effb1-effb2)*ones(nS,1)])*eigmat;    % 2&3-color D

                for ii=1:length(photoncolor)-1;
                    ratemat=diag(exp(diag(ratemateig)*photoninterval(ii)));
                    probonesub=ratemat*Emat(:,:,photoncolor(ii))*probonesub;
                    if rem(ii,30)==0,
                        normprobone=norm(probonesub);
                        logamp=logamp+log(normprobone);
                        probonesub=probonesub/normprobone;
                    end
                end
            probone(k)=-log(pfin*eigmat*Emat(:,:,photoncolor(end))*probonesub*inveigmat*pini)-logamp;
            
        logmlhval=-sum(probone);