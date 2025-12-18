% So far written for calculation of only one step transition.
function logmlhval=mlhrateeigtrpathrevcalvit3cAlexABIparC(params,TPtime,frburstdata,cumindex,indexone,cntrate,stateid,numstep,tpE,pI)

        offbndid=find(params(1:10) <= 0);params(offbndid)=0.01;
        offbndid=find(params(1:10) >= 1);params(offbndid)=0.99;

        eff1F=params(1);
        eff1U=params(2);
        eff2F=params(3);
        eff2U=params(4);
        eff12F=params(5);
        eff12U=params(6);
        if isempty(tpE)
            effsteps1=interp1q([0 1]',[eff1F eff1U]',(1:numstep)'/(numstep+1))';
            effsteps2=interp1q([0 1]',[eff2F eff2U]',(1:numstep)'/(numstep+1))';
            effsteps12=interp1q([0 1]',[eff12F eff12U]',(1:numstep)'/(numstep+1))';
        else
%             effsteps1=repmat(tpE(:,1),[1 numstep]);
%             effsteps2=repmat(tpE(:,2),[1 numstep]);
            effsteps1=tpE(:,1)';
            effsteps2=tpE(:,2)';
            effsteps12=tpE(:,3)';
        end
        nIntm=length(pI);

        effs2F=params(7);
        effs2U=params(8);
        effs3F=params(9);
        effs3U=params(10);
        if isempty(tpE)
            effstepss2=repmat((effs2F+effs2U)/2,[nIntm 1]);
            effstepss3=repmat((effs3F+effs3U)/2,[nIntm 1]);
        else
            effstepss2=effs2U+(effs2F-effs2U)*(sum(tpE,2)-eff1U-eff2U-eff12U)/(eff1F+eff2F+eff12F-eff1U-eff2U-eff12U);
            effstepss3=effs3U+(effs3F-effs3U)*(sum(tpE,2)-eff1U-eff2U-eff12U)/(eff1F+eff2F+eff12F-eff1U-eff2U-eff12U);
            if size(tpE,2) > 3 && tpE(1,4) > 0, effstepss2=tpE(:,4); end
            if size(tpE,2) == 5, effstepss3=tpE(:,5); end
        end
        offbndid=find(effstepss2 <= 0);effstepss2(offbndid)=0.01;
        offbndid=find(effstepss2 >= 1);effstepss2(offbndid)=0.99;
        offbndid=find(effstepss3 <= 0);effstepss3(offbndid)=0.01;
        offbndid=find(effstepss3 >= 1);effstepss3(offbndid)=0.99;
        
        fr1=params(12);
        ratesum=params(11);

%         effb1=params(11);
%         effb2=params(12);
%         effb12=params(13);
%         ktobrt=params(14:15);
%         frnb0=params(16:17);    % for a photon count rate of 100 ms-1

        trsteprate=1./TPtime/2;       
%         if numstep==2, trsteprate=4/3./TPtime; end
%         if numstep==3, trsteprate=5/2./TPtime; end        

        effvec=[eff1F; effsteps1'; eff1U; eff2F; effsteps2'; eff2U; eff12F; effsteps12'; eff12U; effs2F; effstepss2; effs2U; effs3F; effstepss3; effs3U];
        
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
            
%        inputparam=repmat([effvec; params(7:end)],1,size(TPtime,2));
        inputparam=repmat([effvec; params(13:end)],1,2);
        LUboundsdummy=[min(inputparam,[],2)*0.8 max(inputparam,[],2)*1.2]; % dummy LUbounds
        logmlhval=-mlhTPgen3cAlexABIcal_MT(inputparam,LUboundsdummy,frburstdata,cumindex,indexone,cntrate,ratemat0,pinput);