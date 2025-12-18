function logmlhval=mlhrateeigtrpathrevcalvitABIparC(params,TPtime,frburstdata,cumindex,indexone,cntrate,stateid,numstep,varargin)
%
% logmlhval = mlhrateeigtrpathrevcalvitABIparC(params,TPtime,frburstdata,cumindex,indexone,cntrate,stateid,numstep,varargin)
%
% Calculate a likelihood value for transition path data with given
% parameters and transition path time.
%
% inputdata (frburstdata) is a 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% Input parameters (params) contain 7 parameters including 4 2-state transition parameters,
% and additional parameters:
%
% eff1: Apparent FRET efficiency of State 1 (e.g., folded)
% eff2: Apparent FRET efficiency of State 2 (unfolded)
% ratesum: relaxation rate of the two-state kinetic model (between states 1 and 2)
% fr1: fraction of State 1 in the two-state model
% effb: Apparent FRET efficiency of the acceptor dark state
% ktobrt: Rate constant of a transition from acceptor dark to bright state
% frnb0: Acceptor bright state population at a reference photon count rate (100 ms-1)
%
% "TPtime" is the transition path time for which the likelihood is calculated.
% So far developed for the calculation of only one TPtime input.
% The number of TPtime input is two but mlhTPgenABIcal_MT calculates only for the first value.
%
% Likelihood value is calculated for the same type of transitions.
% "stateid" is the beginning state in the transition. 1: transition 1 --> 2, 2: transition 2 --> 1
%
% "numstep" is the number of steps in the transition path, but so far only
% the one-step model has been developed.
%
% The additional input arguments (must be included) are the FRET efficiencies of the transition paths (E_TP)
% and the fractions of the TPs, pI.
%
% Output is a likelihood value.

        eff1=params(1);
        eff2=params(2);
%        effsteps=interp1q([0 1]',[eff1 eff2]',(1:numstep)'/(numstep+1));
        if nargin == 10
            Eintvec=varargin{1};
            pI=varargin{2};
        end
        nIntm=length(Eintvec);
        
        fr1=params(4);
        ratesum=params(3);
        trsteprate=1./TPtime/2;       
%        trsteprate=numstep./TPtime/2;
%        if numstep==2, trsteprate=4/3./TPtime; end
%        if numstep==3, trsteprate=5/2./TPtime; end        

        effb=params(5);
        ktobrt=params(6);
        frnb0=params(7);        % for a photon count rate of 100 ms-1

        effvec=[eff1; Eintvec; eff2; effb];
        
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
        for ii=1:nIntm;
            ratemat0([1 ii+1 end],ii+1)=[1; -2; 1]*trsteprate(ii,1);
        end
        
        peq=[-1/ratemat0(1,1); pI./trsteprate(:,1); -1/ratemat0(end,end)];
        pinput=[peq./sum(peq) pini pfin];
            
        inputparam=repmat([effvec; ktobrt; frnb0],1,size(TPtime,2));
        LUboundsdummy=[min(inputparam,[],2)*0.8 max(inputparam,[],2)*1.2]; % dummy LUbounds
        logmlhval=-mlhTPgenABIcal_MT(inputparam,LUboundsdummy,frburstdata,cumindex,indexone,cntrate,ratemat0,pinput);