function logmlhval=mlhrateeigtrpathrevcalvitABIC(params,TPtime,frburstdata,cumindex,indexone,cntrate,stateid,numstep,varargin)
%
% logmlhval = mlhrateeigtrpathrevcalvitABIC(params,TPtime,frburstdata,cumindex,indexone,cntrate,stateid,numstep,varargin)
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
% "numstep" is the number of steps in the transition path, ranging from 1 to 3.
%
% The additional input argument selects the FRET efficiency of the transition path (E_TP).
% [] or 1: E_TPs interpolates eff1 and eff2 with the number of steps.
% 2: E_TP = (eff1 + eff2)/2 for all TP steps
% x (0 < x < 1): E_TP = x for all TP steps
%
% Output is a likelihood value.

        eff1=params(1);
        eff2=params(2);
        effsteps=interp1q([0 1]',[eff1 eff2]',(1:numstep)'/(numstep+1));
        if nargin == 9
            if varargin{1} == 2
                effsteps=repmat((eff1+eff2)/2,numstep,1);
            elseif varargin{1} ~= 1
                effsteps=repmat(varargin{1},numstep,1);
            end
        end
        
        fr1=params(4);
        ratesum=params(3);
        trsteprate=numstep./TPtime/2;
        if numstep==2, trsteprate=4/3./TPtime; end
        if numstep==3, trsteprate=5/2./TPtime; end        

        effb=params(5);
        ktobrt=params(6);
        frnb0=params(7);        % for a photon count rate of 100 ms-1

        ratesumF=[repmat(2*ratesum*(1-fr1),1,length(TPtime)); repmat(trsteprate,numstep,1)];
        ratesumB=[repmat(trsteprate,numstep,1); repmat(2*ratesum*fr1,1,length(TPtime))]; 

        ratesumNew=ratesumF+ratesumB;
        peqtemp=cumprod([ones(1,size(ratesumF,2)); ratesumF./ratesumB],1);
        peqNew=zeros(size(peqtemp));
        for jj=1:size(peqNew,2)
            peqNew(:,jj)=peqtemp(:,jj)/sum(peqtemp(:,jj));
        end
        frnNew=zeros(size(peqNew));
        frnNew(1,:)=peqNew(1,:);
        for jj=1:numstep
            fdiv=prod(1-frnNew(1:jj,:),1);
            frnNew(jj+1,:)=peqNew(jj+1,:)./fdiv;
        end
        frnNew(end,:)=[];
        
        inputparam=[repmat([eff1; effsteps; eff2; effb],1,length(TPtime)); ratesumNew; frnNew; repmat([ktobrt; frnb0],1,length(TPtime))];
        LUboundsdummy=[min(inputparam,[],2)*0.8 max(inputparam,[],2)*1.2]; % dummy LUbounds
        logmlhval=-mlhTPABIcal_MT(inputparam,LUboundsdummy,frburstdata,cumindex,indexone,cntrate,stateid);