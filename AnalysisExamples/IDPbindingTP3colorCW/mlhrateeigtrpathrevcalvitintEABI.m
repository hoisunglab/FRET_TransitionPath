function logmlhval=mlhrateeigtrpathrevcalvitintEABI(params,frburstdata,cumindex,indexone,stateid,numstep,varargin)
%
% logmlhval = mlhrateeigtrpathrevcalvitintEABI(params,frburstdata,cumindex,indexone,stateid,numstep,varargin)
%
% Calculate a likelihood value for transition path data with given
% parameters and transition path time.
%
% inputdata (frburstdata) is a 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% Input parameters (params) contain 8 parameters including 4 2-state transition parameters,
% and additional parameters:
%
% eff1: Apparent FRET efficiency of State 1 (e.g., folded)
% eff2: Apparent FRET efficiency of State 2 (unfolded)
% ratesum: relaxation rate of the two-state kinetic model (between states 1 and 2)
% fr1: fraction of State 1 in the two-state model
% TP time
% effb: Apparent FRET efficiency of the acceptor dark state
% ktobrt: Rate constant of a transition from acceptor dark to bright state
% frnb0: Acceptor bright state population at a reference photon count rate (100 ms-1)
%
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
        if nargin == 7,
            if varargin{1} == 2,
                effsteps=repmat((eff1+eff2)/2,numstep,1);
            elseif varargin{1} ~= 1,
                effsteps=repmat(varargin{1},numstep,1);
            end
        end
        
        fr1=params(4);
        ratesum=params(3);
        trsteprate=numstep/params(5)/2;        
        if numstep==2, trsteprate=4/3/params(5); end
        if numstep==3, trsteprate=5/2/params(5); end        

        effb=params(6);
        ktobrt=params(7);
        frnb0=params(8);        % for a photon count rate of 100 ms-1

        clear probone
        for k=1:length(indexone);
            if stateid(k)==1,
                pini0=[1; zeros(1+numstep,1)];
                pfin0=[zeros(1,1+numstep) 1];
            else
                pini0=[zeros(1+numstep,1); 1];
                pfin0=[1 zeros(1,1+numstep)];
            end

            oneburst=frburstdata(cumindex(indexone(k))+1:cumindex(indexone(k)+1),:);
            photoncolor=oneburst(:,end);
            photoninterval=diff(oneburst(:,end-2))*1e-4;        % in ms
            probonesub=diag(ones(2*(2+numstep),1));logamp=0;

            cntrate=1/mean(photoninterval);
            ktodark=ktobrt*(1-frnb0)/frnb0*cntrate/100;
%            ktodark=ktobrt*(1-frnb0)/frnb0;
            ratesumb=ktobrt+ktodark;
            frnb=ktobrt/ratesumb;

            peqb=[frnb; 1-frnb];
            pini=[pini0*peqb(1); pini0*peqb(2)];
            pfin=[pfin0 pfin0];

            effmat=diag([eff1 effsteps' eff2 repmat(effb,[1 2+numstep])]);
            ratemat00=zeros(2+numstep);
            ratemat00(1:2,1)=[-1; 1]*2*ratesum*(1-fr1);
            ratemat00(end-1:end,end)=[1; -1]*2*ratesum*fr1;
            for ii=2:length(ratemat00(:,1))-1;
                ratemat00(ii-1:ii+1,ii)=[1; -2; 1]*trsteprate;
            end
            ratemat0=zeros(2*(2+numstep));
            ratemat0(1:(2+numstep),1:(2+numstep))=ratemat00;
            ratemat0(3+numstep:end,3+numstep:end)=ratemat00;
            for ii=1:2+numstep;
                ratemat0([ii ii+2+numstep],[ii ii+2+numstep])=ratemat0([ii ii+2+numstep],[ii ii+2+numstep])+ratesumb*[-peqb(2) peqb(1);peqb(2) -peqb(1)];
            end
        
            [eigmat ratemateig]=eig(ratemat0);
            inveigmat=inv(eigmat);
            Emat(:,:,1)=inveigmat*effmat*eigmat;
            Emat(:,:,2)=inveigmat*(eye(size(effmat))-effmat)*eigmat;
            
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
        end
        
        logmlhval=-sum(probone);



