function [respconv exitflag]=mlhTPABIfit1TPindi(initparams,LUbounds,effb,resparamfix,oneburst,stateid,fitparamid,trcc)
%
% [resparam exitflag] = mlhTPABIfit1TPindi(initparams,LUbounds,effb,resparamfix,oneburst,stateid,fitparamid,trcc)
%
% Optimize transition path parameters of a single 2-color transition with a
% single-path, single-step TP model by maximizing likelihood function calculated by
% mlhrateeigtrpathrevcalvitintEABI.m.
%
% inputdata (onet3rkk) is a 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% Optimization parameters are selected by "fitparamid":
% fretEF: Apparent FRET efficiency of State 1 (e.g., folded)
% fretEU: Apparent FRET efficiency of State 2 (unfolded)
% ktobrt: Rate constant of a transition from acceptor dark to bright state
% frnb0: Acceptor bright state population at a reference photon count rate (100 ms-1)
% trplottime: TP time
% Eint: FRET efficiency of the TP
%
% "effb" is the FRET efficiency of the acceptor dark state.
%
% "resparamfix" has 4 parameters State 1 and State 2 FRET efficiencies,
% which are not used (optimization parameters), relaxation rate, which is
% not so important unless too high because it is scaled by a factor of
% 1000, and fraction of state 1, which does not affect the optimization.
%
% "stateid" is the beginning state.
%
% Note that the photon detection channels in the 5th column of
% the photon trajectory data should be pre-converted.
% [Acceptor donor] = [1 2]
%
% "fitparamid" defines the parameters to be optimized. The rest of
% parameters (e.g., acceptor blinking parameters) are fixed during optimization.
%
% "trcc" is transition type.
% 1: 3-color (not optimized here)
% 2: Donor and acceptor 1
% 3: Donor and acceptor 2 (effb is the same as unbound state E in binding experiment)
%
% Output is a set of optimized parameters.

    function logmlh=mlhTPfitsub(params)
        pconv=initparams;
        pconv(fitparamid)=diff(LUbounds(fitparamid,:),1,2)./(1+params.^2)+LUbounds(fitparamid,1);
        fretEF=pconv(1);
        fretEU=pconv(2);
        ktobrt=pconv(3);
        frnb0=pconv(4);        % for a photon count rate of 100 ms-1
        trplottime=pconv(5);
        Eint=pconv(6);

        resparamone=resparamfix(end-1:end);resparamone(1)=resparamone(1)/1000;
        resparamone=[fretEF; fretEU; resparamone];
        
        if trcc == 3, effb=fretEU; end      % For DA2 segments

%        logmlhfu=mlhrateeigtrpathrevcalvitABIC([resparamone; effb; ktobrt; frnb0],[trplottime trplottime],oneburst,cumindex,indexone,cntrate,stateid,nintm,Eint);
        logmlhfu=mlhrateeigtrpathrevcalvitintEABI([resparamone; trplottime; effb; ktobrt; frnb0],oneburst,cumindex,indexone,stateid,1,Eint);

        logmlh=-logmlhfu;
    end

burstphotons=size(oneburst,1);
indexone=1:length(burstphotons);
cumindex=[0; cumsum(burstphotons)];

options=optimset('MaxFunEval', 100000, 'MaxIter', 3000);

invinitparams=sqrt(diff(LUbounds,1,2)./(initparams-LUbounds(:,1))-1);
[resparams fval exitflag]=fminsearch(@mlhTPfitsub,invinitparams(fitparamid),options);
respconv=initparams;
respconv(fitparamid)=diff(LUbounds(fitparamid,:),1,2)./(1+resparams.^2)+LUbounds(fitparamid,1);

end