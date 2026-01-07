function [respconv exitflag]=mlhTP3cABIfit1TPindi(initparams,LUbounds,fixparam,resparamfix,oneburst,stateid,fitparamid)
%
% [resparam exitflag] = mlhTP3cABIfit1TPindi(initparams,LUbounds,fixparam,resparamfix,oneburst,stateid,fitparamid)
%
% Optimize transition path parameters of a single 3-color transition with a
% single-path, single-step TP model by maximizing likelihood function calculated by
% mlhrateeigtrpathrevcalvit3cintEABI.m.
% This analysis is for a 3-color binding experiment. IDP labeled with donor (D) 
% and acceptor 1 (A1) and binding partner is labeled with acceptor 2 (A2).
%
% inputdata (oneburst) is a 5 column matrix with burst-bin-tagged T3R data.
% Defalt unittime (time resolution) is 100 ns.
%
% Optimization parameters are selected by "fitparamid":
% eff1s1: Acceptor 1 (A1) fractions of State 1 (bound) and State 2 (unbound)
% eff2s1: Acceptor 2 (A2) fractions of State 1 and State 2
% ktobrt: Rate constants of transitions from A1 and A2 dark to bright states
% frnb0: A1 and A2 bright state populations at a reference photon count rate (100 ms-1)
% trplottime: TP time
% tpE: Acceptor fractions of the TP
%
% "fixparam" contains donor leak into A1 and A2 channel, A1 leak into A2
% channel, and predetermined DA1 and DA2 FRET efficiencies for 3-color-only analysis.
% effs2: Apparent FRET efficiency E1 for the data with donor and A1 (DA1)
% effs3: Apparent FRET efficiency E2 for the data with donor and A2 (DA2)
% tpE2c: 2-color FRET efficiencies of the TPs. This can be empty or obtained from fixparam{5}. 
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
% [A2 A1 donor] = [1 2 3]
%
% "fitparamid" defines the parameters to be optimized. The rest of
% parameters (e.g., acceptor blinking parameters) are fixed during optimization.
%
% Output is a set of optimized parameters.

    function logmlh=mlhTPfitsub(params)
        pconv=initparams;
        pconv(fitparamid)=diff(LUbounds(fitparamid,:),1,2)./(1+params.^2)+LUbounds(fitparamid,1);
        darkparams=pconv(5:8);
        ktobrt=darkparams(1:2);  % rate from dark to bright state of A1 and A2
        frnb0=darkparams(3:4);    % A1 and A2 bright population at photon count rate of 100 ms-1
        eff1s1=pconv(1:2);             % 3 color E1app
        eff2s1=pconv(3:4);    % 3 color E2app
        
        resparamone=resparamfix(end-1:end);resparamone(1)=resparamone(1)/1000;
        effb1=fixparam{1}(1);  % E1 of A1,A2 dark state (= A1/(A2 + A1 + D))
        effb2=fixparam{1}(2);  % E2 of A1,A2 dark state (= A2/(A2 + A1 + D)) (depends on the A2-protein concentration)
        effs2=fixparam{2};  % Pre-determined E_DA1 (= (A1+A2)/(D+A1+A2)) will be used unless modified below using fitting parameters
        effs3=fixparam{3};  % Pre-determined E_DA2 (= A2/(D+A1+A2)) will be used unless modified below using fitting parameters
        effb12=fixparam{4};  % E12 when no A2: A1 leak into A2 channel
        tpE2c=[];
        if length(fixparam) == 5, tpE2c=fixparam{5}; end

                trplottime=pconv(9);
                tpE=pconv(10:11)';  % Eint
                tpE=[tpE tpE2c];
        
                resparamtwo=[eff1s1; eff2s1; effs2(1); eff1s1(2)+eff2s1(2); effs3; resparamone];
                effb12=eff2s1(2)/(eff1s1(2)+eff2s1(2));       % E12 of unbound state
                logmlhfu=mlhrateeigtrpathrevcalvit3cintEABI([resparamtwo; effb1; effb2; effb12; ktobrt; frnb0],trplottime,oneburst,cumindex,indexone,stateid,1,tpE);
        
            logmlh=-logmlhfu;
       
    end

burstphotons=size(oneburst,1);
indexone=1:length(burstphotons);
cumindex=[0; cumsum(burstphotons)];

options=optimset('MaxFunEval', 20000, 'MaxIter', 3000);

invinitparams=sqrt(diff(LUbounds,1,2)./(initparams-LUbounds(:,1))-1);
[resparams fval exitflag]=fminsearch(@mlhTPfitsub,invinitparams(fitparamid),options);
respconv=initparams;
respconv(fitparamid)=diff(LUbounds(fitparamid,:),1,2)./(1+resparams.^2)+LUbounds(fitparamid,1);

end