clear all

accchan=1;  % acceptor channel
donchan=2;  % donor channel

E_fold = 0.70; % Simulation FRET efficiency of the folded state
E_unfold = 0.41; % Simulation FRET efficiency of the unfolded state
effb=0.0583;    % FRET efficiency of acceptor dark state (i.e., donor-only)
tTP = 0.0045;        % Simulation TP time in us
effTP=0.58;         % Simulation FRET efficiency of the TP

% Load simulation data
% t3rdata: photon trajectories
% stateiddata: segmentation information
% param2st: 2-state folding kinetic parameters
% effb: FRET efficiency of acceptor dark state (i.e., donor-only)
% tTP: Simulation TP time in us
% effTP: Simulation FRET efficiency of the TP
% mincntrate: Minimum photon countrates in the two states used for trajectory selection
% minN: Mimimum number of photons contained in the two states used for trajectory selection
load ProteinFoldingTP2colorExample.mat

%% *** Likelihood analysis, instantaneous transition ************
confplotlevel=0.95;     % 95% confidence
resparamfix=param2st;

drvintv=0.002;
drvpnt=5;
initparams=[0.7 0.42 500 0.99]';        % [EF, EU, blinking rate (kb), pb (bright state population)], for mlhTPABIfit0
LUbounds=[0.6 1; 0.3 0.5; 50 3000; 0.1 0.99999];
isincE=1;   % 2 for constant E_int = (E_F + E_U)/2
fitfu=3;    % combine folding and unfolding transitions

clear mlhfitparam errorparam logmlh
nintm=1;
resparams=mlhTPABIsimfit0C(initparams,LUbounds,effb,resparamfix,t3rdata,stateiddata,fitfu,nintm,isincE);
mlhfitparam=resparams
[errorparam logmlh]=mlhTPABIsimfit0errorC(mlhfitparam,effb,resparamfix,t3rdata,stateiddata,fitfu,nintm,isincE,drvintv,drvpnt);
save(['TPfit0param.mat'],'mlhfitparam','errorparam','logmlh','-mat');

%% Compute likelihood values for various tTP and E_TP using fitted parameters
nintm=1;
fitfu=3;
confplotlevel=0.95;     % 95% confidence
trplottime0=[1.6 2.5 4 6 7 8 10];trplottime=[1e-5*trplottime0 1e-4*trplottime0 1e-3*trplottime0 0.01*trplottime0]; % in ms
intE=[0.1:0.1:0.4 0.5 0.6 0.65 0.7 0.8];           % *************** intE: list of intermediate FRET efficiencies
resparamfix=param2st;
load(['TPfit0param.mat']);

clear logmlhvalintE
for uu=1:length(intE)
    uu
    logmlhres=mlhTPABIsimcalC(mlhfitparam,trplottime,effb,resparamfix,t3rdata,stateiddata,nintm,intE(uu));
    logmlhvalintE{1}(:,uu)=logmlhres{1};
    logmlhvalintE{2}(:,uu)=logmlhres{2};
end
save(['TPintE.mat'],'logmlhvalintE','-mat');

load(['TPintE.mat']);
figure;
for uu=1:length(intE)
    subplot(3,4,uu);plot(trplottime,logmlhvalintE{1}(:,uu)-logmlhvalintE{1}(1,uu)+logmlhvalintE{2}(:,uu)-logmlhvalintE{2}(1,uu),'o-','markersize',6);set(gca,'xscale','log','xtick',10.^(-4:-1),'ticklength',[0.02 0.02]);
    ylim([-15 15]);xlim([1e-5 3e-2]);axis square
    line([1e-5 3e-2],-[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
    line([1e-5 3e-2],[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
    title(['E = ' num2str(intE(uu))]);        
end

%% Fit tptime together with other parameters
resparamfix=param2st;
drvintv=0.005;
drvpnt=5;
nintm=1;
isincE=1; % 2: constant intermediate E
fitfu=3;
initparams=[0.7 0.42 500 0.95 0.001 0.55]';        % [EF, EU, blinking rate (kb), pb (bright state population), tpt in ms, Eint]
LUbounds=[0.6 0.8; 0.3 0.5; 50 3000; 0.5 0.99999; 1e-4 0.1; 0.45 0.65];

clear mlhfitparam errorparam logmlh dellogmlh;
for nintm=1
    for fitfu=3
        resparams=mlhTPABIsimfitEintNintmEC(initparams,LUbounds,effb,resparamfix,t3rdata,stateiddata,fitfu,nintm);
        mlhfitparam(:,fitfu,nintm)=resparams;
        resparams
        [errorparam(:,fitfu,nintm) logmlh(fitfu,nintm) dellogmlh(fitfu,nintm)]=mlhTPABIsimfitEintNintmEerrorC(mlhfitparam(:,fitfu,nintm),effb,resparamfix,t3rdata,stateiddata,fitfu,nintm,drvintv,drvpnt);
    end
end
save(['TPfitparam.mat'],'mlhfitparam','errorparam','logmlh','-mat');

%% Fit int E with TPT for more than one intermediate (note that this is not consistent with the simulated model of the example)
resparamfix=param2st;
nintm=1;    % not used in this fitting
drvintv=0.005;
drvpnt=5;
fitfu=3;
% % 1 pathway
% jj=1;
% initparams=[0.7 0.42 500 0.95 0.001 0.55]';        % [EF, EU, blinking rate (kb), pb (bright state population), TP time, E_intm]
% LUbounds=[0.6 0.8; 0.3 0.5; 50 3000; 0.5 0.99999; 1e-4 0.1; 0.45 0.65];
% 2 pathways
jj=2;
initparams=[0.7 0.42 500 0.95 ...      % [EF, EU, blinking rate (kb), pb (bright state population),
        0.001 0.005 0.6 0.55 0.5]';        % [TP time, E_intm, fraction pathways f], p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
LUbounds=[0.6 0.8; 0.3 0.5; 50 3000; 0.5 0.99999; ...
        1e-4 0.01; 5e-4 0.01; 0.45 0.68; 0.45 0.68; 0.01 0.99];
% % 3 pathways
% jj=3;
% initparams=[0.65 0.28 200 0.95 ...      % [EF, EU, blinking rate (kb), pb (bright state population),
%         0.1 0.6 0.5 0.3 0.5]';        % [TP time, E_intm, fraction pathways f], p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
% LUbounds=[0.5 1; 0.2 0.5; 100 2000; 0.8 1; ...
%         0.05 0.5; 0.2 0.9; 0.32 0.5; 0.4 0.58; 0.1 0.9];

clear mlhfitparam errorparam logmlh dellogmlh bic;
resparams=mlhTPABIsimfitEintParC(initparams,LUbounds,effb,resparamfix,t3rdata,stateiddata,fitfu,nintm);
nparI=(length(initparams)-3)/3;  % number of parallel pathways
pfactor=cumprod([1; 1-resparams(4+(2*nparI+1:3*nparI-2))]);
mlhfitparam(:,jj)=[resparams(1:4+2*nparI); resparams(5+2*nparI:end).*pfactor]; % Convert f to population p
mlhfitparam(:,jj)
[errorparam(:,jj) logmlh(jj) dellogmlh(jj) bic(jj)]=mlhTPABIsimfitEintParerrorC(mlhfitparam(:,jj),effb,resparamfix,t3rdata,stateiddata,fitfu,nintm,drvintv,drvpnt);

save(['TPfitparamParallelfit.mat'],'mlhfitparam','errorparam','logmlh','dellogmlh','bic','-mat');
