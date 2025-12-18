clear all

chans=[1 3 4];  % [Acceptor 2, Acceptor 1, Donor]
numchans=length(chans);
resparamfix=[0.0649; 0.3472; 0.7195; 0.1010; 0.7347; 0.6188; 0.0301; 0.3634; 0.2429; 0.0846; 0.3916; 0.3593]; % from 3-color 2-state ABI analysis of the mid intensity data
% [eff1B, eff1U, eff2B, eff2U, E1B_DA1, E2B_DA2, E2U_DA1, element 8 - 10 not used, relaxation rate k, bound fraction]

% experimental parameters for recoloring with a double-path model
load TAD_A488U_A594M_NCBD_CF680R_30nM0.2msTP3cABIfitparamParallelfitC_2TP_hs90All.mat
mlhfitparamPar1=mlhfitparamPar;
clear mlhfitparamPar

npath=2;
effD=[0.05 mlhfitparamPar1{2}(8,6)];
resparamfix(1:7)=mlhfitparamPar1{2}([1:5 7 8],6);
Eint=[mlhfitparamPar1{2}(15:16,6) mlhfitparamPar1{2}(17:18,6)];     % intermediate FRET efficiency ([eps1 eps2; ...])
numintstates=1;
trpathrate=1./mlhfitparamPar1{2}(13:14,6); % [parhway1; pathway2; ...]
ppath=mlhfitparamPar1{2}(23,6);ppath(2)=1-ppath(1);     % Fraction of pathways
Eint2c=[mlhfitparamPar1{2}(19:20,6) mlhfitparamPar1{2}(21:22,6)];     % E1 (DA1) and E2 (DA2) of the TP

effb12=resparamfix(4)/sum(resparamfix([2 4]));  % A1 leak into A2 channel
eFUd1=[effD(1)*(1-resparamfix(6)) resparamfix(6); effD(1)*(1-resparamfix(7)) resparamfix(7)];   % A1 dark
eFUd2=[resparamfix(5)*(1-effb12) resparamfix(5)*effb12; resparamfix(2) resparamfix(4)];   % A2 dark
ed3=[effD(1) effD(2)];   % A1 & A2 dark

Eintd1=[effD(1)*(1-Eint2c(:,2)) Eint2c(:,2)];
Eintd2=[Eint2c(:,1)*(1-effb12) Eint2c(:,1)*effb12];

confplotlevel=0.95;

drvintv=0.005;
drvpnt=5;

trplottime0=[1.6 2.5 4 6.5 10];trplottime=[1e-5*trplottime0(1) 1e-3*trplottime0 0.01*trplottime0 0.1*trplottime0]; % in microsecond
load IDPbindingTP2path3colorCWExample.mat

%% *** Likelihood analysis, instantaneous transition ************
drvintv=0.002;
drvpnt=5;
initparams=[0.08 0.35 ...   % eps1_3c
            0.7 0.11 ...    % epa2_3c
            0.73 0.45 ...   % E1_DA1
            0.6 0.04 ...    % E2_DA2
            100 100 0.95 0.95]';        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2] for mlhTP3cABIfit0
LUbounds=[0.03 0.2; 0.25 0.4; ...
        0.6 0.8; 0.05 0.2;
        0.65 0.85; 0.4 0.55;
        0.4 0.8; 0.001 0.2; 
        20 2000; 20 2000; 0.9 0.9999; 0.9 0.9999];
fixparam={effD', []', []', 0.2};   % {[effb1 effb2], effs2, effs3, effb12}
                                    % [effb1 effb2 effb12] are donor and accetor leak into other channels (effb1 = A1/(A1+D), effb2 = A2/(A2+D), effb12 = A2/(A1+A2)
                                    % effs2 and effs3 are FRET efficiencies in state 2 (A2 dark, E = (A1+A2)/(D+A1+A2)) and 3 (A1 dark, E = A2/(D+A2)). These values are fitting parameters in some conditions
nintm=1;
isincE=1;   % 2 for constant E_int = (E_F + E_U)/2
%fitidin=1:6;  % 1: 3-color segments, 2: DA1 segments, 3: DA2 segments, 4: 3-color + DA1, 5: DA2 segments, 6: 3c + DA1 + DA2
fitidin=[2 3 6];
fitfu=3;

resparams=mlhTP3cABIsimfit0C(initparams,LUbounds,fixparam,resparamfix,t3rdata,stateiddata,chans,fitfu,nintm,isincE,fitidin)
mlhfitparam=resparams;
[errorparams logmlhs]=mlhTP3cABIsimfit0errorC(mlhfitparam,fixparam,resparamfix,t3rdata,stateiddata,chans,fitfu,nintm,isincE,fitidin,drvintv,drvpnt);
errorparam=errorparams;
logmlh=logmlhs;
save(['TPfit0param3c.mat'],'mlhfitparam','errorparam','logmlh','-mat');

%% Compute likelihood values using fitted parameters, DA1 and DA2 intE
fitidin=[2 3]; % 1: 3-color segments, 2: DA1 segments, 3: DA2 segments, 4: 3-color + DA1, 5: DA2 segments, 6: 3c + DA1 + DA2
intE1E2{1}=[0.06 0.1:0.1:0.5 0.55:0.05:0.8];
intE1E2{2}=[0.05 0.1:0.1:0.8];
fitfu=3;
load(['TPfit0param3c.mat']);
clear logmlhvalintE
for ff=1:length(fitidin)
    for uu=1:length(intE1E2{ff})
        logmlhres=mlhTP3cABIsimcalC(mlhfitparam,trplottime,fixparam,resparamfix,t3rdata,stateiddata,chans,1,intE1E2{ff}(uu),fitidin(ff));
        logmlhvalintE{ff}{1}(:,uu)=logmlhres{1}(:,fitidin(ff));
        logmlhvalintE{ff}{2}(:,uu)=logmlhres{2}(:,fitidin(ff));
    end
end
save(['TPintE2c.mat'],'logmlhvalintE','-mat');

load(['TPintE2c.mat']);
for ff=1:length(fitidin)
    figure;
    for uu=1:length(intE1E2{ff})
        subplot(3,4,uu);plot(trplottime,logmlhvalintE{ff}{1}(:,uu)-logmlhvalintE{ff}{1}(1,uu)+logmlhvalintE{ff}{2}(:,uu)-logmlhvalintE{ff}{2}(1,uu),'o-','markersize',6);set(gca,'xscale','log','xtick',10.^(-4:-1),'ticklength',[0.02 0.02]);
        ylim([-20 60]);xlim([1e-4 1]);axis square
        line([1e-5 1e0],-[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
        line([1e-5 1e0],[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
        title(['E = ' num2str(intE1E2{ff}(uu))]);        
    end
end

%% Compute likelihood values using fitted parameters, 3-color intE1 intE2
intE3c={[0.025:0.025:0.5],[0.05:0.05:0.8]};    % {E1, E2}
nintm=1;
isincE=1;   % 2 for constant E_int = (E_F + E_U)/2
fitidin=6;  % 1: 3-color segments, 2: DA1 segments, 3: DA2 segments, 4: 3-color + DA1, 5: DA2 segments, 6: 3c + DA1 + DA2
fitfu=3;
load(['TPfit0param3c.mat']);
clear logmlhkk
logmlhres=mlhTP3cABIsimintEcalC(mlhfitparam,trplottime,fixparam,resparamfix,t3rdata,stateiddata,chans,nintm,intE3c,fitidin);
for jj=1:length(trplottime)
    logmlhres3c(:,:,jj)=logmlhres{1}(:,:,jj)-logmlhres{1}(:,:,1)+logmlhres{2}(:,:,jj)-logmlhres{2}(:,:,1);
end
save(['TP3cintEFitid6.mat'],'logmlhres3c','-mat');

load(['TP3cintEFitid6.mat']);
tptoffset=6;
clevel=-1:0.2:1;
intEplotid=[3 13; 15 3];
    logmlhone=logmlhres3c;
    maxmlh=max([10 round(max(max(max(logmlhone))))]);
    figure;
    for jj=1:length(trplottime(tptoffset:end))
        subplot(3,4,jj);
        contourf(intE3c{2},intE3c{1},logmlhone(:,:,jj+tptoffset-1),maxmlh*clevel);
        axis equal;axis([0 intE3c{2}(end) 0 intE3c{1}(end)]);
        title(sprintf('t_{TC} = %.1f \\mus',trplottime(tptoffset+jj-1)*1000));
        colormap(cmap2d(8));
        caxis(maxmlh*[-1 1]);
    end
    subplot(3,4,12);    % dummy plot for colorbar
    contourf(intE3c{2},intE3c{1},logmlhone(:,:,jj+tptoffset-1),maxmlh*clevel);
    axis equal;axis([0 intE3c{2}(end) 0 intE3c{1}(end)]);
    title('Dummy for colorbar plot');xlabel('E_2');ylabel('E_1');
    colormap(cmap2d(8));
    caxis(maxmlh*[-1 1]);
    colorbar('Ticks',maxmlh*(-1:0.5:1));
figure;
    logmlhone=logmlhres3c;
    maxmlh=max([10 round(max(max(max(logmlhone))))]);
    subplot(2,3,1);plot(trplottime,permute(logmlhone(intEplotid(1,1),intEplotid(1,2),:),[3 1 2]),'bo-','markersize',8);hold on
    plot(trplottime,permute(logmlhone(intEplotid(2,1),intEplotid(2,2),:),[3 1 2]),'ro-','markersize',8);hold off
    set(gca,'xscale','log','xtick',10.^(-4:-1),'ticklength',[0.02 0.02]);
    ylim([-1 1]*maxmlh*1.05);xlim([1e-4 1]);axis square
    line([1e-5 1e0],-[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
    line([1e-5 1e0],[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
    title(sprintf('[E1 = %.2f, E2 = %.2f], [E1 = %.2f, E2 = %.2f]',intE3c{1}(intEplotid(1,1)),intE3c{2}(intEplotid(1,2)),intE3c{1}(intEplotid(2,1)),intE3c{2}(intEplotid(2,2))));

%% Compute likelihood values of individual transitions using fitted parameters, 3-color intE1 intE2
tptplotid=[1 13 14]; 
intE3cTPTindi=[intE3c{1}([1; intEplotid(:,1)])' intE3c{2}([1; intEplotid(:,2)])' trplottime(tptplotid)'];    % [E1, E2, tpt], the first parameters for instantaneous transition
fitidin=1; % 1: 3-color segments, 2: DA1 segments, 3: DA2 segments, 4: 3-color + DA1, 5: DA2 segments, 6: 3c + DA1 + DA2
clear logmlhres3cindi
logmlhres3cindi=mlhTP3cABIsimintEcalindiC(mlhfitparam,intE3cTPTindi,fixparam,resparamfix,t3rdata,stateiddata,chans,1,fitidin);
save(['TP3cintEindi_Fitid1.mat'],'logmlhres3cindi','-mat');

load(['TP3cintEindi_Fitid1.mat']);
figure;
logmlhindi=logmlhres3cindi{1}(:,2:end)-repmat(logmlhres3cindi{1}(:,1),[1 size(intE3cTPTindi,1)-1]);
logmlhindi=[logmlhindi; logmlhres3cindi{2}(:,2:end)-repmat(logmlhres3cindi{2}(:,1),[1 size(intE3cTPTindi,1)-1])];

subplot(2,3,1);plot(logmlhindi(:,1),logmlhindi(:,2),'o');xlabel(num2str(intE3cTPTindi(2,:)));ylabel(num2str(intE3cTPTindi(3,:)));axis([-5 5 -5 5]);line([0 0],[-5 5],'color','r');line([-5 5],[0 0],'color','r');axis square
text(3,3,num2str(round(length(find(logmlhindi(:,1) > 0 & logmlhindi(:,2) > 0))/size(logmlhindi,1)*100)));
text(3,-3,num2str(round(length(find(logmlhindi(:,1) > 0 & logmlhindi(:,2) < 0))/size(logmlhindi,1)*100)));
text(-3,3,num2str(round(length(find(logmlhindi(:,1) < 0 & logmlhindi(:,2) > 0))/size(logmlhindi,1)*100)));
text(-3,-3,num2str(round(length(find(logmlhindi(:,1) < 0 & logmlhindi(:,2) < 0))/size(logmlhindi,1)*100)));

%% Fit int E with TPT for more than one intermediate
%isincE=1; % 2: constant intermediate E
resparamfix=[0.0748; 0.337; 0.688; 0.101; 0.715; 0.592; 0.0303; 0.316; 0.226; 0.09628; 0.552; 0.368]; % from 3-color 2-state analysis of the mid intensity data
drvintv=0.002;
drvpnt=5;
nintm=1;    % not used in this fitting
isincE=1;   % not used in this fitting, 2 for constant E_int = (E_F + E_U)/2
clear mlhfitparamPar errorparamPar logmlh dellogmlh bic;
fitfu=3;
jj=1;
  switch jj
    case 1
    % 1 pathway
%    fitidin=1:6;
    fitidin=6; % 1: 3-color segments, 2: DA1 segments, 3: DA2 segments, 4: 3-color + DA1, 5: DA2 segments, 6: 3c + DA1 + DA2
    initparams=[0.08 0.35 ...   % E1_3c
            0.75 0.11 ...    % E2_3c
            0.75 0.45 ...   % E1_DA1
            0.6 0.04 ...    % E2_DA2
            100 100 0.95 0.95 ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2]
            0.5 ...         % TP time
            0.1 0.65 0.6 0.55]';     % E1TP E2TP E1TP_DA1 E1TP_DA2
    LUbounds=[0.03 0.2; 0.25 0.4; ...
        0.6 0.9; 0.05 0.2;
        0.6 1; 0.4 0.6;
        0.4 0.8; 0.001 0.2;
        20 2000; 20 2000; 0.8 0.9999; 0.9 0.9999;
        0.1 1.5;
        0.01 0.25; 0.2 0.75; 0.45 0.75; 0.4 0.65];
%    fixparam={[0.05 0.0473]', [0.7664 0.4614]', [0.7203 0.0473]', 0.2};   % {[effb1 effb2], effs2, effs3, effb12}
    fixparam={effD', []', []', 0.2};   % {[effb1 effb2], effs2, effs3, effb12}
                                    % [effb1 effb2 effb12] are donor and accetor leak into other channels (effb1 = A1/(A1+D), effb2 = A2/(A2+D), effb12 = A2/(A1+A2)
                                    % effs2 and effs3 are FRET efficiencies in state 2 (A2 dark, E = (A1+A2)/(D+A1+A2)) and 3 (A1 dark, E = A2/(D+A2)). These values are fitting parameters in some conditions
    fixparam1=fixparam;
%    fixparam1{5}=[0.5638 0.5620];
    fixparam1{5}=[];

      case 2
    % 2 pathways
%    fitidin=[1 3:6];
    fitidin=6; % 1: 3-color segments, 2: DA1 segments, 3: DA2 segments, 4: 3-color + DA1, 5: DA2 segments, 6: 3c + DA1 + DA2
    initparams=[0.07 0.33 ...   % E1_3c
            0.75 0.10 ...    % E2_3c
            0.75 0.45 ...   % E1_DA1
            0.6 0.04 ...    % E2_DA2
            100 100 0.9728 0.9881 ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2]
            0.1 0.5 ...         % TP time
            0.15 0.08 0.45 0.65 0.55 0.55 0.4 0.55 ...     % E1TP1 E1TP2 E2TP1 E2TP2 E1TP1_DA1 E1TP2_DA1 E2TP1_DA2 E2TP2_DA2
            0.6]';     % fraction pathways f, p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
    LUbounds=[0.03 0.1; 0.25 0.4; ...
        0.6 0.9; 0.05 0.2;
        0.6 0.9; 0.4 0.55;
        0.4 0.8; 0.001 0.2;
        20 2000; 20 2000; 0.9 0.9999; 0.9 0.9999;
        0.01 1.5; 0.05 1.5;
        0.01 0.25; 0.01 0.15; 0.3 0.7; 0.3 0.8; 0.35 0.75; 0.35 0.75; 0.2 0.55; 0.4 0.65;
        0.2 0.8];
%    fixparam={[0.05 0.0473]', [0.7664 0.4614]', [0.7161 0.0473]', 0.2};   % {[effb1 effb2], effs2, effs3, effb12}
    fixparam={effD', []', []', 0.2};   % {[effb1 effb2], effs2, effs3, effb12}
                                    % [effb1 effb2 effb12] are donor and accetor leak into other channels (effb1 = A1/(A1+D), effb2 = A2/(A2+D), effb12 = A2/(A1+A2)
                                    % effs2 and effs3 are FRET efficiencies in state 2 (A2 dark, E = (A1+A2)/(D+A1+A2)) and 3 (A1 dark, E = A2/(D+A2)). These values are fitting parameters in some conditions
    fixparam1=fixparam;
%    fixparam1{5}=[0.5638 0.5052; 0.5638 0.5584];
    fixparam1{5}=[];

      case 3
    % 3 pathways
    fitidin=1;
    initparams=[0.07 0.33 ...   % E1_3c
            0.8 0.11 ...    % E2_3c
            0.75 0.45 ...   % E1_DA1
            0.6 0.04 ...    % E2_DA2
            100 100 0.98 0.98 ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2]
            0.5 0.5 0.2 ...         % TP time
            0.37 0.05 0.3 0.15 0.65 0.2 0.54 0.56 0.55 0.1 0.35 0.55 ...     % E1TP1 E1TP2 E2TP1 E2TP2 E1TP1_DA1 E1TP2_DA1 E2TP1_DA2 E2TP2_DA2
            0.4 0.8]';     % fraction pathways f, p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
    LUbounds=[0.03 0.1; 0.25 0.4; ...
        0.7 0.85; 0.05 0.2;
        0.65 0.85; 0.4 0.55;
        0.4 0.8; 0.001 0.2;
        20 2000; 20 2000; 0.9 0.9999; 0.9 0.9999;
        0.05 2.5; 0.05 1; 0.01 1.2;
        0.3 0.45; 0.03 0.15; 0.2 0.4; 0.1 0.3; 0.5 0.8; 0.1 0.6; 0.4 0.8; 0.4 0.8; 0.4 0.8; 0.001 0.2; 0.22 0.5; 0.4 0.65;
        0.2 0.8; 0.4 0.99];
    fixparam={[0.05 0.0474]', [0.7664 0.4614]', [0.7216 0.0474]', 0.2};   % {[effb1 effb2], effs2, effs3, effb12}
                                    % [effb1 effb2 effb12] are donor and accetor leak into other channels (effb1 = A1/(A1+D), effb2 = A2/(A2+D), effb12 = A2/(A1+A2)
                                    % effs2 and effs3 are FRET efficiencies in state 2 (A2 dark, E = (A1+A2)/(D+A1+A2)) and 3 (A1 dark, E = A2/(D+A2)). These values are fitting parameters in some conditions
    fixparam1=fixparam;
%    fixparam1{5}=[0 0; 0 0; 0 0];
  end
                
resparams=mlhTP3cABIsimfitEintParC(initparams,LUbounds,fixparam1,resparamfix,t3rdata,stateiddata,chans,fitfu,nintm,fitidin)
nparI=(length(initparams)-11)/6;  % number of parallel pathways
if nparI < 3
    mlhfitparamPar{jj}=resparams;
else
    for nfid=1:length(fitidin)
        pfactor=cumprod([1; 1-resparams(12+(5*nparI+1:6*nparI-2),fitidin(nfid))],1);
        mlhfitparamPar{jj}=[resparams(1:12+5*nparI,fitidin(nfid)); resparams(13+5*nparI:end,fitidin(nfid)).*pfactor]; % Convert f to population p
    end
end            
[errorparams logmlhs bics]=mlhTP3cABIsimfitEintParerrorC(mlhfitparamPar{jj},fixparam1,resparamfix,t3rdata,stateiddata,chans,fitfu,nintm,fitidin,drvintv,drvpnt);
errorparamPar{jj}=errorparams;
logmlh(:,jj)=logmlhs;
bic(:,jj)=bics;

save(['TPfitparam3cPar.mat'],'mlhfitparamPar','errorparamPar','logmlh','bic','-mat');

%% Fit E_TP and t_TP of individual transitions
load(['TPfitparam3cPar.mat']);
mlhfitparamPar1=mlhfitparamPar;
drvintv=0.002;
drvpnt=5;

trcc=1;initparams{trcc}(1:6)={[0.08 0.35 ...   % E1_3c
            0.75 0.11 ...    % E2_3c
...%            100 100 0.95 0.95 ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2]
            mlhfitparamPar1{1}(9:12,6)' ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2], fix to the collective fitting values
            0.5 ...         % TP time
            0.1 0.65]'};     % E1TP E2TP
        initparams{trcc}{2}(end-1:end)=[0.1 0.4]';      % For different initial values of E2
        initparams{trcc}{3}(end-1:end)=[0.1 0.15]';
        initparams{trcc}{4}(end-2)=0.05;      % For different initial values of TPT and E2
        initparams{trcc}{5}(end-2:end)=[0.05 0.1 0.4]';      
        initparams{trcc}{6}(end-2:end)=[0.05 0.1 0.15]';
    LUbounds{trcc}=[0.01 0.25; 0.15 0.5; ...
        0.4 0.9; 0.01 0.25;
%        20 2000; 20 2000; 0.8 0.9999; 0.9 0.9999;
        mlhfitparamPar1{1}(9,6)*[0.99999 1.00001]; mlhfitparamPar1{1}(10,6)*[0.99999 1.00001]; mlhfitparamPar1{1}(11,6)*[0.9999 1.0001]; mlhfitparamPar1{1}(12,6)*[0.99999 1.00001];
        0.001 5;
        0.001 0.45; 0.2 0.75];

trcc=2;initparams{trcc}(1:2)={[0.75 0.45 ...    % E1_DA1
...%            100 0.95 ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2]
            mlhfitparamPar1{1}([9 11],6)' ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2]
            0.5 ...         % TP time
            0.6]'};     % E1TP E2TP E1TP_DA1 E2TP_DA2
        initparams{trcc}{2}(end-1)=0.05;
    LUbounds{trcc}=[0.6 1; 0.4 0.6;
%        20 2000; 0.9 0.9999;
        mlhfitparamPar1{1}(9,6)*[0.99999 1.00001]; mlhfitparamPar1{1}(11,6)*[0.9999 1.0001];
        0.001 5;
        0.4 0.75];

trcc=3;initparams{trcc}(1:6)={[0.6 0.04 ...    % E2_DA2
...%            100 0.95 ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2]
            mlhfitparamPar1{1}([10 12],6)' ...        % blinking rate (kb), pb (bright state population), [kb1 kb2 pb1 pb2]
            0.5 ...         % TP time
            0.55]'};     % E1TP E2TP E1TP_DA1 E2TP_DA2
        initparams{trcc}{2}(end)=0.35;
        initparams{trcc}{3}(end)=0.15;
        initparams{trcc}{4}(end-1)=0.05;
        initparams{trcc}{5}(end-1:end)=[0.05 0.35]';
        initparams{trcc}{6}(end-1:end)=[0.05 0.15]';
    LUbounds{trcc}=[0.4 0.8; 0.001 0.2;
%        20 2000; 0.9 0.9999;
        mlhfitparamPar1{1}(10,6)*[0.99999 1.00001]; mlhfitparamPar1{1}(12,6)*[0.99999 1.00001];
        0.001 5;
        0.01 0.7];
    
    
    fixparam={[0.05 mlhfitparamPar1{1}(8,6)]', [mlhfitparamPar1{1}(5,6) sum(mlhfitparamPar1{1}([2 4],6))]', mlhfitparamPar1{1}(7:8,6), 0.2};   % {[effb1 effb2], effs2, effs3, effb12}
                                    % [effb1 effb2 effb12] are donor and accetor leak into other channels (effb1 = A1/(A1+D), effb2 = A2/(A2+D), effb12 = A2/(A1+A2)
                                    % effs2 and effs3 are FRET efficiencies in state 2 (A2 dark, E = (A1+A2)/(D+A1+A2)) and 3 (A1 dark, E = A2/(D+A2)). These values are fitting parameters in some conditions
    fixparam1=fixparam;
    fixparam1{5}=mlhfitparamPar1{1}(16,6);
%    fixparam1{5}=[0.5638 0.5620];

clear mlhfitparam errorparam logmlh dellogmlh bic mlhfitparam0 errorparam0 logmlh0 bic0
clear fretFUall crFUall trjlengthall nphotonsall logmlhvaltrplot1 binfrhistTPall
donchan=chans(3);
ac1chan=chans(2);
ac2chan=chans(1);

maxiterid(1:3)={{[],[]}};maxiterid0=maxiterid;
for trcc=1:3     % 1: 3-color segments, 2: 2-color segments
    trcc
    for zz=1:length(t3rdata{trcc})
        stateidvit=stateiddata{trcc}{zz}(1,3);
        fufu=stateidvit;
        fuplotid=zz;
        onet3r=t3rdata{trcc}{zz};
    
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
                    
            burstbint3r=onet3r;
            burstphotons=size(onet3r,1);
    
            indexone=1:length(burstphotons);
            cumindex=[0; cumsum(burstphotons)];
    
            effb1=fixparam{1}(1);  % E1 of A1,A2 dark state (= A1/(A2 + A1 + D))
            effb2=fixparam{1}(2);  % E2 of A1,A2 dark state (= A2/(A2 + A1 + D)) (depends on the A2-protein concentration)
                    if trcc == 1
                        fitparamid=1:4;
                        initparamone=initparams{trcc}{1};initparamone(9)=1e-6;     % instantaneous transition
                        initparamone(5:8)=mlhfitparamPar{1}(9:12,6);
                        LUboundone=LUbounds{trcc};LUboundone(5:8,:)=[initparamone(5:8)*0.9999 initparamone(5:8)*1.0001];
                        [resparams exitflag]=mlhTP3cABIfit1TPindi(initparamone,LUboundone,fixparam1,resparamfix,burstbint3r,stateidvit,fitparamid);
                        mlhfitparam0{trcc}{fufu}(:,fuplotid)=resparams;
                        [errorparam0{trcc}{fufu}(:,fuplotid) logmlh0{trcc}{fufu}(fuplotid) dellogmlhtemp bic0{trcc}{fufu}(fuplotid)]=mlhTP3cABIfit1TPindierror(resparams,fixparam1,resparamfix,burstbint3r,stateidvit,fitparamid,drvintv,drvpnt);
                        if exitflag == 0, maxiterid0{trcc}{fufu}=[maxiterid0{trcc}{fufu}; mm kk zz]; end
                    elseif trcc == 2
                        fitparamid=1:2;
                        initparamone=initparams{trcc}{1};initparamone(5)=1e-6;     % instantaneous transition
                        initparamone(3:4)=mlhfitparamPar{1}([9 11],6);
                        LUboundone=LUbounds{trcc};LUboundone(3:4,:)=[initparamone(3:4)*0.9999 initparamone(3:4)*1.0001];
                        [resparams exitflag]=mlhTPABIfit1TPindi(initparamone,LUboundone,effb1+effb2,resparamfix,burstbint3r,stateidvit,fitparamid,trcc);
                        mlhfitparam0{trcc}{fufu}(:,fuplotid)=resparams;
                        [errorparam0{trcc}{fufu}(:,fuplotid) logmlh0{trcc}{fufu}(fuplotid) dellogmlhtemp bic0{trcc}{fufu}(fuplotid)]=mlhTPABIfit1TPindierror(resparams,effb2,resparamfix,burstbint3r,stateidvit,fitparamid,trcc,drvintv,drvpnt);
                        if exitflag == 0, maxiterid0{trcc}{fufu}=[maxiterid0{trcc}{fufu}; zz]; end
                    else
                        fitparamid=1:2;
                        initparamone=initparams{trcc}{1};initparamone(5)=1e-6;     % instantaneous transition
                        initparamone(3:4)=mlhfitparamPar{1}([10 12],6);
                        LUboundone=LUbounds{trcc};LUboundone(3:4,:)=[initparamone(3:4)*0.9999 initparamone(3:4)*1.0001];
                        [resparams exitflag]=mlhTPABIfit1TPindi(initparamone,LUboundone,effb2,resparamfix,burstbint3r,stateidvit,fitparamid,trcc);
                        mlhfitparam0{trcc}{fufu}(:,fuplotid)=resparams;
                        [errorparam0{trcc}{fufu}(:,fuplotid) logmlh0{trcc}{fufu}(fuplotid) dellogmlhtemp bic0{trcc}{fufu}(fuplotid)]=mlhTPABIfit1TPindierror(resparams,effb2,resparamfix,burstbint3r,stateidvit,fitparamid,trcc,drvintv,drvpnt);
                        if exitflag == 0, maxiterid0{trcc}{fufu}=[maxiterid0{trcc}{fufu}; zz]; end
                    end                
                for ff=1:length(initparams{trcc})
                    if trcc == 1
                        fitparamid=[1:4 9:11];
                        initparamone=initparams{trcc}{ff};
                        initparamone(5:8)=mlhfitparamPar{1}(9:12,6);
                        LUboundone=LUbounds{trcc};LUboundone(5:8,:)=[initparamone(5:8)*0.9999 initparamone(5:8)*1.0001];
                        [resparams exitflag]=mlhTP3cABIfit1TPindi(initparamone,LUboundone,fixparam1,resparamfix,burstbint3r,stateidvit,fitparamid);
                        mlhfitparam{ff}{trcc}{fufu}(:,fuplotid)=resparams;
                        [errorparam{ff}{trcc}{fufu}(:,fuplotid) logmlh{ff}{trcc}{fufu}(fuplotid) dellogmlh{ff}{trcc}{fufu}(fuplotid) bic{ff}{trcc}{fufu}(fuplotid)]=mlhTP3cABIfit1TPindierror(resparams,fixparam1,resparamfix,burstbint3r,stateidvit,fitparamid,drvintv,drvpnt);
                        if exitflag == 0, maxiterid{trcc}{fufu}=[maxiterid{trcc}{fufu}; ff zz]; end
                    elseif trcc == 2
                        fitparamid=[1:2 5:6];
                        initparamone=initparams{trcc}{ff};
                        initparamone(3:4)=mlhfitparamPar{1}([9 11],6);
                        LUboundone=LUbounds{trcc};LUboundone(3:4,:)=[initparamone(3:4)*0.9999 initparamone(3:4)*1.0001];
                        [resparams exitflag]=mlhTPABIfit1TPindi(initparamone,LUboundone,effb1+effb2,resparamfix,burstbint3r,stateidvit,fitparamid,trcc);
                        mlhfitparam{ff}{trcc}{fufu}(:,fuplotid)=resparams;
                        [errorparam{ff}{trcc}{fufu}(:,fuplotid) logmlh{ff}{trcc}{fufu}(fuplotid) dellogmlh{ff}{trcc}{fufu}(fuplotid) bic{ff}{trcc}{fufu}(fuplotid)]=mlhTPABIfit1TPindierror(resparams,effb2,resparamfix,burstbint3r,stateidvit,fitparamid,trcc,drvintv,drvpnt);
                        if exitflag == 0, maxiterid{trcc}{fufu}=[maxiterid{trcc}{fufu}; ff zz]; end
                    else
                        fitparamid=[1:2 5:6];
                        initparamone=initparams{trcc}{ff};
                        initparamone(3:4)=mlhfitparamPar{1}([10 12],6);
                        LUboundone=LUbounds{trcc};LUboundone(3:4,:)=[initparamone(3:4)*0.9999 initparamone(3:4)*1.0001];
                        [resparams exitflag]=mlhTPABIfit1TPindi(initparamone,LUboundone,effb2,resparamfix,burstbint3r,stateidvit,fitparamid,trcc);
                        mlhfitparam{ff}{trcc}{fufu}(:,fuplotid)=resparams;
                        [errorparam{ff}{trcc}{fufu}(:,fuplotid) logmlh{ff}{trcc}{fufu}(fuplotid) dellogmlh{ff}{trcc}{fufu}(fuplotid) bic{ff}{trcc}{fufu}(fuplotid)]=mlhTPABIfit1TPindierror(resparams,effb2,resparamfix,burstbint3r,stateidvit,fitparamid,trcc,drvintv,drvpnt);
                        if exitflag == 0, maxiterid{trcc}{fufu}=[maxiterid{trcc}{fufu}; ff zz]; end
                    end                
                end
    
    end     % zz
end     % trcc
save(['TP3cABIfitparamindi.mat'],'mlhfitparam','errorparam','logmlh','dellogmlh','bic','maxiterid','mlhfitparam0','errorparam0','logmlh0','bic0','maxiterid0','-mat');

%% Extract and plot simulated tptime and transition window
% Copy and run the definition of initparams and LUbounds in the previous section before running this section
load(['TP3cABIfitparamindi.mat']);
numfit=[6 2 6];   % number of fitting with varied initial parameters
clear plotprmid
for trcc=1:3
    if trcc == 1, fitparamid=[1:4 9:11];
    else fitparamid=[1 2 5 6];
    end
    for fufu=1:2
        plotprmid{trcc}{fufu}=[];
        nonzeroid=[];maxiteridone=[];
        if ~isempty(logmlh{1}{trcc}{fufu}) && size(logmlh{1}{trcc}{fufu},2) >= 1, nonzeroid=find(logmlh{1}{trcc}{fufu} < 0); end
        if ~isempty(maxiterid{trcc}{fufu}), maxiteridone=maxiterid{trcc}{fufu}; end
        for jj=1:length(nonzeroid)
            survidone=ones(numfit(trcc),1);
            clear mlhparamtmp logmlhtemp
            for ff=1:numfit(trcc)
                mlhparamtmp(:,ff)=mlhfitparam{ff}{trcc}{fufu}(:,nonzeroid(jj));
                if any(find(abs(mlhparamtmp(fitparamid,ff) - LUbounds{trcc}(fitparamid,1)) < LUbounds{trcc}(fitparamid,1)*0.01 | abs(mlhparamtmp(fitparamid,ff) - LUbounds{trcc}(fitparamid,2)) < LUbounds{trcc}(fitparamid,2)*0.01))
                    survidone(ff)=0;
                end
                logmlhtemp(ff)=logmlh{ff}{trcc}{fufu}(nonzeroid(jj));
            end
            if ~isempty(maxiteridone), survidone(maxiteridone(find(maxiteridone(:,end) == nonzeroid(jj)),end-1))=0; end
            if trcc == 1
                survidone(find(mlhparamtmp(end,:) > mlhparamtmp(3,:) | mlhparamtmp(end,:) < mlhparamtmp(4,:)))=0;   % consider only E2
            else
                survidone(find(mlhparamtmp(end,:) > mlhparamtmp(1,:) | mlhparamtmp(end,:) < mlhparamtmp(2,:)))=0;
            end
            nzsurvid=find(survidone == 1);
            if ~isempty(nzsurvid)
                [maxval maxid]=max(logmlhtemp(nzsurvid));
                plotprmid{trcc}{fufu}=[plotprmid{trcc}{fufu}; nonzeroid(jj) nzsurvid(maxid) jj];
            end
        end
    end
end

clear TCparam TCerror TPTsim
for trcc=1:3
    if trcc == 1, fitparamid=[1 3 2 4 10 11 9];   % [e1b e2b e1u e2u e1tc e2tc tpt + delta(logmlh)]
    else fitparamid=[1 2 6 5];   % [e2b e2u e2tc tpt + delta(logmlh)]
    end
    
    TCparam{trcc}=[];TCerror{trcc}=[];
    for fufu=1:2
        for jj=1:size(plotprmid{trcc}{fufu},1)
            ptmp=plotprmid{trcc}{fufu}(jj,:);
            TCparam{trcc}=[TCparam{trcc}; mlhfitparam{ptmp(2)}{trcc}{fufu}(fitparamid,ptmp(1))' logmlh{ptmp(2)}{trcc}{fufu}(ptmp(1))-logmlh0{trcc}{fufu}(ptmp(1))];
            TCerror{trcc}=[TCerror{trcc}; errorparam{ptmp(2)}{trcc}{fufu}(fitparamid,ptmp(1))'];
        end
    end
end

trcc=1;
xhist = [-(2/40):(1/40):(41/40)];
tptaxis=0:100:2000;
TCparamone=TCparam{trcc};
figure;
for jj=1:6
    subplot(3,3,jj+floor((jj-1)/2));bar(xhist+diff(xhist(1:2))/2,histc(TCparamone(:,jj),xhist,1),1); xlim([0 1])
end
subplot(3,3,3);bar(tptaxis+diff(tptaxis(1:2))/2,histc(TCparamone(:,7)*1000,tptaxis,1),1);

tptcrit=[0 50 200 1000 10000]/1000;
logmlhcrit=0:0.5:2;
figure;
for jj=1:length(logmlhcrit)
    subplot(2,3,jj);
    plotsubid=find(TCparamone(:,7) > tptcrit(1) & TCparamone(:,7) <= tptcrit(2) & TCparamone(:,8) > logmlhcrit(jj));
    plot(TCparamone(plotsubid,6),TCparamone(plotsubid,5),'bo','markersize',8);
    plotsubid=find(TCparamone(:,7) > tptcrit(2) & TCparamone(:,7) <= tptcrit(3) & TCparamone(:,8) > logmlhcrit(jj));
    hold on; plot(TCparamone(plotsubid,6),TCparamone(plotsubid,5),'co','markersize',8);
    plotsubid=find(TCparamone(:,7) > tptcrit(3) & TCparamone(:,7) <= tptcrit(4) & TCparamone(:,8) > logmlhcrit(jj));
    plot(TCparamone(plotsubid,6),TCparamone(plotsubid,5),'go','markersize',8);
    plotsubid=find(TCparamone(:,7) > tptcrit(4) & TCparamone(:,7) <= tptcrit(5) & TCparamone(:,8) > logmlhcrit(jj));
    plot(TCparamone(plotsubid,6),TCparamone(plotsubid,5),'ro','markersize',8); hold off
    axis equal;xlim([0 0.8]);ylim([0 0.5]);
    title(['\Delta logL > ' num2str(logmlhcrit(jj))]);
end    

trcc=2;
tptaxis=0:100:2000;
TCparamone=TCparam{trcc};
figure;
for jj=1:3
    subplot(3,3,3*jj-1);bar(xhist+diff(xhist(1:2))/2,histc(TCparamone(:,jj),xhist,1),1);xlim([0 1]);
end
subplot(3,3,3);bar(tptaxis+diff(tptaxis(1:2))/2,histc(TCparamone(:,4)*1000,tptaxis,1),1);
figure;
for jj=1:length(logmlhcrit)
    subplot(2,3,jj);
    logLid=find(TCparamone(:,5) > logmlhcrit(jj));
    plot(TCparamone(logLid,3),TCparamone(logLid,4),'o');xlim([0 1]);
    title(['\Delta logL > ' num2str(logmlhcrit(jj))]);
end

trcc=3;
tptaxis=0:100:2000;
TCparamone=TCparam{trcc};
figure;
for jj=1:3
    subplot(3,3,3*jj-1);bar(xhist+diff(xhist(1:2))/2,histc(TCparamone(:,jj),xhist,1),1);xlim([0 1]);
end
subplot(3,3,3);bar(tptaxis+diff(tptaxis(1:2))/2,histc(TCparamone(:,4)*1000,tptaxis,1),1);
figure;
for jj=1:length(logmlhcrit)
    subplot(2,3,jj);
    logLid=find(TCparamone(:,5) > logmlhcrit(jj));
    plot(TCparamone(logLid,3),TCparamone(logLid,4),'o','markersize',8);xlim([0 0.8]);set(gca,'yscale','log');ylim([0.01 10]);
    title(['\Delta logL > ' num2str(logmlhcrit(jj))]);
end
