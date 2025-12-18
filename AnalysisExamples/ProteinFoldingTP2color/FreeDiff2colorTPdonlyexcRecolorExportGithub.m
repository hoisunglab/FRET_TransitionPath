clear all

filepathcommon='D:\HoiSung_Data\hschung\lab\Data\Experimental\';
%filepathcommon='C:\User1\hschung\lab\Data\Experimental\';
Dleak = 0.06;

accchan=1;
donchan=2;
pulsedid=false;
%MLEmodel=1;     % 2-state
MLEmodel=2;     % 2-state with acceptor blinking/bleaching

% Mid intensity, to obtain E and kinetic parameters for viberbi segmentation
filepathM{1}={'2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_1000CW.sptw\GroupMeas_1\'};
filepathH{1}={'2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW.sptw\GroupMeas_1\','2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW.sptw\GroupMeas_2\', ...
                '2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW.sptw\GroupMeas_3\','2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW.sptw\GroupMeas_4\', ...
                '2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW.sptw\GroupMeas_5\','2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW.sptw\GroupMeas_6\'};
filepathH{2}={'2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW_r2.sptw\GroupMeas_1\','2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW_r2.sptw\GroupMeas_2\', ...
                '2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW_r2.sptw\GroupMeas_3\','2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW_r2.sptw\GroupMeas_4\', ...
                '2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW_r2.sptw\GroupMeas_5\','2020\2023\20231223_ZMW_free_diff_WWdomain_CW\WWdomain_HEPES_50mM_Urea_4p5M_TSQ_2mM_Tw20_0p01per_3500CW_r2.sptw\GroupMeas_6\'};
filenameMid = 'WWwaveM';
filenameHid = 'WWwaveH';

switch MLEmodel
    case 1
        fretpeakrange{1}={[0.3 1.1; 0.25 1.1; 0.25 1.1; 0.2 1.1; 0.2 1.1]};
        fretpeakrangeH={[0.3 1.1; 0.25 1.1; 0.25 1.1; 0.2 1.1; 0.2 1.1],[0.3 1.1; 0.25 1.1; 0.25 1.1; 0.2 1.1; 0.2 1.1]};
    case 2
        fretpeakrange{1}={[0.3 1.1; 0.25 1.1; 0.25 1.1; 0.2 1.1; 0.2 1.1]};
        fretpeakrangeH={[0.3 1.1; 0.25 1.1; 0.25 1.1; 0.2 1.1; 0.2 1.1],[0.3 1.1; 0.25 1.1; 0.25 1.1; 0.2 1.1; 0.2 1.1]};
end     
fretpeakrangeHvit=[0.3 0.8; 0.3 0.8];
donlyEcut=[0.2 0.2];
analysisids=[3 1; 3 1];  % [nn k] for high intensity data
effbH=[0.0583 0.0607];    % determined from the donly burst photons below
tTP = 0.0045;        % 2 us
effTP=0.58;

% ------------ remove donly segment and recolor
phintvid=3;     % 3, t_th <= 3 us for I=6000, 5 us for I = 3500
Nthid=1;    % 1, N_th = 60
numrec=5;

load(strcat(filenameHid,'_furatesABIdonlyexc.mat'));
resparamM=permute(mean(resparam(:,1,3,:,:),4),[1 5 2 3 4]);
%load(strcat(filenameMid,'_furatesABI.mat'));
%resparamM(3,2:3)=permute(resparam(3,1,3,:),[1 4 2 3]);

% load(strcat(filenameHid,'ViterbiHMMABI.mat'));      % Viterbi analysis using resparamfix that was obtained from the 2-state maximum likelihood analysis for the mid intensity data
load(strcat(filenameHid,'ViterbiHMMABIhighk.mat'));      % Viterbi analysis using resparamfix that was obtained from the 2-state maximum likelihood analysis for the mid intensity data
drvintv=0.002;
drvpnt=5;
clear burstbint3rrec indexonerec cumindexrec statetrj
for mm=1:length(filepathH)
    mm
    recparam=[resparamM(:,mm); effbH(mm); tTP; effTP];
    for kk=1:length(filepathH{mm})
        kk
        load([filepathcommon filepathH{mm}{kk} 'binburstparams.mat']);
        filerootend = 'free';
        load([filepathcommon filepathH{mm}{kk} 'Sub1and2' filerootend '.mat']);

        % note: nnn and nn replace nn and ii above.
        nn=analysisids(mm,1);     % nn=3, t_th <= 3 us for I=6000, 5 us for I = 3500
        k=analysisids(mm,2);     % k=1, N_th = 60

        datafilenames=[filepathcommon filepathH{mm}{kk} burstbinfiles{nn} 'ThrS' sprintf('%d',thresholds(k))];
        datafiles = textread(datafilenames, '%s', 'delimiter', '\n');
        numfiles = length(datafiles);
        burstbint3r1=[];
        for jj=1:numfiles
            load('-mat', [filepathcommon filepathH{mm}{kk} datafiles{jj}], 'burstbint3r');
            burstbint3r1=[burstbint3r1; burstbint3r];
        end
        burstbint3r=burstbint3r1;
        if any(kk == pulsedid)
            burstbint3r(:,end-2)=burstbint3r(:,end-2)*50/100;     % 50 ns --> 100 ns
        else
            burstbint3r(:,end-2)=burstbint3r(:,end-2)*1e-5;     % 1 ps --> 100 ns
        end
        burstfridone=unique(burstbint3r(:,1));
        histone=histc(burstbint3r(:,1),burstfridone);
        cumindex=[0; cumsum(histone)];
        indexone=1:length(burstfridone);

        burstdelid=find(Totalfretbursts{nn,k}(:,2) < fretpeakrangeHvit(mm,1) | Totalfretbursts{nn,k}(:,2) >= fretpeakrangeHvit(mm,2));
        indexone(burstdelid)=[];

        burstbint3r2=zeros(size(burstbint3r,1),size(burstbint3r,2)+numrec);nphotons=zeros(size(burstbint3r,1),1);nphtot=0;nt3r=1;statetrjone=[];
        for ii=1:length(indexone)
            onet3rb=burstbint3r(cumindex(indexone(ii))+1:cumindex(indexone(ii)+1),:);
            stateidone=[statestrid{mm}{kk}{nn}{ii}; [onet3rb(end,end-2)*1e-7 size(onet3rb,1)+1 -1]];
            FUidone=find(stateidone(1:end-1,end) == 1 | stateidone(1:end-1,end) == 2);
            FUidone=[FUidone FUidone];
            for jj=size(FUidone,1):-1:2
                if FUidone(jj,1)-FUidone(jj-1,2) == 1
                    FUidone(jj-1,2)=FUidone(jj,2);
                    FUidone(jj,:)=[];
                end
            end

            brightidone=[stateidone(FUidone(:,1),2) stateidone(FUidone(:,2)+1,2)-1];
            for jj=1:size(brightidone,1)
                FUphotons=onet3rb(brightidone(jj,1):brightidone(jj,2),:);
                nphotonone=size(FUphotons,1);
                if nphotonone > 1
                    nphotons(nt3r)=nphotonone;
                    recstring=[];
                    for rr=1:numrec
                        [recstringone statetrjtwo]=recolorstrTPABI(recparam,FUphotons);
                        recstring=[recstring recstringone];
                        statetrjone=[statetrjone; statetrjtwo; 0 0 0];
                    end
                    burstbint3r2(nphtot+1:nphtot+nphotonone,:)=[FUphotons recstring];
                    nphtot=nphtot+nphotonone;
                    nt3r=nt3r+1;
                end
            end
        end
        burstbint3r2(nphtot+1:end,:)=[];nphotons(nt3r:end)=[];
        cumindexrec{mm}{kk}=[0; cumsum(nphotons)];
        indexonerec{mm}{kk}=1:length(nphotons);
        burstbint3rrec{mm}{kk}=burstbint3r2;
        statetrj{mm}{kk}=statetrjone;

        save(strcat(filenameHid,'_furatesABIdonlyexcRec',num2str(tTP*1000),'us.mat'),'burstbint3rrec','indexonerec','cumindexrec','statetrj','-mat','-v7.3');
    end
end

% ------------ run maximum likelihood
numrec=5;
load(strcat(filenameHid,'_furatesABIdonlyexcRec',num2str(tTP*1000),'us.mat'));
drvintv=0.002;
drvpnt=5;
clear resparam errorparam logmlh bic
for mm=1:length(filepathH)
    mm
    for kk=1:length(filepathH{mm})
        kk
        
        burstbint3rall=burstbint3rrec{mm}{kk};
        cumindex=cumindexrec{mm}{kk};
        indexone=indexonerec{mm}{kk};
        statetrjone=statetrj{mm}{kk};

        initparams=[0.7 0.45 15 0.5 500 0.99]'; % E1, E2, k, p1, kb, pb
        LUbounds=[0.6 0.9; 0.3 0.55; 0.1 100; 0.1 0.8; 0.1 1000; 0.1 0.999999];      

        for nn=1:numrec
            burstbint3r=burstbint3rall(:,[1:4 5+nn]);
            resparams=mlhrateeiglingenABIC(initparams,LUbounds,burstbint3r,cumindex,indexone);
            nstate=length(initparams)/3;
            pfactor=cumprod([1; 1-resparams(2*nstate:end-3)]);
            resparam(:,nn,kk,mm)=[resparams(1:2*nstate-1); resparams(2*nstate:end-2).*pfactor; resparams(end-1:end)]; % Convert f to population p
            resparam(:,nn,kk,mm)
            % Error is calculated for k (sum of rates between adjacent states) and relative population p.
            [errorparams logmlh(nn,kk,mm) bic(nn,kk,mm)]=mlhrateeiglingenABIerrorC(resparam(:,nn,kk,mm),burstbint3r,cumindex,indexone,drvintv,drvpnt);
            errorparam(:,nn,kk,mm)=errorparams;

            save(strcat(filenameHid,'_furatesABIdonlyexcRecABI',num2str(tTP*1000),'us.mat'),'resparam','errorparam','logmlh','bic','-mat');
        end
    end
end

% Viterbi -------- 1: folded, 2: unfolded ---
numrec=5;
numrows=5;numcols=2;
load(strcat(filenameHid,'_furatesABIdonlyexcRec',num2str(tTP*1000),'us.mat'));
statetrj1=statetrj;
load(strcat(filenameHid,'_furatesABIdonlyexcRecABI',num2str(tTP*1000),'us.mat'));
resparamM=mean(mean(resparam(1:4,:,:,1),3),2);
%resparamM(3)=resparamM(3)/4;
isABI=true;
isplot=false;
maxPlotN=50;
clear logprob statestrid
for mm=1:length(filepathH)
    mm
    for kk=1:length(filepathH{mm})
        kk

        burstbint3rall=burstbint3rrec{mm}{kk};
        cumindex=cumindexrec{mm}{kk};
        indexone=indexonerec{mm}{kk};
        statetrjone=statetrj1{mm}{kk};
        
        logprob{mm}{kk}(1:numrec)={[]};statestrid{mm}{kk}(1:numrec)={[]};plotid=1;
        for ii=1:length(indexone)
            onet3rb0=burstbint3rall(cumindex(indexone(ii))+1:cumindex(indexone(ii)+1),:);
            for nn=1:numrec
                onet3rb=onet3rb0(:,[1:4 5+nn]);
                if ~isABI
%                         [logprobone stateidone]=viterbi(resparam(:,kk,hh,mm),onet3rb);
                    [logprobone stateidone]=viterbiAbk(resparamM,onet3rb);       % viterbiAbk is equivalent to viterbi
                else
                    ratemat00=[[-1; 1]*(1-resparamM(4)) [1; -1]*resparamM(4)]*resparamM(3);
                    Ablinkparam=[effbH(mm); [1; 1]*resparamM(3)/2];
                    [logprobone stateidone]=viterbigenABI(resparamM(1:2),ratemat00,onet3rb,Ablinkparam);
                end

                logprob{mm}{kk}{nn}=[logprob{mm}{kk}{nn}; ii kk indexone(ii) nn logprobone];
                statetrj=[stateidone; [stateidone(2:end,1)-1 stateidone(1:end-1,2)]; [length(onet3rb(:,1)) stateidone(end,2)]];
                statetrj=[onet3rb(statetrj(:,1),end-2)*1e-7 statetrj];
                [junk statetrjsortedid]=sort(statetrj(:,2));
                statetrj=statetrj(statetrjsortedid,:);
                stateidone=[onet3rb(stateidone(:,1),end-2)*1e-7 stateidone];
                statestrid{mm}{kk}{nn}=[statestrid{mm}{kk}{nn} {stateidone}];

                % Used for plot photon trajectory and dividied segments
                if isplot
                    donstring=find(onet3rb(:,end)==2);
                    accstring=find(onet3rb(:,end)==1);
                    if plotid==1 || plotid > numrows*numcols, figure; plotid=1; end
                    subplot(numrows,numcols,plotid);
                    plot(onet3rb(accstring,end-2)*1e-4,ones(size(accstring))*1.15,'ro',onet3rb(donstring,end-2)*1e-4,ones(size(donstring))*.85,'go');
                    if ~isABI
                        hold on;plot(statetrj(:,1)*1e3,2.5-statetrj(:,3),'r-');ylim([0 2]);title([num2str(burstfridone(indexone(ii))) ',' num2str(indexone(ii))]);hold off
                    else
                        hold on;plot(statetrj(:,1)*1e3,2-0.5*statetrj(:,3),'r-');ylim([0 2]);title([num2str(burstfridone(indexone(ii))) ',' num2str(indexone(ii))]);hold off
                    end                        
%                     remtimeone=remtimeall{mm}{hh}{kk}(find(remtimeall{mm}{hh}{kk}(:,1)==fileindexnew{mm}{kk}(nn) & remtimeall{mm}{hh}{kk}(:,2) == nonzerofrid(ii)),7:8);
%                     for xx=1:lengthremid;
%                         line([1 1]*remtimeone(xx,1)*tick/100*1e-4,[0 2]);
%                         line([1 1]*remtimeone(xx,2)*tick/100*1e-4,[0 2]);
%                     end                            
                    plotid=plotid+1;
                    if ii > maxPlotN, break; end
                end
            end
        end
    end
end
if ~isABI
%    save(strcat(filenameHid,'ViterbiHMMdonlyexc.mat'),'logprob','statestrid','-mat');
%     save(strcat(filenameHid,'ViterbiHMMdiv4.mat'),'logprob','statestrid','-mat');
else
%    save(strcat(filenameHid,'ViterbiHMMABI.mat'),'logprob','statestrid','-mat');
    save(strcat(filenameHid,'ViterbiHMMABIhighkdonlyexcRecABI',num2str(tTP*1000),'us.mat'),'logprob','statestrid','-mat');    % Using high intensity rate
%    save(strcat(filenameHid,'ViterbiHMMABIhighkdonlyexcdiv4.mat'),'logprob','statestrid','-mat');    % Using high intensity rate
end

% ----Selection of trajectories with transitions and appropriate parameters
load(strcat(filenameHid,'_furatesABIdonlyexcRec',num2str(tTP*1000),'us.mat'));
phintvid=3;     % 3, t_th <= 3 us for I=6000, 5 us for I = 3500
Nthid=1;    % 1, N_th = 60
numrec=5;
load(strcat(filenameHid,'_furatesABIdonlyexcRecABI',num2str(tTP*1000),'us.mat'));
resparamM=mean(mean(resparam(1:4,:,:,1),3),2);

%load(strcat(filenameHid,'ViterbiHMMABIhighk.mat'));      % Viterbi analysis using resparamfix that was obtained from the 2-state maximum likelihood analysis for the mid intensity data
load(strcat(filenameHid,'ViterbiHMMABIhighkdonlyexcRecABI',num2str(tTP*1000),'us.mat'));      % Viterbi analysis using resparamfix that was obtained from the 2-state maximum likelihood analysis for the mid intensity data
%minsegL=[0.03 0.03];   % minimum F and U segment length in ms
%mincntrate=[2000 2000];   % minimum F and U count rates in ms-1
minsegL=[0.04 0.04];   % minimum F and U segment length in ms      Higher peak
mincntrate=[1800 1800];   % minimum F and U count rates in ms-1
minN=round(minsegL.*mincntrate);
Estd=1.5;
% Erng=1.5*sqrt(resparamM(1:2).*(1-resparamM(1:2))./minN');
% FRETcrit=resparamM(1:2)+Erng*[-1 1];
clear Nphoton fretFU trjLFU crFU onet3rSel logprobSel stateidSel
for mm=1:length(filepathH)
    for kk=1:length(filepathH{mm})
        FUids=ones(numrec,1); % F->U transition
        UFids=ones(numrec,1); % U->F transition
        logprobSel{mm}{kk}(1:5)={{[],[]}};
        
        burstbint3rall=burstbint3rrec{mm}{kk};
        cumindex=cumindexrec{mm}{kk};
        indexone=indexonerec{mm}{kk};
        statetrjone=statetrj{mm}{kk};

        for ii=1:length(indexone)
            onet3rb0=burstbint3rall(cumindex(indexone(ii))+1:cumindex(indexone(ii)+1),:);
            for nn=1:numrec
                onet3rb=onet3rb0(:,[1:4 5+nn]);
                stateidone=[statestrid{mm}{kk}{nn}{ii}; [onet3rb(end,end-2)*1e-7 size(onet3rb,1)+1 -1]];
                FUidone=find(stateidone(1:end-2,end) == 1 & stateidone(2:end-1,end) == 2);
                UFidone=find(stateidone(1:end-2,end) == 2 & stateidone(2:end-1,end) == 1);
                for jj=1:length(FUidone)
                    Fphotons=onet3rb(stateidone(FUidone(jj),2):stateidone(FUidone(jj)+1,2)-1,:);
                    Uphotons=onet3rb(stateidone(FUidone(jj)+1,2):stateidone(FUidone(jj)+2,2)-1,:);
                    Allphotons=onet3rb(stateidone(FUidone(jj),2):stateidone(FUidone(jj)+2,2)-1,:);
                    Nphone=[length(find(Fphotons(:,end) == accchan)) size(Fphotons,1) length(find(Uphotons(:,end) == accchan)) size(Uphotons,1)];
                    fretFUone=[Nphone(1)/Nphone(2) Nphone(3)/Nphone(4)];
    %                trjLone=[diff(Fphotons([1 end],end-2)) diff(Uphotons([1 end],end-2)) diff(Allphotons([1 end],end-2))]*1e-4;
                    trjLone=[diff(stateidone(FUidone(jj)+[0 1],1)) diff(stateidone(FUidone(jj)+[1 2],1)) diff(stateidone(FUidone(jj)+[0 2],1))]*1e3;
                    cntrateone=1e4./[mean(diff(Fphotons(:,end-2))) mean(diff(Uphotons(:,end-2))) mean(diff(Allphotons(:,end-2)))];
    %                cntrateone=[Nphone(2) Nphone(4) sum(Nphone([2 4]))]./trjLone;

                    isNphotonOK=Nphone(2) >= minN(1) & Nphone(4) >= minN(2);
                    Erng=Estd*sqrt(resparamM(1:2).*(1-resparamM(1:2))./Nphone([2 4])');
                    FRETcrit=resparamM(1:2)+Erng*[-1 1];
                    isFRETOK=fretFUone(1) >= FRETcrit(1,1) && fretFUone(1) < FRETcrit(1,2) && fretFUone(2) >= FRETcrit(2,1) && fretFUone(2) < FRETcrit(2,2);
                    isTrjLOK=true;
    %                 isTrjLOK=trjLone(1) >= minsegL(1) & trjLone(2) >= minsegL(2);
    %                isCntrateOK=true;
                     isCntrateOK=cntrateone(1) >= mincntrate(1) & cntrateone(2) >= mincntrate(2);
                    if isNphotonOK && isFRETOK && isCntrateOK && isTrjLOK
                       Nphoton{mm}{kk}{nn}{1}(FUids(nn),:)=Nphone;
                       fretFU{mm}{kk}{nn}{1}(FUids(nn),:)=fretFUone;
                       trjLFU{mm}{kk}{nn}{1}(FUids(nn),:)=trjLone;
                       crFU{mm}{kk}{nn}{1}(FUids(nn),:)=cntrateone;
                       onet3rSel{mm}{kk}{nn}{1}{FUids(nn)}=Allphotons;
                       logprobSel{mm}{kk}{nn}{1}=[logprobSel{mm}{kk}{nn}{1}; logprob{mm}{kk}{nn}(ii,:)];
                       stateidSel{mm}{kk}{nn}{1}{FUids(nn)}=stateidone(FUidone(jj)+[0:2],:);
                       FUids(nn)=FUids(nn)+1;
                    end
                end
                for jj=1:length(UFidone)
                    Uphotons=onet3rb(stateidone(UFidone(jj),2):stateidone(UFidone(jj)+1,2)-1,:);
                    Fphotons=onet3rb(stateidone(UFidone(jj)+1,2):stateidone(UFidone(jj)+2,2)-1,:);
                    Allphotons=onet3rb(stateidone(UFidone(jj),2):stateidone(UFidone(jj)+2,2)-1,:);
                    Nphone=[length(find(Fphotons(:,end) == accchan)) size(Fphotons,1) length(find(Uphotons(:,end) == accchan)) size(Uphotons,1)];
                    fretFUone=[Nphone(1)/Nphone(2) Nphone(3)/Nphone(4)];
    %                trjLone=[diff(Fphotons([1 end],end-2)) diff(Uphotons([1 end],end-2)) diff(Allphotons([1 end],end-2))]*1e-4;
                    trjLone=[diff(stateidone(UFidone(jj)+[1 2],1)) diff(stateidone(UFidone(jj)+[0 1],1)) diff(stateidone(UFidone(jj)+[0 2],1))]*1e3;
                    cntrateone=1e4./[mean(diff(Fphotons(:,end-2))) mean(diff(Uphotons(:,end-2))) mean(diff(Allphotons(:,end-2)))];
    %                cntrateone=[Nphone(2) Nphone(4) sum(Nphone([2 4]))]./trjLone;

                    isNphotonOK=Nphone(2) >= minN(1) & Nphone(4) >= minN(2);
                    Erng=Estd*sqrt(resparamM(1:2).*(1-resparamM(1:2))./Nphone([2 4])');
                    FRETcrit=resparamM(1:2)+Erng*[-1 1];
                    isFRETOK=fretFUone(1) >= FRETcrit(1,1) && fretFUone(1) < FRETcrit(1,2) && fretFUone(2) >= FRETcrit(2,1) && fretFUone(2) < FRETcrit(2,2);
                    isTrjLOK=true;
    %                 isTrjLOK=trjLone(1) >= minsegL(1) & trjLone(2) >= minsegL(2);
    %                isCntrateOK=true;
                     isCntrateOK=cntrateone(1) >= mincntrate(1) & cntrateone(2) >= mincntrate(2);
                    if isNphotonOK && isFRETOK && isCntrateOK && isTrjLOK
                       Nphoton{mm}{kk}{nn}{2}(UFids(nn),:)=Nphone;
                       fretFU{mm}{kk}{nn}{2}(UFids(nn),:)=fretFUone;
                       trjLFU{mm}{kk}{nn}{2}(UFids(nn),:)=trjLone;
                       crFU{mm}{kk}{nn}{2}(UFids(nn),:)=cntrateone;
                       onet3rSel{mm}{kk}{nn}{2}{UFids(nn)}=Allphotons;
                       logprobSel{mm}{kk}{nn}{2}=[logprobSel{mm}{kk}{nn}{2}; logprob{mm}{kk}{nn}(ii,:)];
                       stateidSel{mm}{kk}{nn}{2}{UFids(nn)}=stateidone(UFidone(jj)+[0:2],:);
                       UFids(nn)=UFids(nn)+1;
                    end
                end
            end
        end
    end
end
save([filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'indiEstd' num2str(Estd) 'highkdonlyexcRec',num2str(tTP*1000),'us.mat'], ...
        'Nphoton','fretFU','trjLFU','crFU','onet3rSel','logprobSel','stateidSel','-mat');


% Plot E, cntrate, segment length, number of photons
xhist = [-(2/40):(1/40):(41/40)];
xcnt = 2000:200:6000;
xtrjL = 0:20:2000;
xNph = 0:40:1000;

%load([filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) '.mat']);
%load([filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'highk.mat']);
%load([filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'Estd1.0highk.mat']);
%load([filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'indiEstd' num2str(Estd) 'highk.mat']);
load([filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'indiEstd' num2str(Estd) 'highkdonlyexcRec',num2str(tTP*1000),'us.mat']);
clear meanNph meanE meantrjL meancnt meanNphall meanEall meantrjLall meancntall Nphhistall Ehistall trjLhistall cnthistall
recnum=1;
for mm=1:length(filepathH)
    Nphmm={[],[]};fretmm={[],[]};trjLmm={[],[]};crmm={[],[]};
    Nphhist{mm}=zeros(length(xNph),3,length(filepathH{mm}),2);
    Ehist{mm}=zeros(length(xhist),2,length(filepathH{mm}),2);
    trjLhist{mm}=zeros(length(xtrjL),3,length(filepathH{mm}),2);
    cnthist{mm}=zeros(length(xcnt),3,length(filepathH{mm}),2);
    for fufu=1:2
        for kk=1:length(filepathH{mm})
            Nphmm{fufu}=[Nphmm{fufu}; [Nphoton{mm}{kk}{recnum}{fufu}(:,[2 4]) sum(Nphoton{mm}{kk}{recnum}{fufu}(:,[2 4]),2)]];
            fretmm{fufu}=[fretmm{fufu}; fretFU{mm}{kk}{recnum}{fufu}];
            trjLmm{fufu}=[trjLmm{fufu}; trjLFU{mm}{kk}{recnum}{fufu}];
            crmm{fufu}=[crmm{fufu}; crFU{mm}{kk}{recnum}{fufu}];
            
            Nphhist{mm}(:,:,kk,fufu)=histc([Nphoton{mm}{kk}{recnum}{fufu}(:,[2 4]) sum(Nphoton{mm}{kk}{recnum}{fufu}(:,[2 4]),2)],xNph,1);
            Ehist{mm}(:,:,kk,fufu)=histc(fretFU{mm}{kk}{recnum}{fufu},xhist,1);
            trjLhist{mm}(:,:,kk,fufu)=histc(trjLFU{mm}{kk}{recnum}{fufu}*1e3,xtrjL,1);
            cnthist{mm}(:,:,kk,fufu)=histc(crFU{mm}{kk}{recnum}{fufu},xcnt,1);
        end
        Nphhistall(:,:,fufu,mm)=histc(Nphmm{fufu},xNph,1);
        Ehistall(:,:,fufu,mm)=histc(fretmm{fufu},xhist,1);
        trjLhistall(:,:,fufu,mm)=histc(trjLmm{fufu}*1e3,xtrjL,1);
        cnthistall(:,:,fufu,mm)=histc(crmm{fufu},xcnt,1);    
        
        meanNph(fufu,:,mm)=mean(Nphmm{fufu},1);
        meanE(fufu,:,mm)=mean(fretmm{fufu},1);
        meantrjL(fufu,:,mm)=mean(trjLmm{fufu},1);
        meancnt(fufu,:,mm)=mean(crmm{fufu},1);
    end
end

% titletext={'F -> U','U -> F'};
% for mm=2:length(filepathH)
%     figure;
%     for kk=1:length(filepathH{mm})
%         for fufu=1:2
%             for nn=1:size(Ehist{mm},2)
%                 subplot(length(filepathH{mm}),6,6*(kk-1)+nn+3*(fufu-1));bar(xhist+(xhist(2)-xhist(1))/2,Ehist{mm}(:,nn,kk,fufu),1);set(gca,'xlim',[-0.06 1.06],'ylim',[0 max(Ehist{mm}(:,nn,kk,fufu))*1.15],...
%                         'ticklength',[0.02 0.02],'fontsize',12,'xtick',0:0.2:1);
%                 if kk==length(filepathH{mm}), xlabel('FRET Efficiency'); end
%                 if kk==1 && nn==1, title(titletext{fufu},'fontsize',12); end
%             end
%         end
%     end
% end
% for mm=2:length(filepathH)
%     figure;
%     for kk=1:length(filepathH{mm})
%         for fufu=1:2
%             for nn=1:size(Nphhist{mm},2)
%                 subplot(length(filepathH{mm}),6,6*(kk-1)+nn+3*(fufu-1));bar(xNph,Nphhist{mm}(:,nn,kk,fufu),1);set(gca,'xlim',[0 1000],'ylim',[0 max(Nphhist{mm}(:,nn,kk,fufu))*1.15],...
%                         'ticklength',[0.02 0.02],'fontsize',12,'xtick',0:200:1000);
%                 if kk==length(filepathH{mm}) && nn==1, xlabel('Number of Photons (F, U, F+U)'); end
%                 if kk==1 && nn==1, title(titletext{fufu},'fontsize',12); end
%             end
%         end
%     end
% end
% for mm=2:length(filepathH)
%     figure;
%     for kk=1:length(filepathH{mm})
%         for fufu=1:2
%             for nn=1:size(trjLhist{mm},2)
%                 subplot(length(filepathH{mm}),6,6*(kk-1)+nn+3*(fufu-1));bar(xtrjL+(xtrjL(2)-xtrjL(1))/2,trjLhist{mm}(:,nn,kk,fufu),1);set(gca,'xlim',[0 400],'ylim',[0 max(trjLhist{mm}(:,nn,kk,fufu))*1.15],...
%                         'ticklength',[0.02 0.02],'fontsize',12,'xtick',0:100:400);
%                 if kk==length(filepathH{mm}), xlabel('segment length (\mus)'); end
%                 if kk==1 && nn==1, title(titletext{fufu},'fontsize',12); end
%             end
%         end
%     end
% end
% for mm=2:length(filepathH)
%     figure;
%     for kk=1:length(filepathH{mm})
%         for fufu=1:2
%             for nn=1:size(cnthist{mm},2)
%                 subplot(length(filepathH{mm}),6,6*(kk-1)+nn+3*(fufu-1));bar(xcnt+(xcnt(2)-xcnt(1))/2,cnthist{mm}(:,nn,kk,fufu),1);set(gca,'xlim',[2000 5000],'ylim',[0 max(cnthist{mm}(:,nn,kk,fufu))*1.15],...
%                         'ticklength',[0.02 0.02],'fontsize',12,'xtick',2000:1000:5000);
%                 if kk==length(filepathH{mm}), xlabel('count rate (ms^-^1)'); end
%                 if kk==1 && nn==1, title(titletext{fufu},'fontsize',12); end
%             end
%         end
%     end
% end

titletext={'F -> U','U -> F'};
for mm=1:length(filepathH)
    figure;
    for fufu=1:2
        for nn=1:size(Ehistall,2)
            subplot(4,6,6*(1-1)+nn+3*(fufu-1));bar(xhist+(xhist(2)-xhist(1))/2,Ehistall(:,nn,fufu,mm),1);set(gca,'xlim',[-0.06 1.06],'ylim',[0 max(Ehistall(:,nn,fufu,mm))*1.15],...
                    'ticklength',[0.02 0.02],'fontsize',12,'xtick',0:0.2:1);
            xlabel('FRET Efficiency');
            if nn==1, title(titletext{fufu},'fontsize',12); end
        end
        for nn=1:3
            subplot(4,6,6*(2-1)+nn+3*(fufu-1));bar(xNph,Nphhistall(:,nn,fufu,mm),1);set(gca,'xlim',[0 1000],'ylim',[0 max(Nphhistall(:,nn,fufu,mm))*1.15],...
                    'ticklength',[0.02 0.02],'fontsize',12,'xtick',0:200:1000);
            if nn==1, xlabel('Number of Photons (F, U, F+U)'); end

            subplot(4,6,6*(3-1)+nn+3*(fufu-1));bar(xtrjL+(xtrjL(2)-xtrjL(1))/2,trjLhistall(:,nn,fufu,mm),1);set(gca,'xlim',[0 400],'ylim',[0 max(trjLhistall(:,nn,fufu,mm))*1.15],...
                    'ticklength',[0.02 0.02],'fontsize',12,'xtick',0:100:400);
            if nn==1, xlabel('segment length (\mus)'); end

            subplot(4,6,6*(4-1)+nn+3*(fufu-1));bar(xcnt+(xcnt(2)-xcnt(1))/2,cnthistall(:,nn,fufu,mm),1);set(gca,'xlim',[2000 5000],'ylim',[0 max(cnthistall(:,nn,fufu,mm))*1.15],...
                    'ticklength',[0.02 0.02],'fontsize',12,'xtick',2000:1000:5000);
            if nn==1, xlabel('count rate (ms^-^1)'); end
        end
    end
end

% ----Data export
load(strcat(filenameHid,'_furatesABIdonlyexcRec',num2str(tTP*1000),'us.mat'));
phintvid=3;     % 3, t_th <= 3 us for I=6000, 5 us for I = 3500
Nthid=1;    % 1, N_th = 60
numrec=5;
load(strcat(filenameHid,'_furatesABIdonlyexcRecABI',num2str(tTP*1000),'us.mat'));
resparamM=permute(mean(mean(resparam(1:4,:,:,:),3),2),[1 4 2 3]);

load(strcat(filenameHid,'ViterbiHMMABIhighkdonlyexcRecABI',num2str(tTP*1000),'us.mat'));      % Viterbi analysis using resparamfix that was obtained from the 2-state maximum likelihood analysis for the mid intensity data
minsegL=[0.04 0.04];   % minimum F and U segment length in ms      Higher peak
mincntrate=[1800 1800];   % minimum F and U count rates in ms-1
minN=round(minsegL.*mincntrate);
Estd=1.5;

filenameSel=[filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'indiEstd' num2str(Estd) 'highkdonlyexcRec',num2str(tTP*1000),'us'];
load([filenameSel '.mat']);

exportmm=1;
exportrec=1;   % recolor number

param2st=resparamM(:,exportmm);
effb=effbH(exportmm);    % determined from the donly burst photons below

t3rdata=[];stateiddata=[];
for kk=1:length(filepathH{exportmm})
    for fufu=1:2
        t3rdata=[t3rdata onet3rSel{exportmm}{kk}{exportrec}{fufu}];
        stateiddata=[stateiddata stateidSel{exportmm}{kk}{exportrec}{fufu}];
    end
end
save ProteinFoldingTP2colorExample.mat t3rdata stateiddata param2st mincntrate minN effb tTP effTP -mat

% Likelihood analysis ************
% Include acceptor blinking
confplotlevel=0.95;
load(strcat(filenameHid,'_furatesABIdonlyexcRecABI',num2str(tTP*1000),'us.mat'));
resparamM=permute(mean(mean(resparam(1:4,:,:,:),3),2),[1 4 2 3]);
trplottime0=[1.6 2.5 4 6 7 8 10];trplottime=[1e-5*trplottime0 1e-4*trplottime0 1e-3*trplottime0 0.01*trplottime0]; % in ms
%filenameSel=[filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2))];
%filenameSel=[filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'highk'];
%filenameSel=[filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'Estd1.0highk'];
%filenameSel=[filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'indiEstd' num2str(Estd) 'highk'];
filenameSel=[filenameHid 't3rSelected_phintvid' num2str(phintvid) 'Nthid' num2str(Nthid) 'crF' num2str(mincntrate(1)) 'crU' num2str(mincntrate(2)) 'NF' num2str(minN(1)) 'NU' num2str(minN(2)) 'indiEstd' num2str(Estd) 'highkdonlyexcRec',num2str(tTP*1000),'us'];
load([filenameSel '.mat']);

drvintv=0.002;
drvpnt=5;
initparams=[0.7 0.42 500 0.99]';        % [EF, EU, blinking rate (kb), pb (bright state population)], for mlhTPABIfit0
LUbounds=[0.6 1; 0.3 0.5; 50 3000; 0.1 0.99999];
isincE=1;   % 2 for constant E_int = (E_F + E_U)/2
fitfu=3;
% initparams=[0.88; 0.48; 0.15; 500; 0.95; 5000];    % for mlhTPABI2fastUfit0, [EF, EU, deltaEU, blinking rate (kb), pb (bright state population), ratesum2U]
% LUbounds=[0.7 1; 0.4 0.65; 0.001 0.35; 50 3000; 0.1 1; 2000 20000];

clear mlhfitparam errorparam logmlh
for mm=1:length(filepathH)
    effb=effbH(mm);
    resparamfix=resparamM(:,mm);
    for nn=1:numrec
        burstbint3r=[];stateidkk=[];
        for kk=1:length(filepathH{mm})
            for fufu=1:2
                burstbint3r=[burstbint3r onet3rSel{mm}{kk}{nn}{fufu}];
                stateidkk=[stateidkk stateidSel{mm}{kk}{nn}{fufu}];
            end
        end

        for nintm=1
            resparams=mlhTPABIsimfit0C(initparams,LUbounds,effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm,isincE);
            mlhfitparam(:,nn,mm)=resparams
            [errorparam(:,nn,mm) logmlh(nn,mm)]=mlhTPABIsimfit0errorC(mlhfitparam(:,nn,mm),effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm,isincE,drvintv,drvpnt);

    %         resparams=mlhTPABI2fastUsimfit0C(initparams,LUbounds,effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm,isincE);
    %         mlhfitparam(:,mm)=resparams
    %         [errorparam(:,mm) logmlh(mm)]=mlhTPABI2fastUsimfit0errorC(mlhfitparam(:,mm),effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm,isincE,drvintv,drvpnt);
        end
    end
end
%save([filenameSel 'TPfit0param.mat'],'mlhfitparam','errorparam','logmlh','-mat');
save([filenameSel 'TPfit0paramdonlyexcRec',num2str(tTP*1000),'us.mat'],'mlhfitparam','errorparam','logmlh','-mat');
% save([filenameSel 'TP2fastUfit0param.mat'],'mlhfitparam','errorparam','logmlh','-mat');

% Compute likelihood values using fitted parameters    
nintm=1;
isincE=1;   % 2 for constant E_int = (E_F + E_U)/2
fitfu=3;
%load([filenameSel 'TPfit0param.mat']);
load([filenameSel 'TPfit0paramdonlyexcRec',num2str(tTP*1000),'us.mat']);
% load([filenameSel 'TP2fastUfit0param.mat']);
clear logmlhkk
titletext={'F -> U','U -> F','All'};
for mm=1:length(filepathH)
    figure;
    effb=effbH(mm);
    resparamfix=resparamM(:,mm);
    for nn=1:numrec
        nn
        burstbint3r=[];stateidkk=[];
        for kk=1:length(filepathH{mm})
            for fufu=1:2
                burstbint3r=[burstbint3r onet3rSel{mm}{kk}{nn}{fufu}];
                stateidkk=[stateidkk stateidSel{mm}{kk}{nn}{fufu}];
            end
        end

        logmlhres=mlhTPABIsimcalC(mlhfitparam(:,nn,mm),trplottime,effb,resparamfix,burstbint3r,stateidkk,nintm,isincE);
    %     logmlhres=mlhTPABIsim2fastUcalC(mlhfitparam(:,mm),trplottime,effb,resparamfix,burstbint3r,stateidkk,nintm,isincE);

        logmlhkk(:,1)=logmlhres{1}-logmlhres{1}(1);
        logmlhkk(:,2)=logmlhres{2}-logmlhres{2}(1);
        logmlhkk(:,3)=sum(logmlhkk(:,1:2),2);

        for fitfu=1:3
            subplot(3,numrec,nn+numrec*(fitfu-1));
            plot(trplottime,logmlhkk(:,fitfu),'o-','markersize',8);set(gca,'xscale','log','xtick',10.^(-5:-1),'ticklength',[0.02 0.02]);
            ylim([-15 15]);xlim([1e-5 3e-2]);
            line([1e-5 3e-2],-[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
            line([1e-5 3e-2],[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');axis square
            title(titletext{fitfu});
        end
    end
end

% Compute likelihood values using fitted parameters, intE
intE=[0.1:0.1:0.4 0.5 0.6 0.65 0.7 0.8];           % *************** intE: list of intermediate FRET efficiencies for *intE analysis
%load([filenameSel 'TPfit0param.mat']);
load([filenameSel 'TPfit0paramdonlyexcRec',num2str(tTP*1000),'us.mat']);
% load([filenameSel 'TP2fastUfit0param.mat']);
for mm=1:length(filepathH)
    effb=effbH(mm);
    resparamfix=resparamM(:,mm);
    for nn=1:numrec
        nn
        burstbint3r=[];stateidkk=[];
        for kk=1:length(filepathH{mm})
            for fufu=1:2
                burstbint3r=[burstbint3r onet3rSel{mm}{kk}{nn}{fufu}];
                stateidkk=[stateidkk stateidSel{mm}{kk}{nn}{fufu}];
            end
        end
        for uu=1:length(intE)
            uu
            logmlhres=mlhTPABIsimcalC(mlhfitparam(:,nn,mm),trplottime,effb,resparamfix,burstbint3r,stateidkk,nintm,intE(uu));
        %     logmlhres=mlhTPABIsim2fastUcalC(mlhfitparam(:,mm),trplottime,effb,resparamfix,burstbint3r,stateidkk,nintm,intE(uu));
            logmlhvalintE{1}(:,uu,nn,mm)=logmlhres{1};
            logmlhvalintE{2}(:,uu,nn,mm)=logmlhres{2};
        end
    end
end
%save([filenameSel 'ABIintEC.mat'],'logmlhvalintE','-mat');
save([filenameSel 'ABIintECdonlyexcRec',num2str(tTP*1000),'us.mat'],'logmlhvalintE','-mat');
%save([filenameSel 'ABI2fastUintEC.mat'],'logmlhvalintE','-mat');

%load([filenameSel 'ABIintEC.mat']);
load([filenameSel 'ABIintECdonlyexcRec',num2str(tTP*1000),'us.mat']);
%load([filenameSel 'ABI2fastUintEC.mat']);
recnum=1;
for mm=1:length(filepathH)
    figure;
    for uu=1:length(intE)
        subplot(3,4,uu);plot(trplottime,logmlhvalintE{1}(:,uu,recnum,mm)-logmlhvalintE{1}(1,uu,recnum,mm)+logmlhvalintE{2}(:,uu,recnum,mm)-logmlhvalintE{2}(1,uu,recnum,mm),'o-','markersize',6);set(gca,'xscale','log','xtick',10.^(-4:-1),'ticklength',[0.02 0.02]);
        ylim([-15 15]);xlim([1e-5 3e-2]);axis square
        line([1e-5 3e-2],-[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
        line([1e-5 3e-2],[1 1]*log(confplotlevel/(1-confplotlevel)),'color','r');
        title(['E = ' num2str(intE(uu))]);        
    end
end

% Fit tptime together with other parameters
drvintv=0.005;
drvpnt=5;
nintm=1;
isincE=1; % 2: constant intermediate E
fitfu=3;
%initparams=[0.8 0.45 500 0.95 0.002]';        % [EF, EU, blinking rate (kb), pb (bright state population), tpt in ms]
%LUbounds=[0.7 1; 0.4 0.65; 50 3000; 0.5 1; 1e-4 0.01];
initparams=[0.7 0.42 500 0.95 0.001 0.55]';        % [EF, EU, blinking rate (kb), pb (bright state population), tpt in ms, Eint]
LUbounds=[0.6 0.8; 0.3 0.5; 50 3000; 0.5 0.99999; 1e-4 0.1; 0.45 0.65];

clear mlhfitparam errorparam logmlh dellogmlh;
for mm=1:length(filepathH)
    effb=effbH(mm);
    resparamfix=resparamM(:,mm);
    for nn=1:numrec
        nn
        burstbint3r=[];stateidkk=[];
        for kk=1:length(filepathH{mm})
            for fufu=1:2
                burstbint3r=[burstbint3r onet3rSel{mm}{kk}{nn}{fufu}];
                stateidkk=[stateidkk stateidSel{mm}{kk}{nn}{fufu}];
            end
        end

%        for nintm=1:3
        for nintm=1
%            for fitfu=1:3
            for fitfu=3
    %             resparams=mlhTPABIsimfitNintmEC(initparams,LUbounds,effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm,isincE);
    %             mlhfitparam(:,fitfu,mm)=resparams
    %             [errorparam(:,fitfu,mm) logmlh(fitfu,mm)]=mlhTPABIsimfitNintmEerrorC(mlhfitparam(:,fitfu,mm),effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm,isincE,drvintv,drvpnt);

                resparams=mlhTPABIsimfitEintNintmEC(initparams,LUbounds,effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm);
                mlhfitparam(:,fitfu,nintm,nn,mm)=resparams;
                resparams
                [errorparam(:,fitfu,nintm,nn,mm) logmlh(fitfu,nintm,nn,mm) dellogmlh(fitfu,nintm,nn,mm)]=mlhTPABIsimfitEintNintmEerrorC(mlhfitparam(:,fitfu,nintm,nn,mm),effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm,drvintv,drvpnt);
            end
        end
    end
end
%save([filenameSel 'TPfitparamincEC.mat'],'mlhfitparam','errorparam','logmlh','-mat');
%save([filenameSel 'TPfitparamconstEfitCEint.mat'],'mlhfitparam','errorparam','logmlh','-mat');
save([filenameSel 'TPfitparamconstEfitCEintdonlyexcRec',num2str(tTP*1000),'us.mat'],'mlhfitparam','errorparam','logmlh','-mat');

% Fit int E with TPT for more than one intermediate
%isincE=1; % 2: constant intermediate E
nintm=1;    % not used in this fitting
drvintv=0.005;
drvpnt=5;
fitfu=3;
% % 1 pathway
% jj=1;
% initparams=[0.65 0.28 200 0.95 0.2 0.4]';        % [EF, EU, blinking rate (kb), pb (bright state population), TP time, E_intm]
% LUbounds=[0.5 1; 0.2 0.5; 100 2000; 0.8 1; 0.1 0.8; 0.32 0.5];
% 2 pathways
jj=2;
initparams=[0.8 0.45 500 0.95 ...      % [EF, EU, blinking rate (kb), pb (bright state population),
        0.001 0.01 0.65 0.55 0.5]';        % [TP time, E_intm, fraction pathways f], p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
LUbounds=[0.7 1; 0.4 0.65; 50 3000; 0.5 1; ...
        1e-4 0.05; 5e-4 0.05; 0.45 0.75; 0.45 0.75; 0.01 0.99];
% % 3 pathways
% jj=3;
% initparams=[0.65 0.28 200 0.95 ...      % [EF, EU, blinking rate (kb), pb (bright state population),
%         0.1 0.6 0.5 0.3 0.5]';        % [TP time, E_intm, fraction pathways f], p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
% LUbounds=[0.5 1; 0.2 0.5; 100 2000; 0.8 1; ...
%         0.05 0.5; 0.2 0.9; 0.32 0.5; 0.4 0.58; 0.1 0.9];

clear mlhfitparam errorparam logmlh dellogmlh bic;
for mm=1:length(filepathH)
    effb=effbH(mm);
    resparamfix=resparamM(:,mm);
    for nn=1:numrec
        burstbint3r=[];stateidkk=[];
        for kk=1:length(filepathH{mm})
            for fufu=1:2
                burstbint3r=[burstbint3r onet3rSel{mm}{kk}{nn}{fufu}];
                stateidkk=[stateidkk stateidSel{mm}{kk}{nn}{fufu}];
            end
        end

        resparams=mlhTPABIsimfitEintParC(initparams,LUbounds,effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm);
        nparI=(length(initparams)-3)/3;  % number of parallel pathways
        pfactor=cumprod([1; 1-resparams(4+(2*nparI+1:3*nparI-2))]);
        mlhfitparam(:,jj,nn,mm)=[resparams(1:4+2*nparI); resparams(5+2*nparI:end).*pfactor]; % Convert f to population p
        mlhfitparam(:,jj,nn,mm)
        [errorparam(:,jj,nn,mm) logmlh(jj,nn,mm) dellogmlh(jj,nn,mm) bic(jj,nn,mm)]=mlhTPABIsimfitEintParerrorC(mlhfitparam(:,jj,nn,mm),effb,resparamfix,burstbint3r,stateidkk,fitfu,nintm,drvintv,drvpnt);
    end
end
%save([filenameSel 'TPfitparamParallelfitC.mat'],'mlhfitparam','errorparam','logmlh','dellogmlh','bic','-mat');
save([filenameSel 'TPfitparamParallelfitCdonlyexcRec',num2str(tTP*1000),'us.mat'],'mlhfitparam','errorparam','logmlh','dellogmlh','bic','-mat');
