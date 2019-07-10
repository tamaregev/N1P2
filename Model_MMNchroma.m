%Adapting from Model_N1P2 to MMNchroma experiment
%May 3 2019

%% Definitions
cd('L:\Experiments\MMNchroma\Analysis')
Definitions_MMNchroma;%new idea for the first time! move definitions into a separate script!

order = cs;
N1P2GHfolder = ['L:\Experiments\N1P2\Analysis\N1P2_GH'];

addpath('L:\Z backup\Tamar\fromZ\Documents\MATLAB\MatlabFunctions\mine')
addpath(N1P2GHfolder)
%% compare EDAT to events
% don;t do this for MMNchroma because of recoding markers.
%
% s=13;
% %load expdata and markers
% FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
% load([EDATfolder FileName '_expdata.mat' ])
% %load and read markers
% VMRKfile = [ExportFolder FileName  '_dt_RDI_imported.vmrk'];
% [eventCoInd, ~]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed
% stimCoInd = eventCoInd(~ismember(eventCoInd(:,1),[210,200,100,110,254]),:);
% trials = expdata.trials;
% phaseName = 'Passive';
% block_list = expdata.(phaseName).block_list;
% nBlocks = length(block_list);
% if(~length(trials)==length(stimCoInd))
%     warning('nTrials is not equal in edat and events indexes')
% end
% newbli = find(eventCoInd(:,1)==100);%new block index
% if(~length(newbli)==nBlocks)
%     warning('nBlocks is not equal in edat and events indexes')
% end
% 
% % compare SOAs
% expdata_plannedSOA = nan(length(trials)-1,1);
% expdata_measuredSOA = nan(length(trials)-1,1);
% events_byIndSOA = nan(length(trials)-1,1);
% sample = nan(length(trials)-1,1);
% difs_reals = nan(length(trials)-1,1);
% eventInd = nan(length(trials)-1,1);
% 
% for it = 1:(length(trials)-1)
%     expdata_plannedSOA(it,1) = trials(it).SOA;
%     expdata_measuredSOA(it,1) = trials(it+1).time-trials(it).time;
%     events_byIndSOA(it,1) = (stimCoInd(it+1,2)-stimCoInd(it,2))/srate;
%     sample(it,1) = stimCoInd(it,2);
%     eventInd(it,1) = it;
% end
% difs_reals(:,:) = abs(expdata_measuredSOA-events_byIndSOA);
% compareSOA = table(eventInd,sample,expdata_plannedSOA,expdata_measuredSOA,events_byIndSOA,difs_reals);
% bigDifInd = find(difs_reals>=0.02);
% compareSOA(bigDifInd,:)
% %correct SOA is expdata_measuredSOA, because it is based on measured timing

%
%% calculate RA and plot - one participant
        
%parameters:
R0=0.5;
phaseName = 'stimParamsBlocksMMN';
sigma = 10;%MIDI
tau = 2;%seconds
SOA_threshold = 0.6;%seconds
Mis = nan;%
%Mis=[20:130]';
vmrkfiletag = 'ImportMarkers_ReCoded';
%vmrkfiletag = 'RDI_imported';%for N1P2

for s=23
    %whichSubjects
    %load expdata and markers
    FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
    %cd(EDATfolder)
    if str2num(sessions{s})>1
        %fix for s23:
        load(['D' EDATfolder(2:end) FileName '2_expdata.mat' ])
    else
        load(['D' EDATfolder(2:end) FileName '_expdata.mat' ])
    end
    %load([FileName '_expdata.mat' ])
    %cd(AnalysisFolder)
    %load and read markers
    VMRKfile = [ExportFolder FileName  '_dt_' vmrkfiletag '.vmrk'];
    [eventCoInd, artInd]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed

    %calculate expected activity RA
    [ RA, smpls, stimCodes, seqInd ] = calcRA_MMNchroma( R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold, Mis, bls, order);

% pause
    %plot RA
    if 1
        ERPfigure;
        isp=0;
        for it=1:size(RA,1)
            for ib=1:size(RA,2)
                isp=isp+1;
                subplot(size(RA,1),size(RA,2),isp);
                if isnan(Mis)
                    imagesc(1:size(RA,4),5,squeeze(RA(it,ib,:,:)))
                else
                    imagesc(1:size(RA,4),Mis,squeeze(RA(it,ib,:,:)))
                end
                if ib==1
                    ylabel(['Type: ' num2str(it)])
                end
                if it==1
                    title(['Sequence #: ' num2str(ib)])
                end
                colorbar;
            end
        end
        pause(0.1)
        suptitle(['Subj ' num2str(s) ', Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % average RA as ERP
    phaseName = 'stimParamsBlocksMMN';
    expdata.blocks = expdata.(phaseName);
    ERRA = nan(size(stimCodes,1),size(stimCodes,2),5);
    CIRA = nan(size(stimCodes,1),size(stimCodes,2),5);

    for it = 1:size(stimCodes,1)
        for is = 1:size(stimCodes,2)
            Codes = squeeze(stimCodes(it,is,:));
            Codes = Codes - floor(Codes/100)*100;
            Codes = Codes - floor(Codes/1000)*1000;
            whichCodes = unique(Codes);
            whichCodes = whichCodes([2,3,1,4,5]);
            for ic=1:length(whichCodes)
                MIDIs = [freq2MIDI([expdata.blocks(it).toneSynth.Fstandards{:}]) , freq2MIDI([expdata.blocks(it).toneSynth.Fdeviants{:}])];
                MIDIs = MIDIs(order);%MIDIs of all stimuli
                Mj = MIDIs(whichCodes==whichCodes(ic));%current stimulus MIDI
               
                is5first = squeeze(seqInd(it,is,:))<=5;
                if isnan(Mis)
                    ERRA(it,is,ic) = mean(RA(it,is,ic,Codes==whichCodes(ic)),4);
                    CI = Confidence(RA(it,is,ic,Codes==whichCodes(ic) & ~is5first));
                else
                    ERRA(it,is,ic) = mean(RA(it,is,Mis==Mj,Codes==whichCodes(ic)),4);
                    CI = Confidence(RA(it,is,Mis==Mj,Codes==whichCodes(ic) & ~is5first));
                end
                CIRA(it,is,ic) = ERRA(it,is,ic)-CI(2);
            end
        end
    end
    ERRAm = squeeze(mean(ERRA,2));
    CIRAm = squeeze(mean(CIRA,2));

    ERPfigure;
    %plot(ERRAm','linewidth',2)
    errorbar(ERRAm',CIRAm','linewidth',2)
    set(gca,'fontsize',14)
    legend(blocks(bls))
    suptitle(['Subj ' num2str(s) ', Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)])
end
%% calculate RA - all participants, taus, sigmas

mkdir(modelFolder,date)
saveFolder = [modelFolder date filesep];

%parameters:
R0=0.5;
phaseName = 'stimParamsBlocksMMN';
%as in N1P2 - 19-Dec-2018 : full. sigmas 1:18. taus 0.2:0.2:10
sigmas = [1:15];%MIDI
taus = [0.01:0.05:1];%seconds
SOA_threshold = 0.6;%seconds
%Mis = (20:130)';
Mis = nan; %for having only the 5 Mis of the stimuli
artIndss = cell(size(whichSubjects));

vmrkfiletag = 'ImportMarkers_ReCoded';
%vmrkfiletag = 'RDI_imported';%for N1P2

tictot = tic;
for s=whichSubjects
    ticsubj = tic;
    fprintf(['subj %2.0f loading...'],s)
    %load expdata and markers
    FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
    if length(sessions{s})>1
        %fix for s23:
        load(['D' EDATfolder(2:end) FileName '2_expdata.mat' ])
    else
        load(['D' EDATfolder(2:end) FileName '_expdata.mat' ])
    end
    %load and read markers
    VMRKfile = [ExportFolder FileName  '_dt_' vmrkfiletag '.vmrk'];
    [eventCoInd, artInd]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed
    artIndss{s} = artInd;
    fprintf('Done\n')
    for sigma = sigmas
        ticsig = tic;
        fprintf(['Subj = %2.0f, sigma = %2.2f'],s,sigma)
        for tau = taus
            %calculate expected activity RA
            if sigma == sigmas(1) && tau == taus(1)
                [ RA, smpls, stimCodes, seqInds ] = calcRA_MMNchroma(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold,Mis,bls,order);

                if s==whichSubjects(1)
                    %initialize only once
                    smplss = nan([whichSubjects(end),size(smpls)]); 
                    stimCodess = nan([whichSubjects(end),size(stimCodes)]);
                    seqIndss = nan([whichSubjects(end),size(seqInds)]);        
                    %initialize one big array for everything
                    RAs_SigTau = nan([whichSubjects(end),size(RA),length(sigmas),length(taus)]);
                    RAs_SigTau = single(RAs_SigTau);
                end
            else
                [ RA, ~, ~, ~ ] = calcRA_MMNchroma(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold,Mis,bls,order);
            end
            RAs_SigTau(s,:,:,:,1:size(RA,4),sigmas==sigma,taus==tau) = RA;
            %clear RA
            if sigma == sigmas(1) && tau == taus(1)%once per subject            
                smplss(s,:,:,1:size(smpls,3)) = smpls;
                stimCodess(s,:,:,1:size(stimCodes,3)) = stimCodes;
                seqIndss(s,:,:,1:size(seqInds,3)) = seqInds;
            end
            fclose('all');
%             disp(['Done in ' num2str(toc(ticsigtau)) ' sec.'])
        end
        fprintf(['Done in %2.1f sec.\n'],toc(ticsig))
    end
    disp(['Done Subj ' Subjects{s} ' in ' num2str(toc(ticsubj)) ' sec.'])
end
fprintf('saving RAs_SigTau ...');ticsaving = tic;
save([saveFolder 'RAs_SigTau'],'RAs_SigTau','-v7.3')
fprintf(['done saving in %2.2f sec.\n'],toc(ticsaving))
save([saveFolder 'Metadata'],'smplss','stimCodess','seqIndss','artIndss')
save([saveFolder 'Params'],'sigmas','taus','Mis','R0','SOA_threshold')
disp(['Done all in ' num2str(toc(tictot)) ' sec.'])

%% compare to ERPs -
%% calc. expected peaks
%RAdate = '22-May-2019';
RAdate = '17-Jun-2019';

loadFolder = [modelFolder RAdate filesep];
saveFolder = loadFolder;

tic
fprintf(['Loading...'])
load([loadFolder 'Params'])
load([loadFolder 'RAs_SigTau'])
load([loadFolder 'Metadata'])
fprintf(['Done in %4.1f sec \n'], toc)

nStim = 5;
%initiatize peaks array
RApeaks = nan(whichSubjects(end),length(bls),nStim,nStim,length(sigmas),length(taus));%subjects x types x cond x prevcond x sigmas x taus
conditions = {'note 1','note 2','note 3','note 4','note 5'};
prevConditions = {'note 1','note 2','note 3','note 4','note 5'};

%STIMcodes = [1,2,3,4,5];
STIMcodes = [20,30,11,40,50];
EventCodes=cell(1,length(conditions));
phaseCodes = [1000, 3000];

tictot = tic;
ibl=0;
for bl = bls
    ticblock = tic;
    ibl=ibl+1;
    blockName = blocks{bl};
    disp(blockName)
    
    Bcode = 100*bl;
    EventCodes=cell(1,length(conditions));
    for i=1:length(conditions)
        EventCodes{i} = phaseCodes + Bcode + STIMcodes(i);
    end
    prevEventCodes = EventCodes;

    grandRA = calcGrandRA_CondPrev( ibl, whichSubjects, RAs_SigTau, EventCodes, prevEventCodes, stimCodess, smplss, seqIndss, artIndss);
    RApeaks(:,ibl,:,:,:,:) = squeeze(grandRA);
end
save([saveFolder 'RApeaks'],'RApeaks')
%% Documentation of runs
%%22-May-2019 : full as in - N1P2 19-Dec-2018 
%sigmas 1:18, taus 0.2:0.2:10
%%03-Jun-2019 : for calculating ellipsoid of N1: sigmas 1:9 taus 0.9:0.1:7
%17-Jun-2019 : for checking lower bound on tau of P2 : sigmas = 1:15, taus= 0.01:0.05:1
%% average and plot U-shape
%RAdate = '22-May-2019';
RAdate = '03-Jun-2019';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
load([loadFolder 'Metadata'])

plotIndividual = false;
phaseName = 'Passive';
sigma = 5;
tau = 6;
nStim = 5;

ERRA = nan(whichSubjects(end),nStim,size(stimCodess,2));
CIRA = nan(nStim,size(stimCodess,2));

%for each subject:
for s = whichSubjects
    stimCodes = squeeze(stimCodess(s,:,:,:));
    for it = 1:size(stimCodes,1)
        Codes = squeeze(stimCodes(it,1,:));
        Codes = Codes - floor(Codes/100)*100;
        Codes = Codes - floor(Codes/1000)*1000;
        whichCodes = unique(Codes);
        ERRA(s,:,it) = nanmean(RApeaks(s,it,:,:,sigmas==sigma,taus==tau),4);
    end
    if plotIndividual
        ERPfigure;
        plot(squeeze(ERRA(s,:,:)),'linewidth',2)
        set(gca,'fontsize',14)
        legend(Bnames)
        suptitle(['Subj ' num2str(s) ', Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)])
    end
end

%grand average across subjects:
for it = 1:size(stimCodes,1)
    for ic = 1:nStim
        CI = Confidence(ERRA(:,ic,it));
        CIRA(ic,it) = squeeze(nanmean(ERRA(:,ic,it),1)) - CI(1);
    end
end
%save([saveFolder 'ERRA'],'ERRA','CIRA','sigma','tau')

ERPfigure;
plot(squeeze(nanmean(ERRA(:,:,:),1)),'linewidth',2)
%errorbar(squeeze(nanmean(ERRA(:,:,:),1)),CIRA,'linewidth',2)
set(gca,'fontsize',14)
legend(blocks(bls))
suptitle(['Grand RA, Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)])
%% average and plot prevcond bars
RAdate = '22-May-2019';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
load([loadFolder 'Metadata'])


sigma = 7;
tau = 2;
nStim = 5;

ERRA = nan(whichSubjects(end),nStim,nStim,size(stimCodess,2));
%for each subject:
for s = whichSubjects
    stimCodes = squeeze(stimCodess(s,:,:,:));
    for it = 1:size(stimCodes,1)
        ERRA(s,:,:,it) = RApeaks(s,it,:,:,sigmas==sigma,taus==tau);
    end
end

for it=1:length(bls)
    ERPfigure;
    %set(gcf,'Position',[10+300*it 10 300 700])
    for curr = 1:5 
        subplot(5,1,curr)
        RA = squeeze(nanmean(ERRA(:,curr,:,it),1));
        bar(1:5,RA,'linewidth',2)
        %errorbar(squeeze(nanmean(ERRA(:,:,:),1)),CIRA,'linewidth',2)
        set(gca,'fontsize',14)
        legend(blocks(bls(it)))
        title(['Curr note ' num2str(curr)])
    end
    xlabel('previous note')
    suptitle({['Grand RA prevcond, Block type: ' num2str(it)] [' Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)]})
end

%% fitting parameters to all block types together -
%% find best scaling factor
weightedFlag = false;
bestSFflag = false;
%load Model
RAdate = '22-May-2019';
%RAdate = '19-Dec-2018';
%RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks.'])
load([loadFolder 'Params'])

%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '18-Nov-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));
        
SFs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks));
  
bls = [2 4];
for ipeak = 1:length(whichpeaks)
    ibl=0;
    for bl=bls
        ibl=ibl+1;
        peakAmps(:,ibl,:,:) = allPeak_amps{bl}(whichSubjects,:,:,ipeak);
        %subj x bl x con x prevcon 
%         peakAmps(:,bl,:) = allGrandcon_amps{bl}(whichSubjects,:,ipeak);
    end
    ste = nanstd(peakAmps)/sqrt(size(peakAmps,1));
    ste = ste(:);ste(isnan(ste))=[];
    data = nanmean(peakAmps);
    data = data(:);data(isnan(data))=[];
    for isig = 1:length(sigmas)
        for itau = 1:length(taus)
            model = 1-nanmean(RApeaks(whichSubjects,:,:,:,isig,itau));
            model = model(:);model(isnan(model))=[];
            if weightedFlag
                SFs(isig,itau,ipeak) = model\(data./ste);%performs SS linear regression
            elseif bestSFflag
                SFs(isig,itau,ipeak) = model\data;%performs SS linear regression
            else
                SFs(isig,itau,ipeak) = max(abs(data))*sign(data(1));
                %SFs(isig,itau,ipeak) = max(abs(data))*sign(data(1))/max(model);
            end
        end
    end
end
ERPfigure
for ipeak=1:length(whichpeaks)
    subplot(1,2,ipeak)
    imagesc(SFs(:,:,ipeak))
end
suptitle('scaling factors')
saveFolder = loadFolder;
if bestSFflag
    save([saveFolder 'scalingFactors_bestSF'],'SFs','bestSFflag','sigmas','taus','whichpeaks')                   
else
    save([saveFolder 'scalingFactors_constSF'],'SFs','bestSFflag','sigmas','taus','whichpeaks')
end
%% plot on top of bargraph N1 data
bestSFflag = false;

RAdate = '22-May-2019';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])

if bestSFflag
    load([loadFolder 'scalingFactors_bestSF'])
else
    load([loadFolder 'scalingFactors_constSF'])
end
whichpeaks = {'N1','P2'};
peak = 'N1';
ipeak = find(strcmp(whichpeaks,peak));
tol = 0.0001;

[xlocs,N1_means] = BarPlotN1_MMNchroma(peak);
title([peak ' data (bars) versus model (lines)'])
%notice - peak_means dim is blocks X notes whereas ERR is subj X notes X
%blocks. For historical reasons.
meanV = mean(mean(N1_means));
maxV = max(max(abs(N1_means)));
Predict = squeeze(nanmean((1-RApeaks)));
PredictU = squeeze(nanmean(Predict,3));
Psigs = [1,6,10];
Pts = [1,3,10];
nComb=numel(Psigs)*numel(Pts);
cmap = zeros(numel(Psigs)*numel(Pts),3);
cmap(1:ceil(numel(Psigs)*numel(Pts)/2),1) = linspace(0,1,ceil(numel(Psigs)*numel(Pts)/2));
cmap(ceil(numel(Psigs)*numel(Pts)/2):end,1) = ones(floor(numel(Psigs)*numel(Pts)/2)+1,1);
cmap(ceil(numel(Psigs)*numel(Pts)/2):end,2) = linspace(0,0.8,floor(numel(Psigs)*numel(Pts)/2)+1);
cmap(ceil(numel(Psigs)*numel(Pts)/2):end,3) = linspace(0,0.8,floor(numel(Psigs)*numel(Pts)/2)+1);

hold on
bls=[2,4];
Bnames = blocks(bls);
ibl=0;
for bl=1:length(Bnames)
    ibl=ibl+1;
    xind=((bl-1)*5+1):(bl*5);
    ist=0;cbarTicks=cell(1,nComb);
    for sig = Psigs
        for t = Pts
            ist = ist+1;
            SF = SFs(sig==sigmas,abs(t-taus)<tol,ipeak);
            RAs = squeeze(nanmean(nanmean(RApeaks(whichSubjects,ibl,:,:,sig==sigmas,abs(t-taus)<tol)),4));
            Predict = 1-RAs;
            cbarTicks{ist}=['Sig=' num2str(sig) ', tau=' num2str(t)];
            plot(xlocs(xind),Predict'.*SF,'Color',cmap(ist,:),'linew',2)
        end
    end
end
ax=subplot(10,10,[60,70,80,90,100]);
colormap(ax,cmap);c=colorbar;
set(ax,'visible','off')
set(c,'Ticks',1/nComb:1/nComb:1)
set(c,'TickLabels',cbarTicks)
set(c,'Position',[0.75,0.1,0.02,0.5])
set(c,'YAxisLocation','right')
set(c,'fontsize',14)
%% calc and plot model error
weightedFlag = false;
bestSFflag = true;

%load Model
%RAdate = '10-Apr-2019';
RAdate = '22-May-2019';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
if bestSFflag
    load([loadFolder 'scalingFactors_bestSF'])
else
    load([loadFolder 'scalingFactors_constSF'])
end

%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '18-Nov-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));

Errs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks));
bls=[2,4];
for ipeak = 1:length(whichpeaks)
    ibl=0;
    for bl=bls
        ibl=ibl+1;
        peakAmps(:,ibl,:,:) = allPeak_amps{bl}(whichSubjects,:,:,ipeak);
        %subj x bl x con x prevcon 
%         peakAmps(:,bl,:) = allGrandcon_amps{bl}(whichSubjects,:,ipeak);
    end
    ste = nanstd(peakAmps)/sqrt(size(peakAmps,1));
    ste = ste(:);ste(isnan(ste))=[];
    data = nanmean(peakAmps);
    data = data(:);data(isnan(data))=[];
    for isig = 1:length(sigmas)
        for itau = 1:length(taus)
            model = 1-nanmean(RApeaks(whichSubjects,:,:,:,isig,itau));
            model = model(:);model(isnan(model))=[];
            SF = SFs(isig,itau,ipeak);
            %SF = SFs(1,1,ipeak);
            if weightedFlag
                Errs(isig,itau,ipeak) = rms((data - model*SF)./ste);%performs SS linear regression
            else
                %this was a mistake:
                Errs(isig,itau,ipeak) = rms((data - model*SF))/mean(ste);%performs SS linear regression
                %Errs(isig,itau,ipeak) = rms((data - model*SF));%performs SS linear regression
            end
        end
    end
end
saveFolder = loadFolder;
if bestSFflag
    save([saveFolder 'Errs_bestSF'],'Errs','bestSFflag','sigmas','taus','whichpeaks')                   
else
    save([saveFolder 'Errs_constSF'],'Errs','bestSFflag','sigmas','taus','whichpeaks')                   
end

ERPfigure;
set(gcf,'Position', [50 50 1200 450])
errmin = min(min(min(Errs(:,:,:,:))));
errmax = max(max(max(Errs(:,:,:,:))));
bestaus = nan(2,1);
bestsigmas = nan(2,1);
for ipeak=1:length(whichpeaks)
    ax(ipeak) = subplot(1,2,ipeak);
    [row col]=find(Errs(:,:,ipeak)==min(min(Errs(:,:,ipeak))));
    imagesc(taus,sigmas,Errs(:,:,ipeak));
    %caxis([errmin errmax])
     if exist('bestSFflag','var')
        if ~bestSFflag
            caxis([1 2])
        end
    end
    title(whichpeaks{ipeak})
    xlabel('\tau (sec)','fontsize',20);
    ylabel('\sigma (semitones)','fontsize',20)
    hold on
    disp([whichpeaks{ipeak} ': sigma = ' num2str(sigmas(row)) ', tau = ' num2str(taus(col)) ])
    plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',1)
    plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',1)
    set(gca,'fontsize',20)
    bestaus(ipeak) = taus(col);
    bestsigmas(ipeak) = sigmas(row);
    text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',20);
end
suptitle('normalized error between model and data')
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
end
if bestSFflag
    save([saveFolder filesep 'bestFitParams_bestSF'],'bestaus','bestsigmas')
else
    save([saveFolder filesep 'bestFitParams_constSF'],'bestaus','bestsigmas')
end
% saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'fig')
% saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'jpg')
%% calc and plot model error - contours

%load Model
%RAdate = '08-Apr-2019';
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
load([loadFolder 'scalingFactors'])                   

%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));

Errs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks));

for ipeak = 1:length(whichpeaks)
    for bl=1:length(allPeak_amps)
        peakAmps(:,bl,:,:) = allPeak_amps{bl}(whichSubjects,:,:,ipeak);
        %subj x bl x con x prevcon 
%         peakAmps(:,bl,:) = allGrandcon_amps{bl}(whichSubjects,:,ipeak);
    end
    ste = nanstd(peakAmps)/sqrt(size(peakAmps,1));
    ste = ste(:);ste(isnan(ste))=[];
    data = nanmean(peakAmps);
    data = data(:);data(isnan(data))=[];
    for isig = 1:length(sigmas)
        for itau = 1:length(taus)
            model = 1-nanmean(RApeaks(whichSubjects,:,:,:,isig,itau));
            model = model(:);model(isnan(model))=[];
            SF = SFs(isig,itau,ipeak);
            Errs(isig,itau,ipeak) = rms((data - model*SF)/mean(ste));%performs SS linear regression
        end
    end
end
saveFolder = loadFolder;
save([saveFolder 'Errs'],'Errs','sigmas','taus','whichpeaks')                   


ERPfigure;
lastaupos = 25;
set(gcf,'Position', [50 50 1200 450])
errmin = min(min(min(Errs(:,1:lastaupos,:,:))));
errmax = max(max(max(Errs(:,1:lastaupos,:,:))));
for ipeak=1:length(whichpeaks)
    ax(ipeak) = subplot(1,2,ipeak);
    [row, col]=find(Errs(:,1:lastaupos,ipeak)==min(min(Errs(:,1:lastaupos,ipeak))));
    [X, Y]=meshgrid(taus(1:lastaupos),sigmas);
    contourf(X,Y,Errs(:,1:lastaupos,ipeak),[1.1:0.05:1.5 1.6:0.1:2.3],'linewidth',2);
    caxis([errmin errmax])
    
    title(whichpeaks{ipeak})
    xlabel('\tau (sec)','fontsize',20);
    ylabel('\sigma (semitones)','fontsize',20)
    hold on
    disp([whichpeaks{ipeak} ': sigma = ' num2str(sigmas(row)) ', tau = ' num2str(taus(col)) ])
    plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 1 1],'linewidth',3)
    plot([taus(1) taus(lastaupos)],[sigmas(row) sigmas(row)],'Color',[1 1 1],'linewidth',3)
    set(gca,'fontsize',14)
end
suptitle('normalized error between model and data')
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
end
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'fig')
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'jpg')
%saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'eps')
%set(0, 'DefaultFigureRenderer', 'painters');
print(gcf,'-depsc','-painters',[AnalysisFolder 'Figures\model\norm_err_model_data_' date '.eps']);

epsclean([AnalysisFolder 'Figures\model\norm_err_model_data_' date '.eps'],[AnalysisFolder 'Figures\model\norm_err_model_data_' date '_clean.eps']);
epsclean([AnalysisFolder 'Figures\model\norm_err_model_data_' date '.eps'],[AnalysisFolder 'Figures\model\norm_err_model_data_' date '_clean_combineAreas.eps'],'combineAreas',true)

%% fitting parameters to each block type separately -
%% find best scaling factor (each block type)

%load Model
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])

%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
peakAmps = nan(size(squeeze(RApeaks(whichSubjects,1,:,:,1,1))));
        
SFs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks),size(RApeaks,4));%sigmas, taus, peaks, blockTypes
    
for ipeak = 1:length(whichpeaks)
    for bl=1:length(allPeak_amps)
        peakAmps(:,:,:) = allPeak_amps{bl}(whichSubjects,:,:,ipeak);
        %subj x con x prevcon 

        data = squeeze(nanmean(peakAmps(:,:,:)));
        data = data(:);data(isnan(data))=[];
        for isig = 1:length(sigmas)
            for itau = 1:length(taus)
                model = squeeze(1-nanmean(RApeaks(whichSubjects,bl,:,:,isig,itau)));
                model = model(:);model(isnan(model))=[];
                SFs(isig,itau,ipeak,bl) = model\data;%performs SS linear regression
            end
        end
    end
end
ERPfigure
ii=0;
for ipeak=1:length(whichpeaks)
    for bl=1:length(allPeak_amps)
        ii=ii+1;
        subplot(2,5,ii)
        imagesc(SFs(:,:,ipeak,bl))
    end
end
suptitle('scaling factors')

saveFolder = loadFolder;
save([saveFolder 'scalingFactors_blockTypes'],'SFs','sigmas','taus','whichpeaks')                   

%% calc and plot model error - individual blocks

SF_block = true;

%load Model
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
if SF_block
    load([loadFolder 'scalingFactors_blockTypes'])                   
else
    load([loadFolder 'scalingFactors'])
end
%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
%peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));
peakAmps = nan(size(squeeze(RApeaks(whichSubjects,1,:,:,1,1))));
        
%SFs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks),size(RApeaks,4));%sigmas, taus, peaks, blockTypes

Errs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks),size(RApeaks,4));%sigmas, taus, peaks, blockTypes

for ipeak = 1:length(whichpeaks)
    for bl=1:length(allPeak_amps)
        peakAmps(:,:,:) = allPeak_amps{bl}(whichSubjects,:,:,ipeak);
        %subj x bl x con x prevcon 
%         peakAmps(:,bl,:) = allGrandcon_amps{bl}(whichSubjects,:,ipeak);
        ste = nanstd(peakAmps)/sqrt(size(peakAmps,1));
        ste = ste(:);ste(isnan(ste))=[];
        data = squeeze(nanmean(peakAmps));
        data = data(:);data(isnan(data))=[];
        
        for isig = 1:length(sigmas)
            for itau = 1:length(taus)
                model = squeeze(1-nanmean(RApeaks(whichSubjects,bl,:,:,isig,itau)));
                model = model(:);model(isnan(model))=[];
                if SF_block
                    SF = SFs(isig,itau,ipeak,bl);
                else
                    SF = SFs(isig,itau,ipeak);
                end
                Errs(isig,itau,ipeak,bl) = rms((data - model*SF)/mean(ste));%performs SS linear regression
            end
        end
    end
end
saveFolder = loadFolder;
save([saveFolder 'Errs_blocks'],'Errs','sigmas','taus','whichpeaks')                   

ERPfigure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
ii=0;
for ipeak=1:length(whichpeaks)
   % errmin = min(min(min(Errs(:,1:25,ipeak,:))));
   % errmax = max(max(max(Errs(:,1:25,ipeak,:))));
    for bl=1:length(allPeak_amps)
        ii=ii+1;        
        ax(ii) = subplot(2,5,ii);
        [row col]=find(Errs(:,1:25,ipeak,bl)==min(min(Errs(:,1:25,ipeak,bl))));
        imagesc(taus(1:25),sigmas,Errs(:,1:25,ipeak,bl));
        hold on 
        plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 1 1],'linewidth',3)
        plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 1 1],'linewidth',3)
        title(whichpeaks{ipeak})
        set(gca,'fontsize',12)
        xlabel('\tau (sec)','fontsize',14);
        if bl==1
            ylabel('\sigma (semitones)','fontsize',14)
        end
        %pause(1)
        %c(ii) = colorbar(ax(ii));
    end
    %ylabel(c(ii),'STE units')
end
suptitle('normalized error between model and data')
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_bType_' date],'fig')
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_bType_' date],'jpg')

%% fitting parameters to each ~~condition~~ separately -
%% find best scaling factor (each condition)
spreads = 1:3;%these are the 3 conditions
block2cond = [1,2,2,3,3];

%load Model
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])

%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};

SFs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks),length(spreads));%sigmas, taus, peaks, conditions
    
for ipeak = 1:length(whichpeaks)
    for ispread=1:length(spreads)
        bls = find(block2cond == ispread);
        peakAmps = nan(size((RApeaks(whichSubjects,1:length(bls),:,:,1,1))));

        for ib=1:length(bls)
            peakAmps(:,ib,:,:) = allPeak_amps{bls(ib)}(whichSubjects,:,:,ipeak);
        end

        data = squeeze(nanmean(peakAmps));
        data = data(:);data(isnan(data))=[];
        for isig = 1:length(sigmas)
            for itau = 1:length(taus)
                model = nan(length(bls),size(RApeaks,3),size(RApeaks,4));
                for ib=1:length(bls)
                    model(ib,:,:) = squeeze(1-nanmean(RApeaks(whichSubjects,bls(ib),:,:,isig,itau)));         
                end
                model = model(:);model(isnan(model))=[];
                SFs(isig,itau,ipeak,ispread) = model\data;%performs SS linear regression
            end
        end
    end
end
ERPfigure
ii=0;
for ipeak=1:length(whichpeaks)
    for ispread=1:length(spreads)
        ii=ii+1;
        subplot(2,3,ii)
        imagesc(SFs(:,:,ipeak,ispread))
    end
end
suptitle('scaling factors')

saveFolder = loadFolder;
save([saveFolder 'scalingFactors_spreads'],'SFs','sigmas','taus','whichpeaks')                   

%% calc and plot model error - individual conditions

SF_cond = true;

%load Model
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
if SF_cond
    load([loadFolder 'scalingFactors_spreads'])                   
else
    load([loadFolder 'scalingFactors'])
end
%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
%peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));
peakAmps = nan(size(squeeze(RApeaks(whichSubjects,1,:,:,1,1))));
        
%SFs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks),size(RApeaks,4));%sigmas, taus, peaks, blockTypes

Errs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks),length(spreads));%sigmas, taus, peaks, blockTypes

for ipeak = 1:length(whichpeaks)
    for ispread=1:length(spreads)

        bls = find(block2cond == ispread);
        peakAmps = nan(size((RApeaks(whichSubjects,1:length(bls),:,:,1,1))));

        for ib=1:length(bls)
            peakAmps(:,ib,:,:) = allPeak_amps{bls(ib)}(whichSubjects,:,:,ipeak);
        end

        ste = nanstd(peakAmps)/sqrt(size(peakAmps,1));
        ste = ste(:);ste(isnan(ste))=[];
        data = squeeze(nanmean(peakAmps));
        data = data(:);data(isnan(data))=[];
        
        for isig = 1:length(sigmas)
            for itau = 1:length(taus)
                model = nan(length(bls),size(RApeaks,3),size(RApeaks,4));
                for ib=1:length(bls)
                    model(ib,:,:) = squeeze(1-nanmean(RApeaks(whichSubjects,bls(ib),:,:,isig,itau)));         
                end
                model = model(:);model(isnan(model))=[];
                if SF_cond
                    SF = SFs(isig,itau,ipeak,ispread);
                else
                    SF = SFs(isig,itau,ipeak);
                end
                Errs(isig,itau,ipeak,ispread) = rms((data - model*SF)/mean(ste));%performs SS linear regression
            end
        end
    end
end
saveFolder = loadFolder;
save([saveFolder 'Errs_spreads'],'Errs','sigmas','taus','whichpeaks')                   

bestParams = nan(length(spreads),length(whichpeaks),2);
ERPfigure;
set(gcf,'units','normalized','outerposition',[0 0 0.7 0.7])
ii=0;
errmin = min(min(min(Errs(:,1:25,ipeak,:))));
errmax = max(max(max(Errs(:,1:25,ipeak,:))));
   
for ipeak=1:length(whichpeaks)
   for ispread=1:length(spreads)
        ii=ii+1;        
        ax(ii) = subplot(2,3,ii);
        [row col]=find(Errs(:,1:25,ipeak,ispread)==min(min(Errs(:,1:25,ipeak,ispread))));
        imagesc(taus(1:25),sigmas,Errs(:,1:25,ipeak,ispread));
        caxis([errmin errmax])
        hold on 
        plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 1 1],'linewidth',3)
        plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 1 1],'linewidth',3)
        bestParams(ispread,ipeak,1)=sigmas(row);
        bestParams(ispread,ipeak,2)=taus(col);
        title({whichpeaks{ipeak},['Condition ' num2str(ispread)]})
        set(gca,'fontsize',12)
        if ipeak ==2
            xlabel('\tau (sec)','fontsize',14);
        end
        if ispread==1
            ylabel('\sigma (semitones)','fontsize',14)
        end
        %pause(1)
        %c(ii) = colorbar(ax(ii));
    end
    c(ipeak) = colorbar;
    ylabel(c(ipeak),'STE units','fontsize',14)
    if ipeak ==1
        set(c(ipeak),'Position',[0.92 0.562 0.02 0.323])
    else
        set(c(ipeak),'Position',[0.92 0.108 0.02 0.323])
    end
end
suptitle('Normalized error between model and data - fitting each conditions separately')
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_spreads_' date],'fig')
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_spreads_' date],'jpg')
%%  
freqspreads = [[10,10,9,10];[7 7 7 7];[4,4,3,3]];
meanfreqspreads = mean(freqspreads');
ERPfigure
hold on
plot(meanfreqspreads,bestParams(:,1,1),'b','linewidth',2)
plot(meanfreqspreads,bestParams(:,2,1),'r','linewidth',2)
set(gca,'fontsize',20)
legend({'N1','P2'},'fontsize',20,'Location','nw')
plot(meanfreqspreads,bestParams(:,1,1),'.b','MarkerSize',50)
plot(meanfreqspreads,bestParams(:,2,1),'.r','MarkerSize',50)
title('\sigma as a function of frequency spread','fontsize',20)
ylabel('Best fit \sigma (semitones)','fontsize',20)
xlabel('Mean spread in the sequence (semitones)','fontsize',20)

%% Fit Model for each participant separately
% April 2019
% calc. expected peaks
%RAdate = '19-Dec-2018';
RAdate = '10-Apr-2019';
loadFolder = [modelFolder RAdate filesep];
saveFolder = [loadFolder 'individual' filesep];
mkdir(saveFolder)

tic
fprintf(['Loading...'])
load([loadFolder 'Metadata'])
load([loadFolder 'Params'])
load([loadFolder 'RApeaks'])
fprintf(['Done in %4.1f sec \n'], toc)

%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));
        
SFs = nan(length(whichSubjects),size(RApeaks,5),size(RApeaks,6),length(whichpeaks));
    
for ipeak = 1:length(whichpeaks)
    for bl=1:length(allPeak_amps)
        peakAmps(:,bl,:,:) = allPeak_amps{bl}(whichSubjects,:,:,ipeak);
        %subj x bl x con x prevcon 
%         peakAmps(:,bl,:) = allGrandcon_amps{bl}(whichSubjects,:,ipeak);
    end
    for s=1:length(whichSubjects)
        data = peakAmps(s,:,:,:);
        data = data(:);data(isnan(data))=[];
        for isig = 1:length(sigmas)
            for itau = 1:length(taus)
                model = 1-(RApeaks(whichSubjects(s),:,:,:,isig,itau));
                model = model(:);model(isnan(model))=[];
                SFs(s,isig,itau,ipeak) = model\data;%performs SS linear regression
            end
        end
    end
end
for ipeak=1:length(whichpeaks)
    ERPfigure
    for s=1:length(whichSubjects)
        subplot(6,6,s)
        imagesc(squeeze(SFs(s,:,:,ipeak)))
        title(['Subject ' num2str(whichSubjects(s))])
    end
    suptitle(['scaling factors ' whichpeaks{ipeak}])
end
save([saveFolder 'scalingFactors'],'SFs','sigmas','taus','whichpeaks')                   

%%%%%%%%%%%% calc model error

peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));
Errs = nan(length(whichSubjects),size(RApeaks,5),size(RApeaks,6),length(whichpeaks));
for ipeak = 1:length(whichpeaks)
    for bl=1:length(allPeak_amps)
        peakAmps(:,bl,:,:) = allPeak_amps{bl}(whichSubjects,:,:,ipeak);
        %subj x bl x con x prevcon 
%         peakAmps(:,bl,:) = allGrandcon_amps{bl}(whichSubjects,:,ipeak);
    end
    for s=1:length(whichSubjects)
%         ste = nanstd(peakAmps)/sqrt(size(peakAmps,1));
%         ste = nanstd(peakAmps(s,:,:,:))/sqrt(size(peakAmps(s,:,:,:),1));
%         ste = ste(:);ste(isnan(ste))=[];
        data = peakAmps(s,:,:,:);
        data = data(:);data(isnan(data))=[];
        for isig = 1:length(sigmas)
            for itau = 1:length(taus)
                model = 1-RApeaks(whichSubjects(s),:,:,:,isig,itau);
                model = model(:);model(isnan(model))=[];
                SF = SFs(s,isig,itau,ipeak);
                Errs(s,isig,itau,ipeak) = rms(data - model*SF);
            end
        end
    end
end
save([saveFolder 'Errs'],'Errs','sigmas','taus','whichpeaks')                   
%% plot error spaces
plotflag = false;
ERPfigure;
lastsigpos = length(sigmas);
lastaupos = length(taus);
set(gcf,'Position', [50 50 1200 450])
errmin = min(min(min(min(Errs(:,1:lastsigpos,1:lastaupos,:)))));
errmax = max(max(max(max(Errs(:,1:lastsigpos,1:lastaupos,:)))));
bestsigmas = nan(length(whichSubjects),2);
bestaus = nan(length(whichSubjects),2);
minErrs = nan(length(whichSubjects),2);
meanErrs = nan(length(whichSubjects),2);

for s=1:length(whichSubjects)
    for ipeak=1:length(whichpeaks)
        hold off
        ax(ipeak) = subplot(1,2,ipeak);
        [row col]=find(squeeze(Errs(s,1:lastsigpos,1:lastaupos,ipeak))==min(min(Errs(s,1:lastsigpos,1:lastaupos,ipeak))));
        bestsigmas(s,ipeak) = sigmas(row);
        bestaus(s,ipeak) = taus(col);
        %[X, Y]=meshgrid(taus(1:lastaupos),sigmas);
        %contourf(X,Y,squeeze(Errs(s,1:lastsigpos,1:lastaupos,ipeak)),[1.1:0.05:1.5 1.6:0.1:2.3],'linewidth',2);
        imagesc(taus,sigmas,squeeze(Errs(s,1:lastsigpos,1:lastaupos,ipeak)));
        %caxis([errmin errmax])
        minErrs(s,ipeak) = squeeze(Errs(s,row,col,ipeak));
        meanErrs(s,ipeak) = mean(mean(Errs(s,:,:,ipeak)));
        title(whichpeaks{ipeak})
        xlabel('\tau (sec)','fontsize',20);
        ylabel('\sigma (semitones)','fontsize',20)
        hold on
        disp(['Subj ' num2str(s) ' ' whichpeaks{ipeak} ': sigma = ' num2str(sigmas(row)) ', tau = ' num2str(taus(col)) ])
        plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 1 1],'linewidth',3)
        plot([taus(1) taus(lastaupos)],[sigmas(row) sigmas(row)],'Color',[1 1 1],'linewidth',3)
        set(gca,'fontsize',14)
    end
    suptitle(['Subj ' num2str(whichSubjects(s))])
    for ipeak = 1:length(whichpeaks)
        c(ipeak) = colorbar(ax(ipeak));
        ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
    end
    if plotflag
        pause
    end
end
save([saveFolder 'individual_bestFits'],'bestsigmas','bestaus','minErrs','meanErrs')
disp('done')
%saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'fig')
%saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'jpg')
%saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'eps')
%set(0, 'DefaultFigureRenderer', 'painters');
%print(gcf,'-depsc','-painters',[AnalysisFolder 'Figures\model\norm_err_model_data_' date '.eps']);

%epsclean([AnalysisFolder 'Figures\model\norm_err_model_data_' date '.eps'],[AnalysisFolder 'Figures\model\norm_err_model_data_' date '_clean.eps']);
%epsclean([AnalysisFolder 'Figures\model\norm_err_model_data_' date '.eps'],[AnalysisFolder 'Figures\model\norm_err_model_data_' date '_clean_combineAreas.eps'],'combineAreas',true)
%% plot best fits
ERPfigure
datas = {bestsigmas,bestaus};
strings = {'sigma','tau'};
unitstring = {'semitones','seconds'};
for i=1:2
    subplot(1,2,i)
    data = datas{i};
    h=boxplot(data,'notch','on','Labels',{'N1','P2'});
    set(h,{'linew'},{2});
    ylabel(['\' strings{i} ' (' unitstring{i} ')'])
    set(gca,'fontsize',16)
    hold on
    X=ones(length(whichSubjects),2);X(:,2)=X(:,2)+1;
    Xjit=X+randn(size(X))/10;
    plot(Xjit',data','o','markersize',6,'Color','k')
    title(['best fit \' strings{i}])
    %plot([1,2],[mean(data(:,1)),mean(data(:,2))],'g','linewidth',2)
    [p,H]=signrank(data(:,2),data(:,1));
    text(1,9,['Wilcoxon p=' num2str(p,2)],'fontsize',14,'Color','g')
end
%% investigate extreme participants
XSi=find(data(:,1)>10);
XS=whichSubjects(XSi);
sigmasX=datas{1}(XSi,:);
tausX=datas{2}(XSi,:);
save([saveFolder 'XtremeS'],'XS','XSi')

%% plot individual participants data (extremes)
%comes from:
% plot on top of bargraph N1 data
RAdate = '10-Apr-2019';
loadFolder = [modelFolder RAdate filesep];
saveFolder = [loadFolder 'individual' filesep];
load([saveFolder 'XtremeS']);
%tic
%fprintf(['Loading...'])
%load([loadFolder 'Metadata'])
%load([loadFolder 'Params'])
%load([loadFolder 'RApeaks'])
%fprintf(['Done in %4.1f sec \n'], toc)
whichpeaks = {'N1','P2'};

peak = 'N1';
ipeak = find(strcmp(whichpeaks,peak));

plotflag = true;
[xlocs,N1_means,N1_errs,bpeaks] = BarPlotN1_individuals(plotflag,peak);
title([peak ' grand data (bars) versus individual participants (lines)'])
%notice - peak_means dim is blocks X notes whereas ERR is subj X notes X
%blocks. For historical reasons.
hold on
%XtremeS
for bl=1:length(Bnames)
    xind=((bl-1)*length(Bnames)+1):(bl*length(Bnames));
    plot(xlocs(xind),bpeaks(:,xind),'Color',[0.5 0.5 0.5],'linew',1)
    plot(xlocs(xind),bpeaks(XSi,xind),'Color','r','linew',1)
end
allSi = 1:length(whichSubjects);
nonXSi=allSi(~ismember(allSi,XSi));
diffsInd = bpeaks-repmat(mean(bpeaks),size(bpeaks,1),1);
noise = rms(diffsInd');

ERPfigure
for si=1:length(whichSubjects)
    subplot(6,6,si)
    h = barwitherr(N1_errs, N1_means);% Plot with errorbar
    ylim([-5,2])
    hold on
    for bl=1:length(Bnames)
        xind=((bl-1)*length(Bnames)+1):(bl*length(Bnames));
        if ismember(si,XSi)
            plot(xlocs(xind),bpeaks(si,xind),'Color','r','linew',2)
        else
            plot(xlocs(xind),bpeaks(si,xind),'Color','g','linew',2)
        end
    end    
end

figure
subplot(1,3,1)
data=noise;
boxplot([data(XSi) data(nonXSi)],[zeros(size(data(nonXSi))),ones(size(data(XSi)))])
hold on
plot(ones(size(data(XSi))),data(XSi),'.','Markersize',16)
hold on
plot(2*ones(size(data(nonXSi))),data(nonXSi),'.','Markersize',16)
title(['noise (rms diff from mean)'])

load([saveFolder 'individual_bestFits'])
ipeak = 1;%N1
subplot(1,3,2)
data=minErrs(:,ipeak)';
boxplot([data(XSi) data(nonXSi)],[zeros(size(data(nonXSi))),ones(size(data(XSi)))])
hold on
plot(ones(size(data(XSi))),data(XSi),'.','Markersize',16)
hold on
plot(2*ones(size(data(nonXSi))),data(nonXSi),'.','Markersize',16)
title(['min Errs'])

subplot(1,3,3)
data=meanErrs(:,ipeak)';
boxplot([data(XSi) data(nonXSi)],[zeros(size(data(nonXSi))),ones(size(data(XSi)))])
hold on
plot(ones(size(data(XSi))),data(XSi),'.','Markersize',16)
hold on
plot(2*ones(size(data(nonXSi))),data(nonXSi),'.','Markersize',16)
title(['mean Errs'])

%% %%%%%%%Fit model to each pariticipant and condition separately!
%
%% find best scaling factor (each condition, participant)
spreads = 1:3;%these are the 3 conditions
block2cond = [1,2,2,3,3];

%load Model
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
saveFolder = [loadFolder 'individual' filesep];

load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])

%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};

%SFs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks),length(spreads));%sigmas, taus, peaks, conditions
SFs = nan(length(whichSubjects),size(RApeaks,5),size(RApeaks,6),length(whichpeaks),length(spreads));
%conditions:
for ipeak = 1:length(whichpeaks)
    for ispread=1:length(spreads)
        bls = find(block2cond == ispread);
        peakAmps = nan(size((RApeaks(whichSubjects,1:length(bls),:,:,1,1))));

        for ib=1:length(bls)
            peakAmps(:,ib,:,:) = allPeak_amps{bls(ib)}(whichSubjects,:,:,ipeak);
        end
        for s=1:length(whichSubjects)
            data = peakAmps(s,:,:,:);
            data = data(:);data(isnan(data))=[];
            for isig = 1:length(sigmas)
                for itau = 1:length(taus)
                    model = nan(length(bls),size(RApeaks,3),size(RApeaks,4));
                    for ib=1:length(bls)
                        model(ib,:,:) = 1-RApeaks(whichSubjects(s),bls(ib),:,:,isig,itau);
                    end
                    model = model(:);model(isnan(model))=[];
                    SFs(s,isig,itau,ipeak,ispread) = model\data;%performs SS linear regression
                end
            end
        end
    end
end

for ipeak=1:length(whichpeaks)
    ii=0;
    ERPfigure
    for ispread=1:length(spreads)
        for s=1:length(whichSubjects)
            ii=ii+1;
            subplot(3,length(whichSubjects),ii)
            imagesc(squeeze(SFs(s,:,:,ipeak,ispread)))
        end
    end
    suptitle(['scaling factors ' whichpeaks{ipeak}])
end

save([saveFolder 'scalingFactors_spreads'],'SFs','sigmas','taus','whichpeaks')                   

%% calc and plot model error - individual conditions, participant

SF_cond = true;

%load Model
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
if SF_cond
    load([loadFolder 'individual' filesep 'scalingFactors_spreads'])                   
else
    load([loadFolder 'individual' filesep 'scalingFactors'])
end
%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
%peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));
peakAmps = nan(size(squeeze(RApeaks(whichSubjects,1,:,:,1,1))));
        
%SFs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks),size(RApeaks,4));%sigmas, taus, peaks, blockTypes

Errs = nan(length(whichSubjects),size(RApeaks,5),size(RApeaks,6),length(whichpeaks),length(spreads));%sigmas, taus, peaks, blockTypes

for ipeak = 1:length(whichpeaks)
    for ispread=1:length(spreads)
        bls = find(block2cond == ispread);
        peakAmps = nan(size((RApeaks(whichSubjects,1:length(bls),:,:,1,1))));

        for ib=1:length(bls)
            peakAmps(:,ib,:,:) = allPeak_amps{bls(ib)}(whichSubjects,:,:,ipeak);
        end
        for s=1:length(whichSubjects)
%             ste = nanstd(peakAmps)/sqrt(size(peakAmps,1));
%             ste = ste(:);ste(isnan(ste))=[];
            data = peakAmps(s,:,:,:);
            data = data(:);data(isnan(data))=[];
            for isig = 1:length(sigmas)
                for itau = 1:length(taus)
                    model = nan(length(bls),size(RApeaks,3),size(RApeaks,4));
                    for ib=1:length(bls)
                        model(ib,:,:) = 1-RApeaks(whichSubjects(s),bls(ib),:,:,isig,itau);         
                    end
                    model = model(:);model(isnan(model))=[];
                    if SF_cond
                        SF = SFs(s,isig,itau,ipeak,ispread);
                    else
                        SF = SFs(s,isig,itau,ipeak);
                    end
                    Errs(s,isig,itau,ipeak,ispread) = rms(data - model*SF);%performs SS linear regression
                end
            end
        end
    end
end

save([saveFolder 'Errs_spreads'],'Errs','sigmas','taus','whichpeaks')                   

bestParams = nan(length(whichSubjects),length(spreads),length(whichpeaks),2);
%plot error spaces
for s=1:length(whichSubjects)
    ERPfigure;
    set(gcf,'units','normalized','outerposition',[0 0 0.7 0.7])
    ii=0;
    errmin = min(min(min(min(min(Errs(s,:,:,ipeak,:))))));
    errmax = max(max(max(max(max(Errs(s,:,:,ipeak,:))))));
   
    for ipeak=1:length(whichpeaks)
       for ispread=1:length(spreads)
            ii=ii+1;        
            ax(ii) = subplot(2,3,ii);
            [row col]=find(squeeze(Errs(s,:,:,ipeak,ispread))==min(min(Errs(s,:,:,ipeak,ispread))));
            imagesc(taus(1:25),sigmas,squeeze(Errs(s,:,:,ipeak,ispread)));
%             caxis([errmin errmax])
            hold on 
            plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 1 1],'linewidth',3)
            plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 1 1],'linewidth',3)
            bestParams(s,ispread,ipeak,1)=sigmas(row);
            bestParams(s,ispread,ipeak,2)=taus(col);
            title({whichpeaks{ipeak},['Condition ' num2str(ispread)]})
            set(gca,'fontsize',12)
            if ipeak ==2
                xlabel('\tau (sec)','fontsize',14);
            end
            if ispread==1
                ylabel('\sigma (semitones)','fontsize',14)
            end
            %pause(1)
            %c(ii) = colorbar(ax(ii));
        end
        c(ipeak) = colorbar;
        ylabel(c(ipeak),'STE units','fontsize',14)
        if ipeak ==1
            set(c(ipeak),'Position',[0.92 0.562 0.02 0.323])
        else
            set(c(ipeak),'Position',[0.92 0.108 0.02 0.323])
        end
    end
    suptitle(['Subj ' num2str(s) '. Normalized error - each conditions separately'])
    % saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_spreads_' date],'fig')
    % saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_spreads_' date],'jpg')
end
%%
freqspreads = [[10,10,9,10];[7 7 7 7];[4,4,3,3]];
meanfreqspreads = mean(freqspreads');
%bestParams: s x ispread x ipeak x sigortau
errorbars = nan(length(meanfreqspreads),2);
for ipeak = 1:2
    for ispread = 1:length(meanfreqspreads) 
        %CIs = Confidence(bestParams(:,ispread,ipeak,1));
        %errorbars(ispread,ipeak) = mean(bestParams(:,ispread,ipeak,1)) - CIs(1);
        x = bestParams(:,ispread,ipeak,1);
        errorbars(ispread,ipeak) = nanstd(x)/sqrt(length(x(~isnan(x))));
    end
end
ERPfigure
hold on
plot(meanfreqspreads,mean(bestParams(:,:,1,1)),'b','linewidth',2)
plot(meanfreqspreads,mean(bestParams(:,:,2,1)),'r','linewidth',2)
set(gca,'fontsize',20)
legend({'N1','P2'},'fontsize',20,'Location','nw')
plot(meanfreqspreads,mean(bestParams(:,:,1,1)),'.b','MarkerSize',50)
plot(meanfreqspreads,mean(bestParams(:,:,2,1)),'.r','MarkerSize',50)
errorbar(meanfreqspreads,mean(bestParams(:,:,1,1)),errorbars(:,1),'b')
errorbar(meanfreqspreads,mean(bestParams(:,:,2,1)),errorbars(:,1),'r')

title('\sigma as a function of frequency spread','fontsize',20)
ylabel('Best fit \sigma (semitones)','fontsize',20)
xlabel('Mean spread in the sequence (semitones)','fontsize',20)

