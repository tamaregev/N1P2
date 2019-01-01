%The Model as in Herrmann et al.
% Nov 25 2018, Tamar Regev
%% Definitions
cd('L:\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
sbjcts = 1:37;
mss = [5,14];
whichSubjects = sbjcts(~ismember(sbjcts,[badSubjects,mss]));

grandFolder = [AnalysisFolder 'grandAverage'];
mixedFolder = [AnalysisFolder 'MixedModel\'];
modelFolder = [AnalysisFolder 'Model\'];

srate = 512;
Expinfo.srate = srate;

%codes and names table:
Bnames = {'1','2a','2b','3a','3b'}';
Bcodes = [10,20,30,40,50]';
Bcodes_names = table(Bnames,Bcodes);
bls=1:5;

%% compare EDAT to events

s=4;
%load expdata and markers
FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
load([EDATfolder FileName '_expdata.mat' ])
%load and read markers
VMRKfile = [ExportFolder FileName  '_dt_RDI_imported.vmrk'];
[eventCoInd, ~]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed
stimCoInd = eventCoInd(~ismember(eventCoInd(:,1),[210,200,100,110,254]),:);
trials = expdata.trials;
phaseName = 'Passive';
block_list = expdata.(phaseName).block_list;
nBlocks = length(block_list);
if(~length(trials)==length(stimCoInd))
    warning('nTrials is not equal in edat and events indexes')
end
newbli = find(eventCoInd(:,1)==100);%new block index
if(~length(newbli)==nBlocks)
    warning('nBlocks is not equal in edat and events indexes')
end

% compare SOAs
expdata_plannedSOA = nan(length(trials)-1,1);
expdata_measuredSOA = nan(length(trials)-1,1);
events_byIndSOA = nan(length(trials)-1,1);
sample = nan(length(trials)-1,1);
difs_reals = nan(length(trials)-1,1);
eventInd = nan(length(trials)-1,1);

for it = 1:(length(trials)-1)
    expdata_plannedSOA(it,1) = trials(it).SOA;
    expdata_measuredSOA(it,1) = trials(it+1).time-trials(it).time;
    events_byIndSOA(it,1) = (stimCoInd(it+1,2)-stimCoInd(it,2))/srate;
    sample(it,1) = stimCoInd(it,2);
    eventInd(it,1) = it;
end
difs_reals(:,:) = abs(expdata_measuredSOA-events_byIndSOA);
compareSOA = table(eventInd,sample,expdata_plannedSOA,expdata_measuredSOA,events_byIndSOA,difs_reals);
bigDifInd = find(difs_reals>=0.02);
compareSOA(bigDifInd,:)
%correct SOA is expdata_measuredSOA, because it is based on measured timing
%% calculate RA and plot - one participant
        
%parameters:
R0=0.5;
phaseName = 'Passive';
sigma = 5;%MIDI
tau = 1.8;%seconds
SOA_threshold = 0.6;%seconds
Mis = nan;%[20:130]';

for s=whichSubjects
    
    %load expdata and markers
    FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
    load(['D' EDATfolder(2:end) FileName '_expdata.mat' ])
    %load and read markers
    VMRKfile = [ExportFolder FileName  '_dt_RDI_imported.vmrk'];
    [eventCoInd, artInd]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed

    %calculate expected activity RA
    [ RA, smpls, stimCodes, seqInd ] = calcRA(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold,Mis);
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
    phaseName = 'Passive';
    ERRA = nan(size(stimCodes,1),size(stimCodes,2),5);
    CIRA = nan(size(stimCodes,1),size(stimCodes,2),5);

    for it = 1:size(stimCodes,1)
        for is = 1:size(stimCodes,2)
            Codes = squeeze(stimCodes(it,is,:));
            Codes = Codes - floor(Codes/10)*10;
            whichCodes = unique(Codes);
            for ic=1:length(whichCodes)
                Mj = expdata.(phaseName).blocks(it).MIDIs{expdata.(phaseName).blocks(it).CODE_NOTES==ic};%current stimulus MIDI
                is5first = squeeze(seqInd(it,is,:))<=5;
                if isnan(Mis)
                    ERRA(it,is,ic) = mean(RA(it,is,ic,Codes==ic),4);
                    CI = Confidence(RA(it,is,ic,Codes==ic & ~is5first));
                else
                    ERRA(it,is,ic) = mean(RA(it,is,Mis==Mj,Codes==ic),4);
                    CI = Confidence(RA(it,is,Mis==Mj,Codes==ic & ~is5first));
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
    legend(Bnames)
    suptitle(['Subj ' num2str(s) ', Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)])
end
%% calculate RA - all participants, taus, sigmas

mkdir(modelFolder,date)
saveFolder = [modelFolder date filesep];

%parameters:
R0=0.5;
phaseName = 'Passive';
%sigmas = [1:12,14,16,18];%MIDI
sigmas = [7];%MIDI
taus = [3];%seconds
SOA_threshold = 0.6;%seconds
%Mis = (20:130)';
Mis = nan; %for having only the 5 Mis of the stimuli
artIndss = cell(size(whichSubjects));

tictot = tic;
for s=whichSubjects
    ticsubj = tic;
    fprintf(['subj %2.0f loading...'],s)
    %load expdata and markers
    FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
    load(['D' EDATfolder(2:end) FileName '_expdata.mat' ])
    %load and read markers
    VMRKfile = [ExportFolder FileName  '_dt_RDI_imported.vmrk'];
    [eventCoInd, artInd]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed
    artIndss{s} = artInd;
    fprintf('Done\n')
    for sigma = sigmas
        ticsig = tic;
        fprintf(['Subj = %2.0f, sigma = %2.2f'],s,sigma)
        for tau = taus
            %calculate expected activity RA
            if sigma == sigmas(1) && tau == taus(1)
                [ RA, smpls, stimCodes, seqInds ] = calcRA(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold,Mis);
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
                [ RA, ~, ~, ~ ] = calcRA(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold,Mis);
            end
            RAs_SigTau(s,:,:,:,:,sigmas==sigma,taus==tau) = RA;
            %clear RA
            if sigma == sigmas(1) && tau == taus(1)%once per subject            
                smplss(s,:,:,:) = smpls;
                stimCodess(s,:,:,:) = stimCodes;
                seqIndss(s,:,:,:) = seqInds;
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
RAdate = '23-Dec-2018';
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
RApeaks = nan(whichSubjects(end),length(Bnames),nStim,nStim,length(sigmas),length(taus));%subjects x types x cond x prevcond x sigmas x taus
conditions = {'note 1','note 2','note 3','note 4','note 5'};
prevConditions = {'note 1','note 2','note 3','note 4','note 5'};
STIMcodes = [1,2,3,4,5];

tictot = tic;
for bl = 1:5
    ticblock = tic;
    
    blockName = Bnames{bl};
    disp(blockName)
    Bcode = Bcodes_names.Bcodes(strcmp(Bcodes_names.Bnames,blockName));

    EventCodes=cell(1,length(conditions));
    for i=1:length(conditions)
        EventCodes{i} = Bcode + STIMcodes(i);
    end
    prevEventCodes = EventCodes;

    grandRA = calcGrandRA_CondPrev( bl, whichSubjects, RAs_SigTau, EventCodes, prevEventCodes, stimCodess, smplss, seqIndss, artIndss);
    RApeaks(:,bl,:,:,:,:) = squeeze(grandRA);
end
save([saveFolder 'RApeaks'],'RApeaks')
%% Documentation of runs
%19-Dec-2018 : full
%23-Dec-2018 : full with larger sigmas, half resolution, only odd sigmas
%24-Dec-2018 : smaller taus
%27-Dec-2018 : smaller sigmas
%% average and plot U-shape
RAdate = '23-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
load([loadFolder 'Metadata'])

plotIndividual = false;
phaseName = 'Passive';
sigma = 7;
tau = 2;
nStim = 5;

ERRA = nan(whichSubjects(end),nStim,size(stimCodess,2));
CIRA = nan(nStim,size(stimCodess,2));

%for each subject:
for s = whichSubjects
    stimCodes = squeeze(stimCodess(s,:,:,:));
    for it = 1:size(stimCodes,1)
        Codes = squeeze(stimCodes(it,1,:));
        Codes = Codes - floor(Codes/10)*10;
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
legend(Bnames)
suptitle(['Grand RA, Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)])
%% average and plot prevcond bars
RAdate = '23-Dec-2018';
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

for it=1:5
    ERPfigure;
    set(gcf,'Position',[10+300*it 10 300 700])
    for curr = 1:5 
        subplot(5,1,curr)
        RA = squeeze(nanmean(ERRA(:,curr,:,it),1));
        bar(1:5,RA,'linewidth',2)
        %errorbar(squeeze(nanmean(ERRA(:,:,:),1)),CIRA,'linewidth',2)
        set(gca,'fontsize',14)
        legend(Bnames)
        title(['Curr note ' num2str(curr)])
    end
    xlabel('previous note')
    suptitle({['Grand RA prevcond, Block type: ' num2str(it)] [' Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)]})
end
%% fitting parameters to all block types together -
%% find best scaling factor

%load Model
RAdate = '27-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])

%load data
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
whichpeaks = {'N1','P2'};
peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));
        
SFs = nan(size(RApeaks,5),size(RApeaks,6),length(whichpeaks));
    
for ipeak = 1:length(whichpeaks)
    for bl=1:length(allPeak_amps)
        peakAmps(:,bl,:,:) = allPeak_amps{bl}(whichSubjects,:,:,ipeak);
        %subj x bl x con x prevcon 
%         peakAmps(:,bl,:) = allGrandcon_amps{bl}(whichSubjects,:,ipeak);
    end
    data = nanmean(peakAmps);
    data = data(:);data(isnan(data))=[];
    for isig = 1:length(sigmas)
        for itau = 1:length(taus)
            model = 1-nanmean(RApeaks(whichSubjects,:,:,:,isig,itau));
            model = model(:);model(isnan(model))=[];
            SFs(isig,itau,ipeak) = model\data;%performs SS linear regression
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
save([saveFolder 'scalingFactors'],'SFs','sigmas','taus','whichpeaks')                   
%% plot on top of bargraph N1 data

RAdate = '24-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
save([loadFolder 'scalingFactors'],'SFs')
peak = 'P2';
ipeak = find(strcmp(whichpeaks,peak));

plotflag = true;

[xlocs,N1_means] = BarPlotN1(plotflag,peak);
title([peak ' data (bars) versus model (lines)'])
%notice - peak_means dim is blocks X notes whereas ERR is subj X notes X
%blocks. For historical reasons.
meanV = mean(mean(N1_means));
maxV = max(max(abs(N1_means)));
Predict = squeeze(nanmean((1-RApeaks)));
PredictU = squeeze(nanmean(Predict,3));
Psigs = [1,3,11,17,21,29];
Pts = [1.8];
nComb=numel(Psigs)*numel(Pts);
cmap = zeros(numel(Psigs)*numel(Pts),3);
cmap(1:ceil(numel(Psigs)*numel(Pts)/2),1) = linspace(0,1,ceil(numel(Psigs)*numel(Pts)/2));
cmap(ceil(numel(Psigs)*numel(Pts)/2):end,1) = ones(floor(numel(Psigs)*numel(Pts)/2)+1,1);
cmap(ceil(numel(Psigs)*numel(Pts)/2):end,2) = linspace(0,0.8,floor(numel(Psigs)*numel(Pts)/2)+1);
cmap(ceil(numel(Psigs)*numel(Pts)/2):end,3) = linspace(0,0.8,floor(numel(Psigs)*numel(Pts)/2)+1);

hold on
for bl=1:length(Bnames)
    xind=((bl-1)*length(Bnames)+1):(bl*length(Bnames));
    ist=0;cbarTicks=cell(1,nComb);
    for sig = Psigs
        for t = Pts
            ist = ist+1;
            SF = SFs(sig==sigmas,t==taus,ipeak);
            RAs = squeeze(nanmean(nanmean(RApeaks(whichSubjects,bl,:,:,sig==sigmas,t==taus)),4));
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
set(c,'Position',[0.85,0.1,0.02,0.5])
set(c,'YAxisLocation','right')
set(c,'fontsize',14)
%% calc and plot model error

%load Model
RAdate = '23-Dec-2018';
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
for ipeak=1:length(whichpeaks)
    ax(ipeak) = subplot(1,2,ipeak);
    imagesc(taus,sigmas,Errs(:,:,ipeak));
    title(whichpeaks{ipeak})
    set(gca,'fontsize',12)
    xlabel('\tau (sec)','fontsize',14);
    ylabel('\sigma (semitones)','fontsize',14)
   
end
suptitle('normalized error between model and data')
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'STE units')
end
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'fig')
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_' date],'jpg')

%% fitting parameters to each block type separately -
%% find best scaling factor (each block type)

%load Model
RAdate = '24-Dec-2018';
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

        data = nanmean(peakAmps(:,:,:));
        data = data(:);data(isnan(data))=[];
        for isig = 1:length(sigmas)
            for itau = 1:length(taus)
                model = 1-nanmean(RApeaks(whichSubjects,bl,:,:,isig,itau));
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

%% calc and plot model error

%load Model
RAdate = '23-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
load([loadFolder 'scalingFactors_blockTypes'])                   

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
        data = nanmean(peakAmps);
        data = data(:);data(isnan(data))=[];
        for isig = 1:length(sigmas)
            for itau = 1:length(taus)
                model = 1-nanmean(RApeaks(whichSubjects,bl,:,:,isig,itau));
                model = model(:);model(isnan(model))=[];
                SF = SFs(isig,itau,ipeak,bl);
                Errs(isig,itau,ipeak,bl) = rms((data - model*SF)/mean(ste));%performs SS linear regression
            end
        end
    end
end
saveFolder = loadFolder;
save([saveFolder 'Errs'],'Errs','sigmas','taus','whichpeaks')                   

ERPfigure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
ii=0;
for ipeak=1:length(whichpeaks)
    for bl=1:length(allPeak_amps)
        ii=ii+1;        
        ax(ii) = subplot(2,5,ii);
        imagesc(taus,sigmas,Errs(:,:,ipeak,bl));
        title(whichpeaks{ipeak})
        set(gca,'fontsize',12)
        xlabel('\tau (sec)','fontsize',14);
        ylabel('\sigma (semitones)','fontsize',14)
        c(ii) = colorbar(ax(ii));
    end
    ylabel(c(ii),'STE units')
end
suptitle('normalized error between model and data')
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_bType_' date],'fig')
saveas(gcf,[AnalysisFolder 'Figures\model\norm_err_model_data_bType_' date],'jpg')
