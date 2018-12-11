%The Model as in Herrmann et al.
% Nov 25 2018, Tamar Regev
%% Definitions
cd('L:\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
sbjcts = 1:37;
whichSubjects = sbjcts(~ismember(sbjcts,badSubjects));

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
    
whichSubjects = sbjcts(~ismember(sbjcts,badSubjects));
    
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
sigmas = [1:18];%MIDI
taus = [0.2:0.2:5];%seconds
SOA_threshold = 0.6;%seconds
%Mis = (20:130)';
Mis = nan; %for having only the 5 Mis of the stimuli
artIndss = cell(size(whichSubjects));

tictot = tic;
for s=whichSubjects
    ticsubj = tic;
    fprintf(['subj %2f loading...'],s)
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
        fprintf(['Subj = %2f, sigma = %2f' num2str(sigma)],s,sigma)
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
        fprintf(['Done in %2f sec.\n'],toc(ticsig))
    end
    disp(['Done Subj ' Subjects{s} ' in ' num2str(toc(ticsubj)) ' sec.'])
end
fprintf('saving RAs_SigTau ...');ticsaving = tic;
save([saveFolder 'RAs_SigTau'],'RAs_SigTau','-v7.3')
fprintf(['done saving in %2.2f sec.\n'],toc(ticsaving))
save([saveFolder 'Metadata'],'smplss','stimCodess','seqIndss','artIndss')
save([saveFolder 'Params'],'sigmas','taus','Mis','R0','SOA_threshold')
disp(['Done all in ' num2str(toc(tictot)) ' sec.'])

%% compare to ERPs-
%% calc. expected peaks
RAdate = '10-Dec-2018';
loadFolder = [modelFolder RAdate filesep];

tic
fprintf(['Loading...'])
load([loadFolder 'Params'])
load([loadFolder 'RAs_SigTau'])
load([loadFolder 'Metadata'])
fprintf(['Done in %4.1f sec \n'], toc)

nStim = 5;
%initiatize peaks array
RApeaks = nan(whichSubjects(end),length(Bnames),nStim,nStim,length(sigmas),length(taus));
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
save([modelFolder 'RApeaks'],'RApeaks')
%% average and plot U-shape
load([modelFolder 'RApeaks'])
load([loadFolder 'Metadata'])

phaseName = 'Passive';
sigma = 10;
tau = 1.8;
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
    ERPfigure;
    plot(squeeze(ERRA(s,:,:)),'linewidth',2)
    set(gca,'fontsize',14)
    legend(Bnames)
    suptitle(['Subj ' num2str(s) ', Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)])
end
%grand average across subjects:
for it = 1:size(stimCodes,1)
    for ic = 1:nStim
        CI = Confidence(ERRA(:,ic,it));
        CIRA(ic,it) = squeeze(nanmean(ERRA(:,ic,it),1)) - CI(1);
    end
end
ERPfigure;
plot(squeeze(nanmean(ERRA(:,:,:),1)),'linewidth',2)
%errorbar(squeeze(nanmean(ERRA(:,:,:),1)),CIRA,'linewidth',2)
set(gca,'fontsize',14)
legend(Bnames)
suptitle(['Grand RA, Sigma = ' num2str(sigma) ', Tau = ' num2str(tau) ', R0 = ' num2str(R0)])
                   
%% plot on top of bargraph N1 data

