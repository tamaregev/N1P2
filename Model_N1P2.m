%The Model as in Herrmann et al.
% Nov 25 2018, Tamar Regev
%% Definitions
Definitions_N1P2
sbjcts = 1:37;

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
Mis = [20:130]';

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
    if 0
        ERPfigure;
        isp=0;
        for it=1:size(RA,1)
            for ib=1:size(RA,2)
                isp=isp+1;
                subplot(size(RA,1),size(RA,2),isp);
                imagesc(1:size(RA,4),Mis,squeeze(RA(it,ib,:,:)))

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
                ERRA(it,is,ic) = mean(RA(it,is,Mis==Mj,Codes==ic),4);
                is5first = squeeze(seqInd(it,is,:))<=5;
                CI = Confidence(RA(it,is,Mis==Mj,Codes==ic & ~is5first));
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
saveFolder = [modelFolder date];

%parameters:
R0=0.5;
phaseName = 'Passive';
sigmas = 1:18;%MIDI
taus = 0.2:0.2:5;%seconds
SOA_threshold = 0.6;%seconds
Mis = (20:130)';

whichSubjects = sbjcts(~ismember(sbjcts,badSubjects));

for sigma = sigmas
    for tau = taus
        tictausig = tic;
        for s=whichSubjects
            ticsubj = tic;
            disp(['Subj = ' num2str(s) ', sigma = ' num2str(sigma) ', tau = ' num2str(tau) ])

            %load expdata and markers
            FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
            load(['D' EDATfolder(2:end) FileName '_expdata.mat' ])
            %load and read markers
            VMRKfile = [ExportFolder FileName  '_dt_RDI_imported.vmrk'];
            [eventCoInd, artInd]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed

            %calculate expected activity RA
            if sigma == sigma(1) && tau == tau(1)
                [ RA, smpls, stimCodes, seqInds ] = calcRA(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold,Mis);
            else
                [ RA, ~, ~, ~ ] = calcRA(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold,Mis);
            end
            if s==whichSubjects(1)
                RAs = nan([whichSubjects(end),size(RA)]);
                smplss = nan([whichSubjects(end),size(smpls)]); 
                stimCodess = nan([whichSubjects(end),size(stimCodes)]);
                seqIndss = nan([whichSubjects(end),size(seqInds)]);        
            end
            RAs(s,:,:,:,:) = RA;
        
            if sigma == sigma(1) && tau == tau(1)            
                smplss(s,:,:,:) = smpls;
                stimCodess(s,:,:,:) = stimCodes;
                seqIndss(s,:,:,:) = seqInds;
            end
            fclose('all');
            disp(['Done Subj ' Subjects{s} ' in ' num2str(toc(ticsubj)) ' sec.'])
        end
        disp(['...saving sigma = ' num2str(sigma) ', tau = ' num2str(tau) '...'])
        save([saveFolder filesep 'RAs_Sig' num2str(find(sigmas==sigma)) 'Tau' num2str(find(taus==tau)) '_' date],'RAs')
        disp(['done tausig in ' num2str(toc(tictausig)) 'sec.'])
    end
end
save([saveFolder 'Metadata_' date],'smplss','stimCodess','seqIndss')
save([saveFolder 'Params_' date],'Mis','sigmas','taus')
        
%% compare to ERPs
RAdate = '06-Dec-2018';
load([modelFolder 'RAs_' RAdate])

RApeaks = nan();
