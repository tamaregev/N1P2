%calc model single trials
Folder = 'L:\Experiments\N1P2\Analysis\Model\singleTrials';
mkdir(Folder)
%% compute single trials peak amps 
%no need to run again
ExpN = 3;
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
mkdir(Expdir)

[ singleAmps, isArtefact ] = calcAmpSingleTrials( ExpN );
%Done all in 2520.1113

save([Expdir filesep 'singleAmps'],'singleAmps','isArtefact')
%% prep model structures
%no need to run again
ExpN = 3;
switch ExpN
    case 1
    case 2
    case 3
        Definitions_N1P2
end

%load model - RASigTau
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];

fprintf(['Loading...']);tic
load([loadFolder 'RAs_SigTau'])
load([loadFolder 'Params'])
load([loadFolder 'Metadata'])
fprintf(['Done in %4.1f sec \n'], toc)
%takes about 6.5 minutes to load RASigTau
%RAs_SigTau is %subjs x types x sequences x channels x timepoints

%Keep only RA of relevant stim:
relRA = nan(size(RAs_SigTau,1),size(RAs_SigTau,2),size(RAs_SigTau,3),size(RAs_SigTau,5),size(RAs_SigTau,6),size(RAs_SigTau,7));
Codess = stimCodess - 10*floor(stimCodess/10);
for s=whichSubjects
    for bi = 1:size(RAs_SigTau,2)
        for seq=1:size(RAs_SigTau,3)
            for ti=1:size(RAs_SigTau,5)
                relRA(s,bi,seq,ti,:,:) = RAs_SigTau(s,bi,seq,Codess(s,bi,seq,ti),ti,:,:);
            end
        end
    end
end
save([loadFolder 'relRA'],'relRA','-v7.3');clear RAs_SigTau

%concatenate sequences:
load([loadFolder 'relRA'])
RAcat = [];
for is=1:size(relRA,3)
    RAcat = cat(4,RAcat,relRA(:,:,is,:,:,:));
end
RAcat = squeeze(RAcat);
clear relRA

%concatenate sequences and stimuli:
load([loadFolder 'Metadata'])%artIndss, seqIndss, smplss, stimCodess
seqIndcat = [];stimCodecat = []; smplcat = [];
for is=1:size(seqIndss,3)
    seqIndcat = cat(4,seqIndcat,seqIndss(:,:,is,:));
    stimCodecat = cat(4,stimCodecat,stimCodess(:,:,is,:));
    smplcat = cat(4,smplcat,smplss(:,:,is,:));    
end
seqIndcat = squeeze(seqIndcat);
stimCodecat = squeeze(stimCodecat);
smplcat = squeeze(smplcat); 

save([loadFolder 'RAcat'],'RAcat','-v7.3');
save([loadFolder 'metadatacat'],'seqIndcat','stimCodecat','smplcat')
%reduce size of model :
itaumax = 25;
RAcat = RAcat(:,:,:,:,1:itaumax);
taus = taus(1:itaumax);
save([loadFolder 'RAcat_reduced'],'RAcat','taus','-v7.3')
%% compare data and model
ExpN = 3;
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
RAdate = '19-Dec-2018';
tic
[ SFs, Errs, Evs , sigmas, taus, ws ] = modelSingleTrials( ExpN, Folder, RAdate );
disp(['Done in ' num2str(toc) ' sec.'])
%8.1833 minutes
save([Expdir filesep 'stModel'],'SFs', 'Errs', 'sigmas', 'taus')
save([Expdir filesep 'stModel_Ev_' num2str(ws)],'Evs','ws')
%% plot SFs
ExpN = 3;
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
RAdate = '19-Dec-2018';
load([Expdir filesep 'stModel'])
switch ExpN
    case 1
    case 2
    case 3
        Definitions_N1P2
end

whichpeaks = {'N1','P2'};
loadFolder = [modelFolder RAdate filesep];
biasFlag = true;
ERPfigure
labels = {'intercept','slope'};
isp=0;%for biasFlag
for ipeak=1:length(whichpeaks)
    if biasFlag
        for i=1:2
            isp=isp+1;
            subplot(2,2,isp)
            imagesc(taus,sigmas,SFs(:,:,ipeak,i))
            caxis([min(min(SFs(3:end,2:end,ipeak,i))) max(max(SFs(3:end,2:end,ipeak,i)))])
            title(labels{i})
            colorbar
        end
    else
        subplot(1,2,ipeak)
        imagesc(taus,sigmas,SFs(:,:,ipeak))
    end
end
suptitle('scaling factors')
%% plot error space
EvFlag = false;
if EvFlag
    Val = Evs;
else
    Val = Errs;
end
ExpN = 3;
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
RAdate = '19-Dec-2018';
load([Expdir filesep 'stModel'])
%Errs=Errs./8.882;
whichpeaks = {'N1','P2'};

ERPfigure;
set(gcf,'Position', [50 50 1200 450])
bestaus = nan(2,1);
bestsigmas = nan(2,1);
for ipeak=1:length(whichpeaks)
    errmin = min(min(min(Val(:,:,ipeak))));
    errmax = max(max(max(Val(:,:,ipeak))));
    ax(ipeak) = subplot(1,2,ipeak);
    if EvFlag
        [row col]=find(Val(:,:,ipeak)==max(max(Val(:,:,ipeak))));
    else
        [row col]=find(Val(:,:,ipeak)==min(min(Val(:,:,ipeak))));
    end
    imagesc(taus,sigmas,Val(:,:,ipeak));
    caxis([errmin errmax])
     if exist('bestSFflag','var')
        if ~bestSFflag
            caxis([1.2 2])
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
    text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',22);
end
if EvFlag
    suptitle('explained variance of model in data')
else
    suptitle('normalized error between model and data')
end
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    if EvFlag
        ylabel(c(ipeak),'explained variance %','fontsize',16)
    else
        ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
    end
end
%% calc significance, permutations
num_permutations = 1000;

ExpN = 3;
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
RAdate = '19-Dec-2018';
tic
plotFlag = true;
compareEstimates( ExpN, Folder, RAdate, num_permutations, plotFlag )
disp(['Done in ' num2str(toc) ' sec.'])

%% prep structures for Eli:
ExpN = 3;

%load data single trials
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
load([Expdir filesep 'singleAmps'])
%singleAmps
%isArtefact

RAdate = '19-Dec-2018';
switch ExpN
    case 1
    case 2
    case 3
        Definitions_N1P2
        nBlockTypes = 5;
        PhaseName = 'Passive';
end
loadFolder = [modelFolder RAdate filesep];

load([loadFolder 'Params'])
load([loadFolder 'Metadata'])%artIndss, seqIndss, smplss, stimCodess
load([loadFolder 'metadatacat']);
fprintf(['Loading...']);tic
load([loadFolder 'RAcat_reduced'])
fprintf(['Done in %4.1f sec \n'], toc)

EliFolder = [modelFolder 'Eli' filesep];
mkdir(EliFolder);
save([EliFolder 'taus'],'taus')
save([EliFolder 'sigmas'],'sigmas')
save([EliFolder 'singleAmps'],'singleAmps')
save([EliFolder 'isArtefact'],'isArtefact')
save([EliFolder 'seqIndcat'],'seqIndcat')
RA = RAcat;
save([EliFolder 'RA'],'RA')
save([EliFolder 'whichSubjects'],'whichSubjects')
save([EliFolder 'stimCodecat'],'stimCodecat')

% prep freqs: block types x tone for the frequencies used in the experiment
% prep MIDIs:
addpath('L:\Z backup\Tamar\fromZ\Documents\MATLAB\MatlabFunctions\mine')
s=2;
FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
load(['L' EDATfolder(2:end) FileName '_expdata.mat' ])
MIDIs = nan(nBlockTypes,5);
for ibt = 1:nBlockTypes
     m = expdata.(phaseName).blocks(ibt).MIDIs;%MIDIs of all stimuli
     for im=1:length(m)
        MIDIs(ibt,im) = m{im};
     end
end
freqs = MIDI2freq(MIDIs);

save([EliFolder 'MIDIs'],'MIDIs')
save([EliFolder 'freqs'],'freqs')

%% single subjects
    %% compare data and model
ExpN = 3;
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
RAdate = '19-Dec-2018';

[ SFs, Errs, data, model, sigmas, taus ] = modelSingleTrials_individual( ExpN, Folder, RAdate );
save([Expdir filesep 'stModel_individual'],'SFs', 'Errs', 'data', 'model', 'sigmas', 'taus')
    %% plot SFs
ExpN = 3;
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
RAdate = '19-Dec-2018';
load([Expdir filesep 'stModel_individual'])
switch ExpN
    case 1
    case 2
    case 3
        Definitions_N1P2
end

whichpeaks = {'N1','P2'};
loadFolder = [modelFolder RAdate filesep];
biasFlag = true;
labels = {'intercept','slope'};
ERPfigure
for s=whichSubjects
    isp=0;%for biasFlag
    for ipeak=1:length(whichpeaks)
        if biasFlag
            for i=1:2
                isp=isp+1;
                subplot(2,2,isp)
                imagesc(taus,sigmas,squeeze(SFs(s,:,:,ipeak,i)))
                caxis([min(min(SFs(s,3:end,2:end,ipeak,i))) max(max(SFs(s,3:end,2:end,ipeak,i)))])
                title(labels{i})
                colorbar
            end
        else
            subplot(1,2,ipeak)
            imagesc(taus,sigmas,squeeze(SFs(s,:,:,ipeak)))
        end
    end
    suptitle(['Subj ' num2str(s) ', scaling factors'])
    pause
    hold off
end

    %% plot error space

ExpN = 3;
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
RAdate = '19-Dec-2018';
load([Expdir filesep 'stModel_individual'])
%Errs=Errs./8.882;
whichpeaks = {'N1','P2'};
    
ERPfigure;
set(gcf,'Position', [50 50 1200 450])
bestaus = nan(whichSubjects(end),2,1);
bestsigmas = nan(whichSubjects(end),2,1);
  
for s=whichSubjects
    for ipeak=1:length(whichpeaks)
        errmin = min(min(min(Errs(s,:,:,ipeak))));
        errmax = max(max(max(Errs(s,:,:,ipeak))));
        ax(ipeak) = subplot(1,2,ipeak);
        [row col]=find(squeeze(Errs(s,:,:,ipeak))==min(min(squeeze(Errs(s,:,:,ipeak)))));
        imagesc(taus,sigmas,squeeze(Errs(s,:,:,ipeak)));
        caxis([errmin errmax])
         if exist('bestSFflag','var')
            if ~bestSFflag
%                 caxis([1.2 2])
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
        bestaus(s,ipeak) = taus(col);
        bestsigmas(s,ipeak) = sigmas(row);
        text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',22);
    end
    suptitle(['Subj ' num2str(s) ', normalized error between model and data'])
    for ipeak = 1:length(whichpeaks)
        c(ipeak) = colorbar(ax(ipeak));
        ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
    end
    %pause
    hold off
end
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
    X=ones(whichSubjects(end),2);X(:,2)=X(:,2)+1;
    Xjit=X+randn(size(X))/10;
    plot(Xjit',data','o','markersize',6,'Color','k')
    title(['best fit \' strings{i}])
    %plot([1,2],[mean(data(:,1)),mean(data(:,2))],'g','linewidth',2)
    [p,H]=signrank(data(:,2),data(:,1));
    text(1,4,['Wilcoxon p=' num2str(p,2)],'fontsize',14,'Color','g')
end