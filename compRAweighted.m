%% Compare RA to weighted responses
% Nov 18 2019 Tamar
% compare model predictions based on a single RA channel (corresponding to
% current stimulus) to a calculation based on weighting by the tuning curve

%% Definitions
%For Exp 3 - N1P2
% single participant s102 (this is only for model predictions so doesn't matter which)
%from script Model_N1P2
addpath('S:\Lab-Shared\NewDataArch\CommonResources\Tools\Matlab_Tools')
addpath('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
addpath('S:\Lab-Shared\Experiments\MMNchroma\Analysis')
addpath('S:\Lab-Shared\Experiments\MMNchromaF\Analysis')

Definitions_N1P2
GHfolder = [AnalysisFolder 'N1P2_GH'];
cd(GHfolder)
FigFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\PaperFigures\Model\weightedRA';
%% calculate RA and plot - one participant
        
%parameters:
R0=0.5;
phaseName = 'Passive';
sigma = 8;%MIDI
tau = 2;%seconds
SOA_threshold = 0.6;%seconds
%Mis = nan;%
Mis=[20:130]';

s=2;
%load expdata and markers
FileName = [ExpName '_' Subjects{s} '_' sessions{s}(1)];
cd(EDATfolder)
load([FileName '_expdata.mat' ])
cd(AnalysisFolder)
%load and read markers
VMRKfile = [ExportFolder FileName  '_dt_RDI_imported.vmrk'];
[eventCoInd, artInd]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed

%calculate expected activity RA
[ RA, smpls, stimCodes, seqInd ] = calcRA(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold,Mis);
% trim for all 3 sequences
RAt=nan(size(RA,1),size(RA,3),size(RA,4)*size(RA,2));
smplst = nan(size(smpls,1),size(smpls,2)*size(smpls,3));
stimCodest = nan(size(smplst));seqIndt = nan(size(smplst));
it=1;
for is=1:size(RA,2)
    RAt(:,:,it:(it+size(RA,4)-1))=RA(:,is,:,:);
    smplst(:,it:(it+size(RA,4)-1))=smpls(:,is,:);
    stimCodest(:,it:(it+size(RA,4)-1))=stimCodes(:,is,:);
    seqIndt(:,it:(it+size(RA,4)-1))=seqInd(:,is,:);
    it=it+size(RA,4);
end

% pause
blocks = expdata.(phaseName).blocks;

%%%%%%%% plot RA %%%%%%%%%%%
if 1
    ERPfigure;
    isp=0;
    for it=1:size(RAt,1)
        isp=isp+1;
        subplot(size(RAt,1),1,isp);
        if isnan(Mis)
            imagesc(1:size(RAt,3),5,squeeze(RAt(it,:,:)))
        else
            imagesc(1:size(RAt,3),Mis,squeeze(RAt(it,:,:)))
        end
        title(['Block Type: ' num2str(Bnames{it})])
        colorbar;
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
            Mj = blocks(it).MIDIs{blocks(it).CODE_NOTES==ic};%current stimulus MIDI
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

% save  
mkdir(modelFolder,date)
saveFolder = [modelFolder date filesep];
save([saveFolder 'RA_Subj' num2str(s)],'RA','RAt','ERRA','CIRA','tau','sigma','R0','Mis','Bnames','stimCodes','s','blocks','seqInd','smpls','phaseName','seqIndt','smplst','stimCodest')

%% calculate weighted responses from RAt
RAdate = '18-Nov-2019';
saveFolder = [modelFolder RAdate filesep];
s=2;
load([saveFolder 'RA_Subj' num2str(s)])
resp = nan(size(RAt,1),size(RAt,3));
respw = nan(size(RAt,1),size(RAt,3));
hf=ERPfigure;
set(hf,'Position',[10 -300 400 800]);
colors = {[0.3 0.745 0.93],[0.75 0 0.75],[1 0 0],[0.75 0 0.75],[0.3 0.745 0.93]};

for ibt=1:size(RAt,1)
    subplot(size(RAt,1),1,ibt)
    Codes = squeeze(stimCodes(ibt,:));
    Codes = Codes - floor(Codes/10)*10;
    whichCodes = unique(Codes); 
    didi = zeros(size(whichCodes));
    for t=1:size(RAt,3)
            ct=Codes(t);
            
            Mt = blocks(ibt).MIDIs{blocks(ibt).CODE_NOTES==ct};%current stimulus MIDI
            resp(ibt,t)=RAt(ibt,Mis==Mt,t);
            %weighted sum...:
            for i=1:length(Mis)
                Mi=Mis(i);
                tuning(i) = exp(-0.5*(((Mi-Mt)/sigma).^2));
            end
            if ~didi(ct)
                plot(Mis,tuning,'linewidth',2,'Color',colors{whichCodes==ct})
                hold on
            end
            didi(ct)=1;
            respw(ibt,t)=tuning*RAt(ibt,:,t)';
            xlim([Mis(1) Mis(end)])
    end
    if ibt ~=5
        set(gca,'xticklabels',{},'yticklabels',{})
    end
    title([Bnames{ibt}])
    set(gca,'fontsize',12)
end
xlabel('Frequency channels (Midi)')
ylabel('weights')
    
saveas(gcf,[FigFolder filesep 'TuningCurves'],'pdf')
%% compare response to weighted
figure
r=resp(:);rw=respw(:);
sqi=seqIndt(:);
clf
for ibt=1:size(resp,1)
    plot(resp(ibt,:),respw(ibt,:),'.')
    hold on
end
p=polyfit(resp(:),respw(:),1);
y=polyval(p,resp(:));
plot(resp(:),y,'linewidth',2,'Color','k')
legend(Bnames,'location','nw')
%plot(resp(seqIndt<6),respw(seqIndt<6),'.r')
xlabel('RA')
ylabel('RA-weighted')
title('Single channel vs. weighted by tuning curve')
set(gca,'fontsize',16)
saveas(gcf,[FigFolder filesep 'compRA2weighted'],'pdf')
