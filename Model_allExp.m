%compute the model for data from all Experiments together
GHfolder = 'S:\Lab-Shared\N1P2\Analysis\N1P2_GH';
saveFolder = 'S:\Lab-Shared\N1P2\Analysis\Model\AllExp';
whichpeaks = {'N1','P2'};
%% Exp 1 - load 
Definitions_MMNchroma;
whichSubjectsAll{1} = whichSubjects;
%load Model:
%RAdate = '03-Jun-2019';
RAdate = '22-May-2019';
%RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks.'])
RApeaksAll{1} = RApeaks;
load([loadFolder 'Params'])

%load data:
electrodeName = 'Cz';addtag = '';peakdate = '18-Nov-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate],'allPeak_amps')
allPeak_ampsAll{1} = allPeak_amps;
%rename variables
blsAll{1} = [2,4];
%% Exp 2 - load 
Definitions_MMNchromaF;
whichSubjectsAll{2} = whichSubjects;

%load Model
%RAdate = '03-Jun-2019';
RAdate = '20-May-2019';
%RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks.'])
RApeaksAll{2} = RApeaks;
load([loadFolder 'Params'])

%load data:
peakdate = '16-Nov-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate],'allPeak_amps')
allPeak_ampsAll{2} = allPeak_amps;
blsAll{2} = 2;

%% Exp 3 - load
cd(GHfolder)
Definitions_N1P2
whichSubjectsAll{3} = whichSubjects;

%load Model
%RAdate = '03-Jun-2019';

RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks.'])
RApeaksAll{3} = RApeaks;
load([loadFolder 'Params'])

%load data
peakdate = '30-Oct-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate],'allPeak_amps')
allPeak_ampsAll{3} = allPeak_amps;
blsAll{3} = 1:5;

%% calc SFs and Errs all together:
bestSFflag = true;

weightedFlag = false;

SFs = nan(length(sigmas),length(taus),length(whichpeaks));
Errs = nan(length(sigmas),length(taus),length(whichpeaks));

for ipeak = 1:length(whichpeaks)
    for iexp = 1:3
        whichSubjects = whichSubjectsAll{iexp};
        RApeaks = RApeaksAll{iexp};
        allPeak_amps = allPeak_ampsAll{iexp};
        bls = blsAll{iexp};
        ibl=0;
        peakAmps = nan(size(RApeaks(whichSubjects,:,:,:,1,1)));
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
        steAll{iexp}=ste;dataAll{iexp}=data;
    end
    data = [dataAll{1}; dataAll{2}; dataAll{3}];
    ste = [steAll{1}; steAll{2}; steAll{3}];
    for isig = 1:length(sigmas)
        for itau = 1:length(taus)
            for iexp=1:3
                RApeaks = RApeaksAll{iexp};
                whichSubjects = whichSubjectsAll{iexp};
                model = 1-nanmean(RApeaks(whichSubjects,:,:,:,isig,itau));
                model = model(:);model(isnan(model))=[];
                modelAll{iexp} = model;
            end
            model = [modelAll{1}; modelAll{2}; modelAll{3}];
            if weightedFlag
                SFs(isig,itau,ipeak) = model\(data./ste);%performs SS linear regression
            elseif bestSFflag
                SFs(isig,itau,ipeak) = model\data;%performs SS linear regression
            else
                SFs(isig,itau,ipeak) = max(abs(data))*sign(data(1));
                %SFs(isig,itau,ipeak) = max(abs(data))*sign(data(1))/max(model);

            end
            SF = SFs(isig,itau,ipeak);
            Errs(isig,itau,ipeak) = rms((data - model*SF)/mean(ste));%performs SS linear regression
        end
    end
end
if bestSFflag
    save([saveFolder filesep 'scalingFactors_bestSF'],'SFs','bestSFflag','sigmas','taus','whichpeaks')                   
    save([saveFolder filesep 'Errs_bestSF'],'Errs','bestSFflag','sigmas','taus','whichpeaks')                   
else
    save([saveFolder filesep 'scalingFactors_constSF'],'SFs','bestSFflag','sigmas','taus','whichpeaks')
    save([saveFolder filesep 'Errs_constSF'],'Errs','bestSFflag','sigmas','taus','whichpeaks')                   
end

%% plot Err space
%bestSFflag = false;
if bestSFflag
    load([saveFolder filesep 'scalingFactors_bestSF'])                   
    load([saveFolder filesep 'Errs_bestSF'])                   
else
    load([saveFolder filesep 'scalingFactors_constSF'])
    load([saveFolder filesep 'Errs_constSF'])                   
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
            caxis([1.3 2])
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
    
    text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1 0 0],'fontsize',22);
    
end
suptitle('normalized error between model and data')
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
end
save([saveFolder filesep 'bestFitParams'],'bestaus','bestsigmas')

%% 