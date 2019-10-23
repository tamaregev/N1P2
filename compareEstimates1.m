function [bestauss, bestsigmass  ] = compareEstimates( ExpN, Folder, RAdate, num_permutations, plotFlag )
%COMAREESTIMATES compares the values of tau and sigma estimated for N1 and P2 
% by computing the empirical distribution of difference betwween
%parameter estimates given the null hypothesis - that N1
%and P2 give the same values. We permute between N1 and P2 of the different
%participants and calculate the estimates.

switch ExpN
    case 1
    case 2
    case 3
        Definitions_N1P2
end

whichpeaks = {'N1','P2'};

%load data single trials
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
load([Expdir filesep 'singleAmps'])

%load model - 
loadFolder = [modelFolder RAdate filesep];

load([loadFolder 'Params'])
load([loadFolder 'Metadata'])%artIndss, seqIndss, smplss, stimCodess
load([loadFolder 'metadatacat']);
fprintf(['Loading...']);tic
load([loadFolder 'RAcat_reduced'])
fprintf(['Done in %4.1f sec \n'], toc)
%takes less than 1 minute to load RAcat_reduced

%the REAL permutation -
real_vec = ones(whichSubjects(end),1);

perm_vecs = [ real_vec , round(rand(whichSubjects(end),num_permutations)) + 1];
bestauss = nan(num_permutations,2);
bestsigmass = nan(num_permutations,2);

data=nan(size(singleAmps));
data(whichSubjects,:,:,1)=zscore(singleAmps(whichSubjects,:,:,1),0,1);
data(whichSubjects,:,:,2)=zscore(singleAmps(whichSubjects,:,:,1),0,2);
%permutations:
for ip = 1:num_permutations + 1
    ticperm=tic;
    disp(ip)
    perm_vec = perm_vecs(:,ip);
    dataperm = nan(size(singleAmps));
    dataperm(perm_vec==1,:,:,:) = (data(perm_vec==1,:,:,:));
    dataperm(perm_vec==2,:,:,[1,2]) = (data(perm_vec==2,:,:,[2,1]));
    Errs = calcErrModelData(dataperm,RAcat, sigmas, taus, whichpeaks, whichSubjects, isArtefact, seqIndcat);
    if plotFlag
       [bestaus, bestsigmas] = plotErrs(Errs, whichpeaks, taus, sigmas );
    else
        [bestaus, bestsigmas] = findBestParams(Errs, whichpeaks, taus, sigmas );
    end
    bestauss(ip,:) = bestaus;
    bestsigmass(ip,:) = bestsigmas;
    disp(['Done perm ' num2str(ip) ' in ' num2str(toc(ticperm)) ' sec.'])
end


   
end

%% functions

function Errs = calcErrModelData(datafull,RA, sigmas, taus, whichpeaks, whichSubjects, isArtefact, seqIndcat)
        SFs = nan(length(sigmas),length(taus),length(whichpeaks),2);
        Errs = nan(length(sigmas),length(taus),length(whichpeaks));
        isArt = isArtefact(whichSubjects,:,:);isArt = isArt(:);   
        seqi = seqIndcat(whichSubjects,:,:);seqi=seqi(:);
        for ipeak = 1:length(whichpeaks)
            %disp([whichpeaks(ipeak)])
            for isig = 1:length(sigmas)
                %disp(['sigma = ' num2str(sigmas(isig))])
                for itau = 1:length(taus)
                    %RAcat is %subjs x types x timepoints
                    data = datafull(whichSubjects,:,:,ipeak); % data is subjs x types x timepoints
                    data = data(:);
                    model = 1-RA(whichSubjects,:,:,isig,itau);
                    model = model(:);
                    %exclude artifacts and initial conditions (f for fix):
                    modelf = model(isArt==0 & seqi>10); dataf = data(isArt==0 & seqi>10);
                    model = modelf; data = dataf;
                    newmodel = [ones(size(model)),model];
                    SF = newmodel\data;%performs SS linear regression
                    SFs(isig,itau,ipeak,:) = SF;
                    newmodel = [ones(size(model)), model];
                    predict = newmodel*SF;
                    Errs(isig,itau,ipeak) = rms((data - predict));
                 end
            end
        end
end
    
function [bestaus, bestsigmas] = plotErrs(Errs, whichpeaks, taus, sigmas )
    Val = Errs;
    ERPfigure;
    set(gcf,'Position', [50 50 1200 450])
    bestaus = nan(2,1);
    bestsigmas = nan(2,1);
    for ipeak=1:length(whichpeaks)
        errmin = min(min(min(Val(:,:,ipeak))));
        errmax = max(max(max(Val(:,:,ipeak))));
        ax(ipeak) = subplot(1,2,ipeak);
        [row col]=find(Val(:,:,ipeak)==min(min(Val(:,:,ipeak))));
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
    suptitle('normalized error between model and data')
    for ipeak = 1:length(whichpeaks)
        c(ipeak) = colorbar(ax(ipeak));
        ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
    end
end

function [bestaus, bestsigmas] = findBestParams(Errs,whichpeaks,taus,sigmas)
    Val = Errs;
    bestaus = nan(2,1);
    bestsigmas = nan(2,1);
    for ipeak=1:length(whichpeaks)
        [row col]=find(Val(:,:,ipeak)==min(min(Val(:,:,ipeak))));
        bestaus(ipeak) = taus(col);
        bestsigmas(ipeak) = sigmas(row);
    end
end