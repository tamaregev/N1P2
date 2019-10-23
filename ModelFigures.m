%% ModelFigures
FigFolder = 'L:\Experiments\N1P2\Analysis\N1P2_GH\PaperFigures\Model';
lw=2;
addpath('S:\Lab-Shared\NewDataArch\CommonResources\Tools\Matlab_Tools')
%% Experiment 1
cd('L:\Experiments\MMNchroma\Analysis')
Definitions_MMNchroma;%new idea for the first time! move definitions into a separate script!
bls = [2,4];
cs = [1 2 5 3 4];
Exp = 'Exp1';

%%      model

weightedFlag = false;
bestSFflag = true;
%load Model
%RAdate = '10-Apr-2019';
RAdate = '22-May-2019';
loadFolder = [modelFolder RAdate filesep];
whichpeaks = {'N1','P2'};
if bestSFflag
    load([loadFolder 'Errs_bestSF'],'Errs','sigmas','taus','whichpeaks')                   
else
    load([loadFolder 'Errs_constSF'],'Errs','sigmas','taus','whichpeaks')                   
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
    errmin = (min(min(Errs(:,:,ipeak))));
    if bestSFflag
        %caxis([errmin errmax])
    else
        caxis([errmin 2])
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
    text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0, 0],'fontsize',40);
end
suptitle('normalized error between model and data')
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
end
if bestSFflag
    saveas(gcf,[FigFolder filesep 'ErrSpace_bestSF_' Exp],'pdf')
else
    saveas(gcf,[FigFolder filesep 'ErrSpace_constSF_' Exp],'pdf')
end
%% Experiment 2
cd('L:\Experiments\MMNchromaF\Analysis')
Definitions_MMNchromaF;
N1P2GHfolder = ['L:\Experiments\N1P2\Analysis\N1P2_GH'];
addpath(N1P2GHfolder)
bls = 2;
order = [1 2 5 3 4];
addpath('L:\Z backup\Tamar\fromZ\Documents\MATLAB\MatlabFunctions\mine')
Exp = 'Exp2';

%%      model
bestSFflag = true;
weightedFlag = false;
%RA and Errs are calculated in script
%  L:\Experiments\N1P2\Analysis\N1P2_GH\Model_MMNchromaF.m
% calc and plot model error

%load Model
%RAdate = '10-Apr-2019';
RAdate = '20-May-2019';
loadFolder = [modelFolder RAdate filesep];
if bestSFflag
    load([loadFolder 'Errs_bestSF'],'Errs','sigmas','taus','whichpeaks')                   
else
    load([loadFolder 'Errs_constSF'],'Errs','sigmas','taus','whichpeaks')                   
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
    errmin = (min(min(Errs(:,:,ipeak))));
    if bestSFflag
        %caxis([errmin errmax])
    else
        caxis([errmin 2]);
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
    text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',40);
end
suptitle('normalized error between model and data')
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
end
if bestSFflag
    saveas(gcf,[FigFolder filesep 'ErrSpace_bestSF_' Exp],'pdf')
else
    saveas(gcf,[FigFolder filesep 'ErrSpace_constSF_' Exp],'pdf')
end

%% Experiment 3
cd('L:\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
cd(GHfolder)
Exp = 'Exp3';
%%      model
weightedFlag = false;
besSFflag = false;
%RA and Errs are calculated in script
%   L:\Experiments\N1P2\Analysis\N1P2_GH\Model_N1P2.m
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
if bestSFflag
    load([loadFolder 'Errs_bestSF'])                   
else
    load([loadFolder 'Errs_constSF'])                   
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
    errmin = (min(min(Errs(:,:,ipeak))));
    if bestSFflag
        %caxis([errmin errmax])
    else
        caxis([errmin 2]);
    end
    
    title(whichpeaks{ipeak})
    xlabel('\tau (sec)','fontsize',20);
    ylabel('\sigma (semitones)','fontsize',20)
    hold on
    disp([whichpeaks{ipeak} ': sigma = ' num2str(sigmas(row)) ', tau = ' num2str(taus(col)) ])
    plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',lw)
    plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',lw)
    set(gca,'fontsize',20)
    bestaus(ipeak) = taus(col);
    bestsigmas(ipeak) = sigmas(row);
    text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',40);
end
suptitle('normalized error between model and data')
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
end

if bestSFflag
    saveas(gcf,[FigFolder filesep 'ErrSpace_bestSF_' Exp],'pdf')
else
    saveas(gcf,[FigFolder filesep 'ErrSpace_constSF_' Exp],'pdf')
end

%% All Exp together!
%model calculated at script:
%L:\Experiments\N1P2\Analysis\N1P2_GH\Model_allExp.m

loadFolder = 'L:\Experiments\N1P2\Analysis\Model\AllExp';
Exp = 'All';
%% model
bestSFflag = true;
% loads Errs which were calculated in :
% L:\Experiments\N1P2\Analysis\N1P2_GH\Model_allExp.m
if bestSFflag
    load([loadFolder filesep 'Errs_bestSF'])                   
else
    load([loadFolder filesep 'Errs_constSF'])                   
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
    errmin = (min(min(Errs(:,:,ipeak))));
    if bestSFflag
        %caxis([errmin errmax])
    else
        caxis([errmin 2]);
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
    text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',40);
end
suptitle('normalized error between model and data')
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
end
if bestSFflag
    saveas(gcf,[FigFolder filesep 'ErrSpace_bestSF_' Exp],'pdf')
else
    saveas(gcf,[FigFolder filesep 'ErrSpace_constSF_' Exp],'pdf')
end

%% ellipsoid
% Run this after running any 'Model' part above
bestSFflag = false;
if bestSFflag
    load([loadFolder filesep 'Errs_bestSF'])                   
else
    load([loadFolder filesep 'Errs_constSF'])                   
end

tol = 0.1;
ERPfigure
for ipeak=1:length(whichpeaks)
    ax(ipeak) = subplot(1,2,ipeak);
    [row, col]=find(Errs(:,:,ipeak)==min(min(Errs(:,:,ipeak))));
    ih(ipeak) = imagesc(taus,sigmas,Errs(:,:,ipeak));
    errmin = (min(min(Errs(:,:,ipeak))));
    if bestSFflag
        %caxis([errmin errmax])
    else
        caxis([errmin 3]);
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
    text(4,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',40);
    
    %ellipsoid
    [rows, cols] = find(abs(Errs(:,:,ipeak) - errmin - 1) < tol);
    plot(taus(cols),sigmas(rows),'.','markersize',8,'Color',[1 0 0])
end
axes(ax(1))
plot(taus(cols),sigmas(rows),'.','markersize',12,'Color',[0 1 0]);

suptitle(Exp)
for ipeak = 1:length(whichpeaks)
    c(ipeak) = colorbar(ax(ipeak));
    ylabel(c(ipeak),'normalized error (STE units)','fontsize',16)
end