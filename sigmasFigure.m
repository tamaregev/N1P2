%% sigmasFigure
FigFolder = 'L:\Experiments\N1P2\Analysis\N1P2_GH\PaperFigures\adaptiveWidth';
lw=2;

%% Experiment 3
cd('L:\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
cd(GHfolder)
Exp = 'Exp3';
%%      model
bestSFflag = false;
weightedFlag = false;
givenTauflag = true;

if bestSFflag
    givenTaus = {2.8,0.4};
else
    givenTaus = {2,0.6};
end
tol=0.001;
%RA and Errs are calculated in script
%   L:\Experiments\N1P2\Analysis\N1P2_GH\Model_N1P2.m
%under %% fitting parameters to each ~~condition~~ separately
spreads = 1:3;%these are the 3 conditions
block2cond = [1,2,2,3,3];

RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
if bestSFflag
    load([loadFolder 'Errs_spreads_bestSF'])                   
else
    load([loadFolder 'Errs_spreads_constSF'])
end
bestParams = nan(length(spreads),length(whichpeaks),2);
ERPfigure;
set(gcf,'units','normalized','outerposition',[0 0 0.7 0.7])
ii=0;
errmin = min(min(min(Errs(:,:,:,:))));
errmax = max(max(max(Errs(:,:,:,:))));

for ipeak=1:length(whichpeaks)
   for ispread=1:length(spreads)
        ii=ii+1;        
        ax(ii) = subplot(2,3,ii);
        if givenTauflag
            col = find(abs(taus - givenTaus{ipeak}) < tol);
            row = find(Errs(:,col,ipeak,ispread)== min(Errs(:,col,ipeak,ispread)));
        else
            [row, col] = find(Errs(:,:,ipeak,ispread)==min(min(Errs(:,:,ipeak,ispread))));
        end
        imagesc(taus(1:25),sigmas,Errs(:,1:25,ipeak,ispread));
        if bestSFflag
            %caxis([errmin errmax])
        else
            errmin = min(min(min(Errs(:,:,ipeak,:))));
            caxis([errmin 2])
        end
        hold on 
        plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',2)
        plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',2)
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

if bestSFflag
    saveas(gcf,[FigFolder filesep 'ErrSpace_spreads_bestSF_' Exp],'pdf')
else
    saveas(gcf,[FigFolder filesep 'ErrSpace_spreads_constSF_' Exp],'pdf')
end

%%% plot best fit sigma
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

if bestSFflag
    saveas(gcf,[FigFolder filesep 'bestSig_spreads_bestSF_' Exp],'pdf')
else
    saveas(gcf,[FigFolder filesep 'bestSig_spreads_constSF_' Exp],'pdf')
end

