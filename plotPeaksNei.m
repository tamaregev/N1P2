function plotPeaksNei(peaksFolder,ExpN,includeSubjects,ibls,pos)
%this function is used in P2 modulation Figure
%transformation matrix:
load([peaksFolder filesep 'peaks_allExp'])
readfrom = NPallexp_peak_meanz;
errsfrom = NPallexp_peakz;
%peaks_allExp are currently saved in the following way -
%run each experiments definitions and then the P2 Bargraph section. Saved
%at end of Exp 3 Bargraph section

nCond=5;
whichpeaks={'N1','P2'};
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};
peakColors = {[0 0 1],[1 0 0]};
lw=2;lwe=1;
fromynei=-0.22;toynei=0.5;

h=ERPfigure;
set(h,'Position',pos)
nnei=2;%plot until this neighbor
cmap=parula(4);
for isp=1:length(ibls)
    subplot(1,length(ibls),isp)
    %calc only for current ibl
    ibl=ibls(isp);
    neiPeak = nan(length(legendlabels),length(ibl)*8,length(whichpeaks)); 
    neiPeakall = nan(length(includeSubjects),length(legendlabels),length(ibl)*8,length(whichpeaks));
    
    count = zeros(length(legendlabels),1);
    for con=1:nCond
        for prevcon = 1:nCond
            nei=transM(con,prevcon);
            if nei
                for ipeak=1:length(whichpeaks)
                    neiPeak(nei,count(nei)+1,ipeak) = readfrom{ipeak,ExpN}(isp,con,prevcon);
                    neiPeakall(:,nei,count(nei)+1,ipeak) = errsfrom{ipeak,ExpN}(:,isp,con,prevcon);
                 end
                count(nei)=count(nei)+1;
            end
        end
    end

    %plot
    clear hp
    ddatas=nan(1,2);derrdatas=nan(size(neiPeakall,1),2);
    for ipeak = 1:length(whichpeaks)
        data = squeeze(nanmean(neiPeak(1:nnei,:,ipeak),2));
        ddata=data(2)-data(1);ddatas(ipeak)=ddata;
%         %calc CI across participants:
         errdata = squeeze(nanmean(neiPeakall(:,1:nnei,:,ipeak),3));
         derrdata = errdata(:,2) - errdata(:,1);derrdatas(:,ipeak)=derrdata;
%         errs = squeeze(nanstd(errdata)/sqrt(length(includeSubjects)));
%         errs=nanmean(errs,2);
        %calc error across conditions:
%        errdata = neiPeak(1:nnei,:,ipeak);
%         %ste across conditions:
%         errs = nanstd(errdata,0,2)./sqrt(sum(~isnan(errdata),2));
        %CI across conditions:
%         CI=nan(nnei,2);
%         for in=1:nnei
%             CI(in,:) = Confidence(derrdata);
%         end
%         errs = CI(:,2) - data;
        %CI across participants:
        CI = Confidence(derrdata);
        err = CI(2) - ddata;errss(ipeak)=err;
        hold on
        hp(ipeak) = bar(ipeak,ddata,'FaceColor',peakColors{ipeak},'linewidth',lw);% Plot with errorbar
        
        %hp(ipeak) = barwitherr(errs(2),ipeak, data(2)-data(1),'FaceColor',peakColors{ipeak});% Plot with errorbar
        %These would work for presenting more neighbors:
        %hp(ipeak)=plot([1 2],data-data(1),'Color',peakColors{ipeak},'linewidth',lw);
        %plot([1 2],data-data(1),'.','Color',peakColors{ipeak},'MarkerSize',26,'HandleVisibility','off')
        errorbar(ipeak,ddata,err,':','Color','k','linewidth',lwe,'HandleVisibility','off')
        ylim([fromynei toynei])
        xlim([0 nnei+1])
        if isp==1
           ylabel({'Z-scored Voltage (\muV)', 'relative to 1 neighbor'})
       end
        set(gca,'xtick',[1 2],'xticklabels',{'1', '2'})
        xlabel('Neighbor')
        set(gca,'fontsize',16)
    end
%      hp(3) = bar(3,ddatas(2)-ddatas(1),'FaceColor',[0.5 0.5 0.5],'linewidth',lw);% Plot with errorbar
%      CI=Confidence(derrdatas(:,2)-derrdatas(:,1));
%      err=CI(2)-(ddatas(2)-ddatas(1));
%      errorbar(3,ddatas(2)-ddatas(1),err,':','Color','k','linewidth',lwe,'HandleVisibility','off')
%        
    
end
suptitle(['Exp ' num2str(ExpN)])        
legend(hp,whichpeaks)
end