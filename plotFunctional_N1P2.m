function hf = plotFunctional_N1P2( allTables, relCols, blocks, ibs )
% Tamar, Nov 13 2018
% a function of the plot within functionalPlots_GH

% relCols - relevant comuns in allTables

%params
whichpeaks = {'N1','P2'};
markers = {'o','+','d','*','s'};
cmap = [1 0 1;...
        1 0 0;...
        1 0.6 0;...
        0 0.7 0;...
        0 0 1 ];


hf=ERPfigure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])

xvars = {'currMIDI','dist_mean','size_jump'};
xnames = {'curr. freq.','dist. from mean freq.','dist. from prev. freq.'};
nx=length(xvars);

for ipeak = 1:2
    peak = whichpeaks{ipeak};
    for ix=1:length(xvars)
        subplot(2,3,(ipeak-1)*nx+ix)
        x=nan(1,5*4*length(ibs));y=nan(1,5*4*length(ibs));
        iib=0;R=cell(length(ibs));p=cell(length(ibs));
        for ib=ibs
            iib=iib+1;
            T = allTables.(blocks{ib});
            t = T(1,relCols);
            nlines = height(T(T.subject==T{1,1},:));
            meanPeaks = nan(nlines,1);
            errPeaks = nan(nlines,1);
            for line = 1:nlines
                tl = T(line,:);
                tlallS = T(T.currMIDI==tl.currMIDI & T.prevMIDI==tl.prevMIDI,:);
                meanPeaks(line) = mean(tlallS.(peak));
                CI = Confidence(tlallS.(peak));
                errPeaks(line) = abs(meanPeaks(line)-CI(1));
                t(line,:) = T(line,relCols);
            end
            from = (iib-1)*nlines+1;
            to = iib*nlines;
            x(from:to) = t.(xvars{ix});
            
            switch peak
                case 'N1'
                    y(from:to) = -meanPeaks;
                case 'P2'
                    y(from:to) = meanPeaks;
            end
            plot(x(from:to),y(from:to),markers{ib},'MarkerSize',8,'MarkerEdgeColor',cmap(ib,:),'linewidth',1);
            %plot(x(from:to),y(from:to),'o','MarkerSize',4,'MarkerFaceColor',cmap(ib,:),'MarkerEdgeColor',cmap(ib,:));
            hold on
            p1=polyfit(x(from:to),y(from:to),1);
            f1=polyval(p1,x(from:to));
            plot(x(from:to),f1,'Color',cmap(ib,:),'HandleVisibility','off','linewidth',1)
            [R{iib},p{iib}]=corrcoef(x(from:to),y(from:to));
        end
        ylabel([peak ' [\muV]'])
        xlabel([xnames{ix} ' [MIDI]'])
        set(gca,'fontsize',16)

        hold on
        p1=polyfit(x,y,1);
        f1=polyval(p1,x);
        plot(x,f1,'k','linewidth',3)
        
        xlim = get(gca,'Xlim');
        dx=(xlim(2)-xlim(1))*0.45;
        ylim = get(gca,'Ylim');
        dy=abs((ylim(2)-ylim(1)))*0.07;
        iib=0;
        for ib=ibs
            iib=iib+1;
            text(xlim(2)-dx,ylim(2)-(iib)*dy,['R=' num2str(R{iib}(1,2),2) ', p=' num2str(p{iib}(1,2),2)],'Color',cmap(ib,:),'fontsize',18);
        end
        %[xs,I] = sort(x);
        %y=y(I);
        [Rtot,ptot]=corrcoef(x,y);
        text(xlim(2)-dx,ylim(2)-(iib+1)*dy,['R=' num2str(Rtot(1,2),2) ', p=' num2str(ptot(1,2),2)],'Color','k','fontsize',18);
%         p2=polyfit(x,y,2);
%         f2=polyval(p2,x);
%         plot(x,f2,'.b','MarkerSize',20)
    end
end
legendstring = {};
for iib=1:length(ibs)
    legendstring(iib)={['b' num2str(ibs(iib))]};
end
    
lh=legend(legendstring,'linear fit');
set(lh,'position',[0.93 0.33 0.05 0.1],'fontsize',18)

end

