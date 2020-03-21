%N1 modulation
addpath('S:\Lab-Shared\NewDataArch\CommonResources\Tools\Matlab_Tools')
FigFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\PaperFigures\N1modulation';
MixedFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\MixedModel';
HandyFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\HandyStructures';
addpath(FigFolder)
addpath('S:\Lab-Shared\Experiments\MMNchroma\Analysis')
addpath('S:\Lab-Shared\Experiments\MMNchromaF\Analysis')
addpath('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
addpath('S:\Lab-Shared\Z backup\Tamar\fromZ\Documents\MATLAB\MatlabFunctions\downloaded');
colors = {[0.3 0.745 0.93],[0.75 0 0.75],[1 0 0],[0.75 0 0.75],[0.3 0.745 0.93]};
linewidths = [2,2,2,2,2];
linestyles = {':',':','-','-','-'};
fromxERP=-100;toxERP=300;fromyERP=-2.5;toyERP=3.2;
fromx=-0.1;tox=0.3;
srate = 512;
nBlockTypes=8;
nCond=5;
megaERP = nan(nBlockTypes,nCond,round((toxERP-fromxERP)/1000*srate));
fromymega=-2.3;toymega=3.2;

mb=0;%global index to count blocks combining all 3 exp
%% Experiment 1
cd('S:\Lab-Shared\Experiments\MMNchroma\Analysis')
Definitions_MMNchroma;%new idea for the first time! move definitions into a separate script!
ibls = [1 2];%relative to blocks in all experiments
ExpN=1;
%%      ERPs

if 1 %parameters
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP' 
    load('S:\Lab-Shared\Experiments\MMNchroma\Analysis\Properties')
    [ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );
    mode = 'lme';
    matrixdate = '06-Mar-2016';
    bslwin = [-0.1 0];
    cs = [1 2 5 3 4];
   
end

hf=ERPfigure;
set(hf,'Position',[100 100 460 300])

isp=0;
for bl = bls
    mb=mb+1;
    %bl
    isp=isp+1;
    subplot(1,length(bls),isp)
    load([mixedFolder filesep 'grandMatrix_condComb_' mode '_' blocks{bl} '_' matrixdate])
    bls=[2,4];
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;

%    r=rectangle('Position',[80, fromy, 40, toy-fromy],'facecolor',[0.7 0.7 0.7],'edgecolor','none');
    hold on
    conditions_ordered=conditions(cs);
    for c=cs
        i=find(cs==c);
        data = nanmean(grandMatrix(:,includeSubjects,c,:,electrode),4);
        meandata=nanmean(data,2);
        megaERP(mb,i,:)=meandata(1:size(megaERP,3));
        t=0:1:size(data,1)-1;t=t/srate;t=t+winbegin;t=t*1000;
        plot(t,nanmean(data,2),'DisplayName',conditions{c},'Color',colors{i},'linew',linewidths(i),'linestyle',linestyles{i})
        %plot(t,meandata,'LineWidth',3,'DisplayName',conditions{c},'Color',Colors{i})
        hold all
    end
    axis([fromxERP, toxERP, fromyERP, toyERP])
    set(gca,'fontsize',20)
    line([0,0],[fromyERP, toyERP],'linestyle','-.','col','k');
    line([fromxERP, toxERP],[0, 0],'linestyle','-.','col','k');
    box on
end

FigName = 'Exp1ERPs';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
%% plot mega ERP per Exp
load(['S:\Lab-Shared\Experiments\N1P2\Analysis\grandAverage' filesep 'megaERP'])


legendlabels={'tone 1','tone 2','tone 3','tone 4','tone 5'};
h=ERPfigure;
set(h,'Position',[100 200 400 300])
for con=1:nCond
    data = megaERP(ibls,con,:);
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    data = bsxfun(@minus,data,mean(data(:,:,newbsl),3));
    t=0:1/srate:(length(data)-1)/srate;t=t+winbegin;
    if con==1
        r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
    end
    hold on
    data=squeeze(data);
    hp(con)=plot(t,mean(data),'color',colors{con},'linew',3,'linestyle',linestyles{con});
    
    set(gca,'fontsize',20)
    line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
    line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
    ylim([fromymega toymega])
    set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
end
llidx=5:-1:1;
legend(hp(llidx),legendlabels(llidx),'location','nw','fontsize',16)
   
FigName = ['megaERP_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')

plot2svg([FigFolder filesep FigName],gcf)

%%      Bargraphs
% plot bargraph 5 cond N1
if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    addtag = '';
    peakdate = '14-Nov-2018';
    whichpeaks = {'N1','P2'};
    include = false;%include subjects for which there was a manual inspection
    order = [1,2,5,3,4];
end
%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
bls = [2, 4]; 

nCond = size(allGrandcon_times{1},2);

for ipeak = 1:2
    %:length(whichpeaks)
    %exclude subjects for which grandcon peak was calculated as median
    if include
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
        mss = [];
        for bl=1:bls
            allGrandcon_as_median{bl}(:,:,ipeak);
            [r, c]=find(allGrandcon_as_median{bl}(:,:,ipeak));
            mss = [mss, r'];
        end
        mss = unique(mss);
        includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
    end    
    peak_means = nan(length(bls),nCond);
    peak_errs = nan(length(bls),nCond);
    bpeaks = nan(length(includeSubjects),length(bls)*nCond);
    bxlocs = nan(length(includeSubjects),length(bls)*nCond);
    ib=0;
    for bl=bls
        ib=ib+1;
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(includeSubjects,:,ipeak);
        %cbpeaks(any(isnan(cbpeaks), 2), :) = [];%commented this out for
        cbpeaks = cbpeaks(:,order);
        %the smaller number of participants in classicctrl
        
        %into all block peaks
        bpeaks(:,((ib-1)*nCond+1):((ib-1)*nCond+nCond))=cbpeaks;
        bxlocs(:,((ib-1)*nCond+1):((ib-1)*nCond+nCond))=repmat(ib-1,[size(cbpeaks,1),nCond]);
        for con = 1:nCond
            %for scatter plot -
            bxlocs(:,((ib-1)*nCond+1)+con-1)=bxlocs(:,((ib-1)*nCond+1)+con-1)+0.15*(con-3);
            %for bargraph -
            peaks = cbpeaks(:,con);
            peak_means(ib,con) = nanmean(peaks);
            CI = Confidence(peaks);
            peak_errs(ib,con) = abs(nanmean(peaks)-CI(1));
        end        
    end
    bxlocs_noisy = bxlocs + ones(size(bxlocs)).*(1+(rand(size(bpeaks))-0.5)/100) ;
    
    hf = figure;
    set(hf,'Position',[100 100 450 300])

    hb = barwitherr(peak_errs, peak_means);% Plot with errorbar
    for con = 1:5
        hb(con).FaceColor = colors{con};
    end
    hold on
 
    set(gca,'xticklabels',blocks(bl))
    if ipeak==1
        ylim([-3.03 0])
    end
    xlim([0.5 2.5])
%      boxplot(bpeaks,'notch','on','Position',bxlocs(1,:)+1)
    
    set(gca,'fontsize',20)
    %title([whichpeaks{ipeak} ' of all conditions'],'fontsize',16)
    %legend({'tone 1','tone 2','tone 3','tone 4','tone 5'},'Location','northeastoutside')
    %ylabel(['mean amplitude ? CI 5-95%, \muV'])

FigName = 'Exp1bars';
saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'fig')
saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'pdf')

%save handy peaks:
allPeaks=reshape(bpeaks,[size(bpeaks,1),size(bpeaks,2)/length(bls),length(bls)]);
save([HandyFolder filesep 'Exp' num2str(ExpN) '_' whichpeaks{ipeak}],'allPeaks')
end

%%      LME model
savedate = '13-May-2019';%peaks calculated at the Amplitudes_MMNchroma script at:
% S:\Lab-Shared\Experiments\MMNchroma\Analysis\Amplitudes_MMNchroma.m
include = false;%include subjects for which time of grandcon peak was determined due to median

load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
formulas = {'N1~dist_mean+(1|subject)','P2~dist_mean+(1|subject)'};

lmes = cell(length(blocks),length(formulas));

if include
    %allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
    mss = [];%known
%     for bl=1:length(allGrandcon_as_median)
%         for ipeak = 1:2%here I exclude from all analysis any subject that either had a problem in the N1 or in the P2s
%             [r, c]=find(allGrandcon_as_median{bl}(:,:,ipeak));
%             mss = [mss, r'];
%         end
%     end
%     mss = unique(mss);
    %allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end

for ib=bls
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
    for i=1:length(formulas)
        lmes{ib,i} = fitlme(allTables.(blocks{ib}),formulas{i});
        disp(lmes{ib,i}.Formula)
        disp(num2str(lmes{ib,i}.coefTest))
    end
end

%% Experiment 2
cd('S:\Lab-Shared\Experiments\MMNchromaF\Analysis')
Definitions_MMNchromaF;

addpath('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
bls = 2;
cs = [1 2 5 3 4];
ibls=3;
ExpN=2;
%%      ERPs
% plot waves on top of each other
tag = 'BP1-20';
blocks = {'MMN_F_pure','MMN_Fctrl_pure','MMN_classic_pure','MMN_F_shep','MMN_Fctrl_shep','MMN_classic_shep'};
matrixdate = '26-Apr-2017';
bslwin = [-0.1,0];
pwin = [85,150];%ms
electrodeName = 'Cz';
[ electrode ] = chName2n( electrodeName );

hf=ERPfigure;
set(hf,'Position',[100 100 200 300])

cs=[2,3,1,4,5];
for bl=2
    mb=mb+1;
    load([grandFolder filesep 'grandMatrix_' blocks{bl} '_' tag])
 %   load([mixedFolder filesep 'grandMatrix_condComb_' mode '_' blocks{bl} '_' matrixdate])
    matrixdate = '26-Apr-2017';
    %load([grandFolder filesep 'grandMatrixAll_' blocks{bl} '_notes_12-Feb-2016'])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
    %plot
    for con=cs
        i=find(con==cs);
        plot(t,squeeze(nanmean(grandMatrix(:,includeSubjects,con,electrode),2)),'Color',colors{i},'linew',linewidths(i),'linestyle',linestyles{con==cs})
               
         meandata=squeeze(nanmean(grandMatrix(:,:,con,electrode),2));
         megaERP(mb,i,:)=meandata(1:size(megaERP,3));
        
        hold on
    end  
  %  legend(conditions{cs})
   % title(blocks{bl})
    ax(bl)=gca;
    line([0,0],[fromyERP, toyERP],'linestyle','-.','col','k');
    line([fromxERP, toxERP],[0, 0],'linestyle','-.','col','k');
end
linkaxes(ax)
axis([fromxERP, toxERP, fromyERP, toyERP])
set(gca,'fontsize',20)

FigName = 'Exp2ERPs';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
%% plot mega ERP per Exp
load(['S:\Lab-Shared\Experiments\N1P2\Analysis\grandAverage' filesep 'megaERP'])

legendlabels={'tone 1','tone 2','tone 3','tone 4','tone 5'};
h=ERPfigure;
set(h,'Position',[100 200 400 300])
for con=1:nCond
    data = megaERP(ibls,con,:);
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    data = bsxfun(@minus,data,mean(data(:,:,newbsl),3));
    t=0:1/srate:(length(data)-1)/srate;t=t+winbegin;
    if con==1
        r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
    end
    hold on
    data=squeeze(data);
    hp(con)=plot(t,(data),'color',colors{con},'linew',3,'linestyle',linestyles{con});
    
    set(gca,'fontsize',20)
    line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
    line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
    ylim([fromymega toymega])
    set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
end
llidx=5:-1:1;
legend(hp(llidx),legendlabels(llidx),'location','nw','fontsize',16)
   
FigName = ['megaERP_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')

plot2svg([FigFolder filesep FigName],gcf)

%%      Bargraphs
% plot bargraph 5 cond N1
if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    addtag = '';
    peakdate = '16-Nov-2018';
    whichpeaks = {'N1','P2'};
    include = false;%include subjects for which there was a manual inspection
    order = [1,2,5,3,4];
end
%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
bls=1:2;
nCond = size(allGrandcon_times{1},2);

for ipeak = 1:2
    %:length(whichpeaks)
    %exclude subjects for which grandcon peak was calculated as median
    if include
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
       mss=31;
        includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
    end    
    peak_means = nan(length(bls),nCond);
    peak_errs = nan(length(bls),nCond);
    bpeaks = nan(length(includeSubjects),length(bls)*nCond);
    bxlocs = nan(length(includeSubjects),length(bls)*nCond);
    ib=0;
    for bl=bls
        ib=ib+1;
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(includeSubjects,:,ipeak);
        %cbpeaks(any(isnan(cbpeaks), 2), :) = [];%commented this out for
        cbpeaks = cbpeaks(:,order);
        %the smaller number of participants in classicctrl
        
        %into all block peaks
        bpeaks(:,((ib-1)*nCond+1):((ib-1)*nCond+nCond))=cbpeaks;
        bxlocs(:,((ib-1)*nCond+1):((ib-1)*nCond+nCond))=repmat(ib-1,[size(cbpeaks,1),nCond]);
        for con = 1:nCond
            %for scatter plot -
            bxlocs(:,((ib-1)*nCond+1)+con-1)=bxlocs(:,((ib-1)*nCond+1)+con-1)+0.15*(con-3);
            %for bargraph -
            peaks = cbpeaks(:,con);
            peak_means(ib,con) = nanmean(peaks);
            CI = Confidence(peaks);
            peak_errs(ib,con) = abs(nanmean(peaks)-CI(1));
        end        
    end
    bxlocs_noisy = bxlocs + ones(size(bxlocs)).*(1+(rand(size(bpeaks))-0.5)/100) ;
    
    hf = figure;
    set(hf,'Position',[100 100 250 300])
    hb = barwitherr(peak_errs, peak_means);% Plot with errorbar
    for con = 1:5
        hb(con).FaceColor = colors{con};
    end
    hold on
    %add all participants:
%     for ii=1:size(bxlocs,1)
%         scatter(bxlocs_noisy(ii,:),bpeaks(ii,:),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
%     end
    set(gca,'xticklabels',{'chroma F','ctrl'})
    if ipeak==1
        ylim([-3.03 0])
    end
    xlim([1.5 2.5])
%      boxplot(bpeaks,'notch','on','Position',bxlocs(1,:)+1)
    
    set(gca,'fontsize',20)
    %title([whichpeaks{ipeak} ' of all conditions'],'fontsize',16)
    %legend({'tone 1','tone 2','tone 3','tone 4','tone 5'},'Location','northeastoutside')
    %ylabel(['mean amplitude ? CI 5-95%, \muV'])
    %save handy peaks:
allPeaks=reshape(bpeaks,[size(bpeaks,1),size(bpeaks,2)/length(bls),length(bls)]);
save([HandyFolder filesep 'Exp' num2str(ExpN) '_' whichpeaks{ipeak}],'allPeaks')

end

FigName = 'Exp2bars';
saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'fig')
saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'pdf')

%%      LME model
savedate = '13-May-2019';
include = false;%include subjects for which time of grandcon peak was determined due to median

load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')

formulas = {'N1~dist_mean+(1|subject)','P2~dist_mean+(1|subject)'};

lmes = cell(length(blocks),length(formulas));

if include
    %allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
%     mss = [];
%     for bl=1:length(allGrandcon_as_median)
%         for ipeak = 1:2%here I exclude from all analysis any subject that either had a problem in the N1 or in the P2s
%             [r, c]=find(allGrandcon_as_median{bl}(:,:,ipeak));
%             mss = [mss, r'];
%         end
%     end
%     mss = unique(mss);

    mss = 31;%known
    %allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end

for ib=bls
    
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
    for i=1:length(formulas)
        lmes{ib,i} = fitlme(allTables.(blocks{ib}),formulas{i});
        disp(lmes{ib,i}.Formula)
        disp(num2str(lmes{ib,i}.coefTest))
   
    end
end

%% Experiment 3
cd('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
cd(GHfolder)
ibls=[4:8];
iblss={[4],[5 6],[7 8]};
ExpN=3;
%%      ERPs
% load and plot grand 5 conds -
matrixdate='29-Oct-2018';
tag = 'Bp1-20';
bslwin = [-0.1,0];
electrodeName = 'Cz';
[ electrode ] = chName2n( electrodeName );

h=ERPfigure;
set(h,'Position',[100 100 1200 300])

for bl = 1:5
    mb=mb+1;
    subplot(1,5,bl)
    load([grandFolder filesep 'grandMatrix_condComb_' tag '_Bl' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
     %plot
     r=rectangle('Position',[0.08, fromyERP, 0.04, toyERP-fromyERP],'facecolor',[0.7 0.7 0.7],'edgecolor','none');

    for con=1:size(grandMatrix,3)
        data = nanmean(grandMatrix(:,includeSubjects,con,:,electrode),4);
        meandata=nanmean(data,2);
        ph(con)=plot(t,meandata,'Color',colors{con},'linew',linewidths(con),'linestyle',linestyles{con});
        megaERP(mb,con,:)=meandata(1:size(megaERP,3));
         
        hold on
    end  
    title(blocks{bl})
    ax(bl)=gca;
    set(gca,'fontsize',20)
    line([0,0],[fromyERP, toyERP],'linestyle','-.','col','k');
    line([fromxERP, toxERP],[0, 0],'linestyle','-.','col','k');
    box on
end
hl=legend([ph(5),ph(4),ph(3),ph(2),ph(1)],{' tone 5',' tone 4',' tone 3',' tone 2',' tone 1'});
set(hl,'Position',[0.92 0.4 0.07 0.4])
suptitle(['N = ' num2str(length(whichSubjects)) '. Electrode ' electrodeName ])
linkaxes(ax)
axis([fromxERP, toxERP, fromyERP, toyERP])

FigName = 'Exp3ERPs';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
%% plot mega ERP per Exp
load(['S:\Lab-Shared\Experiments\N1P2\Analysis\grandAverage' filesep 'megaERP'])

legendlabels={'tone 1','tone 2','tone 3','tone 4','tone 5'};
h=ERPfigure;
set(h,'Position',[100 200 1200 300])
for isp=1:length(iblss)
    subplot(1,length(iblss),isp)
    for con=1:nCond
        ibls=iblss{isp};
        data = megaERP(ibls,con,:);
        newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
        data = bsxfun(@minus,data,mean(data(:,:,newbsl),3));
        t=0:1/srate:(length(data)-1)/srate;t=t+winbegin;
        if con==1
            r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
            r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        end
        hold on
        hp(con)=plot(t,squeeze(mean(data,1)),'color',colors{con},'linew',3,'linestyle',linestyles{con});

        set(gca,'fontsize',20)
        line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
        line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
        ylim([fromymega toymega])
        xlim([fromx tox])
        set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
    end
end
llidx=5:-1:1;
legend(hp(llidx),legendlabels(llidx),'location','nw','fontsize',16)
   
FigName = ['megaERP_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')

plot2svg([FigFolder filesep FigName],gcf)

%%      Bargraphs
% plot bargraph 5 cond N1
if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    addtag = '';
    peakdate = '30-Oct-2018';
    whichpeaks = {'N1','P2'};
    include = false;%include subjects for which there was a manual inspection
end
%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
nCond = size(allGrandcon_times{1},2);

for ipeak = 1:length(whichpeaks)
    %:length(whichpeaks)
    %exclude subjects for which grandcon peak was calculated as median
    if include
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
       
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
    end    
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(length(includeSubjects),length(blocks)*nCond);
    bxlocs = nan(length(includeSubjects),length(blocks)*nCond);
    for bl=1:length(blocks)
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(includeSubjects,:,ipeak);
        cbpeaks(any(isnan(cbpeaks), 2), :) = [];
        %into all block peaks
        bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=cbpeaks;
        bxlocs(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=repmat(bl-1,[size(cbpeaks,1),nCond]);
        for con = 1:nCond
            %for scatter plot -
            bxlocs(:,((bl-1)*nCond+1)+con-1)=bxlocs(:,((bl-1)*nCond+1)+con-1)+0.15*(con-3);
            %for bargraph -
            peaks = cbpeaks(:,con);
            peak_means(bl,con) = mean(peaks);
            CI = Confidence(peaks);
            peak_errs(bl,con) = abs(nanmean(peaks)-CI(1));
        end        
            
    end
    bxlocs_noisy = bxlocs + ones(size(bxlocs)).*(1+(rand(size(bpeaks))-0.5)/100) ;
    
    hf = figure;
    set(hf,'Position',[100 100 1200 300])
    hb = barwitherr(peak_errs, peak_means);% Plot with errorbar
    for con = 1:5
        hb(con).FaceColor = colors{con};
    end
    hold on
    %add all participants:
%     for ii=1:size(bxlocs,1)
%         scatter(bxlocs_noisy(ii,:),bpeaks(ii,:),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
%     end

   % set(gca,'xticklabels',blocks)
    set(gca,'xticklabels',{'1','2a','2b','3a','3b'})
%     boxplot(bpeaks,'Position',bxlocs(1,:)+1)
    
    set(gca,'fontsize',20)
    title([whichpeaks{ipeak} ' of all conditions. Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))],'fontsize',16)
    hl=legend([hb(5),hb(4),hb(3),hb(2),hb(1)],{' tone 5',' tone 4',' tone 3',' tone 2',' tone 1'});
    set(hl,'Position',[0.92 0.4 0.06 0.4])

    ylabel(['mean amplitude ? CI 5-95%, \muV'])
    
FigName = 'Exp3bars';
saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'fig')
saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'pdf')

%save handy peaks:
allPeaks=reshape(bpeaks,[size(bpeaks,1),size(bpeaks,2)/length(bls),length(bls)]);
save([HandyFolder filesep 'Exp' num2str(ExpN) '_' whichpeaks{ipeak}],'allPeaks')

end

%%      LME model
savedate = '31-Oct-2018';
include = false;%include subjects for which time of grandcon peak was determined due to median

load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
blocks = {'b1','b2a','b2b','b3a','b3b'};
formulas = {'N1~dist_mean+(1|subject)','P2~dist_mean+(1|subject)'};

lmes = cell(length(blocks),length(formulas));

if include
    allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
%     mss = [];
%     for bl=1:length(allGrandcon_as_median)
%         for ipeak = 1:2%here I exclude from all analysis any subject that either had a problem in the N1 or in the P2s
%             [r, c]=find(allGrandcon_as_median{bl}(:,:,ipeak));
%             mss = [mss, r'];
%         end
%     end
%     mss = unique(mss);
    mss=[5 14];%known
    allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end

for ib=[1,2,3,4,5]
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
    for i=1:length(formulas)
        lmes{ib,i} = fitlme(allTables.(blocks{ib}),formulas{i});
         disp(lmes{ib,i}.Formula)
        disp(num2str(lmes{ib,i}.coefTest))
   
    end
end

%% plot MEGA ERP
%Note: must run all bargraphs plots orderly in order to calc megaERP
%once saved, can be loaded
%save([grandFolder filesep 'megaERP'],'megaERP')
load([grandFolder filesep 'megaERP'])

legendlabels={'tone 1','tone 2','tone 3','tone 4','tone 5'};
h=ERPfigure;
fromymega=-2;toymega=3;
for con=1:nCond
    data = squeeze(megaERP(:,con,:));
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    data = bsxfun(@minus,data,mean(data(:,newbsl),2));
    t=0:1/srate:(length(data)-1)/srate;t=t+winbegin;
    if con==1
        r1=rectangle('Position',[0.08, fromyERP, 0.04, toyERP-fromyERP],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        r2=rectangle('Position',[0.15, fromyERP, 0.04, toyERP-fromyERP],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
    end
    hold on
    hp(con)=varplot(t,data','ci',0.9,'color',colors{con},'linew',3,'linestyle',linestyles{con});
    
    title('mega ERP')
    set(gca,'fontsize',20)
    line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
    line([fromxERP, toxERP],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
    ylim([fromymega toymega])
    set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
end
llidx=5:-1:1;
legend(hp(llidx),legendlabels(llidx),'location','nw','fontsize',16)
   
FigName = 'megaERP';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
saveas(gcf,[FigFolder filesep FigName],'svg')

plot2svg([FigFolder filesep FigName],gcf)

%% LMEmodels using BIG Table
load([MixedFolder filesep 'theTable'])
formulas = {'Voltage ~ size_jump*Potential + (1|subject) + (1|ExpN)',...
            'Voltage ~ Potential*spread*dist_mean + (1|subject) + (1|ExpN)',...     
            'Voltage ~ dist_mean + (1|subject) + (1|ExpN)',...
            'Voltage ~ Potential*size_jump + (1|subject) + (1|ExpN)',...
            'Voltage ~ size_jump + (1|subject) + (1|ExpN)',...
            'Voltage ~ Potential*dist_mean + Potential*size_jump + (1|subject) + (1|ExpN)',...
            'Voltage ~ dist_mean + Potential*size_jump + (1|subject) + (1|ExpN)',...
            'Voltage ~ Potential*dist_mean + size_jump + (1|subject) +(1|ExpN)',...
            'Voltage ~ dist_mean + size_jump + (1|subject) + (1|ExpN)',...
            };

lmes = cell(length(formulas),1);

for i=1:length(formulas)
    %lmes{i} = fitlme(theTable(theTable.ExpN==categorical(3),:),formulas{i});
    lmes{i} = fitlme(theTable,formulas{i});
    lmes{i}.anova;
end
disp('Compare')
compare(lmes{3},lmes{1});
