%P2 modelation
FigFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\PaperFigures\P2modulation';
MixedFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\MixedModel';
peaksFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\PeaksFolder';
%mkdir(FigFolder)
cd(FigFolder)
colors = {[0.3 0.745 0.93],[0.75 0 0.75],[1 0 0],[0.75 0 0.75],[0.3 0.745 0.93]};
peakColors = {[0 0 1],[1 0 0]};
linewidths = [2,2,2,2,2];
linestyles = {'-.','-.','-','-','-'};
fromx=-100;tox=300;fromy=-2.5;toy=3.2;
fromynei=-0.2;toynei=0.5;
fromymega=-2;toymega=3.2;
NPallexp_peak_means=cell(2,3);
NPallexp_peak_meanz=cell(2,3);
lw=3;lwe=1;
nBlockTypes = 8;%all exp together
nCond=5;
srate=512;
megacondERP = nan(nBlockTypes,nCond,nCond,round((tox-fromx)/1000*srate));
mb=0;%global index to count blocks combining all 3 exp
whichpeaks={'N1','P2'};

%% Experiment 1
ExpN=1;
cd('S:\Lab-Shared\Experiments\MMNchroma\Analysis')
addpath('S:\Lab-Shared\Experiments\MMNchroma\Analysis')
Definitions_MMNchroma;%new idea for the first time! move definitions into a separate script!
bls=[2 4];
ibls=[1 2];
%%   P2 Bargraphs
% plot bargraph P2 prevcond
if 1%params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    addtag = '';
    peakdate = '14-Nov-2018';
    whichpeaks = {'N1','P2'};
    include = false;%include subjects for which there was a manual inspection
    order = [1,2,5,3,4];
end
%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
nCond = size(allGrandcon_times{1},2);
nPrevCond = size(allPeak_amps{1},3);

bls = [2, 4];
for ipeak = 1:2 
    %1:length(whichpeaks)
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
    nsub = length(includeSubjects);
    peaks = nan(nsub,length(bls),nCond,nPrevCond);
    peak_means = nan(length(bls),nCond,nCond);
    peak_meanz = nan(length(bls),nCond,nCond);
    
    peak_errs = nan(length(bls),nCond,nCond);
    h=ERPfigure;
    set(h,'Position',[10 100 400 600])
     
        
    for bl= bls  
        %bls
        ib=find(bls==bl);
        for con = 1:size(allPeak_amps{1},2)
            for prevcon = 1:size(allPeak_amps{1},3)
                peaks(:,ib,con,prevcon) = allPeak_amps{bl}(includeSubjects,con,prevcon,ipeak);
                peak_means(ib,con,prevcon) = nanmean(peaks(:,ib,con,prevcon),1);
                CI = Confidence(peaks(:,ib,con,prevcon));
                peak_errs(ib,con,prevcon) = abs(nanmean(peaks(:,ib,con,prevcon))-CI(1));
                %peak_errs(ib,con,prevcon) = nanstd(peaks(:,ib,con,prevcon))/sqrt(size(peaks,1));
            end
        end
        peak_means(ib,:,:) = peak_means(ib,order,order);
        peak_errs(ib,:,:) = peak_errs(ib,order,order);
        
        for i=1:size(allPeak_amps{1},2)
           
            subplot(size(allPeak_amps{1},2),length(bls),(5-i)*length(bls)+ib)
            bar(squeeze(peak_means(ib,i,:)),'FaceColor',colors{i})
           % pause
            barwitherr(squeeze(peak_errs(ib,i,:)),squeeze(peak_means(ib,i,:)),'FaceColor',colors{i})
            hold all
            if ib==1
               % ylabel(['tone ' num2str(i)])
               
            end
            if ib~=1 || i~=1
                set(gca,'xticklabels',{})
                set(gca,'yticklabels',{})
            end
            set(gca,'fontsize',16)
            if ipeak==1
                axis([0 6 -3.5 0])
                %xlim([0 6])
            elseif ipeak==2
                %xlim([0 6])
                axis([0 6 0 4.5])
            end
        end
       
        %suptitle([whichpeaks{ipeak} '. Electrode: ' electrodeName '. Block: ' num2str(bl)])
    
    end
    %     
    FigName = 'Exp1P2bars';
    saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'fig')
    saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'pdf')
    
    %save for all exp together:
    pp = peaks(:);
    if ipeak==1
        pp=-1*pp;
        pp(~isnan(pp))=zscore(pp(~isnan(pp)));
        peakz = reshape(pp,size(peaks));
    else
        pp(~isnan(pp))=zscore(pp(~isnan(pp)));
        peakz = reshape(pp,size(peaks));
    end
    ib=0;
    for bl= bls
        ib=ib+1;
        %bls
        for con = 1:size(allPeak_amps{1},2)
            for prevcon = 1:size(allPeak_amps{1},3)
                peak_meanz(ib,con,prevcon) = nanmean(peakz(:,ib,con,prevcon),1);
            end
        end
    end
    NPallexp_peak_means{ipeak,ExpN}=peak_means; 
    NPallexp_peak_meanz{ipeak,ExpN}=peak_meanz; 
    NPallexp_peaks{ipeak,ExpN}=peaks;
    NPallexp_peakz{ipeak,ExpN}=peakz;

end
% % calculated: 
% % NPallexp_peak_means;%{ipeak x iExp}(ibls x cond x prevcond)
% % NPallexp_peak_meanz;%{ipeak x iExp}(ibls x cond x prevcond)
%% plot peaks of neighbors
pos=[100 200 500 300];
plotPeaksNei(peaksFolder,ExpN,includeSubjects,ibls,pos)

FigName = ['peaksNei_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)

%% calc megacondERP

if 1 %parameters
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP' 
    load('S:\Lab-Shared\Experiments\MMNchroma\Analysis\Properties')
    [ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );
    mode = 'lme';
    matrixdate = '06-Mar-2016';
    bslwin = [-0.1 0];
    order = [1 2 5 3 4];
end

for bl = bls
    mb=mb+1;
    
    ib=find(bl==bls);
    load([mixedFolder filesep 'grandMatrix_condComb_' mode '_' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
    conditions_ordered=conditions(order);
    ic=0;
    for c=1:5
        ic=ic+1;
        ipc=0;
        for pc=1:5
            ipc=ipc+1;
            data = grandMatrix(:,includeSubjects,order(c),order(pc),electrode);
            meandata=nanmean(data,2);
            megacondERP(mb,ic,ipc,:)=meandata(1:size(megacondERP,4));
        end
    end
end
%% plot megacondERP Exp 1 - per block
%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
      
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};


h=ERPfigure;
set(h,'Position',[100 200 500 300])

nnei=2;%plot until this neighbor
cmap=parula(4);
for isp=1:length(ibls)
    subplot(1,length(ibls),isp)
    %calc only for current ibl
    ibl=ibls(isp);
    neiERP = nan(length(legendlabels),length(ibl)*8,size(megacondERP,4)); 
    count = zeros(length(legendlabels),1);
    for con=1:nCond
        for prevcon = 1:nCond
            nei=transM(con,prevcon);
            if nei
                data = squeeze(megacondERP(ibl,con,prevcon,:));
                newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
                data = bsxfun(@minus,data,mean(data(newbsl,:),1));
                neiERP(nei,count(nei)+1,:)=data;
                count(nei)=count(nei)+1;
            end
        end
    end

    %plot
    clear hp
    for nei=1:nnei
        if nei==1
            r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
            r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        end
        hold on
        data = squeeze(neiERP(nei,:,:));
        t=0:1/srate:(length(data)-1)/srate;t=t-0.1;
        %hp(nei)=varplot(t,data','ci',0.9,'color',cmap(nei,:),'linew',3);
        hp(nei)=plot(t,nanmean(data),'color',cmap(nei,:),'linew',3);
        title('mega ERP')
        set(gca,'fontsize',20)
        line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
        line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
        ylim([fromymega toymega])
        xlim([-0.1 0.3])
        set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
    end
end
    legend(hp,legendlabels(1:nnei),'location','nw','fontsize',16)

FigName = ['megacondERP_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)
%% plot megacondERP Exp 1 - all blocks
load('S:\Lab-Shared\Experiments\N1P2\Analysis\grandAverage\megacondERP')
%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
      
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};

nnei=2;%plot until this neighbor
cmap=parula(4);
    
neiERP = nan(length(legendlabels),length(ibls)*8,size(megacondERP,4)); 
count = zeros(length(legendlabels),1);
for con=1:nCond
    for prevcon = 1:nCond
        nei=transM(con,prevcon);
        if nei
            data = squeeze(megacondERP(ibls,con,prevcon,:));
            newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
            data = bsxfun(@minus,data,mean(data(:,newbsl),2));
            neiERP(nei,count(nei)+1:count(nei)+length(ibls),:)=data;
            count(nei)=count(nei)+length(ibls);
        end
    end
end
       
%plot
h=ERPfigure;
set(h,'Position',[100 200 320 300])

clear hp
for nei=1:nnei
    if nei==1
        r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
    end
    hold on
    data = squeeze(neiERP(nei,:,:));
    t=0:1/srate:(length(data)-1)/srate;t=t-0.1;
    %hp(nei)=varplot(t,data','ci',0.95,'color',cmap(nei,:),'linew',3);
    hp(nei)=plot(t,nanmean(data),'Color',cmap(nei,:),'linew',4);
    title('mega ERP')
    set(gca,'fontsize',20)
    line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
    line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
    ylim([fromymega toymega])
    xlim([-0.1 0.3])
    set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
end

    legend(hp,legendlabels(1:nnei),'location','nw','fontsize',16)

FigName = ['megacondERP_Exp' num2str(ExpN) 'all'];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)

%%     cond ERPs
% bls = [2];
% tones = [2];
% condcolors={[0.9 0 0.9],[1 1 1],[0.9 0 0.9],[0.6 0 0.6],[0.1 0 0.1]};
% legendstring=cell(5,length(tones));
% for ir=1:size(legendstring,1)
%     for ic=1:size(legendstring,2)
%         legendstring{ir,ic} = [num2str(ir) ' -> ' num2str(tones(ic))];
%     end
% end
% 
% if 1 %parameters
%     electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP' 
%     load('S:\Lab-Shared\Experiments\MMNchroma\Analysis\Properties')
%     [ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );
%     mode = 'lme';
%     matrixdate = '06-Mar-2016';
%     bslwin = [-0.1 0];
%     order = [1 2 5 3 4];
% end
% 
% for bl = bls
%     ib=find(bl==bls);
%    
%     hf=ERPfigure;
%     set(hf,'Position',[100 100 400 300])
%     load([mixedFolder filesep 'grandMatrix_condComb_' mode '_' blocks{bl} '_' matrixdate])
%     bls=[2];
%     
%     %baseline
%     newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
%     grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
%     t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
% 
% %    r=rectangle('Position',[80, fromy, 40, toy-fromy],'facecolor',[0.7 0.7 0.7],'edgecolor','none');
%     hold on
%     conditions_ordered=conditions(order);
%     for c=1:5
%         if c~=tones(ib)
%             data = grandMatrix(:,:,order(tones(ib)),order(c),electrode);
%             t=0:1:size(data,1)-1;t=t/srate;t=t+winbegin;t=t*1000;
%             plot(t,nanmean(data,2),'DisplayName',legendstring{c},'Color',condcolors{c},'linew',linewidths(c))
%             hold all
%         end
%     end
%     axis([fromx, tox, fromy, toy])
%     set(gca,'fontsize',20)
%     line([0,0],[fromy, toy],'linestyle','-.','col','k','HandleVisibility','off');
%     line([fromx, tox],[0, 0],'linestyle','-.','col','k','HandleVisibility','off');
%     box on
%     hl=legend('Location','nw');
%     set(hl,'fontsize',18)
% end
% 
% FigName = 'Exp1P2ERPs';
% saveas(gcf,[FigFolder filesep FigName],'fig')
% saveas(gcf,[FigFolder filesep FigName],'pdf')
%%      LME model
savedate = '13-May-2019';%peaks calculated at the Amplitudes_MMNchroma script at:
% S:\Lab-Shared\Experiments\MMNchroma\Analysis\Amplitudes_MMNchroma.m
include = false;%include subjects for which time of grandcon peak was determined due to median

load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
formulas = {'N1~dist_mean+(1|subject)','P2~dist_mean+(1|subject)','N1~dist_mean + size_jump+(1|subject)','P2~dist_mean + size_jump+(1|subject)'};

lmes = cell(length(blocks),length(formulas));

if include
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
    mss = [];%known
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end

for ib=bls
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
    for i=1:length(formulas)
         lmes{ib,i} = fitlme(allTables.(blocks{ib}),formulas{i});
        disp(num2str(ib))
        disp(lmes{ib,i}.Formula)
        %disp(num2str(lmes{ib,i}.coefTest))
        disp((lmes{ib,i}.anova))
    end
     disp(num2str(ib))
    disp('Compare for N1')
    compare(lmes{ib,1},lmes{ib,3})
    disp('Compare for P2')
    compare(lmes{ib,2},lmes{ib,4})
end

%% Experiment 2
ExpN=2;
cd('S:\Lab-Shared\Experiments\MMNchromaF\Analysis')
addpath('S:\Lab-Shared\Experiments\MMNchromaF\Analysis')
Definitions_MMNchromaF;
addpath('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
bls=2;%block nums within exp
ibls=3;%relative to all 8 blocks in all 3 exp
%%   P2 Bargraphs
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
bls=2;
nCond = size(allGrandcon_times{1},2);
nPrevCond = size(allPeak_amps{1},3);

% p2 bars
nCond = size(allGrandcon_times{1},2);
for ipeak = 1:2 
    %1:length(whichpeaks)
    %exclude subjects for which grandcon peak was calculated as median
    if include
        allSubjects = 3:33;
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else

        mss=31;
        allSubjects = 3:33;
        includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
    end
    nsub = length(includeSubjects);
    peaks = nan(nsub,length(bls),nCond,nPrevCond);
    peak_means = nan(length(bls),nCond,nPrevCond);
    peak_meanz = nan(length(bls),nCond,nPrevCond);

    peak_errs = nan(length(bls),nCond,nPrevCond);
    h=ERPfigure;
    set(h,'Position',[10 100 200 600])
    bls = [2]; 
        
    for bl= bls  
        %bls
        ib=find(bls==bl);
        for con = 1:size(allPeak_amps{1},2)
            for prevcon = 1:size(allPeak_amps{1},3)
                peaks(:,ib,con,prevcon) = allPeak_amps{bl}(includeSubjects,con,prevcon,ipeak);
                peak_means(ib,con,prevcon) = nanmean(peaks(:,ib,con,prevcon),1);
                CI = Confidence(peaks(:,ib,con,prevcon));
                peak_errs(ib,con,prevcon) = abs(nanmean(peaks(:,ib,con,prevcon))-CI(1));
            end
        end
%         peak_means = peak_means(:,order,order);
%         peak_meanz = peak_meanz(:,order,order);
%         peak_errs = peak_errs(:,order,order);

        peak_means(ib,:,:) = peak_means(ib,order,order);
        peak_meanz(ib,:,:) = peak_meanz(ib,order,order);
        peak_errs(ib,:,:) = peak_errs(ib,order,order);
  
        for i=1:size(allPeak_amps{1},2)
           
            subplot(size(allPeak_amps{1},2),length(bls),(5-i)*length(bls)+ib)
            bar(squeeze(peak_means(ib,i,:)),'FaceColor',colors{i})
            barwitherr(squeeze(peak_errs(ib,i,:)),squeeze(peak_means(ib,i,:)),'FaceColor',colors{i})
            hold all
            if ib==1
               % ylabel(['tone ' num2str(i)])
               
            end
            if ib~=1 || i~=1
                set(gca,'xticklabels',{})
                set(gca,'yticklabels',{})
            end
            set(gca,'fontsize',16)
            if ipeak==1
                axis([0 6 -3.5 0])
            elseif ipeak==2
                axis([0 6 0 4.5])
            end
        end
        %suptitle([whichpeaks{ipeak} '. Electrode: ' electrodeName '. Block: ' num2str(bl)])
    end
    FigName = ['Exp2P2bars'];
saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'fig')
saveas(gcf,[FigFolder filesep FigName '_' whichpeaks{ipeak}],'pdf')
%save for all exp together:
pp = peaks(:);
if ipeak==1
    pp=-1*pp;
    pp(~isnan(pp))=zscore(pp(~isnan(pp)));
    peakz = reshape(pp,size(peaks));
else
    pp(~isnan(pp))=zscore(pp(~isnan(pp)));
    peakz = reshape(pp,size(peaks));
end

ib=0;
for bl= bls
    ib=ib+1;
    %bls
    for con = 1:size(allPeak_amps{1},2)
        for prevcon = 1:size(allPeak_amps{1},3)
            peak_meanz(ib,con,prevcon) = nanmean(peakz(:,ib,con,prevcon),1);
        end
    end
end

NPallexp_peak_means{ipeak,ExpN}=peak_means; 
NPallexp_peak_meanz{ipeak,ExpN}=peak_meanz; 
NPallexp_peaks{ipeak,ExpN}=peaks;
NPallexp_peakz{ipeak,ExpN}=peakz;

end
%% plot peaks of neighbors
% must run P2 Bargraphs first to calculate NPallexp_peak_meanz
pos=[100 200 200 300];
plotPeaksNei(peaksFolder,ExpN,includeSubjects,ibls,pos)

suptitle(['Exp ' num2str(ExpN)])
FigName = ['peaksNei_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)

%% calc megacondERP

mode = 'lme-BP1-20';
matrixdate = '15-Nov-2018';

if 1 %parameters
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP' 
    load('S:\Lab-Shared\Experiments\MMNchroma\Analysis\Properties')
    [ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );
    bslwin = [-0.1 0];
    order = [1 2 5 3 4];
end


for bl = bls
    mb=mb+1;
    
    ib=find(bl==bls);
    load([mixedFolder filesep 'grandMatrix_condComb_' mode '_' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
    conditions_ordered=conditions(order);
    ic=0;
    for c=1:5
        ic=ic+1;
        ipc=0;
        for pc=1:5
            ipc=ipc+1;
            data = grandMatrix(:,includeSubjects,order(c),order(pc),electrode);
            meandata=nanmean(data,2);
            megacondERP(mb,ic,ipc,:)=meandata(1:size(megacondERP,4));
        end
    end
end
%% plot megacondERP Exp 2
ibls = [3];%relative to all 8 blocks across all experiments
%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
      
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};


h=ERPfigure;
set(h,'Position',[100 200 200 300])

nnei=2;%plot until this neighbor
cmap=parula(4);
for isp=1:length(ibls)
    subplot(1,length(ibls),isp)
    %calc only for current ibl
    ibl=ibls(isp);
    neiERP = nan(length(legendlabels),length(ibl)*8,size(megacondERP,4)); 
    count = zeros(length(legendlabels),1);
    for con=1:nCond
        for prevcon = 1:nCond
            nei=transM(con,prevcon);
            if nei
                data = squeeze(megacondERP(ibl,con,prevcon,:));
                newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
                data = bsxfun(@minus,data,mean(data(newbsl,:),1));
                neiERP(nei,count(nei)+1,:)=data;
                count(nei)=count(nei)+1;
            end
        end
    end

    %plot
    clear hp
    for nei=1:nnei
        if nei==1
            r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
            r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        end
        hold on
        data = squeeze(neiERP(nei,:,:));
        t=0:1/srate:(length(data)-1)/srate;t=t-0.1;
        hp(nei)=varplot(t,data','ci',0.9,'color',cmap(nei,:),'linew',3);

        title('mega ERP')
        set(gca,'fontsize',20)
        line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
        line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
        ylim([fromymega toymega])
        set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
    end
    legend(hp,legendlabels(1:nnei),'location','nw','fontsize',16)
end
FigName = ['megacondERP_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)
%% plot megacondERP Exp 2 - all blocks
load('S:\Lab-Shared\Experiments\N1P2\Analysis\grandAverage\megacondERP')
%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
      
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};

nnei=2;%plot until this neighbor
cmap=parula(4);
    
neiERP = nan(length(legendlabels),length(ibls)*8,size(megacondERP,4)); 
count = zeros(length(legendlabels),1);
for con=1:nCond
    for prevcon = 1:nCond
        nei=transM(con,prevcon);
        if nei
            data = squeeze(megacondERP(ibls,con,prevcon,:));
            newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
            data = bsxfun(@minus,data,mean(data(newbsl,:),1));
            neiERP(nei,count(nei)+1:count(nei)+length(ibls),:)=data;
            count(nei)=count(nei)+length(ibls);
        end
    end
end
       
%plot
h=ERPfigure;
set(h,'Position',[100 200 320 300])

clear hp
for nei=1:nnei
    if nei==1
        r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
    end
    hold on
    data = squeeze(neiERP(nei,:,:));
    t=0:1/srate:(length(data)-1)/srate;t=t-0.1;
  %  hp(nei)=varplot(t,data','ci',0.9,'color',cmap(nei,:),'linew',3);
    hp(nei)=plot(t,nanmean(data),'color',cmap(nei,:),'linew',4);

    title('mega ERP')
    set(gca,'fontsize',20)
    line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
    line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
    ylim([fromymega toymega])
    xlim([-0.1 0.3])
    set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
end

   % legend(hp,legendlabels(1:nnei),'location','nw','fontsize',16)

FigName = ['megacondERP_Exp' num2str(ExpN) 'all'];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)

%%     cond ERPs
% bls = [2];
% tones = [4];
% condcolors={[0.1 0 0.1],[0.6 0 0.6],[0.9 0 0.9],[1 1 1],[0.9 0 0.9]};
% legendstring=cell(5,length(tones));
% for ir=1:size(legendstring,1)
%     for ic=1:size(legendstring,2)
%         legendstring{ir,ic} = [num2str(ir) ' -> ' num2str(tones(ic))];
%     end
% end
% 
% mode = 'lme-BP1-20';
% matrixdate = '15-Nov-2018';
% 
% if 1 %parameters
%     electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP' 
%     load('S:\Lab-Shared\Experiments\MMNchroma\Analysis\Properties')
%     [ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );
%     bslwin = [-0.1 0];
%     order = [1 2 5 3 4];
% end
% 
% isp=0;
% for bl = bls
%     ib=find(bl==bls);
%    
%     hf=ERPfigure;
%     set(hf,'Position',[100 100 400 300])
% 
%     load([mixedFolder filesep 'grandMatrix_condComb_' mode '_' blocks{bl} '_' matrixdate])
%     
%     bls=[2];
%     
%     %baseline
%     newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
%     grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
%     t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
% 
% %    r=rectangle('Position',[80, fromy, 40, toy-fromy],'facecolor',[0.7 0.7 0.7],'edgecolor','none');
%     hold on
%     conditions_ordered=conditions(order);
%     for c=1:5
%         if c~=tones(ib)
%             data = grandMatrix(:,:,order(tones(ib)),order(c),electrode);
%             t=0:1:size(data,1)-1;t=t/srate;t=t+winbegin;t=t*1000;
%             plot(t,nanmean(data,2),'DisplayName',legendstring{c},'Color',condcolors{c},'linew',linewidths(c))
%             hold all
%         end
%     end
%     axis([fromx, tox, fromy, toy])
%     set(gca,'fontsize',20)
%     line([0,0],[fromy, toy],'linestyle','-.','col','k','HandleVisibility','off');
%     line([fromx, tox],[0, 0],'linestyle','-.','col','k','HandleVisibility','off');
%     box on
%     hl=legend('Location','nw');
%     set(hl,'fontsize',18)
% 
% end
% 
% FigName = 'Exp2P2ERPs';
% saveas(gcf,[FigFolder filesep FigName],'fig')
% saveas(gcf,[FigFolder filesep FigName],'pdf')
%%      LME model
savedate = '13-May-2019';%peaks calculated at the Amplitudes_MMNchroma script at:
% S:\Lab-Shared\Experiments\MMNchroma\Analysis\Amplitudes_MMNchroma.m
include = false;%include subjects for which time of grandcon peak was determined due to median

load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
formulas = {'N1~dist_mean+(1|subject)','P2~dist_mean+(1|subject)','N1~dist_mean + size_jump+(1|subject)','P2~dist_mean + size_jump+(1|subject)'};

lmes = cell(length(blocks),length(formulas));

if include
    %allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
    mss = 31;%known
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
        disp(num2str(ib))
        disp(lmes{ib,i}.Formula)
        %disp(num2str(lmes{ib,i}.coefTest))
        disp((lmes{ib,i}.anova))
    end
     disp(num2str(ib))
    disp('Compare for N1')
    compare(lmes{ib,1},lmes{ib,3})
    disp('Compare for P2')
    compare(lmes{ib,2},lmes{ib,4})
end

%% Experiment 3
ExpN=3;
cd('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
cd(GHfolder)
ibls = [4:8];%relative to all 8 blocks across all experiments

%%  bargraph P2 
if 1%params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    addtag = '';
    peakdate = '30-Oct-2018';
    whichpeaks = {'N1','P2'};
    include = false;%include subjects for which there was a manual inspection
end
%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
nCond = size(allGrandcon_times{1},2);
nPrevCond = size(allPeak_amps{1},3);

% h=ERPfigure;
%set(h,'Position',[10 100 1000 500])
    
for ipeak = 1:2
    h=ERPfigure;
    set(h,'Position',[10 100 1000 600])
%     
    %1:length(whichpeaks)
    %exclude subjects for which grandcon peak was calculated as median
    if include
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
        mss=[5,14];
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
    end
    nsub = length(includeSubjects);
    peaks = nan(nsub,length(bls),nCond,nPrevCond);
    
    peak_means = nan(length(blocks),nCond);
    peak_meanz = nan(length(blocks),nCond);
    
    peak_errs = nan(length(blocks),nCond);
   
    for bl= bls
        %bls
        for con = 1:size(allPeak_amps{1},2)
            for prevcon = 1:size(allPeak_amps{1},3)
                peaks(:,bl,con,prevcon) = allPeak_amps{bl}(includeSubjects,con,prevcon,ipeak);
                peak_means(bl,con,prevcon) = nanmean(peaks(:,bl,con,prevcon),1);
                CI = Confidence(peaks(:,bl,con,prevcon));
                peak_errs(bl,con,prevcon) = abs(nanmean(peaks(:,bl,con,prevcon))-CI(1));
            end
        end
        
        for i=1:size(allPeak_amps{1},2)
            subplot(size(allPeak_amps{1},2),length(bls),(5-i)*length(bls)+bl)

            bar(squeeze(peak_means(bl,i,:)),'FaceColor',colors{i})
            %ylim([-4 4])
           barwitherr(squeeze(peak_errs(bl,i,:)),squeeze(peak_means(bl,i,:)),'FaceColor',colors{i})
          hold all
            if bl==1
               % ylabel(['tone ' num2str(i)])
               
            end
            if bl~=1 || i~=1
                set(gca,'xticklabels',{})
                set(gca,'yticklabels',{})
            end
            set(gca,'fontsize',16)
            if ipeak==1
                axis([0 6 -3.5 0])
            elseif ipeak==2
                axis([0 6 0 4.5])
            end
            
        end
        
      %  suptitle([whichpeaks{ipeak} '. Electrode: ' electrodeName '. Block: ' num2str(bl)])
%         saveas(gcf,['S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\Geffen figures\prevcond_b' num2str(bl) '_' whichpeaks{ipeak} '.fig'])
%         saveas(gcf,['S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\Geffen figures\prevcond_b' num2str(bl) '_' whichpeaks{ipeak} '.jpg'])

    end

FigName = ['Exp3P2bars_' whichpeaks{ipeak}];
 saveas(gcf,[FigFolder filesep FigName],'fig')
 saveas(gcf,[FigFolder filesep FigName],'pdf')  
%save for all exp together:
pp = peaks(:);
if ipeak==1
    pp=-1*pp;
    pp(~isnan(pp))=zscore(pp(~isnan(pp)));
    peakz = reshape(pp,size(peaks));
else
    pp(~isnan(pp))=zscore(pp(~isnan(pp)));
    peakz = reshape(pp,size(peaks));
end

for bl= bls
    %bls
    for con = 1:size(allPeak_amps{1},2)
        for prevcon = 1:size(allPeak_amps{1},3)
            peak_meanz(bl,con,prevcon) = nanmean(peakz(:,bl,con,prevcon),1);
        end
    end
end
NPallexp_peak_means{ipeak,ExpN}=peak_means;
NPallexp_peak_meanz{ipeak,ExpN}=peak_meanz; 
NPallexp_peaks{ipeak,ExpN}=peaks;
NPallexp_peakz{ipeak,ExpN}=peakz;

end
save([peaksFolder filesep 'peaks_allExp'],'NPallexp_peak_means','NPallexp_peak_meanz','NPallexp_peaks','NPallexp_peakz')
%% plot peaks of neighbors
pos=[100 200 1200 300];
plotPeaksNei(peaksFolder,ExpN,includeSubjects,ibls,pos)
FigName = ['peaksNei_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)

%% calc megacondERP

matrixdate='29-Oct-2018';
tag = 'Bp1-20';

if 1 %parameters
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP' 
    load('S:\Lab-Shared\Experiments\MMNchroma\Analysis\Properties')
    [ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );
    bslwin = [-0.1 0];
    order = [1 2 3 4 5];
end

for bl = bls
    mb=mb+1;
    
    ib=find(bl==bls);
    load([grandFolder filesep 'grandMatrix_condComb_' tag '_Bl' blocks{bl} '_' matrixdate])    
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
    conditions_ordered=conditions(order);
    ic=0;
    for c=1:5
        ic=ic+1;
        ipc=0;
        for pc=1:5
            ipc=ipc+1;
            data = grandMatrix(:,includeSubjects,order(c),order(pc),electrode);
            meandata=nanmean(data,2);
            megacondERP(mb,ic,ipc,:)=meandata(1:size(megacondERP,4));
        end
    end
end
save([grandFolder filesep 'megacondERP'],'megacondERP')
%% plot megacondERP Exp 3
%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
      
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};


h=ERPfigure;
set(h,'Position',[100 200 1200 300])
nnei=2;%plot until this neighbor
cmap=parula(4);
for isp=1:length(ibls)
    subplot(1,length(ibls),isp)
    %calc only for current ibl
    ibl=ibls(isp);
    neiERP = nan(length(legendlabels),length(ibl)*8,size(megacondERP,4)); 
    count = zeros(length(legendlabels),1);
    for con=1:nCond
        for prevcon = 1:nCond
            nei=transM(con,prevcon);
            if nei
                data = squeeze(megacondERP(ibl,con,prevcon,:));
                newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
                data = bsxfun(@minus,data,mean(data(newbsl,:),1));
                neiERP(nei,count(nei)+1,:)=data;
                count(nei)=count(nei)+1;
            end
        end
    end

    %plot
    clear hp
    for nei=1:nnei
        if nei==1
            r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
            r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        end
        hold on
        data = squeeze(neiERP(nei,:,:));
        t=0:1/srate:(length(data)-1)/srate;t=t-0.1;
        hp(nei)=varplot(t,data','ci',0.9,'color',cmap(nei,:),'linew',3);

        title('mega ERP')
        set(gca,'fontsize',20)
        line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
        line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
        ylim([fromymega toymega])
        set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
    end
end
legend(hp,legendlabels(1:nnei),'location','nw','fontsize',16)
FigName = ['megacondERP_Exp' num2str(ExpN)];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)
%% plot megacondERP Exp 3 - all blocks
load('S:\Lab-Shared\Experiments\N1P2\Analysis\grandAverage\megacondERP')
%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
      
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};

nnei=2;%plot until this neighbor
cmap=parula(4);
    
neiERP = nan(length(legendlabels),length(ibls)*8,size(megacondERP,4)); 
count = zeros(length(legendlabels),1);
for con=1:nCond
    for prevcon = 1:nCond
        nei=transM(con,prevcon);
        if nei
            data = squeeze(megacondERP(ibls,con,prevcon,:));
            newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
            data = bsxfun(@minus,data,mean(data(:,newbsl),2));
            neiERP(nei,count(nei)+1:count(nei)+length(ibls),:)=data;
            count(nei)=count(nei)+length(ibls);
        end
    end
end
       
%plot
h=ERPfigure;
set(h,'Position',[100 200 500 300])

clear hp
for nei=1:nnei
    if nei==1
        r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
    end
    hold on
    data = squeeze(neiERP(nei,:,:));
    t=0:1/srate:(length(data)-1)/srate;t=t-0.1;
    hp(nei)=varplot(t,data','ci',0.9,'color',cmap(nei,:),'linew',3);

    title('mega ERP')
    set(gca,'fontsize',20)
    line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
    line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
    ylim([fromymega toymega])
    set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
end

    legend(hp,legendlabels(1:nnei),'location','nw','fontsize',16)

FigName = ['megacondERP_Exp' num2str(ExpN) 'all'];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)
%% plot megacondERP Exp 3 - by spreads
iblss = {4,[5 6],[7 8]};%from all 8 blocks...
load('S:\Lab-Shared\Experiments\N1P2\Analysis\grandAverage\megacondERP')
%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
      
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};

h=ERPfigure;
set(h,'Position',[100 200 1200 300])
nnei=2;%plot until this neighbor
cmap=parula(4);
for isp=1:length(iblss)
    subplot(1,length(iblss),isp)
    %calc only for current ibl
    ibls=iblss{isp};
    neiERP = nan(length(legendlabels),length(ibls)*8,size(megacondERP,4)); 
    count = zeros(length(legendlabels),1);
    for con=1:nCond
        for prevcon = 1:nCond
            nei=transM(con,prevcon);
            if nei
                data = (megacondERP(ibls,con,prevcon,:));
                newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
                data = bsxfun(@minus,data,mean(data(:,:,:,newbsl),4));
                neiERP(nei,count(nei)+1:count(nei)+length(ibls),:)=squeeze(data);
                count(nei)=count(nei)+length(ibls);
            end
        end
    end

    %plot
    clear hp
    for nei=1:nnei
        if nei==1
            r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
            r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        end
        hold on
        data = squeeze(neiERP(nei,:,:));
        t=0:1/srate:(length(data)-1)/srate;t=t-0.1;
        %hp(nei)=varplot(t,data','ci',0.9,'color',cmap(nei,:),'linew',3);
        hp(nei)=plot(t,nanmean(data),'color',cmap(nei,:),'linew',4);
        title('mega ERP')
        set(gca,'fontsize',20)
        line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
        line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
        ylim([fromymega toymega])
        xlim([-0.1 0.3])
        set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
    end
    
end
    legend(hp,legendlabels(1:nnei),'location','nw','fontsize',16)

FigName = ['megacondERP_Exp' num2str(ExpN) 'spreads'];
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)

%%     cond ERPs
% %fromy=[-3, -2.5];toy=[4,3.2];
% % bls = [1,2];
% % tones = [5,5];
% % condcolors={[0 0 0.3],[0 0 0.8],[0 0.1 1],[0.1 0.5 1],[1 1 1]};
% fromy=-1.5;toy = 2.8;
% bls=4;
% tones=3;
% colors = {[0.3 0 0],[1 0 0],[1 1 1],[1 0 0],[0.3 0 0]};
% legendstring=cell(5,length(tones));
% for ir=1:size(legendstring,1)
%     for ic=1:size(legendstring,2)
%         legendstring{ir,ic} = [num2str(ir) ' -> ' num2str(tones(ic))];
%     end
% end
% 
% matrixdate='29-Oct-2018';
% tag = 'Bp1-20';
% 
% if 1 %parameters
%     electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP' 
%     load('S:\Lab-Shared\Experiments\MMNchroma\Analysis\Properties')
%     [ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );
%     bslwin = [-0.1 0];
%     order = [1 2 3 4 5];
% end
% 
% for bl = bls
%     ib=find(bl==bls);
%    
%     hf=ERPfigure;
%     set(hf,'Position',[100 100 400 300])
%     load([grandFolder filesep 'grandMatrix_condComb_' tag '_Bl' blocks{bl} '_' matrixdate])    
%     
%     bls=[1,2];
%     
%     %baseline
%     newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
%     grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
%     t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
%      hold on
%     for c=1:5
%         if c~=tones(ib)
%             data = grandMatrix(:,:,order(tones(ib)),order(c),electrode);
%             t=0:1:size(data,1)-1;t=t/srate;t=t+winbegin;t=t*1000;
%             plot(t,nanmean(data,2),'DisplayName',legendstring{c},'Color',condcolors{c},'linew',linewidths(c))
%             hold all
%         end
%     end
%     axis([fromx, tox, fromy(ib), toy(ib)])
%     set(gca,'fontsize',20)
%     line([0,0],[fromy(ib), toy(ib)],'linestyle','-.','col','k','HandleVisibility','off');
%     line([fromx, tox],[0, 0],'linestyle','-.','col','k','HandleVisibility','off');
%     box on
%     hl=legend('Location','nw');
%     set(hl,'fontsize',18)
% 
% FigName = ['Exp3P2ERPs_bl' num2str(bl)];
% saveas(gcf,[FigFolder filesep FigName],'fig')
% saveas(gcf,[FigFolder filesep FigName],'pdf')
% 
% end
%%      LME model
savedate = '31-Oct-2018';

include = false;%include subjects for which time of grandcon peak was determined due to median

load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
formulas = {'N1~dist_mean+(1|subject)','P2~dist_mean+(1|subject)','N1~dist_mean + size_jump+(1|subject)','P2~dist_mean + size_jump+(1|subject)'};

lmes = cell(length(blocks),length(formulas));

if include
    allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else

    mss=[5 14];%known
    allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end

for ib=bls
    allTables.(['b' blocks{ib}]) = allTables.(['b' blocks{ib}])(ismember(allTables.(['b' blocks{ib}]).subject,includeSubjects),:);
    for i=1:length(formulas)
        lmes{ib,i} = fitlme(allTables.(['b' blocks{ib}]),formulas{i});
        disp(num2str(ib))
        disp(lmes{ib,i}.Formula)
        %disp(num2str(lmes{ib,i}.coefTest))
        disp((lmes{ib,i}.anova))
        
    end
     disp(num2str(ib))
    disp('Compare for N1')
    compare(lmes{ib,1},lmes{ib,3})
    disp('Compare for P2')
    compare(lmes{ib,2},lmes{ib,4})

end

%% All experiments together - prep big table
% One BIG model!!
%moved to new script - bigTable_N1P2 Oct 24
%%      LME models for all experiments!
load('S:\Lab-Shared\Experiments\N1P2\Analysis\MixedModel\theTable')

formulasExp3 = {'Voltage ~ Potential*dist_mean*spread + Potential*size_jump*spread + dist_mean*size_jump + (1|subject)',...
                'Voltage ~ Potential*dist_mean*spread + Potential*size_jump*spread + (1|subject)',...%n.s
                'Voltage ~ Potential*dist_mean*spread + Potential*size_jump*spread - Potential:size_jump:spread + (1|subject)',...%n.s
                };
            
formulas = {'Voltage ~ Potential*dist_mean + Potential*size_jump + dist_mean*size_jump + (1|subject) + (1|ExpN)',...
            'Voltage ~ Potential*dist_mean + Potential*size_jump + (1|subject) + (1|ExpN)',...
            'Voltage ~ Potential*dist_mean + size_jump + (1|subject) + (1|ExpN)',...
            'Voltage ~ Potential*dist_mean + (1|subject) + (1|ExpN)',...
            'Voltage ~ Potential*size_jump + (1|subject) + (1|ExpN)',...
            'Voltage ~ dist_mean:isN1 + dist_mean:isP2 + (1|subject) + (1|ExpN)',...
            'Voltage ~ size_jump:isN1 + size_jump:isP2 + (1|subject) + (1|ExpN)',...
            'Voltage ~ dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + (1|subject) + (1|ExpN)',...
            };

lmes = cell(length(formulas),1);
for i=1:length(formulas)
    disp(i)
    %lmes{i} = fitlme(theTable(theTable.ExpN==categorical(3),:),formulas{i});
    lmes{i} = fitlme(theTable,formulas{i});
end
lmess = cell(length(formulasExp3),1);
for i=1:length(formulasExp3)
    disp(i)
    lmess{i} = fitlme(theTable,formulasExp3{i},'Exclude',theTable.ExpN~=3);
end

%% each potential, all experiments together

formulas = {'Voltage ~ dist_mean + (1|subject)','Voltage ~ dist_mean + size_jump + (1|subject)'};
potentials = {'N1','P2'};
lmes = cell(size(potentials,2),size(formulas,2));
for fi=1:length(formulas)
    for ip=1:length(potentials)
        lmes{ip,fi} = fitlme(theTable(theTable.Potential==categorical(potentials(ip)),:),formulas{fi});
    end
end
%% each potential, each experiment, all conditions
% exp=3;
% formulas = {'Voltage ~ dist_mean + (1|subject)','Voltage ~ dist_mean + size_jump + (1|subject)'};
% potentials = {'N1','P2'};
% lmes = cell(size(potentials,2),size(formulas,2));
% for fi=1:length(formulas)
%     for ip=1:length(potentials)
%         lmes{ip,fi} = fitlme(theTable(theTable.Potential==categorical(potentials(ip)) | theTable.ExpN==categorical((exp)),:),formulas{fi});
%     end
% end

%% P2 figure - all exp together!
%organize per blocks:
load([peaksFolder filesep 'peaks_allExp'])
NPall = cell(2,1);
readfrom = NPallexp_peak_means;

for ipeak=1:2
    NPall{ipeak} = nan(8,5,5);
    ibl=0;
    for iexp=1:3
        ibl=ibl+1;
        for ibe=1:size(readfrom{ipeak,iexp},1)
            NPall{ipeak}(ibl,:,:) = readfrom{ipeak,iexp}(ibe,:,:);
        end
    end
end
meanmean=cell(2,1);
for ipeak=1:2
    meanmean{ipeak} = nan(5,5);
    for ic=1:5
        for ipc=1:5
            meanmean{ipeak}(ic,ipc) = nanmean(NPall{ipeak}(:,ic,ipc));
        end
    end
end   

%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
Nneighbors = cell(2,4);%ipeak, N neighbors
for ipeak=1:2
    for i=1:4
       Nneighbors{ipeak,i}=[]; 
    end
    for ic=1:5
        for ipc=1:5
            if ic==ipc
            else
                Nneighbors{ipeak,transM(ic,ipc)} = [Nneighbors{ipeak,transM(ic,ipc)}, meanmean{ipeak}(ic,ipc)];
            end
        end
    end
end
%add nans to complete 8
for ipeak=1:2
    for in=1:4
        len=length(Nneighbors{ipeak,in});
        if len<8
            for add = (len+1):8
                Nneighbors{ipeak,in} = [Nneighbors{ipeak,in}, nan];
            end
        end
    end
end
%boxplots:
ERPfigure
subplot 211
boxplot([Nneighbors{1,1}',Nneighbors{1,2}',Nneighbors{1,3}'])
ylim([-0.5 0.3])
subplot 212
boxplot([Nneighbors{2,1}',Nneighbors{2,2}',Nneighbors{2,3}'])
ylim([-0.5 0.3])
%P2 bargraph
for ipeak= 1:2
    h=ERPfigure;
    set(h,'Position',[10 100 200 500])
    
    data=meanmean{ipeak};
    %bls
    for ic=1:5
        subplot(5,1,ic)

        bar(data(ic,:),'FaceColor',colors{ic})
        %ylim([-4 4])
  %     barwitherr(squeeze(peak_means(bl,i,:)),squeeze(peak_errs(bl,i,:)),'FaceColor',Colors{i})
        hold all
        ylabel(['tone ' num2str(ic)])
        
    end 
end

%% plot mega cond ERP
%Note: must run all bargraphs plots orderly in order to calc megaERP
%once saved, can be loaded
%save([grandFolder filesep 'megacondERP'],'megacondERP')
load([grandFolder filesep 'megacondERP'])

%transformation matrix:
transM = [0 1 2 3 4; ...
          1 0 1 2 3;...
          2 1 0 1 2;...
          3 2 1 0 1;...
          4 3 2 1 0];
      
legendlabels={'neighbor 1','neighbor 2','neighbor 3','neighbor 4'};
neiERP = nan(length(legendlabels),8*8,size(megacondERP,4)); 

count = zeros(length(legendlabels),1);
for con=1:nCond
    for prevcon = 1:nCond
        nei=transM(con,prevcon);
        if nei
            data = squeeze(megacondERP(:,con,prevcon,:));
            newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
            data = bsxfun(@minus,data,mean(data(:,newbsl),2));
            neiERP(nei,count(nei)+1:count(nei)+size(data,1),:)=data;
            count(nei)=count(nei)+size(data,1);
        end
    end
end

h=ERPfigure;
fromymega=-2;toymega=3;
nnei=2;%plot until this neighbor
cmap=parula(3);

for nei=1:nnei
    if nei==1
        r1=rectangle('Position',[0.08, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
        r2=rectangle('Position',[0.15, fromymega, 0.04, toymega-fromymega],'facecolor',[0.9 0.9 0.9],'edgecolor','none');
    end
    hold on
    data = squeeze(neiERP(nei,:,:));
    t=0:1/srate:(length(data)-1)/srate;t=t-0.1;
    hp(nei)=varplot(t,data','ci',0.9,'color',cmap(nei,:),'linew',3);

    title('mega ERP')
    set(gca,'fontsize',20)
    line([0,0],[fromymega, toymega],'linestyle','-.','col','k','handleVisibility','off');
    line([fromx, tox],[0, 0],'linestyle','-.','col','k','handleVisibility','off');
    ylim([fromymega toymega])
    set(gca,'xtick',[-0.1:0.1:0.3],'xticklabels',[-0.1 0 0.1 0.2 0.3])
end
legend(hp,legendlabels(1:nnei),'location','nw','fontsize',16)

FigName = 'megacondERP';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
plot2svg([FigFolder filesep FigName],gcf)
