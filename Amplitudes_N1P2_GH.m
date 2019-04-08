%Amplitudes_N1P2
% March 15, Tamar - adapting from MixedModelF.m

%% definitions
Definitions_N1P2

grandFolder = [AnalysisFolder 'grandAverage'];
mixedFolder = [AnalysisFolder 'MixedModel\'];

srate = 512;
Expinfo.srate = srate;

%codes and names table:
Bnames = {'1','2a','2b','3a','3b'}';
Bcodes = [10,20,30,40,50]';
Bcodes_names = table(Bnames,Bcodes);
blocks = Bnames;

%% peak detection - new 3
% find peaks of the grand cond ERP, average voltage of prevcons accordingly
matrixdate = '29-Oct-2018';
if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    bls = [1:5];
    whichpeaks = {'N1','P2'};
    disp(whichpeaks)
    pwins = {[50,150],[130,250]};    
    positive_peaks = [false true];
    addtag = '';
    dt = 5;%ms before and after peak - window for calculating N1-P2 amp
    select_largest = true;
    plotflag = true;%plot automatically selected peaks.
    pauseflag = true;
    mode = 'Bp1-20';
    srate = 512;
    bslwin = [-0.1,0];
end
allPeak_smpls = cell(1,length(blocks));
allPeak_amps = cell(1,length(blocks));
allPeak_times = cell(1,length(blocks));
allGrandcon_as_median = cell(1,length(blocks));
allGrandcon_times = cell(1,length(blocks));
allGrandcon_amps = cell(1,length(blocks));

for bl= bls
    disp(blocks{bl}) 
    load([grandFolder filesep 'grandMatrix_condComb_' mode '_bl' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
    %peak detection
    [peak_smpls,  peak_amps, peak_times, grandcon_as_median,  grandcon_peak_times, grandcon_peak_amps ] = PeakDetection_N1P2_GH(grandMatrix, whichSubjects, t, electrodeName, dt, pwins, positive_peaks, select_largest, plotflag, pauseflag, blocks{bl});
    allPeak_smpls{:,bl} = peak_smpls;
    allPeak_amps{:,bl} = peak_amps;
    allPeak_times{:,bl} = peak_times;
    allGrandcon_as_median{:,bl} = grandcon_as_median;
    allGrandcon_times{:,bl} = grandcon_peak_times;
    allGrandcon_amps{:,bl} = grandcon_peak_amps;
end
disp(['saving ' mixedFolder 'N1P2_' electrodeName '_' addtag date ' ...'] )
save([mixedFolder 'N1P2_' electrodeName '_' addtag date],'allPeak_smpls','allPeak_amps','allPeak_times','allGrandcon_as_median','allGrandcon_amps','allGrandcon_times','t','dt','mode','blocks','bls','matrixdate','bslwin','pwins','electrodeName')
disp('done')

%% check that N1 is before P2
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

for bl=1:length(allGrandcon_times)
    tN1 = allGrandcon_times{bl}(:,:,1);
    tP2 = allGrandcon_times{bl}(:,:,2);
    [s, c] = find((tP2 - tN1)<=0);
    disp(['Bl = ' num2str(bl)])
    if ~isempty(s)
        disp(['Problem with s ' num2str(s) ', cond ' num2str(c)])
    end
end

%% compare grandcond peaks and mean prevcons
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
%TODO

ERPfigure
for ipeak = 1:2
    for bl=1:length(allPeak_amps)
        grandcon = allGrandcon_amps{bl}(:,:,ipeak);
        mprevcon = nanmean(allPeak_amps{bl}(:,:,:,ipeak),3);
        subplot(2,5,(ipeak-1)*length(allPeak_amps)+bl)
        hold on
        for is=1:size(grandcon,1)
            for ic=1:size(grandcon,2)
                if ~isnan(grandcon(is,ic))
                    plot(grandcon(is,ic),mprevcon(is,ic),'.','Markersize',10)
    %                 text(grandcon(is,ic),mprevcon(is,ic),[num2str(is) ', ' num2str(ic)])
                end
            end
        end
        xlabel('grandcon peaks');ylabel('mean of prevcons')
        hold on
        xlim([min(min(grandcon)) max(max(grandcon))]);
        ylim([min(min(mprevcon)) max(max(mprevcon))]);
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
        extrem(1) = min([xlims(1) ylims(1)]);
        extrem(2) = max([xlims(2) ylims(2)]);
        plot(extrem,extrem,'k-.')
        title(whichpeaks{ipeak})
    end
end
%% plot bargraph 5 cond N1
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

for ipeak = 1:2
    %:length(whichpeaks)
    %exclude subjects for which grandcon peak was calculated as median
    if include
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
        mss = [];
        for bl=1:length(allGrandcon_as_median)
            allGrandcon_as_median{bl}(:,:,ipeak);
            [r, c]=find(allGrandcon_as_median{bl}(:,:,ipeak));
            mss = [mss, r'];
        end
        mss = unique(mss);
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
    set(hf,'Position',[100 100 1000 500])
    h = barwitherr(peak_errs, peak_means);% Plot with errorbar
    hold on
    %add all participants:
%     for ii=1:size(bxlocs,1)
%         scatter(bxlocs_noisy(ii,:),bpeaks(ii,:),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
%     end

   % set(gca,'xticklabels',blocks)
   set(gca,'xticklabels',{'1','2a','2b','3a','3b'})

%     boxplot(bpeaks,'Position',bxlocs(1,:)+1)
    
    set(gca,'fontsize',14)
    title([whichpeaks{ipeak} ' of all conditions. Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))],'fontsize',16)
    legend({'tone 1','tone 2','tone 3','tone 4','tone 5'},'Location','northeastoutside')
    ylabel(['mean amplitude ± CI 5-95%, \muV'])
end

%% boxplot
%depends on previous
pos = zeros(1,25);
for i=1:5
    for j=1:5
        pos((i-1)*5+j) = (i-1)*5+j + (i-1)*2; 
    end
end 
x=[1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 ];

for ipeak = 1
    %:length(whichpeaks)
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(length(includeSubjects),length(blocks)*nCond);
    bxlocs = repmat(pos,[length(includeSubjects),1]);

    for bl=1:length(blocks)
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(includeSubjects,:,ipeak);
        cbpeaks(any(isnan(cbpeaks), 2), :) = [];
        %into all block peaks
        bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=cbpeaks;
    end
    bxlocs = bxlocs.*(1+(rand(size(bpeaks))-0.5)/100);
    
    figure
    boxplot(bpeaks,'Position',pos)
    set(gca,'xticklabels',{'1','2','3','4','5'})
    set(gca,'fontsize',14)
    title([whichpeaks{ipeak} ' of all conditions'],'fontsize',16)    
    ylabel(['mean amplitude ± CI 5-95%, \muV'])
    hold on
    for ii=1:size(bxlocs,1)
        scatter(bxlocs(ii,:),bpeaks(ii,:),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
    end
end
%% violin
%depends on previous
pos = zeros(1,25);
for i=1:5
    for j=1:5
        pos((i-1)*5+j) = (i-1)*5+j + (i-1)*2; 
    end
end 

for ipeak = 1
    %:length(whichpeaks)
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(length(includeSubjects),length(blocks)*nCond);
    bxlocs = nan(length(includeSubjects),length(blocks)*nCond);

    hf=ERPfigure;
    set(hf,'Position',[100,100,1200,400])
    for bl=1:length(blocks)
        subplot(1,5,bl)
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(includeSubjects,:,ipeak);
        cbpeaks(any(isnan(cbpeaks), 2), :) = [];
        %into all block peaks
        bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=cbpeaks;
        bxlocs(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=repmat(1:5,[length(includeSubjects),1]).*(1+(rand(length(includeSubjects),5)-0.5)/10);
        violin(bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond)))
        ylim([-8 4])
        for ii=1:size(bxlocs,1)
            scatter(bxlocs(ii,((bl-1)*nCond+1):((bl-1)*nCond+nCond)),bpeaks(ii,((bl-1)*nCond+1):((bl-1)*nCond+nCond)),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
        end
        if bl==5
        else
            legend off
        end
        title(['Condition ' num2str(bl)])
    end
    %set(gca,'xticklabels',{'1','2','3','4','5'})
    %set(gca,'fontsize',14)
    suptitle([whichpeaks{ipeak} ' of all conditions'])
   
    %ylabel(['mean amplitude ± CI 5-95%, \muV'])
    
end
%% Violin
%depends on previous
for ipeak = 1
    %:length(whichpeaks)
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(length(includeSubjects),length(blocks)*nCond);
    bxlocs = nan(length(includeSubjects),length(blocks)*nCond);

    hf=ERPfigure;
    set(hf,'Position',[100,100,1200,400])
    for bl=1:length(blocks)
        subplot(1,5,bl)
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(includeSubjects,:,ipeak);
        cbpeaks(any(isnan(cbpeaks), 2), :) = [];
        %into all block peaks
        bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=cbpeaks;
        bxlocs(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=repmat(1:5,[length(includeSubjects),1]).*(1+(rand(length(includeSubjects),5)-0.5)/10);
        violinplot(bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond)),{'1','2','3','4','5'})
        hold on
        ylim([-5.5 1.5])
        plot(1:5,mean(cbpeaks),'+','MarkerSize',16,'MarkerEdgeColor','k')
%         for ii=1:size(bxlocs,1)
%             scatter(bxlocs(ii,((bl-1)*nCond+1):((bl-1)*nCond+nCond)),bpeaks(ii,((bl-1)*nCond+1):((bl-1)*nCond+nCond)),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
%         end
        if bl==5
        else
            legend off
        end
        title(['Condition ' num2str(bl)])
    end

    suptitle([whichpeaks{ipeak} ' of all conditions. Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))])

end

%% plot bargraph P2 prevcond
if 1%params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    addtag = '';
    peakdate = '30-Oct-2018';
    whichpeaks = {'N1','P2'};
    include = true;%include subjects for which there was a manual inspection
end
%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
nCond = size(allGrandcon_times{1},2);
for ipeak = 1:2
    %1:length(whichpeaks)
    %exclude subjects for which grandcon peak was calculated as median
    if include
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
        mss = [];
        for bl=1:length(allGrandcon_as_median)
            allGrandcon_as_median{bl}(:,:,ipeak);
            [r, c]=find(allGrandcon_as_median{bl}(:,:,ipeak));
            mss = [mss, r'];
        end
        mss = unique(mss);
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
    end
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    for bl= 1:5
        %bls
        for con = 1:size(allPeak_amps{1},2)
            for prevcon = 1:size(allPeak_amps{1},3)
                peaks = allPeak_amps{bl}(includeSubjects,con,prevcon,ipeak);
                peak_means(bl,con,prevcon) = nanmean(peaks,1);
                CI = Confidence(peaks);
                peak_errs(bl,con,prevcon) = abs(nanmean(peaks)-CI(1));
            end
        end
        h=ERPfigure;
        set(h,'Position',[10 100 400 500])
        %Colors = {[0 0 0.5],[0.2 0.5 1],[1 0 0],[0.2 0.5 1],[0 0 0.5]};
        Colors = parula(5);
        for i=1:size(allPeak_amps{1},2)
            subplot(size(allPeak_amps{1},2),1,i)
            bar(squeeze(peak_means(bl,i,:)),'FaceColor',Colors(i,:))
    %     barwitherr(squeeze(peak_means(bl,i,:)),squeeze(peak_errs(bl,i,:)),'FaceColor',Colors{i})
            hold all
            
        end
        xlabel('previous note')
        suptitle([whichpeaks{ipeak} '. Electrode: ' electrodeName '. Block: ' num2str(bl)])
%         saveas(gcf,['L:\Experiments\N1P2\Analysis\Figures\Geffen figures\prevcond_b' num2str(bl) '_' whichpeaks{ipeak} '.fig'])
%         saveas(gcf,['L:\Experiments\N1P2\Analysis\Figures\Geffen figures\prevcond_b' num2str(bl) '_' whichpeaks{ipeak} '.jpg'])

    end
end

%% prepare table ... and save

peakdate = '30-Oct-2018';
whichpeaks = {'N1','P2'};
electrodeName = 'Cz';
blockstrings = {'b1','b2a','b2b','b3a','b3b'};
include = false;

%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' peakdate])
blocks = blockstrings;
conditions = {'tone 1','tone 2','tone 3','tone 4','tone 5'};
prevConditions = {'tone 1','tone 2','tone 3','tone 4','tone 5'};
mode = 'lme';

for ii=1:length(whichpeaks)
    whichpeak = whichpeaks{ii};
    allTables = struct;
    
    for bl = 1:length(blocks)
        nSubjs = sum(~isnan(allPeak_amps{bl}(:,1,5,1)));
        numLines = nSubjs*numel(conditions)*(numel(conditions)-1);
        subject = nan(numLines,1);
        y = nan(numLines,1);
        currNote = cell(numLines,1);
        prevNote = cell(numLines,1);
        line = 0;
        for s = 1:size(allPeak_amps{1},1)
            for con = 1:size(allPeak_amps{1},2)
                for prevcon = 1:size(allPeak_amps{1},3)
                    if ~isnan(allPeak_amps{bl}(s,con,prevcon,ii))
                        line = line+1;
                        disp(['bl=' num2str(bl) ', subj=' num2str(s) ', con=' num2str(con) ', prevcon=' num2str(prevcon) ' line=' num2str(line)])
                        subject(line) = s;
                        currNote{line} = conditions{con};
                        prevNote{line} = prevConditions{prevcon};
                        y(line) = allPeak_amps{bl}(s,con,prevcon,ii);
                    end
                end
            end
        end
        
        LMEtable = table(subject,currNote,prevNote,y);
        
        allTables.(blockstrings{bl}) = LMEtable;
        clear LMEtable
    end
    
    % add MIDI numbers to notes:
    MIDI = struct;
    s=2;
    load([EDATfolder ExpName '_' Subjects{s} '_' sessions{s} '_' 'expdata.mat'])
    for b=1:length(blocks)
        for im=1:length(expdata.Passive.blocks(b).MIDIs)
            MIDI.(blocks{b})(im) = expdata.Passive.blocks(b).MIDIs{im};
        end
    end
    
    for bl = 1:length(blocks)
        nSubjs = sum(~isnan(allPeak_amps{bl}(:,1,5,1)));
        numLines = nSubjs*numel(conditions)*(numel(conditions)-1);
        currMIDI = nan(numLines,1);
        prevMIDI = nan(numLines,1);
        dist_mean = nan(numLines,1);
        size_jump = nan(numLines,1);
        line = 0;
        for line = 1:size(allTables.(blocks{bl}),1)
            currStr = allTables.(blocks{bl}).currNote(line);
            currMIDI(line) = MIDI.(blocks{bl})(strcmp(currStr,conditions));
            dist_mean(line) = abs(currMIDI(line) - mean(MIDI.(blocks{bl})));
            prevStr = allTables.(blocks{bl}).prevNote(line);
            prevMIDI(line) = MIDI.(blocks{bl})(strcmp(prevStr,prevConditions));
            size_jump(line) = abs(currMIDI(line) - prevMIDI(line));
        end
        MIDIs = table(currMIDI,prevMIDI,dist_mean,size_jump);
        allTables.(blocks{bl}) = [allTables.(blocks{bl})(:,1:3) , MIDIs, allTables.(blocks{bl})(:,4)];
    end
    % save tables:
    save([mixedFolder 'LMEtables_' whichpeak '_' date ], 'allTables')
end

%% load tables and merge:
savedate = '31-Oct-2018';
whichpeaks = {'N1','P2'};
mode = 'lme';
blocks = {'b1','b2a','b2b','b3a','b3b'};

whichpeak = whichpeaks{1};
load([mixedFolder 'LMEtables_' whichpeak '_' savedate]);%load table: allTables
allTables_N1 = allTables;
whichpeak = whichpeaks{2};
load([mixedFolder 'LMEtables_' whichpeak '_' savedate]);%load table: allTables
allTables_P2 = allTables;
clear allTables
for bl = 1:length(blocks)
    allTables.(blocks{bl}) = allTables_N1.(blocks{bl})(:,1:end-1);
    allTables.(blocks{bl}).N1 = allTables_N1.(blocks{bl}).y;
    allTables.(blocks{bl}).P2 = allTables_P2.(blocks{bl}).y;
    allTables.(blocks{bl}).p2p = allTables_P2.(blocks{bl}).y-allTables_N1.(blocks{bl}).y;
end
save([mixedFolder 'LMEtables_combined_' date ], 'allTables')
%% LME:
savedate = '31-Oct-2018';
include = false;%include subjects for which time of grandcon peak was determined due to median

load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
blocks = {'b1','b2a','b2b','b3a','b3b'};
formulas = {'N1~dist_mean+size_jump+currMIDI+(1|subject)','N1~dist_mean+currMIDI+(1|subject)','N1~dist_mean+(1|subject)',...
            'P2~dist_mean+size_jump+currMIDI+(1|subject)','P2~dist_mean+currMIDI+(1|subject)','P2~size_jump+(1|subject)'};

lmes = cell(length(blocks),length(formulas));

if include
    allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
    mss = [];
    for bl=1:length(allGrandcon_as_median)
        for ipeak = 1:2%here I exclude from all analysis any subject that either had a problem in the N1 or in the P2s
            [r, c]=find(allGrandcon_as_median{bl}(:,:,ipeak));
            mss = [mss, r'];
        end
    end
    mss = unique(mss);
    allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end

for ib=[1,2,3,4,5]
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
    for i=1:length(formulas)
        lmes{ib,i} = fitlme(allTables.(blocks{ib}),formulas{i});
    end
end
