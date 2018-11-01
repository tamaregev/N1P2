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

%% peak detection - old
if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    bls = [1:5];
    whichpeak = 'N1';
    disp(whichpeak)
    switch whichpeak
        case 'N1'
            pwin = [50,150];%100 +- 50%ms    
            positive_peak = false;
            addtag = '';
        case 'P2'
            pwin = [130,250];%190 + - 60%ms
            positive_peak = true;
            addtag = '';
    end
    select_largest = true;
    plotflag = false;%plot automatically selected peaks.
    mode = 'Bp1-20';
    matrixdate = '17-Oct-2018';
    srate = 512;
    bslwin = [-0.1,0];
end
allPeak_smpls = cell(1,length(blocks));
allPeak_amps = cell(1,length(blocks));
for bl=bls
    disp(blocks{bl})
    load([grandFolder filesep 'grandMatrix_condComb_' mode '_bl' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
    figure
    [peak_smpls,  peak_amps, manuals ] = PeakDetection(grandMatrix, t, electrodeName, pwin, positive_peak, select_largest, plotflag );
    allPeak_smpls{bl} = peak_smpls;
    allPeak_amps{bl} = peak_amps; 
end
save([mixedFolder whichpeak '_' electrodeName '_' addtag date],'allPeak_smpls','allPeak_amps','t','blocks','bls','matrixdate','bslwin','pwin','electrodeName','manuals')
% saved:
% peakdates = {'17-Oct-2018','17-Oct-2018'};
% whichpeaks = {'N1','P2'};
%% peak detection - new
% plot all conditions together for each subject (per block)
% arrange all conditions randomly

if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    bls = [1:5];
    whichpeak = 'N1';
    disp(whichpeak)
    switch whichpeak
        case 'N1'
            pwin = [50,150];%100 +- 50%ms    
            positive_peak = false;
            addtag = '';
        case 'P2'
            pwin = [130,250];%190 + - 60%ms
            positive_peak = true;
            addtag = '';
    end
    select_largest = true;
    plotflag = true;%plot automatically selected peaks.
    mode = 'Bp1-20';
    matrixdate = '17-Oct-2018';
    srate = 512;
    bslwin = [-0.1,0];
end
allPeak_smpls = cell(1,length(blocks));
allPeak_amps = cell(1,length(blocks));
allPeak_times = cell(1,length(blocks));
allManuals = cell(1,length(blocks));
for bl= bls
    disp(blocks{bl}) 
    load([grandFolder filesep 'grandMatrix_condComb_' mode '_bl' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
    %peak detection
    [peak_smpls,  peak_amps, peak_times, manuals ] = PeakDetection_N1P2(grandMatrix, whichSubjects, t, electrodeName, pwin, positive_peak, select_largest, plotflag, blocks{bl});
    allPeak_smpls{bl} = peak_smpls;
    allPeak_amps{bl} = peak_amps; 
    allPeak_times{bl} = peak_times;
    allManuals{bl} = manuals;
end
disp('saving...')
save([mixedFolder whichpeak '_' electrodeName '_' addtag date],'allPeak_smpls','allPeak_amps','allPeak_times','t','blocks','bls','matrixdate','bslwin','pwin','electrodeName','allManuals')
disp('done')
%% peak detection - new 2
% as before but add the ERP of all prevcon together and all 4 prevcons arranged under, randomly within the column 
if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    bls = [1:5];
    whichpeak = 'P2';
    disp(whichpeak)
    switch whichpeak
        case 'N1'
            pwin = [50,150];%100 +- 50%ms    
            positive_peak = false;
            addtag = '';
        case 'P2'
            pwin = [130,250];%190 + - 60%ms
            positive_peak = true;
            addtag = '';
    end
    select_largest = true;
    plotflag = true;%plot automatically selected peaks.
    mode = 'Bp1-20';
    matrixdate = '17-Oct-2018';
    srate = 512;
    bslwin = [-0.1,0];
end
allPeak_smpls = cell(1,length(blocks));
allPeak_amps = cell(1,length(blocks));
allPeak_times = cell(1,length(blocks));
allManuals = cell(1,length(blocks));

for bl= bls
    disp(blocks{bl}) 
    load([grandFolder filesep 'grandMatrix_condComb_' mode '_bl' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));    
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
    %peak detection
    [peak_smpls,  peak_amps, peak_times, manuals ] = PeakDetection_N1P2_2(grandMatrix, whichSubjects, t, electrodeName, pwin, positive_peak, select_largest, plotflag, blocks{bl});
    allPeak_smpls{bl} = peak_smpls;
    allPeak_amps{bl} = peak_amps;
    allPeak_times{bl} = peak_times;
    allManuals{bl} = manuals;
end
disp('saving...')
save([mixedFolder whichpeak '_' electrodeName '_' addtag date],'allPeak_smpls','allPeak_amps','allPeak_times','t','blocks','bls','matrixdate','bslwin','pwin','electrodeName','allManuals')
disp('done')
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
    plotflag = false;%plot automatically selected peaks.
    pauseflag = false;
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
    [peak_smpls,  peak_amps, peak_times, grandcon_as_median,  grandcon_peak_times, grandcon_peak_amps ] = PeakDetection_N1P2_3(grandMatrix, whichSubjects, t, electrodeName, dt, pwins, positive_peaks, select_largest, plotflag, pauseflag, blocks{bl});
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

%% plot bargraph 5 cond N1
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
whichpeaks = {'N1','P2'};
%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
nCond = size(allGrandcon_times{1},2);
for ipeak = 1:length(whichpeaks)
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(33,length(blocks)*nCond);
    bxlocs = nan(33,length(blocks)*nCond);
    for bl=1:length(blocks)
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(:,:,ipeak);
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
    bxlocs = bxlocs + ones(size(bxlocs)).*(1+(rand(size(bpeaks))-0.5)/100) ;
    
    hf = figure;
    set(hf,'Position',[100 100 1000 700])
    h = barwitherr(peak_errs, peak_means);% Plot with errorbar
    hold on
    for ii=1:size(bxlocs,1)
        scatter(bxlocs(ii,:),bpeaks(ii,:),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
    end
    set(gca,'xticklabels',blocks)
    set(gca,'fontsize',14)
    title([whichpeaks{ipeak} ' of all conditions'],'fontsize',16)
    legend({'tone 1','tone 2','tone 3','tone 4','tone 5'},'Location','northeastoutside')
    ylabel(['mean amplitude ± CI 5-95%, \muV'])
end
%% boxplot
pos = zeros(1,25);
for i=1:5
    for j=1:5
        pos((i-1)*5+j) = (i-1)*5+j + (i-1)*2; 
    end
end 
x=[1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 ];

for ipeak = 1:length(whichpeaks)
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(33,length(blocks)*nCond);
    bxlocs = repmat(pos,[33,1]);

    for bl=1:length(blocks)
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(:,:,ipeak);
        cbpeaks(any(isnan(cbpeaks), 2), :) = [];
        %into all block peaks
        bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=cbpeaks;
            
    end
    bxlocs = bxlocs.*(1+(rand(size(bpeaks))-0.5)/100) ;
    
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
pos = zeros(1,25);
for i=1:5
    for j=1:5
        pos((i-1)*5+j) = (i-1)*5+j + (i-1)*2; 
    end
end 

for ipeak = 1:length(whichpeaks)
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(33,length(blocks)*nCond);
    bxlocs = nan(33,length(blocks)*nCond);

    hf=ERPfigure;
    set(hf,'Position',[100,100,1200,400])
    for bl=1:length(blocks)
        subplot(1,5,bl)
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(:,:,ipeak);
        cbpeaks(any(isnan(cbpeaks), 2), :) = [];
        %into all block peaks
        bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=cbpeaks;
        bxlocs(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=repmat(1:5,[33,1]).*(1+(rand(33,5)-0.5)/10);
        violin(bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond)))
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

for ipeak = 1:length(whichpeaks)
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(33,length(blocks)*nCond);
    bxlocs = nan(33,length(blocks)*nCond);

    hf=ERPfigure;
    set(hf,'Position',[100,100,1200,400])
    for bl=1:length(blocks)
        subplot(1,5,bl)
        %current block peaks
        cbpeaks = allGrandcon_amps{bl}(:,:,ipeak);
        cbpeaks(any(isnan(cbpeaks), 2), :) = [];
        %into all block peaks
        bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=cbpeaks;
        bxlocs(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond))=repmat(1:5,[33,1]).*(1+(rand(33,5)-0.5)/10);
        violinplot(bpeaks(:,((bl-1)*nCond+1):((bl-1)*nCond+nCond)),{'1','2','3','4','5'})
        hold on
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
    suptitle([whichpeaks{ipeak} ' of all conditions'])
    
end

%% plot bargraph P2 prevcond
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
addtag = '';
peakdate = '30-Oct-2018';
whichpeaks = {'N1','P2'};
%load peaks
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
nCond = size(allGrandcon_times{1},2);
for ipeak = 2 
    %1:length(whichpeaks)
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    for bl= 5
        %bls
        for con = 1:size(allPeak_amps{1},2)
            for prevcon = 1:size(allPeak_amps{1},3)
                peaks = allPeak_amps{bl}(:,con,prevcon,ipeak);
                peak_means(bl,con,prevcon) = nanmean(peaks,1);
                CI = Confidence(peaks);
                peak_errs(bl,con,prevcon) = abs(nanmean(peaks)-CI(1));
            end
        end
        h=ERPfigure;
        set(h,'Position',[10 100 400 500])
        Colors = {[0 0 0.5],[0.2 0.5 1],[1 0 0],[0.2 0.5 1],[0 0 0.5]};
        for i=1:size(allPeak_amps{1},2)
            subplot(size(allPeak_amps{1},2),1,i)
            bar(squeeze(peak_means(bl,i,:)),'FaceColor',Colors{i})
    %     barwitherr(squeeze(peak_means(bl,i,:)),squeeze(peak_errs(bl,i,:)),'FaceColor',Colors{i})
            hold all
        end
        suptitle([whichpeaks{ipeak} '. Electrode: ' electrodeName '. Block: ' num2str(bl)])
    
    end
end

%% prepare table ... and save

peakdate = '30-Oct-2018';
whichpeaks = {'N1','P2'};
electrodeName = 'Cz';
blockstrings = {'b1','b2a','b2b','b3a','b3b'};
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
load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
blocks = {'b1','b2a','b2b','b3a','b3b'};
savedate = '17-Oct-2018';
load([mixedFolder 'LMEtables_combined_' savedate ])
formulas = {'N1~dist_mean+size_jump+currMIDI+(1|subject)','N1~dist_mean+currMIDI+(1|subject)','N1~dist_mean+(1|subject)',...
            'P2~dist_mean+size_jump+currMIDI+(1|subject)','P2~dist_mean+currMIDI+(1|subject)','P2~size_jump+(1|subject)'};
lmes = cell(length(blocks),length(formulas));
for ib=[1,2,3,4,5]
    for i=1:length(formulas)
        lmes{ib,i} = fitlme(allTables.(blocks{ib}),formulas{i});
    end
end

