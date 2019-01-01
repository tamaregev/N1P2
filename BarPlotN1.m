function [xlocs, peak_means] = BarPlotN1(plotflag, whichpeak)
%BARPLOTN1 Summary of this function goes here
%   as done in L:\Experiments\N1P2\Analysis\N1P2_GH\Amplitudes_N1P2_GH.m
% subsection %%plot bargraph 5 cond N1
%output:

%xlocs (1,25) xlocations of the bars, for adding plots on top
%peak_means (5,5) : blocks x conds 

cd('L:\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
mixedFolder = [AnalysisFolder 'MixedModel\'];

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

switch whichpeak
    case 'N1'
        iwhichpeak = 1;
    case 'P2'
        iwhichpeak = 2;
end
    
for ipeak = iwhichpeak
    %:length(whichpeaks)
    %exclude subjects for which grandcon peak was calculated as median
    if include
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
        mss = [];
        for bl=1:length(allGrandcon_as_median)
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

    if plotflag
        hf = figure;
        set(hf,'Position',[100 100 1000 700])
        h = barwitherr(peak_errs, peak_means);% Plot with errorbar
        hold on
        %add all participants:
    %     for ii=1:size(bxlocs,1)
    %         scatter(bxlocs_noisy(ii,:),bpeaks(ii,:),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
    %     end
        set(gca,'xticklabels',blocks)

    %     boxplot(bpeaks,'Position',bxlocs(1,:)+1)

        set(gca,'fontsize',14)
        title([whichpeaks{ipeak} ' of all conditions. Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))],'fontsize',16)
        legend({'tone 1','tone 2','tone 3','tone 4','tone 5'},'Location','northeastoutside')
        ylabel(['mean amplitude ± CI 5-95%, \muV'])
        xlabel('Block type')
    end
end
xlocs = bxlocs(1,:)+ones(size(bxlocs(1,:)));
end

