function [xlocs, peak_means] = BarPlotN1_MMNchromaF(plotflag, whichpeak)
%BARPLOTN1 Summary of this function goes here
%   as done in L:\Experiments\N1P2\Analysis\N1P2_GH\Amplitudes_N1P2_GH.m
% subsection %%plot bargraph 5 cond N1
%output:

%xlocs (1,25) xlocations of the bars, for adding plots on top
%peak_means (5,5) : blocks x conds 

cd('L:\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_MMNchromaF

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
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
       mss=31;
        includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
    end    
    peak_means = nan(length(blocks),nCond);
    peak_errs = nan(length(blocks),nCond);
    bpeaks = nan(length(includeSubjects),length(blocks)*nCond);
    bxlocs = nan(length(includeSubjects),length(blocks)*nCond);
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
    
    if plotflag
        hf = figure;
        set(hf,'Position',[100 100 600 500])
        hb = barwitherr(peak_errs, peak_means);% Plot with errorbar
    %     for con = 1:5
    %         hb(con).FaceColor = colors{con};
    %     end
        hold on
        %add all participants:
    %     for ii=1:size(bxlocs,1)
    %         scatter(bxlocs_noisy(ii,:),bpeaks(ii,:),10,'filled','MarkerFaceColor',[0.5 0.5 0.5])
    %     end
        set(gca,'xticklabels',{'chroma F','ctrl'})
        ylim([-3.03 0])
        xlim([1.5 2.5])
    %      boxplot(bpeaks,'notch','on','Position',bxlocs(1,:)+1)

        set(gca,'fontsize',14)
        title([whichpeaks{ipeak} ' of all conditions. Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))],'fontsize',16)
        legend({'tone 1','tone 2','tone 3','tone 4','tone 5'},'Location','northeastoutside')
        ylabel(['mean amplitude ± CI 5-95%, \muV'])
        xlabel('Block type')
    end
end
xlocs = bxlocs(1,:)+ones(size(bxlocs(1,:)));
end

