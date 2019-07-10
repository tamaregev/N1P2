function [xlocs, peak_means] = BarPlotN1_MMNchroma(whichpeak)
%BARPLOTN1 Summary of this function goes here
%   as done in L:\Experiments\N1P2\Analysis\N1P2_GH\Amplitudes_N1P2_GH.m
% subsection %%plot bargraph 5 cond N1
%output:

%xlocs (1,25) xlocations of the bars, for adding plots on top
%peak_means (5,5) : blocks x conds 

cd('L:\Experiments\MMNchroma\Analysis')
Definitions_MMNchroma
order=cs;
electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP' 
addtag = '';
peakdate = '18-Nov-2018';
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])

bls=[2,4];
nCond = size(allGrandcon_times{1},2);

switch whichpeak
    case 'N1'
        iwhichpeak = 1;
    case 'P2'
        iwhichpeak = 2;
end

for ipeak = iwhichpeak
    peak_means = nan(length(blocks), nCond);
    peak_errs = nan(length(blocks), nCond);
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    bxlocs = nan(length(includeSubjects),length(blocks)*nCond);
    ib=0;
    for bl=bls
        ib=ib+1;
         cbpeaks = allGrandcon_amps{bl}(includeSubjects,:,ipeak);
        cbpeaks = cbpeaks(:,order);
        %the smaller number of participants in classicctrl
        
        %into all block peaks
        bpeaks(:,((ib-1)*nCond+1):((ib-1)*nCond+nCond))=cbpeaks;
       
        bxlocs(:,((ib-1)*nCond+1):((ib-1)*nCond+nCond))=repmat(ib-1,[size(cbpeaks,1),nCond]);
       
        for con = 1:nCond
            bxlocs(:,((ib-1)*nCond+1)+con-1)=bxlocs(:,((ib-1)*nCond+1)+con-1)+0.15*(con-3);
            peaks = cbpeaks(:,con);
            peak_means(ib,con) = nanmean(peaks);
            CI = Confidence(peaks);
            peak_errs(ib,con) = abs(nanmean(peaks)-CI(1));
       end
    end
    hf=ERPfigure;

    %h = bar(peak_means(bls,:));% Plot with errorbars
    hb = barwitherr(peak_errs, peak_means);% Plot with errorbars


    set(gca,'xticklabels',blocks(bls))
    set(gca,'fontsize',14)

    title(['N1 peaks - ' electrodeName],'fontsize',16)
            legend({'tone 1','tone 2','tone 3','tone 4','tone 5'},'Location','northeastoutside')

    ylim([-3.4, 0]);xlim([0.5 2.5])
end
xlocs = bxlocs(1,:)+ones(size(bxlocs(1,:)));