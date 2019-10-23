%functionalAll
FigFolder = 'L:\Experiments\N1P2\Analysis\N1P2_GH\PaperFigures\functionalPlots';
%% Experiment 1
cd('L:\Experiments\MMNchroma\Analysis')
Definitions_MMNchroma;
grandFolder = [AnalysisFolder 'grandAverage\'];
mixedFolder = [AnalysisFolder 'MixedModel\'];
blocks = {'chroma','chromactrl','classic','classicctrl'};
bls = [2,4];
%% plot functional
%savedate = '04-Jan-2017';%'Fz'
if 1 %params
    addpath('L:\Experiments\N1P2\Analysis\N1P2_GH')
    savedate = '14-Nov-2018';%'Cz'
    include = false;%include subjects for which there was a manual inspection
end
load([mixedFolder 'LMEtables_combined_' savedate ])
relCols = 2:7;
if include
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
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end
for ib=bls
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
end
hf = plotFunctional_N1P2( allTables, relCols, blocks, bls );
legendstring = {'1','2'};
lh=legend(legendstring,'all');
set(lh,'position',[0.93 0.33 0.05 0.1],'fontsize',16)

suptitle(['Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))])

FigName = 'Exp1functional';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')

%% Experiment 2
cd('L:\Experiments\MMNchromaF\Analysis')
Definitions_MMNchromaF;
grandFolder = [AnalysisFolder 'grandAverage\'];
mixedFolder = [AnalysisFolder 'MixedModel\'];

allSubjects = 3:33;
badSubjects = [5,7,26];
%whichSubjects = allSubjects(~ismember(allSubjects,badSubjects));
%5 - effexor neurologic drug
%7 - speech deafness in the past
%18 - absolute pitch
%26 - lied that he was a professional musician, very bad performance

blocks = {'MMN_F_pure','MMN_Fctrl_pure','MMN_classic_pure','MMN_F_shep','MMN_Fctrl_shep','MMN_classic_shep'};
addpath('L:\Experiments\N1P2\Analysis\N1P2_GH')
bls = 2;
%% plot all cond, prevcond
if 1 %params
    addpath('L:\Experiments\N1P2\Analysis\N1P2_GH')
    savedate = '16-Nov-2018';%'Cz'
    include = false;%include subjects for which there was a manual inspection
end
load([mixedFolder 'LMEtables_combined_' savedate ])
relCols = 2:7;
%exclude subjects with no clear N1 or P2:
if include
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
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end
for ib=bls
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
end
hf = plotFunctional_N1P2( allTables, relCols, blocks, bls );
legendstring = {'1'};
lh=legend(legendstring,'all');
set(lh,'position',[0.93 0.33 0.05 0.1],'fontsize',16)

suptitle(['Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))])

FigName = 'Exp2functional';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')

%% Experiment 3
cd('L:\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
cd(GHfolder)
%% plot functional
if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    addtag = '';
    peakdate = '30-Oct-2018';
    include = false;%include subjects for which there was a manual inspection
end
%load tables
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])%loaded for allGrandcon_as_median
nCond = size(allGrandcon_times{1},2);
savedate = '31-Oct-2018';
load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
blocks = {'b1','b2a','b2b','b3a','b3b'};
ibs = 1:5;
relCols = 2:7;
%exclude subjects with no clear N1 or P2:
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
    allSubjects = 1:length(Subjects);
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end
for ib=ibs
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
end
%plot here:
hf = plotFunctional_N1P2( allTables, relCols, blocks, ibs );
legendstring = {'1','2a','2b','3a','3b'};
lh=legend(legendstring,'linear fit');
set(lh,'position',[0.93 0.33 0.05 0.1],'fontsize',16)

suptitle(['Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))])

FigName = 'Exp3functional';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
