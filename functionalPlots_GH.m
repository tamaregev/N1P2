% functional plots Nov 1 2018
% plotting N1P2 amplitudes as a function of the explaining variables
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
bls=1:5;
%% plot all subj, cond, prevcond
if 1 %params
    electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
    addtag = '';
    peakdate = '30-Oct-2018';
    whichpeaks = {'N1','P2'};
    include = false;%include subjects for which there was a manual inspection
end
%load tables
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])
nCond = size(allGrandcon_times{1},2);
savedate = '31-Oct-2018';
load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
blocks = {'b1','b2a','b2b','b3a','b3b'};

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
end

ERPfigure
xvars = {'currMIDI','dist_mean','size_jump'};
xnames = {'curr freq','mean freq','prev freq'};
nx=length(xvars);
for ipeak = 1:2
    peak = whichpeaks{ipeak};
    for ix=1:length(xvars)
        subplot(2,3,(ipeak-1)*nx+ix)
        x=nan(1,length(includeSubjects)*5*4*5);y=nan(1,length(includeSubjects)*5*4*5);
        for ib=1:length(blocks)
            T = allTables.(blocks{ib});
            from = ((ib-1)*length(includeSubjects)*5*4+1);
            to = (ib*length(includeSubjects)*5*4);
            x(from:to) = T.(xvars{ix});
            y(from:to) = T.(peak);
        end
        %[xs,I] = sort(x);
        %y=y(I);
        scatter(x,y,10,'filled','MarkerFaceColor',[0.5 0.5 0.5]);
        title([peak ' (' xnames{ix} ')'])
        hold on
        p1=polyfit(x,y,1);
        f1=polyval(p1,x);
        plot(x,f1,'r','linewidth',3)
        p2=polyfit(x,y,2);
        f2=polyval(p2,x);
        
        plot(x,f2,'.b','MarkerSize',20)
        legend('data','linear fit','parabolic fit')
    end
end
%% plot all cond, prevcond
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
for ib=ibs
    allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
end

hf = plotFunctional_N1P2( allTables, relCols, blocks, ibs );

suptitle(['Exp: ''' ExpName '''.  N=' num2str(length(includeSubjects))])

%% plot frequency spread of all blocks
s=10;
load([EDATfolder ExpName '_' Subjects{s} '_' sessions{s} '_expdata.mat'])
frequencies = cell(size(bls));
for ib=1:length(bls)
    bl=bls(ib);
    frequencies{bl} = expdata.Passive.blocks(bl).toneSynth.frequencies;
end

hf=ERPfigure;
set(hf,'Position',[100 100 500 400])
linethick = 10;
fromx = 0; tox=30;
fromy = 200;toy = 2500;

for ib=1:length(bls)
    subplot(1,length(bls),ib); 
    ax(ib) = gca;
    bl=bls(ib);
    for i=1:length(frequencies{bl})
        line([10, 10+linethick],[frequencies{bl}{i}, frequencies{bl}{i}],'linew',4)
        hold on
    end
    title(strrep(blocks{bl},'_',' '))
    if ib==1
        ylabel('log frequencies [Hz]','fontsize',20);
    else
        set(gca,'ytick',[])
    end                       
    set(gca,'yscale','log')
    set(gca,'xtick',[])
end

linkaxes(ax)

axis([fromx, tox, fromy, toy])

