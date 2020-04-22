FigFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\PaperFigures\Stimuli';
MixedFolder = 'L:\Experiments\N1P2\Analysis\MixedModel';
addpath(FigFolder)
addpath('S:\Lab-Shared\Experiments\MMNchroma\Analysis')
addpath('S:\Lab-Shared\Experiments\MMNchromaF\Analysis')
lw = 7;%linewidth
linelen = 10;
fromx = 0; tox=30; 
fromy = 200;toy = 2500;

colors = {[0.3 0.745 0.93],[0.75 0 0.75],[1 0 0],[0.75 0 0.75],[0.3 0.745 0.93]};
%linewidths = [2,2,2,2,2];
linestyles = {'-','-','-','-','-'};
fromxERP=-100;toxERP=300;fromyERP=-2.5;toyERP=3.2;
%% Experiment 1
cd('S:\Lab-Shared\Experiments\MMNchroma\Analysis')
Definitions_MMNchroma;%new idea for the first time! move definitions into a separate script!
%%      stimuli
s=10;
load([EDATfolder ExpName '_' Subjects{s} '_' sessions{s} '_expdata.mat'])
Fstandards = cell(size(bls));
Fdeviants = cell(size(bls));
for ib=1:length(bls)
    bl=bls(ib);
    Fstandards{bl} = expdata.stimParamsBlocksMMN(bl).toneSynth.Fstandards;
    Fdeviants{bl} = expdata.stimParamsBlocksMMN(bl).toneSynth.Fdeviants;
end

hf=ERPfigure;
set(hf,'Position',[100 100 350 400])

for ib=1:length(bls)
    subplot(1,length(bls),ib); 
    ax(ib) = gca;
    bl=bls(ib);
    for i=1:length(Fstandards{bl})
        line([10, 10+linelen],[Fstandards{bl}{i}, Fstandards{bl}{i}],'linew',lw,'Color',colors{cs==i})
        hold on
    end
    line([10, 10+linelen],[Fdeviants{bl}{1}, Fdeviants{bl}{1}],'linew',lw,'Color',colors{3})
    %title(strrep(blocks{bl},'_',' '))
    if ib==1
        ylabel('Frequencies [Hz]','fontsize',20);
    else
        set(gca,'ytick',[])
    end                       
    set(gca,'yscale','log')
    set(gca,'xtick',[])
    set(gca,'ytick',[400 1000 2000],'yticklabels',[400 1000 2000])
    set(gca,'fontsize',20)
end

linkaxes(ax)

axis([fromx, tox, fromy, toy])

FigName = 'Exp1Stim';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
%% Experiment 2
cd('S:\Lab-Shared\Experiments\MMNchromaF\Analysis')
Definitions_MMNchromaF;

addpath('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
bls = 4;%Dec. 27 2019 - this was changed from 2 to 4 to account for the order of blocks in the EDAT. 
cs = [1 2 5 3 4];
%%      stimuli
% plot frequency spread of all blocks
s=10;
load([EDATfolder ExpName '_' Subjects{s} '_' sessions{s} '_expdata.mat'])
Fstandards = cell(size(bls));
Fdeviants = cell(size(bls));
for ib=1:length(bls)
    bl=bls(ib);
    Fstandards{bl} = expdata.MMN.blocks(bl).toneSynth.Fstandards;
    Fdeviants{bl} = expdata.MMN.blocks(bl).toneSynth.Fdeviants;
end

hf=ERPfigure;
set(hf,'Position',[100 100 180 400])

for ib=1:length(bls)
    subplot(1,length(bls),ib); 
    ax(ib) = gca;
    bl=bls(ib);
    for i=1:length(Fstandards{bl})
        line([10, 10+linelen],[Fstandards{bl}{i}, Fstandards{bl}{i}],'linew',lw,'Color',colors{cs==i})
        hold on
    end
    line([10, 10+linelen],[Fdeviants{bl}{1}, Fdeviants{bl}{1}],'linew',lw,'Color',colors{3})
    title(strrep(expdata.MMN.blocks(bl).block_types,'_',' '))
    if ib==1 
        ylabel('Frequencies [Hz]','fontsize',20); 
    else
        set(gca,'ytick',[])
    end                       
    set(gca,'yscale','log')
    set(gca,'ytick',[400 1000 2000],'yticklabels',[400 1000 2000])
    set(gca,'fontsize',20)
    set(gca,'Xtick',[])
end

linkaxes(ax)

axis([fromx, tox, fromy, toy])

FigName = 'Exp2stim';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
%% Experiment 3
cd('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
cd(GHfolder)
%%      stimuli
s=10;
load([EDATfolder ExpName '_' Subjects{s} '_' sessions{s} '_expdata.mat'])
frequencies = cell(size(bls));
for ib=1:length(bls)
    bl=bls(ib);
    frequencies{bl} = expdata.Passive.blocks(bl).toneSynth.frequencies;
end

hf=ERPfigure;
set(hf,'Position',[100 100 800 400])

for ib=1:length(bls)
    subplot(1,length(bls),ib); 
    ax(ib) = gca;
    bl=bls(ib);
    for i=1:length(frequencies{bl})
        line([10, 10+linelen],[frequencies{bl}{i}, frequencies{bl}{i}],'linew',lw,'Color',colors{i})
        hold on
    end
    title(strrep(blocks{bl},'_',' '))
    if ib==1
        ylabel('Frequencies [Hz]','fontsize',20);
    else
        set(gca,'ytick',[])
    end                       
    set(gca,'yscale','log')
    set(gca,'ytick',[400 1000 2000],'yticklabels',[400 1000 2000])
    set(gca,'fontsize',20)
    set(gca,'Xtick',[])
end

linkaxes(ax)

axis([fromx, tox, fromy, toy])

FigName = 'Exp3stim';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
%% Example sequence
%from Exp3
dur = 0.1;%stim duration
bl=3;
bcs=[expdata.trials(:).block_code];
ntones=10;startidx=52;
bcidx=find(bcs==bl,1,'first');
sequence = [expdata.trials(bcidx+startidx-1:bcidx+ntones+startidx-2).trial_code]-10*bl;
SOAs = [expdata.trials(bcidx+startidx-1:bcidx+ntones+startidx-2).SOA]./1000;
h=figure;
set(h,'Position',[100 100 1000 300])
cumSOA = 0;
for ii=1:length(sequence)
    cumSOA = cumSOA + SOAs(ii);
    line([cumSOA, cumSOA+dur],[frequencies{bl}{sequence(ii)}, frequencies{bl}{sequence(ii)}],'linew',lw,'Color',colors{sequence(ii)},'linestyle',linestyles{sequence(ii)});
    hold on
end
df1=frequencies{bl}{2}-frequencies{bl}{1};
df4=frequencies{bl}{5}-frequencies{bl}{4};
ylim([frequencies{bl}{1}-df1/2 frequencies{bl}{5}+df4])
set(gca,'yscale','log')
set(gca,'ytick',[400 1000 2000],'yticklabels',[400 1000 2000])
xlabel('Time (seconds)');ylabel('Frequency (Hz)')
set(gca,'fontsize',24)

FigName = 'Example_sequence';
saveas(gcf,[FigFolder filesep FigName],'fig')
saveas(gcf,[FigFolder filesep FigName],'pdf')
    