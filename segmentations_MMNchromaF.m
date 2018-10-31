%segmentations_MMNchromaF
%% Definitions:
Definitions_MMNchromaF;
%% copy files to D:
% copy data files
for s = 23:25
    Subj = Subjects{s};
    disp(['Copying ' Subj '...'])
    sourceFolder = 'L:\Experiments\MMNchromaF\Experiment\Results\EEG\Export\beforeSegmentations\';
    destFolder = 'D:\Experiments\MMNchromaF\Experiment\Results\EEG\Export\beforeSegmentations\';
    mkdir(destFolder);
    fileName = ['MMNchromaF_' Subj '_' sessions{s}(1) '_dt_ImportMarkers_ReCoded'];
    source ={[sourceFolder fileName '.dat'],[sourceFolder fileName '.vhdr'],[sourceFolder fileName '.vmrk']};
    destination = {[destFolder fileName '.dat'],[destFolder fileName '.vhdr'],[destFolder fileName '.vmrk']};
    for i=1:length(source)
        copyfile(source{i},destination{i})
    end
    % EDAT
    sourceFolder = 'L:\Experiments\MMNchromaF\Experiment\Results\EDAT\';
    destFolder = 'D:\Experiments\MMNchromaF\Experiment\Results\EDAT\';
    mkdir(destFolder);
    fileName = ['MMNchromaF_' Subj '_' sessions{s}(1) '_expdata'];
    source =[sourceFolder fileName '.mat'];
    destination = [destFolder fileName '.mat'];
    copyfile(source,destination)
    disp(['Done Subj ' Subj])
end
%Done until 213
%% read data
%do this in D folder after previous section done
for s = 25
    tic
    Subj = Subjects{s};
    disp(['Loading ' Subj '...'])
    if strcmp(sessions{s},'')
        FileName = [ExpName '_' Subj '_dt_ImportMarkers_ReCoded'];
    else
        FileName = [ExpName '_' Subj '_' sessions{s}(1) '_dt_ImportMarkers_ReCoded'];
    end
    datFile = [ExportFolder FileName '.dat'];
    [allData] = read_analyzer_UnsegmentedData(datFile);
    %save([ExportFolder FileName],'allData')
    save([AnalysisFolder Subj filesep FileName],'allData')
    disp(['Done ' Subj ' in ' num2str(toc) 'sec.'])
end
% L - Done s201 in 117.9241sec.
% D - Done s201 in 118.8348sec.
%Done s202 in 94.7282sec.
%Done until 225
%% individual subject ERPs
%codes and names table:
names = {'Training_pure','Training_shep','MMN_F_pure','MMN_Fctrl_pure','MMN_classic_pure','MMN_F_shep','MMN_Fctrl_shep','MMN_classic_shep','Perception_pure','Perception_shep'}';
PBcodes = [1200, 1400, 2300, 2400, 2500, 2800, 2900, 21000, 3200,3400]';
codes_names = table(names,PBcodes);
%% segment and plot and save?:
%params:
name = 'MMN_classic_pure';
saveFlag = 0;
plotFlag = 1;

%mode = '';%adding to the saving name, e.g. LP30
PBcode = codes_names.PBcodes(strcmp(codes_names.names,name));
conditions = {'dev','stan1','stan2','stan3','stan4'};
STIMcodes = [11,20,30,40,50];
EventCodes=cell(1,length(conditions));
for i=1:length(conditions)
    EventCodes{i} = PBcode + STIMcodes(i);
end
electrodes = 1:74;
srate = 512;
winbegin = -0.1; winend = 0.4;%in sec
win = round(winbegin*srate:winend*srate);%in samples
bsl = 1:floor(0.1*srate);
plot_bsl = 1:floor(0.1*srate);
fromx= -0.1; tox=winend;fromy=-3;toy=3;
lp=50;

elecName = 'Fz';%for plot
for s = 3
    tic
    Subj = Subjects{s};
    subjFolder = [AnalysisFolder Subj '\segmentations\'];
    [ allSEGS, ntrials ] = segmentERPs(s, Expinfo, conditions, EventCodes, electrodes, win, bsl, lp, srate);
    t=0:1:size(allSEGS{1},1)-1;t=t/srate;t=t-0.1;
    if saveFlag
        savedate = date;
        disp('saving...')
        save([subjFolder filesep Subj '_SEGS_' name],'allSEGS','ntrials','conditions','EventCodes','winbegin','winend','bsl','electrodes','srate','lp','t','savedate');
    end
    if plotFlag
        [ electrode ] = chName2n( elecName );
        h = ERPfigure;
        set(h,'position',[100 100 500 800]);
        for i=1:length(allSEGS)
            subplot(length(allSEGS)+1,1,i)
            data = squeeze(allSEGS{i}(:,electrode,:));
            data = bsxfun(@minus,data,mean(data(plot_bsl,:),1));
            hold on
            t=0:1:size(data,1)-1;t=t/srate;t=t-0.1;
            varplot(t,data);
            ax(i)=gca;
            line([0,0],[fromy, toy],'linestyle','-.','col','k');
            line([fromx, tox],[0, 0],'linestyle','-.','col','k');
            title(conditions{i})
        end
        subplot(length(allSEGS)+1,1,length(allSEGS)+1)
        title('overlay')
        hold all
        ax(i+1)=gca;
        for i=1:length(allSEGS)
            data = squeeze(allSEGS{i}(:,electrode,:));
            data = bsxfun(@minus,data,mean(data(plot_bsl,:),1));
            t=0:1:size(data,1)-1;t=t/srate;t=t-0.1;
            varplot(t,data);
        end
        legend(conditions)
        line([0,0],[fromy, toy],'linestyle','-.','col','k');
        line([fromx, tox],[0, 0],'linestyle','-.','col','k');
        linkaxes(ax);
        axis([fromx,tox,fromy,toy]);
        suptitle([Subj ' ' strrep(name,'_',' ') ' ' elecName ' LP' num2str(lp)]);
        figName = [Subj '_' name '_' elecName];
        set(gcf,'name',figName,'numbertitle','off')
    end
    disp(['Done ' name ' in ' num2str(toc)])
end
%% (load and plot:)
%faster if saved already
name = 'Perception_shep';
mode = '';%adding to the saving name, e.g. LP30
%matrixdate = '01-May-2016';    
srate=512;
new_bsl = 1:floor(0.1*srate);
fromx=-0.1;tox=0.8;fromy=-20;toy=20;
elecName = 'Pz';
[ electrode ] = chName2n( elecName );
lp=50;%check that this is the saved case
for s = 3
    Subj = Subjects{s};
    subjFolder = [AnalysisFolder Subj '\segmentations\'];
    %load([subjFolder filesep 'SEGS_' name '_LP' num2str(lp) '_' matrixdate]);
    load([subjFolder filesep Subj '_SEGS_' name]);
    h=ERPfigure;
    set(h,'position',[100 100 500 800])
    for i=1:length(allSEGS)
        subplot(length(allSEGS)+1,1,i)
        data = squeeze(allSEGS{i}(:,electrode,:));
        data = bsxfun(@minus,data,mean(data(new_bsl,:),1));
        hold on
        t=0:1:size(data,1)-1;t=t/srate;t=t-0.1;
        varplot(t,data)
        ax(i)=gca;
        line([0,0],[fromy, toy],'linestyle','-.','col','k');
        line([fromx, tox],[0, 0],'linestyle','-.','col','k');
        title(conditions{i})
    end
    subplot(length(allSEGS)+1,1,length(allSEGS)+1)
    title('overlay')
    hold all
    ax(i+1)=gca;
    for i=1:length(allSEGS)
        data = squeeze(allSEGS{i}(:,electrode,:));
        data = bsxfun(@minus,data,mean(data(new_bsl,:),1));
        t=0:1:size(data,1)-1;t=t/srate;t=t-0.1;
        varplot(t,data)
    end
    legend(conditions)
    line([0,0],[fromy, toy],'linestyle','-.','col','k');
    line([fromx, tox],[0, 0],'linestyle','-.','col','k');
    linkaxes(ax)
    axis([fromx,tox,fromy,toy])
    suptitle([Subj ' ' strrep(name,'_',' ') ' ' elecName ' LP' num2str(lp)])
    figName = [Subj '_' name '_' elecName];
    %saveas(gcf,[AnalysisFolder Subj filesep 'figures\' figName ],'fig')
end
%% topoplot - should fit the last plot
headFile = 'K:\CommonResources\Tools\Matlab_Tools\head64_mastoids.locs';
con = 1;
figure
title([Subj ' ' strrep(name,'_',' ') '. Time: ' num2str(time(1)) ' - ' num2str(time(2)) ' sec.'])
time = [0.08 0.15];
samp(1) = find(abs(t-time(1)) == min(abs(t-time(1))));
samp(2) = find(abs(t-time(2)) == min(abs(t-time(2))));
SEGS = allSEGS{con}(samp(1):samp(2),:,:);
v = mean(mean(SEGS,3),1);
v=v([1:64,70,71]);
topoplot(v,headFile,'electrodes','on','style','map','shading','interp');colorbar
%topoplot(v,headFile,'electrodes','on','style','map','shading','interp','maplimits',[-2 2]);
%% stacked
con = 1;
elecName = 'Fz';
[ electrode ] = chName2n( elecName );
data = squeeze(allSEGS{con}(:,electrode,:));
imagesc(t,1:size(data,2),data')
hold on
line([0 0],[1,size(data,2)],'color','k')
colorbar
%%
