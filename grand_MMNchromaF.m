%grand_MMNchromaF

%% definitions
Definitions_MMNchromaF

%for grand -
grandFolder = [AnalysisFolder 'grandAverage'];
srate = 512;
Expinfo.srate = srate;

%codes and names table:
names = {'Training_pure','Training_shep','MMN_F_pure','MMN_Fctrl_pure','MMN_classic_pure','MMN_F_shep','MMN_Fctrl_shep','MMN_classic_shep','Perception_pure','Perception_shep'}';
PBcodes = [1200, 1400, 2300, 2400, 2500, 2800, 2900, 21000, 3200,3400]';
codes_names = table(names,PBcodes);

%% calculate and save grand matrix
for ii=9:10
    ticcond = tic;
    name = codes_names.names{ii};
    tag = 'BP';

    % function paramettters
    whichSubjects=3:25;
    %winbegin = -0.1; winend = 0.5; 
    winbegin = -0.1; winend = 1; 
    win = round(winbegin*srate:winend*srate);
    bsl = 1:floor(-winbegin*srate);
    nelectrodes = 74;
    Bp = [0.1 50];
    
    PBcode = codes_names.PBcodes(strcmp(codes_names.names,name));
    conditions = {'dev','stan1','stan2','stan3','stan4'};
    STIMcodes = {11,20,30,40,50};
    EventCodes=cell(1,length(conditions));
    for i=1:length(conditions)
        EventCodes{i} = PBcode + STIMcodes{i};
    end

    %run function
    [ grandMatrix, ntrials ] = calcGrandMatrix_MMNchromaF(Expinfo, whichSubjects, conditions, EventCodes, nelectrodes, win, bsl, Bp  );
    tic;disp('saving...');matrixdate = date;
    save([grandFolder filesep 'grandMatrix_' name '_' tag],'grandMatrix','ntrials','conditions','EventCodes','whichSubjects','ntrials','win','bsl','matrixdate','Bp');
    disp(['Done saving ' name ' in ' num2str(toc) ' sec']);
    %Done subject s213 in 14.7718 sec.
    %Done all in 168.8481 sec.
    disp(['Done ' name ' in ' num2str(toc(ticcond)) ' sec']);

end

%% load and plot one
name = 'MMN_classic_pure';
tag = 'BP';
 load([grandFolder filesep 'grandMatrix_' name '_' tag])
%load([grandFolder filesep 'grandMatrix_' name])
plotSubjects = whichSubjects; 
whichCon =1:length(conditions);
colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};

tictot = tic;
srate = 512;nSamps = size(grandMatrix,1);
fromx = -100; tox = 400;
t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms
fromy=-5;toy=5;   
elecName = 'Fz';
[ electrode ] = chName2n( elecName );

axisfontsize = 14;
labelsfontsize=16;
titlefontsize = 14;
legendfontsize = 12;

%plot
h = ERPfigure;
set(h,'Position',[100 100 600 400])
set(h,'Name',['grandAvg_' name],'NumberTitle','off')

for con = whichCon 
    cond = conditions{con};
    grandData = squeeze(grandMatrix(:,plotSubjects,con,electrode));
    varplot(t,grandData,'ci',0.95,'Color',colors{con},'LineWidth',2)
    %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
    hold on
end
axis([fromx tox, fromy toy])
set(gca,'FontSize',axisfontsize)
legend(conditions(whichCon),'FontSize',legendfontsize,'Location','northeast')
ylabel('voltage [\muV]','FontSize',labelsfontsize)
line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
set(gca,'FontSize',axisfontsize)     
xlabel('time [ms]','FontSize',labelsfontsize);
title(['Grand avg. '  strrep(name,'_',' ') ' N = ' num2str(length(plotSubjects)) ', ' elecName ])
%    
disp(['Done  in ' num2str(toc(tictot)) 'sec.'])
%% load and plot all MMNs

%plot
h = ERPfigure;
set(h,'Position',[100 100 800 600])
set(h,'Name','grandAvg_MMNs','NumberTitle','off')
axisfontsize = 12;
labelsfontsize=12;
titlefontsize = 12;
legendfontsize = 10;
    
srate = 512;
fromx = -100; tox = 500;
fromy=-5;toy=6;   
elecName = 'Fz';
[ electrode ] = chName2n( elecName );

i=0;
for n=3:8
    tictot = tic;
    i=i+1;
    subplot(2,3,i)
    name = names{n};
    tag = 'BP';
    load([grandFolder filesep 'grandMatrix_' name '_' tag])
    %load([grandFolder filesep 'grandMatrix_' name])
    %plotSubjects = whichSubjects; 
    plotSubjects =good;
    %plotSubjects = good;
    whichCon =1:length(conditions);
    colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};
    nSamps = size(grandMatrix,1);
    t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms

    for con = whichCon 
        cond = conditions{con};
        grandData = squeeze(grandMatrix(:,plotSubjects,con,electrode));
        varplot(t,grandData,'ci',0.95,'Color',colors{con},'LineWidth',2)
        %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
        hold on
    end
    axis([fromx tox, fromy toy])
    set(gca,'FontSize',axisfontsize)
    ylabel('voltage [microV]','FontSize',labelsfontsize)
    line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
    set(gca,'FontSize',axisfontsize)     
    xlabel('time [ms]','FontSize',labelsfontsize);
    title(strrep(name,'_',' '))
    %    
    disp(['Done  in ' num2str(toc(tictot)) 'sec.'])
end
suptitle(['Grand avg. N = ' num2str(length(plotSubjects)) ', ' elecName ])
    legend(conditions(whichCon),'FontSize',legendfontsize,'Location','northeast')

%% separate according to performance
%% load behaviors
% Training, dprims: from script calcBehaviour_training -
load([AnalysisFolder 'diprime_training'],'diprime_training')
diprimeDiff_training = diprime_training{2} - diprime_training{1};

% Trainig, nerrs: from script behavioral_results -
load([AnalysisFolder 'behaviors_training'],'behaviors_training')
%whichSubjects = [3:4,6:20];
whichSubjects = [3:25];
errOfTotalShep = nan(whichSubjects(end),1);
errOfTotalPure = nan(whichSubjects(end),1);
for i=whichSubjects
    behav = behaviors_training{i};
    errOfTotalShep(i) = (sum(behav.miss(behav.block=='shep'))+sum(behav.no_reaction(behav.block=='shep')))/(sum(behav.no_reaction(behav.block=='shep'))+sum(behav.miss(behav.block=='shep'))+sum(behav.hit(behav.block=='shep')));
    errOfTotalPure(i) = (sum(behav.miss(behav.block=='pure'))+sum(behav.no_reaction(behav.block=='pure')))/(sum(behav.no_reaction(behav.block=='pure'))+sum(behav.miss(behav.block=='pure'))+sum(behav.hit(behav.block=='pure'))); 
end
difShepVSPure = errOfTotalPure - errOfTotalShep;
errTot = mean([errOfTotalPure errOfTotalShep],2);
% 
% figure
% hist(errTot);
diprime_train_pure = nan(whichSubjects(end),1);
diprime_train_shep = nan(whichSubjects(end),1);
diprime_train = nan(whichSubjects(end),1);
    diprime_train_pure = diprime_training{1}(1:whichSubjects(end));
    diprime_train_shep = diprime_training{2}(1:whichSubjects(end));
    diprime_train = mean([diprime_train_pure,diprime_train_shep],2);
measure = diprime_train;
%measure(5) = nan;
%measure = errTot;
mid=median(measure(~isnan(measure)));
good = find(measure>mid);
bad = find(measure<=mid);
%% plot due to behavior
groups = {good,bad}; titles = {'good','bad'};
whichCon = 1:length(conditions);
colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};

tictot = tic;
srate = 512;nSamps = size(grandMatrix,1);
fromx = -100; tox = 500;
t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms
fromy=-5;toy=5;   
elecName = 'Fz';
[ electrode ] = chName2n( elecName );

axisfontsize = 14;
labelsfontsize=16;
titlefontsize = 14;
legendfontsize = 12;

%plot
h = ERPfigure;
set(h,'Position',[100 100 600 400])
set(h,'Name',['ClassicMMNs_behavior'],'NumberTitle','off')

names = {'MMN_classic_pure','MMN_classic_shep'};
tag = 'BP';

for n=1:length(names)
    for ss=1:length(groups)
        subplot(2,2,(n-1)*length(names)+ss)
        plotSubjects = groups{ss};
        name=names{n};
        load([grandFolder filesep 'grandMatrix_' name '_' tag])
        for con = whichCon 
            cond = conditions{con};
            grandData = squeeze(grandMatrix(:,plotSubjects,con,electrode));
            varplot(t,grandData,'ci',0.95,'Color',colors{con},'LineWidth',2)
            %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
            hold on
        end
        axis([fromx tox, fromy toy])
        set(gca,'FontSize',axisfontsize)
        ylabel('voltage [microV]','FontSize',labelsfontsize)
        line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
        set(gca,'FontSize',axisfontsize)     
        xlabel('time [ms]','FontSize',labelsfontsize);
        title([strrep(name,'_',' ') ' ' titles{ss} ])
    end    
end
suptitle(['Grand avg. ' elecName ])
legend(conditions(whichCon),'FontSize',legendfontsize,'Location','northeast')
disp(['Done  in ' num2str(toc(tictot)) 'sec.'])
   %% plot relevant
% measure = diprimeTot_attend(1:20);
% mid=median(measure(~isnan(measure))); 
% good = find(measure>=mid);
% bad = find(measure<mid);

   groups = {good,bad}; titles = {'good','bad'};
whichCon ={1,2:5,6};%6 is control
Conditions = {'dev','stan','dev ctrl'};
colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};
tictot = tic;
srate = 512;nSamps = size(grandMatrix,1);
fromx = -100; tox = 400;
t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms
fromy=-5;toy=5;   
elecName = 'Fz';
[ electrode ] = chName2n( elecName );

axisfontsize = 12;
labelsfontsize=12;
titlefontsize = 14;
legendfontsize = 12;

%plot
h = ERPfigure;
set(h,'Position',[100 100 600 400])
set(h,'Name',['ClassicMMNs_behavior'],'NumberTitle','off')

names = {'MMN_classic_pure','MMN_classic_shep'};
ctrls = {'MMN_Fctrl_pure','MMN_Fctrl_shep'};
tag = 'BP';
 
for n=1:length(names)
    for ss=1:length(groups)
        subplot(2,2,(n-1)*length(names)+ss)
        plotSubjects = groups{ss};
        name=names{n};
        load([grandFolder filesep 'grandMatrix_' name '_' tag])
        for con = 1:length(whichCon)
            cond = Conditions{con};
            if whichCon{con}==6
                load([grandFolder filesep 'grandMatrix_' ctrls{n} '_' tag])
                grandData = squeeze(grandMatrix(:,plotSubjects,3,electrode));
                varplot(t,grandData,'ci',0.95,'Color',colors{con},'LineWidth',2)
               
            elseif length(whichCon{con})>1
                grandData = mean(squeeze(grandMatrix(:,plotSubjects,whichCon{con},electrode)),3);
                varplot(t,grandData,'ci',0.95,'Color',colors{con},'LineWidth',2)
                %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
            else
                grandData = squeeze(grandMatrix(:,plotSubjects,whichCon{con},electrode));
                varplot(t,grandData,'ci',0.95,'Color',colors{con},'LineWidth',2)
                %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
            end
            hold on
        end
        axis([fromx tox, fromy toy])
        set(gca,'FontSize',axisfontsize)
        ylabel('voltage [microV]','FontSize',labelsfontsize)
        line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
        set(gca,'FontSize',axisfontsize)     
        xlabel('time [ms]','FontSize',labelsfontsize);
        title([strrep(name,'_',' ') ' ' titles{ss} ])
    end    
end
suptitle(['Grand avg. ' elecName ])
legend(Conditions,'FontSize',legendfontsize,'Location','northeast')
disp(['Done  in ' num2str(toc(tictot)) 'sec.'])
   %% plot individuals
        
    [Y,I] = sort(errTot);
    subjectList = I(~ismember(I,[1 2]));

    whichCon ={1,2:5,6};%6 is control
Conditions = {'dev','stan','dev ctrl'};
colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};
tictot = tic;
srate = 512;nSamps = size(grandMatrix,1);
fromx = -100; tox = 500;
t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms
fromy=-5;toy=5;   
elecName = 'Fz';
[ electrode ] = chName2n( elecName );

axisfontsize = 14;
labelsfontsize=16;
titlefontsize = 14;
legendfontsize = 12;

%plot
h = ERPfigure;
set(h,'Position',[100 100 600 400])
set(h,'Name',['ClassicMMNs_behavior'],'NumberTitle','off')

names = {'MMN_classic_pure','MMN_classic_shep'};
ctrls = {'MMN_Fctrl_pure','MMN_Fctrl_shep'};
tag = 'BP';
 
for n=1:length(names)
    for ss=1:length(subjectList)
        sub=subjectList(ss);
        subplot(length(subjectList),2,(ss-1)*length(names)+n)
        plotSubjects = sub;
        name=names{n};
        load([grandFolder filesep 'grandMatrix_' name '_' tag])
        for con = 1:length(whichCon)
            cond = Conditions{con};
            if whichCon{con}==6
                load([grandFolder filesep 'grandMatrix_' ctrls{n} '_' tag])
                grandData = squeeze(grandMatrix(:,plotSubjects,3,electrode));
                plot(t,grandData,'Color',colors{con},'LineWidth',2)
                
            elseif length(whichCon{con})>1
                grandData = mean(squeeze(grandMatrix(:,plotSubjects,whichCon{con},electrode)),2);
                plot(t,grandData,'Color',colors{con},'LineWidth',2)
                %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
            else
                grandData = squeeze(grandMatrix(:,plotSubjects,whichCon{con},electrode));
                plot(t,grandData,'Color',colors{con},'LineWidth',2)
                %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
            end
            hold on
        end
        axis([fromx tox, fromy toy])
        set(gca,'FontSize',axisfontsize)
        line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
        set(gca,'FontSize',axisfontsize)     
        if ss==length(subjectList)
            xlabel('time [ms]','FontSize',labelsfontsize);
            ylabel('voltage [microV]','FontSize',labelsfontsize)
        end
        
        if ss==1
            switch name
                case 'MMN_classic_pure'
                    title({[strrep(name,'_',' ')],['s' num2str(sub) ' err = ' num2str(errOfTotalPure(sub))]})
                case 'MMN_classic_shep'
                    title({[strrep(name,'_',' ')],['s' num2str(sub) ' err = ' num2str(errOfTotalShep(sub))]})
            end
        else
            switch name
                case 'MMN_classic_pure'
                    title(['s' num2str(sub) ' err = ' num2str(errOfTotalPure(sub))])
                case 'MMN_classic_shep'
                    title(['s' num2str(sub) ' err = ' num2str(errOfTotalShep(sub))])
            end
        end
        
    end    
end
suptitle(['Grand avg. ' elecName ])
legend(Conditions,'FontSize',legendfontsize,'Location','northeast')
disp(['Done  in ' num2str(toc(tictot)) 'sec.'])
   %% plot individuals - only diff
        
    [Y,I] = sort(errTot);
    %subjectList = I(~ismember(I,[1 2 5]));
    subjectList = I(~ismember(I,[1 2]));

    whichCon ={7};%6 is difference
Conditions = {'dev - dev ctrl'};
colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};
tictot = tic;
srate = 512;nSamps = size(grandMatrix,1);
fromx = -100; tox = 500;
t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms
fromy=-5;toy=5;   
elecName = 'Fz';
[ electrode ] = chName2n( elecName );

axisfontsize = 14;
labelsfontsize=16;
titlefontsize = 14;
legendfontsize = 12;

%plot
h = ERPfigure;
set(h,'Position',[100 100 600 400])
set(h,'Name',['ClassicMMNs_behavior'],'NumberTitle','off')

names = {'MMN_classic_pure','MMN_classic_shep'};
ctrls = {'MMN_Fctrl_pure','MMN_Fctrl_shep'};
tag = 'BP';
 
for n=1:length(names)
    for ss=1:length(subjectList)
        sub=subjectList(ss);
        subplot(length(subjectList),2,(ss-1)*length(names)+n)
        plotSubjects = sub;
        name=names{n};
        load([grandFolder filesep 'grandMatrix_' name '_' tag])
        Dev = squeeze(grandMatrix(:,plotSubjects,1,electrode));
        load([grandFolder filesep 'grandMatrix_' ctrls{n} '_' tag])
        Dev_ctrl = squeeze(grandMatrix(:,plotSubjects,3,electrode));
        plot(t,Dev-Dev_ctrl,'Color',colors{con},'LineWidth',2)
     
        axis([fromx tox, fromy toy])
        set(gca,'FontSize',axisfontsize)
        line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
        set(gca,'FontSize',axisfontsize)     
        if ss==length(subjectList)
            xlabel('time [ms]','FontSize',labelsfontsize);
            ylabel('voltage [microV]','FontSize',labelsfontsize)
        end
        
        if ss==1
            switch name
                case 'MMN_classic_pure'
                    title({[strrep(name,'_',' ')],['s' num2str(sub) ' err = ' num2str(errOfTotalPure(sub))]})
                case 'MMN_classic_shep'
                    title({[strrep(name,'_',' ')],['s' num2str(sub) ' err = ' num2str(errOfTotalShep(sub))]})
            end
        else
            switch name
                case 'MMN_classic_pure'
                    title(['s' num2str(sub) ' err = ' num2str(errOfTotalPure(sub))])
                case 'MMN_classic_shep'
                    title(['s' num2str(sub) ' err = ' num2str(errOfTotalShep(sub))])
            end
        end
        
    end    
end
suptitle(['Grand avg. ' elecName ])
legend(Conditions,'FontSize',legendfontsize,'Location','northeast')
disp(['Done  in ' num2str(toc(tictot)) 'sec.'])

%% load and plot perception


colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};

tictot = tic;
srate = 512;nSamps = size(grandMatrix,1);
fromx = -100; tox = 800;
t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms
fromy=-10;toy=20;   
elecName = 'Pz';
[ electrode ] = chName2n( elecName );

axisfontsize = 14;
labelsfontsize=16;
titlefontsize = 14;
legendfontsize = 12;

%plot
h = ERPfigure;
set(h,'Position',[100 100 600 400])
set(h,'Name',['grandAvg_Perceptions'],'NumberTitle','off')

names = {'Perception_pure','Perception_shep'};
tag = 'BP';
 
for n=1:length(names)
    subplot(1,2,n)
    name=names{n};
    load([grandFolder filesep 'grandMatrix_' name '_' tag])
    plotSubjects = whichSubjects; 
    whichCon =1:length(conditions);
    for con = whichCon 
        cond = conditions{con};
        grandData = squeeze(grandMatrix(:,plotSubjects,con,electrode));
        varplot(t,grandData,'ci',0.95,'Color',colors{con},'LineWidth',2)
        %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
        hold on
    end
    axis([fromx tox, fromy toy])
    set(gca,'FontSize',axisfontsize)
    ylabel('voltage [microV]','FontSize',labelsfontsize)
    line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
    set(gca,'FontSize',axisfontsize)     
    xlabel('time [ms]','FontSize',labelsfontsize);
    title(strrep(name,'_',' '))
    %    
end
suptitle(['Grand avg. N = ' num2str(length(plotSubjects)) ', ' elecName ])
legend(conditions(whichCon),'FontSize',legendfontsize,'Location','northeast')
disp(['Done  in ' num2str(toc(tictot)) 'sec.'])
    %% ... due to behavior
    % Load dprime attend: from script calcBehaviour_attend
load([AnalysisFolder 'behaviour_attend'])
diprimeDiff_attend = behaviour_attend{2}.diprime - behaviour_attend{1}.diprime;
diprimeTot_attend = mean([behaviour_attend{2}.diprime, behaviour_attend{1}.diprime],2);

measure = diprimeTot_attend(1:20);
mid=median(measure(~isnan(measure))); 
good = find(measure>=mid);
bad = find(measure<mid);

groups = {good,bad}; titles = {'good','bad'};
whichCon =1:length(conditions);
colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};

tictot = tic;
srate = 512;nSamps = size(grandMatrix,1);
fromx = -100; tox = 800;
t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms
fromy=-10;toy=20;   
elecName = 'Pz';
[ electrode ] = chName2n( elecName );

axisfontsize = 14;
labelsfontsize=16;
titlefontsize = 14;
legendfontsize = 12;

%plot
h = ERPfigure;
set(h,'Position',[100 100 600 400])
set(h,'Name',['ClassicMMNs_behavior'],'NumberTitle','off')

names = {'Perception_pure','Perception_shep'};
tag = 'BP';
 
for n=1:length(names)
    for ss=1:length(groups)
        subplot(2,2,(n-1)*length(names)+ss)
        plotSubjects = groups{ss};
        name=names{n};
        load([grandFolder filesep 'grandMatrix_' name '_' tag])
        for con = whichCon 
            cond = conditions{con};
            grandData = squeeze(grandMatrix(:,plotSubjects,con,electrode));
            varplot(t,grandData,'ci',0.95,'Color',colors{con},'LineWidth',2)
            %varplot(t,grandData,'Color',colors{con},'LineWidth',2)
            hold on
        end
        axis([fromx tox, fromy toy])
        set(gca,'FontSize',axisfontsize)
        ylabel('voltage [microV]','FontSize',labelsfontsize)
        line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
        set(gca,'FontSize',axisfontsize)     
        xlabel('time [ms]','FontSize',labelsfontsize);
        title([strrep(name,'_',' ') ' ' titles{ss} ])
    end    
end
suptitle(['Grand avg. ' elecName ])
legend(conditions(whichCon),'FontSize',legendfontsize,'Location','northeast')
disp(['Done  in ' num2str(toc(tictot)) 'sec.'])
    %% individuals
   
load([AnalysisFolder 'behaviour_attend'])
diprimeDiff_attend = behaviour_attend{2}.diprime - behaviour_attend{1}.diprime;
diprimeTot_attend = mean([behaviour_attend{2}.diprime, behaviour_attend{1}.diprime],2);
measures = {behaviour_attend{1}.diprime,behaviour_attend{2}.diprime};

whichCon =1:length(conditions);
colors = {[0 0 0],[1 0 0],[0 0 1],[1 0.5 0],[0 0.5 0.1]};
tictot = tic;
srate = 512;nSamps = size(grandMatrix,1);
fromx = -100; tox = 800;
t = (0:(1/srate):(nSamps-1)/srate)*1000+fromx;%in ms
fromy=-10;toy=40;   
elecName = 'Pz';
[ electrode ] = chName2n( elecName );

axisfontsize = 14;
labelsfontsize=16;
titlefontsize = 14;
legendfontsize = 12;

%plot
h = ERPfigure;
set(h,'Position',[100 100 600 400])
set(h,'Name',['ClassicMMNs_behavior'],'NumberTitle','off')

names = {'Perception_pure','Perception_shep'};
tag = 'BP';
 
for n=1:length(names)
    measure = measures{n};

    [Y,I] = sort(measure);
    subjectList = I(~ismember(I,[1 2]) & I<=20);

    for ss=1:length(subjectList)
        sub=subjectList(ss);
        subplot(length(subjectList),2,(ss-1)*length(names)+n)
        plotSubjects = sub;
        name=names{n};
        load([grandFolder filesep 'grandMatrix_' name '_' tag])
        for con = whichCon 
            cond = conditions{con};
            grandData = squeeze(grandMatrix(:,plotSubjects,con,electrode));
            plot(t,grandData,'Color',colors{con},'LineWidth',2)
            hold on
        end
        axis([fromx tox, fromy toy])
        set(gca,'FontSize',axisfontsize)
        line([0 0],[fromy toy],'Color','k','LineStyle','-.');line([fromx, tox],[0 0],'Color','k','LineStyle','-.')
        set(gca,'FontSize',axisfontsize)     
        if ss==length(subjectList)
            xlabel('time [ms]','FontSize',labelsfontsize);
            ylabel('voltage [microV]','FontSize',labelsfontsize)
        end
        
        if ss==1
            switch name
                case 'Perception_pure'
                    title({[strrep(name,'_',' ')],['s' num2str(sub) ' dprime = ' num2str(behaviour_attend{1}.diprime(sub))]})
                case 'Perception_shep'
                    title({[strrep(name,'_',' ')],['s' num2str(sub) ' dprime = ' num2str(behaviour_attend{2}.diprime(sub))]})
            end
        else
            switch name
                case 'Perception_pure'
                    title(['s' num2str(sub) ' dprime = ' num2str(behaviour_attend{1}.diprime(sub))])
                case 'Perception_shep'
                    title(['s' num2str(sub) ' dprime = ' num2str(behaviour_attend{2}.diprime(sub))])
            end
        end
        
    end    
end
suptitle(['Grand avg. ' elecName ])
legend(Conditions,'FontSize',legendfontsize,'Location','northeast')
disp(['Done  in ' num2str(toc(tictot)) 'sec.'])
%% topoplot - should fit the last plot
headFile = 'K:\CommonResources\Tools\Matlab_Tools\head64_mastoids.locs';
con = 1;
time = [0.08 0.15];
figure
title([Subj ' ' strrep(name,'_',' ') '. Time: ' num2str(time(1)) ' - ' num2str(time(2)) ' sec.'])
samp(1) = find(abs(t-time(1)) == min(abs(t-time(1))));
samp(2) = find(abs(t-time(2)) == min(abs(t-time(2))));
SEGS = allSEGS{con}(samp(1):samp(2),:,:);
v = mean(mean(SEGS,3),1);
v=v([1:64,70,71]);
topoplot(v,headFile,'electrodes','on','style','map','shading','interp');colorbar
%topoplot(v,headFile,'electrodes','on','style','map','shading','interp','maplimits',[-2 2]);
