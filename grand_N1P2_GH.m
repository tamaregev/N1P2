%grand_N1P2

%% definitions
Definitions_N1P2
Subjects = 1:37;

%for grand -
grandFolder = [AnalysisFolder 'grandAverage'];
srate = 512;
Expinfo.srate = srate;

%codes and names table:
Bnames = {'1','2a','2b','3a','3b'}';
Bcodes = [10,20,30,40,50]';
Bcodes_names = table(Bnames,Bcodes);

%% calculate and save grand matrix - cond + prev cond
% function parameters

whichSubjects = Subjects(~ismember(Subjects,badSubjects));
%whichSubjects=longSubjects;

winbegin = -0.1; winend = 0.5; 
win = round(winbegin*srate:winend*srate);
bsl = [];
nelectrodes = 74;
tag = 'Bp1-20';
Bp = [1 20];
blocks = {'1','2a','2b','3a','3b'};

tictot = tic;
for bl = 1:5
    ticblock = tic;
    
    blockName = blocks{bl};
    disp(blockName)
    Bcode = Bcodes_names.Bcodes(strcmp(Bcodes_names.Bnames,blockName));

    conditions = {'note 1','note 2','note 3','note 4','note 5'};
    prevConditions = {'note 1','note 2','note 3','note 4','note 5'};
    STIMcodes = [1,2,3,4,5];
    EventCodes=cell(1,length(conditions));
    for i=1:length(conditions)
        EventCodes{i} = Bcode + STIMcodes(i);
    end

    prevEventCodes = EventCodes;

    % run function and save
    [ grandMatrix, ntrials ] = calcGrandMatrix_condPrev_N1P2(Expinfo, whichSubjects, conditions, prevConditions, EventCodes, prevEventCodes, nelectrodes, win, bsl, Bp  );
    disp('saving...'); matrixdate = date;
    save([grandFolder filesep 'grandMatrix_condComb_' tag '_Bl' blockName '_' date],'grandMatrix','ntrials','conditions','prevConditions','EventCodes','whichSubjects','ntrials','win','bsl','winbegin','winend','matrixdate');
    clear grandMatrix ntrials
    disp(['done ' blocks{bl} ' in ' num2str(toc(ticblock)) ' sec.'])
end
disp(['done all in ' num2str(toc(tictot)) ' sec.'])

%
%% load and plot grand 5 conds -
matrixdate='29-Oct-2018';
tag = 'Bp1-20';
blocks = {'1','2a','2b','3a','3b'};
bslwin = [-0.1,0];
electrodeName = 'Cz';
[ electrode ] = chName2n( electrodeName );

colors = {[0 0 1],[0 0.5 1],[0 1 0.5],[0 0.5 0],[1 0 0]};
fromx=-100;tox=500;fromy=-4;toy=4;
h=ERPfigure;
set(h,'Position',[20 1 1000 800])

for bl = 1:5
    switch bl
        case 1
            subplot(3,2,1)
        otherwise
            subplot(3,2,bl+1)
    end
    load([grandFolder filesep 'grandMatrix_condComb_' tag '_Bl' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
     %plot
    for con=1:size(grandMatrix,3)
        varplot(t,nanmean(grandMatrix(:,whichSubjects,con,:,electrode),4),'Color',colors{con},'linew',3)
        hold on
    end  
    title(blocks{bl})
    ax(bl)=gca;
    line([0,0],[fromy, toy],'linestyle','-.','col','k');
    line([fromx, tox],[0, 0],'linestyle','-.','col','k');
end
legend({'note 1','note 2','note 3','note 4','note 5'})
suptitle(['N = ' num2str(length(whichSubjects)) '. Electrode ' electrodeName ])
linkaxes(ax)
axis([fromx, tox, fromy, toy])
%% load and plot grand prevCond
matrixdate='29-Oct-2018';
tag = 'Bp1-20';
blocks = {'1','2a','2b','3a','3b'};
bslwin = [-0.1,0];
electrodeName = 'Cz';
[ electrode ] = chName2n( electrodeName );
ch=electrode;

if 1 %params
    fos = 18;
    fromx=-0.1;tox=0.4;fromy=-5;toy=5;
    cs = [1 2 3 4 5];
    notes = {'tone 1','tone 2','tone 3','tone 4','tone 5'};
end

for bl = 1:5
    h=ERPfigure;
    set(h,'Position',[80,80,900,600]);
    load([grandFolder filesep 'grandMatrix_condComb_' tag '_Bl' blocks{bl} '_' matrixdate])
    %baseline
    newbsl = floor((bslwin(1)-winbegin)*srate)+1:floor((bslwin(2)-winbegin)*srate);
    grandMatrix = bsxfun(@minus,grandMatrix,mean(grandMatrix(newbsl,:,:,:,:),1));
    t=0:1/srate:(size(grandMatrix,1)-1)/srate;t=t+winbegin;t=t*1000;
     
    ncs=length(cs);
    i=0;
    for n=cs
        for pn=cs
            i=i+1;
            data = grandMatrix(:,whichSubjects,n,pn,ch);
            subplot(ncs,ncs,i)
            t=0:1:size(data,1)-1;t=t/srate;t=t+winbegin;
            %varplot(t,data,'LineWidth',2)
            plot(t,nanmean(data,2),'LineWidth',3)
            axis([fromx, tox, fromy, toy])
            %set(gca,'fontsize',fos)
            line([0,0],[fromy, toy],'linestyle','-.','col','k');
            line([fromx,tox],[0, 0],'linestyle','-.','col','k');
            pos=get(gca,'position');
            pos(3)=pos(3)*1.2;
            pos(4)=pos(4)*1.1;
            set(gca,'position',pos);
            if i==25
                xlabel('seconds');ylabel('\muV')
            else
                set(gca,'xticklabels',{});
                set(gca,'yticklabels',{})
            end
            if n==1
                title(['Previous note: ' notes(pn)])
            end
            if pn==1
                ylabel(['Current note: ' notes(n)])
            end
            if n~=pn
            %    text(0.1,2,['mean Ntr = ' num2str(Ntr(n,pn),3)])
            end
        end
    end
    linkaxes
    suptitle(['Block ' num2str(bl)])
end

%% for peak detection Go to -
% Amplitudes_N1P2.m
