addpath('S:\Lab-Shared\NewDataArch\CommonResources\Tools\Matlab_Tools')
addpath('S:\Lab-Shared\Experiments\MMNchroma\Analysis')

load('S:\Lab-Shared\Experiments\N1P2\Analysis\MixedModel\theTable')
            
theTable.subject=categorical(theTable.subject);
T=theTable;
T.range=T.spread*4;

formulas = {'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + dist_mean:size_jump:isN1 + dist_mean:size_jump:isP2 + (isN1-1|subject) + (isP2-1|subject) + (dist_mean:isN1-1|subject) + (dist_mean:isP2-1|subject) + (size_jump:isN1-1|subject) + (size_jump:isP2-1|subject) + (size_jump:dist_mean:isN1-1|subject) + (size_jump:dist_mean:isP2-1|subject)',...%1
            'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + dist_mean:size_jump:isN1 + dist_mean:size_jump:isP2 + (isN1-1|subject) + (isP2-1|subject) + (dist_mean:isN1-1|subject) + (dist_mean:isP2-1|subject) + (size_jump:isN1-1|subject) + (size_jump:isP2-1|subject)',...%2
            'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + (isN1-1|subject) + (isP2-1|subject) + (dist_mean:isN1-1|subject) + (dist_mean:isP2-1|subject) + (size_jump:isN1-1|subject) + (size_jump:isP2-1|subject)',...%3
            };

lmes = cell(length(formulas),1);
for i=1:length(formulas)
    disp(i)
    %lmes{i} = fitlme(theTable(T.ExpN==categorical(3),:),formulas{i});
    lmes{i} = fitlme(T,formulas{i});
end
SaveResultsFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\MixedModel\ResultsFiles';
chooses = [2 3];
for ich=chooses
    lme=lmes{ich};
    char = evalc('disp(lme)');
    fid = fopen([SaveResultsFolder filesep 'allExplmeResults_' num2str(ich)],'wt');
    fprintf(fid,'%s',char); fclose(fid);
end
choose = 2;
compare(lmes{3},lmes{2})
rFormula='Voltage ~ dist_mean + size_jump + dist_mean:size_jump + (1|subject) + (dist_mean-1|subject) + (size_jump-1|subject)';
N1lme = fitlme(T,rFormula,'Exclude',T.isN1==0);
P2lme = fitlme(T,rFormula,'Exclude',T.isN1==1);
%% lme estimates
colorsnp={[0.9608 0.5098 0.1882],[0.2353 0.7059 0.2941]};
Colors={{colorsnp{1},colorsnp{2}},{colorsnp{1},colorsnp{2},colorsnp{1},colorsnp{2},colorsnp{1},colorsnp{2}}};


choose=3;
lme=lmes{choose};
subjs=unique(T.subject);
Names = lme.Coefficients.Name;

Estimates = nan(length(Names),length(subjs));
beta = lme.Coefficients.Estimate;
[re, ne]=lme.randomEffects;
mapne = cell(length(Names),1);
for i=1:size(ne,1)/length(subjs)
    switch i
        case {1,2}
            disp(i)
            mapne{i} = i:2:length(subjs)*2;
            Namesrand{i} = Names{i};
        otherwise
            idx = length(subjs)*(i-1)+1:length(subjs)*i;
            disp(i)
            disp(length(idx))
            Namesrand{i} = ne.Name{idx(1),:}; 
            ind = find(strcmp(Namesrand{i},Names));
            mapne{ind} = idx;
    end
end
%check and assign into Estimates:
for i=1:length(Names)
    if ~isempty(mapne{i})
        disp(['Fixed ' Names{i} ' random ' ne.Name{mapne{i}(1)}]);disp('adding random')
        Estimates(i,:) = beta(i)+re(mapne{i});
    else
        disp(['Fixed ' Names{i}]);disp('no random')
        Estimates(i,:) = beta(i);
    end
end

lmeEstimates = Estimates;

%Effect sizes:
ds=lme.Coefficients.tStat/sqrt(78);
Names = lme.CoefficientNames;
for i=1:length(Names)
    disp([Names{i} ', Estimate = ' num2str(beta(i)) ' d = ' num2str(ds(i)) ' p = ' num2str(lme.Coefficients.pValue(i)) ])
end

%% plot
hf=ERPfigure;
set(hf,'Position',[100 100 400 400])
for in=1:length(Names)
    Namesr{in}=strrep(Names{in},'_',' ');
    Namesr{in}=strrep(Namesr{in},'isN1','N1');
    Namesr{in}=strrep(Namesr{in},'isP2','P2');
end
Namesr{1} = 'Intercept N1';Namesr{2}='Intercept P2';
splocs={[1:2 9:10],[4:8,12:16]};
iNames = {{'Intercept N1','Intercept P2'},{'dist mean:N1','dist mean:P2','size jump:N1','size jump:P2','dist mean:size jump:N1','dist mean:size jump:P2'}};
iest=cell(size(iNames));
for i=1:length(iest)
    iest{i} = nan(size(iNames{i}));
    for ii=1:length(iest{i})
        iest{i}(ii) = find(strcmp(Namesr,iNames{i}{ii}));
    end
end
xlocs={[1 2],[1 2 4 5 7 8]};
for ii=1:2
    subplot(3,8,splocs{ii})
    for ie=iest{ii}
        hold on
        disp([Names{ie}])
        est=Estimates(ie,:);
        x=xlocs{ii}(ie==iest{ii});
        b=bar(x,mean(est),'FaceColor',Colors{ii}{ie==iest{ii}});
        %CI=Confidence(est);err=CI(2)-mean(est);
        tin=tinv([0.025 0.975],78);
        err=lme.Coefficients.SE(ie)*tin(2);
        if ~isempty(mapne{ie})
            errorbar(x,mean(est),err,'.k','linewidth',2)
        end
    end
    set(gca,'xtick',xlocs{ii},'xticklabels',Namesr(iest{ii}))
    xtickangle(45)
    set(gca,'fontsize',12)
    if ii==1;ylabel('Z-scored Voltage (\muV)');xlim([0 3]);ylim([-0.81 0.1]);end
    rectangle('Position',[6.5 -0.001 2 0.001])
end
xlim([0 9]);ylim([-0.014 0.051]);
suptitle('LME estimates')
if 1 %inset
    axes('Position',[.8 .6 .1 .2])
    box on
    
    iii=0;
    for ie=iest{ii}(end-1:end)
        iii=iii+1;
        hold on
        disp([Names{ie}])
        est=Estimates(ie,:);
        b=bar(iii,mean(est),'FaceColor',Colors{ii}{ie==iest{ii}});
     %   CI=Confidence(est);err=CI(2)-mean(est);
       tin=tinv([0.025 0.975],78);
        err=lme.Coefficients.SE(ie)*tin(2);
     
        %if ~isempty(mapne{ie})
            errorbar(iii,mean(est),err,'.k','linewidth',2)
        %end
    end
    set(gca,'xtick',[]);xlim([0 3]);ylim([-0.002 0.001])
    xtickangle(45)
    set(gca,'fontsize',12)
end
saveFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\PaperFigures\Regression';
saveas(gcf,[saveFolder filesep 'LMEestimates_' num2str(choose)],'pdf')     
saveas(gcf,[saveFolder filesep 'LMEestimates'],'jpg')     

%% contrasts:
addpath('S:\Lab-Shared\Z backup\Tamar\fromZ\Documents\MATLAB\MatlabFunctions\downloaded\Violinplot-Matlab-master')
%plot

pairNames = {{'dist mean:N1','size jump:N1'};{'dist mean:P2','size jump:P2'};{'dist mean:N1','dist mean:P2'};{'size jump:N1','size jump:P2'};{'dist mean:size jump:N1','dist mean:size jump:P2'}};
pairs=nan(size(pairNames,1),2);    
for i=1:size(pairNames,1)
    for ii=1:2
        pairs(i,ii) = find(strcmp(Namesr,pairNames{i}{ii}));
    end
end
h=ERPfigure;
nsp=5;
set(h,'Position',[10 -500 800 700])
for estP=1:size(pairs,1)-1
    disp([Namesr{pairs(estP,1)} ' vs. ' Namesr{pairs(estP,2)}])
    est1=Estimates(pairs(estP,1),:);
    est2=Estimates(pairs(estP,2),:);
    dest = est1-est2;
    disp(['Mean = ' num2str(mean(dest))])
    contvec = zeros(size(Names));
    contvec(pairs(estP,1))=1;contvec(pairs(estP,2))=-1;        
    [p, f, df1, df2]=coefTest(lme,contvec);d=sqrt(f)/sqrt(78);disp(['F= ' num2str(f) ' df1= ' num2str(df1) ' df2= ' num2str(df2) ' p= ' num2str(p) ' d= ' num2str(d)])
    disp(' ')
    subplot(1,nsp,estP)
    hold on     
    v(estP)=violinplot(dest',1,'BoxColor',[0 0 0],'ShowMean',true,'ViolinAlpha',0.8,'ViolinColor',[0.9 0.9 0.9],'ShowData',false);

%    disp(['MeanPlot = ' num2str(v(estP).MeanPlot.YData(1))])
 %   disp(' ')
    
    line([0.5 1.5],[0 0],'linestyle','-.','Color','k','linewidth',1)
    set(gca,'xtick',1,'xticklabels',[strrep(Names{pairs(estP,1)},'_',' ') 'vs.' strrep(Names{pairs(estP,2)},'_',' ')])
    set(gca,'fontsize',12)
    ylim([-0.08 0.1])
    if estP~=1;set(gca,'Yticklabels',{});end
end
disp('diff of diffs')
switch choose
    case 2
        contvec = [0 0 1 -1 -1 1 0 0];
    case 3
        contvec = [0 0 1 -1 -1 1];
end
    subplot(1,nsp,5)
    [p, f, df1, df2]=coefTest(lme,contvec);d=sqrt(f)/sqrt(78);disp(['F= ' num2str(f) ' df1= ' num2str(df1) ' df2= ' num2str(df2) ' p= ' num2str(p) ' d= ' num2str(d)])
    ddest = Estimates(3,:)-Estimates(4,:)-Estimates(5,:)+Estimates(6,:);
    hold on
    v(nsp)=violinplot(ddest',1,'BoxColor',[0 0 0],'ShowMean',true,'ViolinAlpha',0.8,'ViolinColor',[0.9 0.9 0.9],'ShowData',false);
    line([0.5 1.5],[0 0],'linestyle','-.','Color','k','linewidth',1)
     set(gca,'xtick',1,'xticklabels',['diff N1 ' 'vs.' ' diff P2'])
     set(gca,'yticklabels',{})
   % xtickangle(30)
    set(gca,'fontsize',12)
%     xlim([0.9 1.1]);
    ylim([-0.08 0.1])
    
    for vi=1:nsp
        v(vi).ShowData=true;
        v(vi).MeanPlot.Color=[0 0 0];
        v(vi).ScatterPlot.MarkerFaceColor=[0.7 0.7 0.7];
        v(vi).ScatterPlot.SizeData=12;
    end
   
%         %compare last pair
    if 0
    %         estP = size(pairs,1);
    %         disp([Names{pairs(estP,1)} ' vs. ' Names{pairs(estP,2)}])
    %         est1 = Estimates(pairs(estP,1),:);
    %         est2 = Estimates(pairs(estP,2),:);
    %         dest = est1-est2;
    %         contvec = zeros(size(Names));
    %         contvec(pairs(estP,1)) = 1;contvec(pairs(estP,2))=-1;        
    %         [p, f] = coefTest(lme,contvec);d=sqrt(f)/sqrt(78);disp(['p= ' num2str(p) ' d = ' num2str(d)])
    %         disp(' ')
    %         subplot(1,6,6)
    %         hold on     
    % %         hb=boxplot(dest,'Notch','on');
    % %         set(hb,'linewidth',2)
    % %         plot(ones(size(dest))+randn(size(dest))/50,dest,'.g')
    %         violinplot(dest)
    %         %bar(mean(dest))
    %         %[h,L,MX,MED,bw]=violin(dest');
    %         line([0.5 1.5],[0 0],'linestyle','-.','Color','k','linewidth',1)
    %      %   CI=Confidence(dest);errs=CI(2)-mean(dest);
    %      %   errorbar(mean(dest),errs,'.m','linewidth',2)
    %         set(gca,'xtick',1,'xticklabels',[strrep(Names{pairs(estP,1)},'_',' ') 'vs.' strrep(Names{pairs(estP,2)},'_',' ')])
    %        % xtickangle(30)
    %         set(gca,'fontsize',12)
    %         xlim([0.9 1.1]);
    %         ylim([-0.0016 0.002])
    end
   saveas(gcf,[saveFolder filesep 'compareLMEestimates'],'jpg')
  saveas(gcf,[saveFolder filesep 'compareLMEestimates'],'pdf')

  %% fit lm for each participant separately
%formula = 'Voltage ~ Potential + dist_mean:isN1 + dist_mean:isP2';
subjs=unique(T.subject);
formula = 'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + dist_mean:size_jump:isN1 + dist_mean:size_jump:isP2';
numEst=8;
lm=cell(length(subjs),1);
Estimates=nan(numEst,length(subjs));SEs=nan(numEst,length(subjs));tStats=nan(numEst,length(subjs));
for s=subjs'
    lm{s}=fitlm(T(T.subject==s,:),formula);
    Estimates(:,s)=lm{s}.Coefficients.Estimate;
    SEs(:,s)=lm{s}.Coefficients.SE;
    tStats(:,s)=lm{s}.Coefficients.tStat;
end
    %% ttest for estimates
    Names=lm{s}.CoefficientNames;
    for estN=1:numEst
        disp(Names{estN})
        [H,p,ci,stats] = ttest(Estimates(estN,:));disp(['Estimate ' num2str(num2str(mean(Estimates(estN,:)))) ' p=' num2str(p) ' tstat=' num2str(stats.tstat) ' se=' num2str(stats.sd/sqrt(78)) ' df=' num2str(stats.df)])
        d=stats.tstat/sqrt(stats.df);disp(['d=' num2str(d)])
        disp(' ')
    end
    %% compare to lmeEstimates
    hf=ERPfigure;
    set(hf,'position',[100 100 400 600])
    for i=1:length(Names)-2
        subplot(3,2,i)
        plot(Estimates(i,:),lmeEstimates(i,:),'.')
        title(Namesr{i})
        hold on
        xdataX = [min(Estimates(i,:)) max(Estimates(i,:))];
        xdataY = [min(lmeEstimates(i,:)) max(lmeEstimates(i,:))];
        line([mean(Estimates(i,:)) mean(Estimates(i,:))],[min([xdataX xdataY]),max([xdataX xdataY])],'Color',[0.7 0.7 0.7])
        line([min([xdataX xdataY]),max([xdataX xdataY])],[mean(lmeEstimates(i,:)) mean(lmeEstimates(i,:))],'Color',[0.7 0.7 0.7])
        if xdataY(1)~=xdataY(2)
            line([min([xdataX xdataY]) max([xdataX xdataY])],[min([xdataX xdataY]) max([xdataX xdataY])])
            d=0;
%            xlim([xdataX(1)-d xdataX(2)+d]);ylim([xdataY(1)-d xdataY(2)+d])
            xlim([min([xdataX xdataY]),max([xdataX xdataY])]);
            ylim([min([xdataX xdataY]),max([xdataX xdataY])]);

        else
            d=0.0001;
            line([mean(Estimates(i,:)) mean(Estimates(i,:))],[xdataY(1)-d xdataY(2)+d],'Color',[0.7 0.7 0.7])
            xlim([xdataX(1)-d xdataX(2)+d]);ylim([xdataY(1)-d xdataY(2)+d])
        end
        set(gca,'fontsize',16)
        if i==5;xlabel('Regressions');ylabel('LME');end
    end
    saveFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\PaperFigures\Regression';
    saveas(gcf,[saveFolder filesep 'compareRegLME'],'pdf')
    saveas(gcf,[saveFolder filesep 'compareRegLME'],'jpg')  
%% plot regression estimates
colorsnp={[0.9608 0.5098 0.1882],[0.2353 0.7059 0.2941]};
Colors={{colorsnp{1},colorsnp{2}},{colorsnp{1},colorsnp{2},colorsnp{1},colorsnp{2},colorsnp{1},colorsnp{2}}};

% plot
hf=ERPfigure;
set(hf,'Position',[100 100 400 400])
for in=1:length(Names)
    Namesr{in}=strrep(Names{in},'_',' ');
    Namesr{in}=strrep(Namesr{in},'isN1','N1');
    Namesr{in}=strrep(Namesr{in},'isP2','P2');
end
Namesr{1} = 'Intercept N1';Namesr{2}='Intercept P2';
splocs={[1:2 9:10],[4:8,12:16]};
iNames = {{'Intercept N1','Intercept P2'},{'dist mean:N1','dist mean:P2','size jump:N1','size jump:P2','dist mean:size jump:N1','dist mean:size jump:P2'}};
iest=cell(size(iNames));
for i=1:length(iest)
    iest{i} = nan(size(iNames{i}));
    for ii=1:length(iest{i})
        iest{i}(ii) = find(strcmp(Namesr,iNames{i}{ii}));
    end
end
xlocs={[1 2],[1 2 4 5 7 8]};
for ii=1:2
    subplot(3,8,splocs{ii})
    for ie=iest{ii}
        hold on
        disp([Names{ie}])
        est=Estimates(ie,:);
        x=xlocs{ii}(ie==iest{ii});
        b=bar(x,mean(est),'FaceColor',Colors{ii}{ie==iest{ii}});
        CI=Confidence(est);err=CI(2)-mean(est);
        %tin=tinv([0.025 0.975],78);
        %err=lme.Coefficients.SE(ie)*tin(2);
        if ~isempty(mapne{ie})
            errorbar(x,mean(est),err,'.k','linewidth',2)
        end
    end
    set(gca,'xtick',xlocs{ii},'xticklabels',Namesr(iest{ii}))
    xtickangle(45)
    set(gca,'fontsize',12)
    if ii==1;ylabel('Z-scored Voltage (\muV)');xlim([0 3]);ylim([-0.81 0.1]);end
    rectangle('Position',[6.5 -0.001 2 0.001])
end
xlim([0 9]);ylim([-0.014 0.051]);
suptitle('LME estimates')
if 1 %inset
    axes('Position',[.8 .6 .1 .2])
    box on
    
    iii=0;
    for ie=iest{ii}(end-1:end)
        iii=iii+1;
        hold on
        disp([Names{ie}])
        est=Estimates(ie,:);
        b=bar(iii,mean(est),'FaceColor',Colors{ii}{ie==iest{ii}});
       CI=Confidence(est);err=CI(2)-mean(est);
       %tin=tinv([0.025 0.975],78);
        %err=lme.Coefficients.SE(ie)*tin(2);
     
        %if ~isempty(mapne{ie})
            errorbar(iii,mean(est),err,'.k','linewidth',2)
        %end
    end
    set(gca,'xtick',[]);xlim([0 3]);ylim([-0.002 0.001])
    xtickangle(45)
    set(gca,'fontsize',12)
end
suptitle('Regression estimates')

saveFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\Figures\PaperFigures\Regression';
saveas(gcf,[saveFolder filesep 'REGestimates'],'pdf')
saveas(gcf,[saveFolder filesep 'REGestimates'],'jpg')
       
%% contrasts Regression :
addpath('S:\Lab-Shared\Z backup\Tamar\fromZ\Documents\MATLAB\MatlabFunctions\downloaded\Violinplot-Matlab-master')
%plot
pairNames = {{'dist mean:N1','size jump:N1'};{'dist mean:P2','size jump:P2'};{'dist mean:N1','dist mean:P2'};{'size jump:N1','size jump:P2'};{'dist mean:size jump:N1','dist mean:size jump:P2'}};
pairs=nan(size(pairNames,1),2);    
for i=1:size(pairNames,1)
    for ii=1:2
        pairs(i,ii) = find(strcmp(Namesr,pairNames{i}{ii}));
    end
end
h=ERPfigure;
set(h,'Position',[10 -500 800 700])
nsp=5;
for estP=1:size(pairs,1)-1
    disp([Namesr{pairs(estP,1)} ' vs. ' Namesr{pairs(estP,2)}])
    est1=Estimates(pairs(estP,1),:);
    est2=Estimates(pairs(estP,2),:);
    dest = est1-est2;
        [H,p,ci,stats] = ttest(dest);disp(['p=' num2str(p) ' tstat=' num2str(stats.tstat) ' df=' num2str(stats.df) ' d=' num2str(mean(dest)/std(dest)) ])
    disp(' ')
    subplot(1,nsp,estP)
    hold on     
%     hb=boxplot(dest,'Notch','on');
%     set(hb,'linewidth',2)
%     plot(ones(size(dest))+randn(size(dest))/50,dest,'.g')
%    
  %bar(mean(dest))
    v(estP)=violinplot(dest',1,'BoxColor',[0 0 0],'ShowMean',true,'ViolinAlpha',0.8,'ViolinColor',[0.9 0.9 0.9],'ShowData',false);
       
    line([0.5 1.5],[0 0],'linestyle','-.','Color','k','linewidth',1)
    %CI=Confidence(dest);errs=CI(2)-mean(dest);
    %errorbar(mean(dest),errs,'.m','linewidth',2)
    set(gca,'xtick',1,'xticklabels',[strrep(Names{pairs(estP,1)},'_',' ') 'vs.' strrep(Names{pairs(estP,2)},'_',' ')])
   % xtickangle(30)
    set(gca,'fontsize',12)
    %xlim([0.9 1.1]);
    ylim([-0.25 0.25])
    if estP~=1;set(gca,'Yticklabels',{});end
    
end
disp('diff of diffs')

    subplot(1,nsp,5)
    ddest = Estimates(3,:)-Estimates(4,:)-Estimates(5,:)+Estimates(6,:);
    [H,p,ci,stats] = ttest(ddest);disp(['p=' num2str(p) ' tstat=' num2str(stats.tstat) ' df=' num2str(stats.df) ' d=' num2str(mean(ddest)/std(ddest))])

    hold on
%     hb=boxplot(ddest,'Notch','on');
%     set(hb,'linewidth',2)
%     plot(ones(size(ddest))+randn(size(ddest))/50,ddest,'.g')
    v(nsp)=violinplot(ddest',1,'BoxColor',[0 0 0],'ShowMean',true,'ViolinAlpha',0.8,'ViolinColor',[0.9 0.9 0.9],'ShowData',false);
    
    for vi=1:nsp
        v(vi).ShowData=true;
        v(vi).MeanPlot.Color=[0 0 0];
        v(vi).ScatterPlot.MarkerFaceColor=[0.7 0.7 0.7];
        v(vi).ScatterPlot.SizeData=12;
    end
    line([0.5 1.5],[0 0],'linestyle','-.','Color','k','linewidth',1)
    %CI=Confidence(ddest);errs=CI(2)-mean(ddest);
    %errorbar(mean(ddest),errs,'.m','linewidth',2)
     set(gca,'xtick',1,'xticklabels',['diff N1 ' 'vs.' ' diff P2'])
     set(gca,'yticklabels',{})
     
   % xtickangle(30)
    set(gca,'fontsize',12)
%     xlim([0.9 1.1]);
    ylim([-0.25 0.25])
%         %compare last pair
    if 0
    %         estP = size(pairs,1);
    %         disp([Names{pairs(estP,1)} ' vs. ' Names{pairs(estP,2)}])
    %         est1 = Estimates(pairs(estP,1),:);
    %         est2 = Estimates(pairs(estP,2),:);
    %         dest = est1-est2;
    %         contvec = zeros(size(Names));
    %         contvec(pairs(estP,1)) = 1;contvec(pairs(estP,2))=-1;        
    %         [p, f] = coefTest(lme,contvec);d=sqrt(f)/sqrt(78);disp(['p= ' num2str(p) ' d = ' num2str(d)])
    %         disp(' ')
    %         subplot(1,6,6)
    %         hold on     
    % %         hb=boxplot(dest,'Notch','on');
    % %         set(hb,'linewidth',2)
    % %         plot(ones(size(dest))+randn(size(dest))/50,dest,'.g')
    %         violinplot(dest)
    %         %bar(mean(dest))
    %         %[h,L,MX,MED,bw]=violin(dest');
    %         line([0.5 1.5],[0 0],'linestyle','-.','Color','k','linewidth',1)
    %      %   CI=Confidence(dest);errs=CI(2)-mean(dest);
    %      %   errorbar(mean(dest),errs,'.m','linewidth',2)
    %         set(gca,'xtick',1,'xticklabels',[strrep(Names{pairs(estP,1)},'_',' ') 'vs.' strrep(Names{pairs(estP,2)},'_',' ')])
    %        % xtickangle(30)
    %         set(gca,'fontsize',12)
    %         xlim([0.9 1.1]);
    %         ylim([-0.0016 0.002])
    end
   saveas(gcf,[saveFolder filesep 'compareRegEstimates'],'jpg')
  saveas(gcf,[saveFolder filesep 'compareRegEstimates'],'pdf')
%% radomize T and run on part of data to test effect size:

rT=T(randperm(height(T)),:);
fractions=[0.01,0.05,0.09,0.1,0.25,0.3,0.5,0.7,0.85,1];
lmer = cell(length(formulas),1);
ii=6;
clear d
for i=1:length(fractions)
    disp(i)
    %lmes{i} = fitlme(T(T.ExpN==categorical(3),:),formulas{i});
    lmer{i} = fitlme(rT(1:fractions(i)*height(rT),:),formulas{ii});
    %Effect size:
    d{i}=lmer{i}.Coefficients.tStat/sqrt(lmer{i}.Coefficients.DF(1));
end

%% average across previous tones:

whichpeaks = {'N1','P2'};
nT = T(1,[1,2,4,6,8:14]);
%loop over participants:
subjs=unique(T.subject);
idx=0;%line in nT
for is=subjs'
    sT=T(T.subject==is,:);
    bNs = unique(sT.blockN);
    for ib = bNs'
        sbT=sT(sT.blockN==ib,:);
        tones = unique(sbT.currMIDI);
        for it = tones'
            for ipeak = 1:2
                idx=idx+1;
                Pot=whichpeaks{ipeak};
                ii=find(sbT.currMIDI==it & sbT.Potential==Pot);
                nT(idx,:)=sbT(ii(1),[1,2,4,6,8:14]);
                nT.Voltage(idx) = mean(sbT.Voltage(ii));
            end
        end
    end
end

%% fit lme
i=6;
formulas{i};
subjs=unique(T.subject);

formula = 'Voltage ~ Potential + dist_mean:isN1 + dist_mean:isP2 + (Potential|subject) + (Potential|ExpN)';
lme = fitlme(nT,formula);
%Effect size:
d=lme.Coefficients.tStat/sqrt(lme.Coefficients.DF(1)); disp(d)
%contrasts:
[p, f ,df1, df2]=coefTest(lme,[0 0 1 -1]);d=sqrt(f)/sqrt(df2);disp(p);disp(d)
%% fit lm for all subjects together 
formula = 'Voltage ~ Potential*subject + dist_mean:isN1 + dist_mean:isP2';
lm=fitlm(nT,formula);

%% for only Exp 3
nSubj=30;
formula =       'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + dist_mean:size_jump:isN1 + dist_mean:size_jump:isP2 + (isN1-1|subject) + (isP2-1|subject) + (dist_mean:isN1-1|subject) + (dist_mean:isP2-1|subject) + (size_jump:isN1-1|subject) + (size_jump:isP2-1|subject)';%8 

formulasExp3 = {'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + dist_mean:size_jump:isN1 + dist_mean:size_jump:isP2 + (isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + dist_mean:size_jump:isN1 + dist_mean:size_jump:isP2):range + (isN1-1|subject) + (isP2-1|subject) + (dist_mean:isN1-1|subject) + (dist_mean:isP2-1|subject) + (size_jump:isN1-1|subject) + (size_jump:isP2-1|subject) + (dist_mean:size_jump:isN1-1|subject) + (dist_mean:size_jump:isP2-1|subject) + ((isN1:range-1)|subject) + ((isP2:range-1)|subject) + ((dist_mean:isN1:range-1)|subject) + ((dist_mean:isP2:range-1)|subject) + ((size_jump:isN1:range-1)|subject) + ((size_jump:isP2:range-1)|subject) + ((dist_mean:size_jump:isN1:range-1)|subject) + ((dist_mean:size_jump:isP2:range-1)|subject)',...
                'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + dist_mean:size_jump:isN1 + dist_mean:size_jump:isP2 + (isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2):range + (isN1-1|subject) + (isP2-1|subject) + (dist_mean:isN1-1|subject) + (dist_mean:isP2-1|subject) + (size_jump:isN1-1|subject) + (size_jump:isP2-1|subject)',...%without triple interactions of mean jump and range, without random factors of range interactions
                'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + (isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2):range + (isN1-1|subject) + (isP2-1|subject) + (dist_mean:isN1-1|subject) + (dist_mean:isP2-1|subject) + (size_jump:isN1-1|subject) + (size_jump:isP2-1|subject)',...%without double interactions of dist_mean, size_jump
                };

lmess = cell(length(formulasExp3),1);
for i=1:length(formulasExp3)
    disp(i)
    lmess{i} = fitlme(T,formulasExp3{i},'Exclude',T.ExpN~=3);
    lme=lmess{i};
    char = evalc('disp(lme)');
    fid = fopen([SaveResultsFolder filesep 'Exp3lmeResults_f' num2str(i)],'wt');
    fprintf(fid,'%s',char); fclose(fid);
end
SaveResultsFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\MixedModel\ResultsFiles';
choose = 3;

stats=compare(lmess{3},lmess{1});
%% comaprisons
choose =2;
lme=lmess{choose};
Names = lme.CoefficientNames;
Estimates = lme.Coefficients.Estimate;
clear Namesr
for in=1:length(Names)
    Namesr{in}=strrep(Names{in},'_',' ');
    Namesr{in}=strrep(Namesr{in},'isN1','N1');
    Namesr{in}=strrep(Namesr{in},'isP2','P2');
end

pairNames = {{'N1:range','P2:range'};{'dist mean:N1:range','dist mean:P2:range'};{'size jump:N1:range','size jump:P2:range'};{'dist mean:N1:range','size jump:N1:range'};{'dist mean:P2:range','size jump:P2:range'}};
%pairNames = {{'spread:N1','spread:P2'};{'dist mean:spread:N1','dist mean:spread:P2'};{'size jump:spread:N1','size jump:spread:P2'}};

pairs=nan(size(pairNames,1),2);    
for i=1:size(pairNames,1)
    for ii=1:2
        pairs(i,ii) = find(strcmp(Namesr,pairNames{i}{ii}));
    end
end
for estP=1:size(pairs,1)
    disp([Namesr{pairs(estP,1)} ' vs. ' Namesr{pairs(estP,2)}])
    est1=Estimates(pairs(estP,1),:);
    est2=Estimates(pairs(estP,2),:);
    dest = est1-est2;
    contvec = zeros(size(Names));
    contvec(pairs(estP,1))=1;contvec(pairs(estP,2))=-1;        
    [p, f, df1, df2]=coefTest(lme,contvec);d=sqrt(f)/sqrt(nSubj);disp(['F(' num2str(df1) ',' num2str(df2) ')= ' num2str(f) ' p= ' num2str(p) ' d= ' num2str(d)])
end

%diff of diffs
disp('diff of diffs')
switch choose
    case 2
        contvec = [0 0 0 0 0 0 0 0 0 0 1 -1 -1 1];
    case 3
         contvec = [0 0 0 0 0 0 0 0 1 -1 -1 1];
end
 [p, f, df1, df2]=coefTest(lme,contvec);d=sqrt(f)/sqrt(nSubj);disp(['F(' num2str(df1) ',' num2str(df2) ')= ' num2str(f) ' p= ' num2str(p) ' d= ' num2str(d)])
%% for P2 in wide condition experiment 3:
formula =       'Voltage ~ isN1-1 + isP2-1 + dist_mean:isN1 + dist_mean:isP2 + size_jump:isN1 + size_jump:isP2 + dist_mean:size_jump:isN1 + dist_mean:size_jump:isP2 + (isN1-1|subject) + (isP2-1|subject) + (dist_mean:isN1-1|subject) + (dist_mean:isP2-1|subject) + (size_jump:isN1-1|subject) + (size_jump:isP2-1|subject)';%8 

lme = fitlme(T,formula,'Exclude',T.ExpN~=3 | T.blockN~='3');
lme.Coefficients
SaveResultsFolder = 'S:\Lab-Shared\Experiments\N1P2\Analysis\MixedModel\ResultsFiles';
