
%% Section 1: generate a data table instead of the matrices
%run this only once then start from Section 2
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];

for ExpN = 1:3
    dataFolder = [Folder filesep 'Exp' num2str(ExpN) filesep];
    cd(dataFolder)
    load RA
    load singleAmps
    load seqIndcat.mat
    load isArtefact.mat
    load whichSubjects.mat
    % Partial tables, so that I can use parfor
    tic
    whichSubjects=whichSubjects(:)';
    for ip=1:size(singleAmps,4)
        parfor it=1:size(singleAmps,3)
            T{it,ip}=table;
            disp([it ip]);
            for ib=1:size(singleAmps,2)
                for is=whichSubjects
                    tmp=struct;
                    tmp.Subject=is;
                    tmp.Block=ib;
                    tmp.Trial=it;
                    tmp.Pot=ip;
                    tmp.Art=isArtefact(is,ib,it);
                    tmp.seqInd=seqIndcat(is,ib,it);
                    tmp.val=singleAmps(is,ib,it,ip);
                    T{it,ip}=[T{it,ip}; struct2table(tmp)];
                end
            end
        end
    end
    toc
    % merge all subtables together
    aT=table;
    for ip=1:size(singleAmps,4)
        for it=1:size(singleAmps,3)
            aT=[aT; T{it,ip}];
        end
    end
    save dataTable T aT
    clear T aT
end
        %% merge data structures

%bT: the big table
% Adding a variable to denote experiment number
% coding subject number by adding 100*ExpN
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
dataFolder = [Folder filesep 'allExp'];
mkdir(dataFolder)
bT = table;
for ExpN = 1:3
    cd([Folder filesep 'Exp' num2str(ExpN) filesep])
    load('dataTable','aT')
    ExpNs = ones(height(aT),1).*ExpN;
    aT.ExpN = ExpNs;
    aT.Subject = aT.Subject+100*ExpN;
    bT = [bT;aT];
end
bT(isnan(bT.val),:)=[];%removed NaNs from bT - because otherwise can't use
%art and seqInd as indexes
cd(dataFolder)
save dataTable bT

%pack some data structures
for ExpN = 1:3
    cd([Folder filesep 'Exp' num2str(ExpN) filesep])
    load RA
    load whichSubjects
    bRA{ExpN} = RA;
    bwhichSubjects{ExpN} = whichSubjects;
end
load taus
load sigmas
cd(dataFolder)
save taus taus
save sigmas sigmas
save([dataFolder filesep 'bRA'],'bRA','-v7.3')
save bwhichSubjects bwhichSubjects
%% Section 2: If dataTable exists already, can start here
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [ dataFolder, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );
%% Section 3: Calculate the lme model for all taus and sigmas - N1 and P2

%whichExp = 'all';%can be 1 2 3 or 'all'
whichExp=3;
modelName = 'lme';
numt = 5;
dataFolder = [Folder filesep 'allExp' filesep];
cd(dataFolder)
infoFlag=false;

for ipeak = 1:2
    ll=nan(length(sigmas),length(taus));
    if infoFlag
        infos=cell(length(sigmas),length(taus));
    end
    switch whichExp
        case {1,2,3}
            disp(['Exp' num2str(whichExp)])
            sel=bT.Pot==ipeak & bT.ExpN==whichExp;
            uRA=bRA{whichExp}(bwhichSubjects{whichExp},:,:,:,:);
        case 'all'
            disp(whichExp)
            sel=bT.Pot==ipeak;
            clear uRA
            for iexp=1:3
                uRA{iexp}=bRA{iexp}(bwhichSubjects{iexp},:,:,:,:);
            end
    end
    %clear bRA
    vv=bT.val(sel);
    art=bT.Art(sel);
    seqInd=bT.seqInd(sel);
    subj=bT.Subject(sel);
    vv=vv(~art & seqInd>numt);
    subj=subj(~art & seqInd>numt);
    csubj=categorical(subj);

    tictot = tic;
    parfor it=1:length(taus)
        [ llst,info ] = calcLLsigmas( it, taus, sigmas, whichExp, uRA, art, seqInd, numt, vv, subj, csubj, modelName, infoFlag);
        ll(:,it)=llst;
        infos(:,it)=info;
    end
    disp(['Done ' whichpeaks{ipeak} ' in ' num2str(toc(tictot))])

    
%     figure
%     clf
%     imagesc(taus,sigmas,ll-ll(1,1));
%     title([whichpeaks{ipeak} ' - loglikelihood relative to 1,1 Exp' num2str(whichExp)])
%     colorbar
%     xlabel('taus')
%     ylabel('sigmas')
%     
    save(['ll' whichpeaks{ipeak} '_Exp' num2str(whichExp) '_' modelName],'ll')
    if infoFlag
        save(['ll' whichpeaks{ipeak} '_Exp' num2str(whichExp) '_' modelName '_info'],'infos')
    end
end
%% Section 4: generate and plot the ll maps relative to extremum
%whichExp = 'all';%can be 1 2 3 or 'all'
whichExp = 'all';%can be 1 2 3 or 'all'
ca=[0 100];
modelName = 'lme';
dataFolder = [Folder filesep 'allExp' filesep];
cd(dataFolder)
%figure
rows=nan(1,2);
cols=nan(1,2);
bestaus=nan(1,2);
bestsigmas=nan(1,2);
h=ERPfigure;
set(h,'Position',[100 100 700 300])
for ipeak = 1:2
    hsp{ipeak}=subplot(1,2,ipeak);
    load(['ll' whichpeaks{ipeak} '_Exp' num2str(whichExp) '_' modelName])
    space = max(max(ll))-ll;
    imagesc(taus,sigmas,2*space)
    [row col]=find(space==min(min(space)));    
    rows(1,ipeak)=row;cols(1,ipeak)=col;
    bestaus(1,ipeak)=taus(col);bestsigmas(1,ipeak)=sigmas(row);
    hold on
    plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',1)
    plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',1)
    text(1,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',18);

    title([ whichpeaks{ipeak}])
    colorbar
    xlabel('taus')
    ylabel('sigmas')
    %caxis([0 100]); % This is the 0.01 critical value of the chi2, df=2 distribution.
    set(gca,'fontsize',16)

end
suptitle(['-2*LLR relative to max. Exp ' num2str(whichExp)])
saveas(gcf,['Exp ' num2str(whichExp) '_LLR'],'png')
saveas(gcf,['Exp ' num2str(whichExp) '_LLR'],'pdf')

% one colorscale for all
for ipeak=1:2
    caxis(hsp{ipeak},ca)
end
saveas(gcf,['Exp ' num2str(whichExp) '_LLR_caxis'],'png')
saveas(gcf,['Exp ' num2str(whichExp) '_LLR_caxis'],'pdf')
% 
% otherpeak=[2 1];
% for ipeak=1:2
%     caxis(hsp{ipeak},[0 6])
%     plot(hsp{otherpeak(ipeak)},[bestaus(1,ipeak) bestaus(1,ipeak)],[sigmas(1) sigmas(end)],'Color',[0 1 0],'linewidth',1)
%     plot(hsp{otherpeak(ipeak)},[taus(1) taus(end)],[bestsigmas(1,ipeak) bestsigmas(1,ipeak)],'Color',[0 1 0],'linewidth',1)       
%     plot(hsp{otherpeak(ipeak)},[bestaus(1,otherpeak(ipeak)) bestaus(1,otherpeak(ipeak))],[sigmas(1) sigmas(end)],'Color',[ 1 0 0],'linewidth',1)
%     plot(hsp{otherpeak(ipeak)},[taus(1) taus(end)],[bestsigmas(1,otherpeak(ipeak)) bestsigmas(1,otherpeak(ipeak))],'Color',[ 1 0 0],'linewidth',1)       
% end
% saveas(gcf,['Exp ' num2str(whichExp) '_LLR95CIchi'],'png')
% saveas(gcf,['Exp ' num2str(whichExp) '_LLR95CIchi'],'pdf')

%% Section 5: Compute distribution under the null hypothesis
 Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [~ , bT, bRA, bwhichSubjects, taus, sigmas, ~ ] = loadVarslmes( Folder );
 dataFolder = [Folder filesep 'allExp' filesep];
     
%whichExp = 'all';%can be 1 2 3 or 'all'
whichExp='all';
models ={'lme','lm'};
im=1;
modelName = models{im};
infoFlag=false;
whichpeaks = {'N1','P2'};
numt = 5;   
np=250;

scrambles=cell(np+1);
switch whichExp
    case {1,2,3}
        disp(['Exp' num2str(whichExp)])
        sel=bT.Pot==1 & bT.ExpN==whichExp;
    case 'all'
        disp(whichExp)
        sel=bT.Pot==1;       
end
ht=height(bT(sel,:));
scrambles{1}=1:ht;
for ip=2:np+1
    scrambles{ip}=randperm(ht);
end

ll=nan(np+1,length(sigmas),length(taus),2);
parfor ip=1:np+1
    disp(ip)
    ticperm = tic;%
    [ llst, ~ ] = calcLL_scrambled( ip, scrambles, taus, sigmas, whichExp, bwhichSubjects, bRA, bT, numt, modelName, infoFlag);
    ll(ip,:,:,:) = llst;
    disp(['Done in ' num2str(toc(ticperm))])
end
save([dataFolder 'testNull_' modelName],'ll','scrambles')
%% plot
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 dataFolder = [Folder filesep 'allExp' filesep];
  
modelName = models{1};
load([dataFolder 'testNull_' modelName])
llmaxs=nan(251,2);
llmins=nan(251,2);
for ip=1:251
    for ipeak=1:2
        llmaxs(ip,ipeak)=max(max(ll(ip,:,:,ipeak)));
        llmins(ip,ipeak)=min(min(ll(ip,:,:,ipeak)));
    end
end
%plot histograms of min and max null LL
insetflag=true;
hf=ERPfigure;
if insetflag
    set(hf,'Position',[100 100 600 300])
else
    set(hf,'Position',[100 100 1200 500])
end
clf
insetlims=[-1123826 -1123810;...
           -1123467 -1123447];
for ipeak =1:2
    subplot(2,1,ipeak)
    histogram(llmins(2:251,ipeak))
    hold on
    histogram(llmaxs(2:251,ipeak))
    plot([llmins(1,ipeak) llmins(1,ipeak)],[0 160],'g','linewidth',3)
    plot([llmaxs(1,ipeak) llmaxs(1,ipeak)],[0 160],'r','linewidth',3)
    ylim([0 160])
    if insetflag
        xlim(insetlims(ipeak,:))
    end
    title(whichpeaks{ipeak})
    
%     plot(llmins(2:251,ipeak),'.b')
%     hold on
%     plot(llmaxs(2:251,ipeak),'.r')
%     plot([0 250],[llmins(1,ipeak) llmins(1,ipeak)],'m')
%     plot([0 250],[llmaxs(1,ipeak) llmaxs(1,ipeak)],'g')
    set(gca,'fontsize',14)
end
if ~insetflag
    legend({'min. LL null','max. LL null','min. LL real','max. LL real'})
end
if insetflag
    saveas(gcf,[dataFolder 'nullLL_inset'],'pdf')
else
    saveas(gcf,[dataFolder 'nullLL'],'pdf')    
end

%%compare to chi square
ERPfigure
clf
theta=[7,6;...
        9, 4];%isig=1: 1 semitone itau=3: 0.6 sec. 
isp=0;
for ii=1:2
    for ipeak = 1:2
        isp=isp+1;
        subplot(2,2,isp)
        llreal=squeeze(ll(1,:,:,ipeak));
        [r, c]=find(llreal==max(max(llreal)));
        %imagesc(squeeze(ll(1,:,:,ipeak)));colorbar        
        llr=abs(-2*(ll(2:251,1,1,ipeak)-ll(2:251,theta(ii,1),theta(ii,2),ipeak)));
        x=0:1:15;
        y=chi2pdf(x,2);
        y=y/sum(y);
        h{isp}=histogram(llr,10,'normalization','pdf');
        %val = h{isp}.Values;
        %binE = h{isp}.BinEdges;
        %hold on
        %plot(binE(1:end-1),val/250,'.')
        %disp(sum(val))
        hold on
        plot(x,y,'r','linewidth',2)        
        ylabel('probability')
        xlabel('LLR value')
        if ii==1
            title(whichpeaks{ipeak})
        end
        xlim([0 15])
        legend({['\sigma = ' num2str(sigmas(theta(ii,1))) ' \tau = ' num2str(taus(theta(ii,2)))],'\chi2 2df pdf'})
        set(gca,'fontsize',14)
    end
end
suptitle(['Null D dist. relative to \chi2 with 2 df'])
saveas(gcf,[dataFolder 'NullDdist_chi'],'pdf')
%% Section 6: Permute N1 and P2 to compare their best estimates
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [~ , bT, bRA, bwhichSubjects, taus, sigmas, ~ ] = loadVarslmes( Folder );

whichExp = 'all';

modelName = 'lme';
savetag = 'allExp_zscoreFlip';
zsFlag = true;
infoFlag = false;
whichpeaks = {'N1','P2'};
otherpeak = [2 1];
numt = 5;     
np=250;

%prep data
for ipeak = 1:length(whichpeaks)
    switch whichExp
        case {1,2,3}
            disp(['Exp' num2str(whichExp)])
            sel=bT.Pot==ipeak & bT.ExpN==whichExp;
            uRA=bRA{whichExp}(bwhichSubjects{whichExp},:,:,:,:);
        case 'all'
            disp(whichExp)
            sel=bT.Pot==ipeak;
            clear uRA
            for iexp=1:3
                uRA{iexp}=bRA{iexp}(bwhichSubjects{iexp},:,:,:,:);
            end
    end
    %clear bRA
    vv=bT.val(sel);
    art=bT.Art(sel);
    seqInd=bT.seqInd(sel);
    subj=bT.Subject(sel);
    vv=vv(~art & seqInd>numt);
    subj=subj(~art & seqInd>numt);

    %zscore N1 and P2 data per participant
    %also flip N1 sign (to positive)
    if zsFlag
        if ipeak==1
            %flip
            vv=-1*vv;
        end
        switch whichExp
            case {1,2,3}
                for is=1:length(bwhichSubjects{whichExp})
                    subjN=bwhichSubjects{whichExp}(is)+100*whichExp;
                    vv(subj==subjN) = zscore(vv(subj==subjN));    
                end
            case 'all'
                for iexp=1:3
                    for is=1:length(bwhichSubjects{iexp})                
                        subjN=bwhichSubjects{iexp}(is)+100*iexp;
                        vv(subj==subjN) = zscore(vv(subj==subjN));
                    end
                end
        end
    end
    
    csubj=categorical(subj);
    tt{ipeak} = table(subj,csubj,vv);
end  
clear aT

ll=nan(np+1,length(sigmas),length(taus),2);
%infos=cell(np+1,length(sigmas),length(taus),2);
real_vec = zeros(height(tt{1}),1);%first permutation is the real data
perm_vecs = [ real_vec , round(rand(height(tt{1}),np))];
tictotal = tic;
parfor ip=1:np+1
    disp(ip)
    ticperm = tic;
    perm = perm_vecs(:,ip);
    %[llstp, infostp] = N1P2permutation_allExp(whichExp,perm,uRA,sigmas,taus,tt,art,seqInd,numt,otherpeak,modelName);
    [llstp, ~] = N1P2permutation_allExp(whichExp,perm,uRA,sigmas,taus,tt,art,seqInd,numt,otherpeak,modelName,infoFlag);
%    infos(ip,:,:,:) = infostp;
    ll(ip,:,:,:) = llstp;
    disp(['Done in ' num2str(toc(ticperm))])
end
disp(['Done all in ' num2str(toc(tictotal))])
save(['compareN1P2_' modelName '_' savetag '_' num2str(np)],'ll','sigmas','taus','perm_vecs')
%save(['compareN1P2_' modelName '_' savetag '_' num2str(np) '_infos'],'infos','-v7.3')
        %% plot results
np=250;
modelName = 'lme';
savetag = 'allExp_zscoreFlip';
whichpeaks = {'N1','P2'};
otherpeak = [2 1];

load(['compareN1P2_' modelName '_' savetag '_' num2str(np)],'ll','sigmas','taus')

bestaus = nan(2,size(ll,1));
bestsigmas = nan(2,size(ll,1));
rows = nan(2,size(ll,1));
cols = nan(2,size(ll,1));
dll = nan(2,size(ll,1));
hf=ERPfigure;
set(gcf,'Position',[100 100 1000 500])
for ip=1:size(ll,1)
    clf
    for ipeak=1:2
        llt=squeeze(ll(ip,:,:,ipeak));
        if any(any(any(isnan(llt))))
        else
            h{ipeak}=subplot(1,2,ipeak);
            space = max(max(llt))-llt;
            imagesc(taus,sigmas,space)
            [row col]=find(space==min(min(space))); 
            bestaus(ipeak,ip)=taus(col);
            bestsigmas(ipeak,ip)=sigmas(row);
            rows(ipeak,ip)=row; cols(ipeak,ip)=col;
            hold on
            plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',1)
            plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',1)
            text(2,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',18);
            title(whichpeaks{ipeak})
            xlabel('\tau (seconds)')
            hc{ipeak}=colorbar;
            if ipeak==1 
                ylabel('\sigma (semitone)') 
            else
                ylabel(hc{ipeak},'LLR')
            end
            caxis([0 6])
            set(gca,'fontsize',12)
        end
    end
    for ipeak = 1:2
        if any(any(any(isnan(llt))))
            
        else
            llt=squeeze(ll(ip,:,:,ipeak));
            llo=squeeze(ll(ip,:,:,otherpeak(ipeak)));
            st = max(max(llt))-llt;
            so = max(max(llo))-llo;
            plot(h{otherpeak(ipeak)},[bestaus(ipeak,ip) bestaus(ipeak,ip)],[sigmas(1) sigmas(end)],'Color',[0 1 0],'linewidth',1)
            plot(h{otherpeak(ipeak)},[taus(1) taus(end)],[bestsigmas(ipeak,ip) bestsigmas(ipeak,ip)],'Color',[0 1 0],'linewidth',1)
            dll(ipeak,ip) = abs(st(rows(otherpeak(ipeak),ip),cols(otherpeak(ipeak),ip))-st(rows(ipeak,ip),cols(ipeak,ip)));disp(dll(ipeak,ip))
        end
    end
    suptitle(['log likelihood ratio relative to extremum'])
%pause
    if ip==1
        saveas(gcf,'llrlme_zscore_chi95percent','png')
    end
end

%plot results
ERPfigure;
pvals=nan(size(dll));
for ipeak = 1:length(whichpeaks)
    subplot(2,2,ipeak)
    h{ipeak}=histogram(dll(ipeak,:),20,'normalization','count');
    hold on
    plot([dll(ipeak,1) dll(ipeak,1)],[0 100],'r','linewidth',2)
    %cs=cumsum(h{ipeak}.BinCounts)/250;
    %plot(h{ipeak}.BinEdges(1:end-1),cs)
    title(whichpeaks{ipeak})
    set(gca,'fontsize',14)
    ylim([0 60])

    subplot(2,2,ipeak+2)
    
    for ip=1:length(dll)
        data=dll(ipeak,1:ip);
        ii=sum(data>=data(1));
        pvals(ipeak,ip) = ii/length(data);
    end
    plot(pvals(ipeak,:));
    hold on
    plot([1 length(dll)],[pvals(ipeak,end) pvals(ipeak,end)],'Color',[0 1 0],'linewidth',1)
    text(length(dll),pvals(ipeak,end),num2str(pvals(ipeak,end)),'fontsize',14)    
    xlabel('Number of permutations')
    ylabel('P-value')
    set(gca,'fontsize',14)
    ylim([0 0.5])
end
saveas(gcf,'permutation_results','pdf')
%% Section 7a: fit individual subjects
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [ ~, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );
dataFolder = [Folder filesep 'individual'];
mkdir(dataFolder);
plotFlag = false;

%whichExp = 'all';%can be 1 2 3 or 'all'
numt = 5;
bestaus=cell(3,1);
bestsigmas=cell(3,1);
maxlls=cell(3,1);
minlls=cell(3,1);
meanlls=cell(3,1);
stdvv=cell(3,1);
if plotFlag
    ERPfigure    
end
for whichExp=1:3
    for su = bwhichSubjects{whichExp}
        
        for ipeak = 1:2
            disp([whichExp su ipeak])
            ll=nan(length(sigmas),length(taus));
            infos=cell(length(sigmas),length(taus));
            %prep data
            sel=bT.Pot==ipeak & bT.Subject==whichExp*100+su;
            vv=bT.val(sel);
            art=bT.Art(sel);
            seqInd=bT.seqInd(sel);
            subj=bT.Subject(sel);
            vv=vv(~art & seqInd>numt);
            subj=subj(~art & seqInd>numt);
            csubj=categorical(subj);
            %prep model
            uRA=bRA{whichExp}(su,:,:,:,:);
            
            tictot = tic;
            for it=1:length(taus)
                for isigma=1:length(sigmas)
                    suRA=squeeze(uRA(:,:,:,isigma,it));
                    suRA=suRA(:); % this assumes that the table is arranged as above
                    suRA(isnan(suRA))=[];
                    if 0
                        %another way to calculate seqInd: (verifies that suRA(:) is correct
                      load([Folder filesep 'Exp' num2str(whichExp) filesep 'seqIndcat'])
                      useqIndcat=seqIndcat(su,:,:);
                      useqIndcat=useqIndcat(:);
                      useqIndcat(isnan(useqIndcat))=[];
                      figure;plot(seqInd);hold all;plot(useqIndcat)
                      isequal(seqInd,useqIndcat)
                    end
                    suRA=suRA(~art & seqInd > numt);
                    tt=table(subj,csubj,suRA,vv);
                    %run linear model
                    formula = 'vv~suRA';
                    tmplm = fitlm(tt,formula);
                    ll(isigma,it)=tmplm.LogLikelihood;
                    %infos{isigma,it}=extractInfolme(tmplm,true);
                end
            end
            save([dataFolder filesep 'll' whichpeaks{ipeak} '_Exp' num2str(whichExp) 's' num2str(su)],'ll','formula')
            %save([dataFolder filesep 'll' whichpeaks{ipeak} '_Exp' num2str(whichExp) 's' num2str(su) '_infos'],'infos')
            
            disp(['Done ' whichpeaks{ipeak} ' in ' num2str(toc(tictot))])
            
            [row col]=find(ll==max(max(ll)));    
            bestaus{whichExp}(ipeak,su)=taus(col);
            bestsigmas{whichExp}(ipeak,su)=sigmas(row);
            maxlls{whichExp}(ipeak,su)=max(max(ll));
            minlls{whichExp}(ipeak,su)=min(min(ll));
            meanlls{whichExp}(ipeak,su)=mean(mean(ll));
            stdvv{whichExp}(ipeak,su)=std(vv);
            if plotFlag
                subplot(1,2,ipeak)
                imagesc(taus,sigmas,ll);colorbar
                title([whichpeaks{ipeak} ])
                colorbar
                xlabel('taus')
                ylabel('sigmas')
                hold on
                plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',1)
                plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',1)
                text(2,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',18);
                drawnow
            end
        end
        if plotFlag
            suptitle(['loglikelihood. Exp' num2str(whichExp) ', subj ' num2str(su)])
            pause
        end
    end
end
save([dataFolder filesep 'summary'],'bestaus','bestsigmas','maxlls','minlls','meanlls','stdvv')
%% Section 7b: fit individuals averaging trials due to RA
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [ ~, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );
dataFolder = [Folder filesep 'individual_avg'];
mkdir(dataFolder);
plotFlag = true;
avgN = 10;%how many trials to average

%whichExp = 'all';%can be 1 2 3 or 'all'
numt = 5;
bestaus=cell(3,1);
bestsigmas=cell(3,1);
maxlls=cell(3,1);
minlls=cell(3,1);
meanlls=cell(3,1);
stdvv=cell(3,1);
if plotFlag
    ERPfigure    
end
for whichExp=2:3
    for su = bwhichSubjects{whichExp}
        for ipeak = 1:2
            disp([whichExp su ipeak])
            ll=nan(length(sigmas),length(taus));
            infos=cell(length(sigmas),length(taus));
            %prep data
            sel=bT.Pot==ipeak & bT.Subject==whichExp*100+su;
            vv=bT.val(sel);
            art=bT.Art(sel);
            seqInd=bT.seqInd(sel);
            subj=bT.Subject(sel);
            vv=vv(~art & seqInd>numt);
            subj=subj(~art & seqInd>numt);
            csubj=categorical(subj);
            %prep model
            uRA=bRA{whichExp}(su,:,:,:,:);
            
            tictot = tic;
            for it=1:length(taus)
                for isigma=1:length(sigmas)
                    suRA=squeeze(uRA(:,:,:,isigma,it));
                    suRA=suRA(:); % this assumes that the table is arranged as above
                    suRA(isnan(suRA))=[];
                  
                    suRA=suRA(~art & seqInd > numt);
                    [sRA,I]=sort(suRA);
                    svv=vv(I);
                    %average N trials:
                    len=length(svv);
                    m = len - mod(len, avgN);
                    y = reshape(svv(1:m), avgN, []);     % Reshape x to a [n, m/n] matrix
                    Avgsvv = transpose(sum(y, 1) / avgN);
                    y = reshape(sRA(1:m), avgN, []);     % Reshape x to a [n, m/n] matrix
                    AvgsRA = transpose(sum(y, 1) / avgN);
                    
                    tt=table(AvgsRA,Avgsvv);
                    %run linear model
                    formula = 'msvv~AvgsRA';
                    tmplm = fitlm(tt,formula);
                    ll(isigma,it)=tmplm.LogLikelihood;
                    %infos{isigma,it}=extractInfolme(tmplm,true);
                end
            end
            save([dataFolder filesep 'll' whichpeaks{ipeak} '_Exp' num2str(whichExp) 's' num2str(su)],'ll','formula')
            %save([dataFolder filesep 'll' whichpeaks{ipeak} '_Exp' num2str(whichExp) 's' num2str(su) '_infos'],'infos')
            
            disp(['Done ' whichpeaks{ipeak} ' in ' num2str(toc(tictot))])
            
            [row col]=find(ll==max(max(ll)));    
            bestaus{whichExp}(ipeak,su)=taus(col);
            bestsigmas{whichExp}(ipeak,su)=sigmas(row);
            maxlls{whichExp}(ipeak,su)=max(max(ll));
            minlls{whichExp}(ipeak,su)=min(min(ll));
            meanlls{whichExp}(ipeak,su)=mean(mean(ll));
            stdvv{whichExp}(ipeak,su)=std(Avgsvv);
            if plotFlag
                subplot(1,2,ipeak)
                imagesc(taus,sigmas,ll);colorbar
                title([whichpeaks{ipeak} ])
                colorbar
                xlabel('taus')
                ylabel('sigmas')
                hold on
                plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',1)
                plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',1)
                text(2,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',18);
                drawnow
            end
        end
        if plotFlag
            suptitle(['loglikelihood. Exp' num2str(whichExp) ', subj ' num2str(su)])
            pause
        end
    end
end
save([dataFolder filesep 'summary'],'bestaus','bestsigmas','maxlls','minlls','meanlls','stdvv')
        %% plot results
figure
for ie=1:3
    ws=bwhichSubjects{ie};
    clear wss
    for is=1:length(ws)
        wss{is}=num2str(ws(is));
    end
    difftaus{ie}=bestaus{ie}(1,:)-bestaus{ie}(2,:);
    subplot(1,3,ie)
    boxplot(difftaus{ie}(ws))
    hold on
    plot(mean(stdvv{ie}(:,ws)),difftaus{ie}(ws),'o','markersize',12)
    text(mean(stdvv{ie}(:,ws)),difftaus{ie}(ws),wss,'fontsize',16)
    title(['Exp' num2str(ie)])
    if ie==1
        ylabel(['N1 \tau - P2 \tau'])
    end
    xlabel('std of data')
end
%% Section 8a: fit sub-groups of subjects
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [ ~, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );

dataFolder = [Folder filesep 'subgroups'];
mkdir(dataFolder);
addpath('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')

plotFlag = false;
whichExp = 'all';%can be 1 2 3 or 'all'
numt = 5;
modelName = 'lme';

percent = 10;
np=100;
wsp = cell(np,3);
for ie=1:3
    ns=length(bwhichSubjects{ie});
    nsp=round(percent/100*ns);
    for ip=1:np
        wsp{ip,ie}=sort(bwhichSubjects{ie}(randperm(ns,nsp)));
    end
end
ll=nan(np,length(sigmas),length(taus),2);
bestaus=nan(np,2);
bestsigmas=nan(np,2);
cols=nan(np,2);
rows=nan(np,2);
infos=cell(np,length(sigmas),length(taus),2);
%ERPfigure;
tictot=tic;
parfor ip=1:np
    %disp([ip])
    ticperm=tic;
    [ llp, bestausp, bestsigmasp, rowsp, colsp ] = subgroupsN1P2(ip, whichExp, bT, wsp, bRA, numt, taus, sigmas, plotFlag, modelName );
    ll(ip,:,:,:)=llp;
    bestaus(ip,:)=bestausp;bestsigmas(ip,:)=bestsigmasp;
    rows(ip,:)=rowsp;cols(ip,:)=colsp;
    disp(['Done ' num2str(ip) ' in ' num2str(toc(ticperm))])
end
disp(['Done all in ' num2str(toc(tictot))])
save([dataFolder filesep 'subgroups_' num2str(percent) 'percent_' num2str(np)],'ll','wsp','bestaus','bestsigmas','rows','cols')
    %save([dataFolder filesep 'subgroups_' num2str(percent) 'percent_' num2str(np) '_infos'])
%% Section 8b - bootstrap to compare N1 and P2
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [ ~, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );

dataFolder = [Folder filesep 'subgroups'];
mkdir(dataFolder);
addpath('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')

plotFlag = false;
whichExp = 'all';%can be 1 2 3 or 'all'
numt = 5;
modelName = 'lme';

np=100;
wsp = cell(np+1,3);
for ie=1:3
    ns=length(bwhichSubjects{ie});
    for ip=1:np+1
        if ip==1
            wsp{ip,ie}=bwhichSubjects{ie};
        else
            wsp{ip,ie}=sort(bwhichSubjects{ie}(randi(ns,size(bwhichSubjects{ie}))));
        end
    end
end
ll=nan(np+1,length(sigmas),length(taus),2);
bestaus=nan(np+1,2);
bestsigmas=nan(np+1,2);
cols=nan(np+1,2);
rows=nan(np+1,2);
infos=cell(np+1,length(sigmas),length(taus),2);
%ERPfigure;
tictot=tic;
parfor ip=1:np+1
    disp([ip])
    ticperm=tic;
    [ llp, bestausp, bestsigmasp, rowsp, colsp ] = bootstrapN1P2(ip, whichExp, bT, wsp, bRA, numt, taus, sigmas, plotFlag, modelName );
    ll(ip,:,:,:)=llp;
    bestaus(ip,:)=bestausp;bestsigmas(ip,:)=bestsigmasp;
    rows(ip,:)=rowsp;cols(ip,:)=colsp;
    disp(['Done ' num2str(ip) ' in ' num2str(toc(ticperm))])
end
disp(['Done all in ' num2str(toc(tictot))])
save([dataFolder filesep 'bootstrap_' num2str(np)],'ll','wsp','bestaus','bestsigmas','rows','cols')
    %% plot
    percent = 100;
    np=100;
    dataFolder = [Folder filesep 'subgroups'];
    if percent==100
        load([dataFolder filesep 'bootstrap_' num2str(np)])
    else
        load([dataFolder filesep 'subgroups_' num2str(percent) 'percent_' num2str(np)])
    end
    ratioFlag = true;
    hf=figure;
    set(hf,'Position',[10 10 400 500])
    if ratioFlag 
        datataus = log(bestaus(:,1)./bestaus(:,2));
        datasigmas = log(bestsigmas(:,1)./bestsigmas(:,2));
    else
        datataus = bestaus(:,1)-bestaus(:,2);
        datasigmas = bestsigmas(:,1)-bestsigmas(:,2);
    end
    hb=boxplot([datasigmas datataus],'Notch','on','Labels',{'\sigma','\tau'});
    set(hb,'linewidth',2)
    title('\theta max N1 versus \theta max P2 in 100 bootstrap runs')
    if ratioFlag
        ylabel('log (\theta max N1 / \theta max P2)')
    else
        ylabel('Differences (\sigma: semitones, \tau: seconds)')
    end
    hold on
    plot(1+ones(size(datataus))+0.02*randn(size(datataus)),datataus,'og','markersize',4)    
    plot(ones(size(datasigmas))+0.02*randn(size(datasigmas)),datasigmas,'og','markersize',4)    
    plot([0 3],[0 0],'k:' )
    set(gca,'fontsize',14)
    if ratioFlag
        saveas(gcf,[dataFolder filesep 'bootstrap_params_rats'],'png')
        saveas(gcf,[dataFolder filesep 'bootstrap_params_rats'],'pdf')
    else
        saveas(gcf,[dataFolder filesep 'bootstrap_params_diffs'],'png')
        saveas(gcf,[dataFolder filesep 'bootstrap_params_diffs'],'pdf')
    end
    %    boxplot([bestaus(:,1),bestaus(:,2)])
    [Htau,ptau]=ttest(difftaus);
    [Hsigma,psigma]=ttest(diffsigmas)
    
    figure
    CIp=90;%confidence interval percent
    tails=(100-CIp)/2;
    CI=nan(2,2);
    realParams=[8,8;3.2,1];
    diffs={bestsigmas,bestaus};
    spi=0;
    for isigtau=1:2
        for ipeak=1:2
            spi=spi+1;
            subplot(2,2,spi)
            bootdata=diffs{isigtau}(:,ipeak);
            histogram(bootdata,10)
            hold on
            %xlim([0 5])
            %ylim([0 60])
            title(whichpeaks{ipeak})
            plot([mean(bootdata) mean(bootdata)],[0 60],'linew',2)
            plot([realParams(isigtau,ipeak) realParams(isigtau,ipeak)], [0 60],'linew',2)
            %CI
            sampsize=length(bootdata);
            leftsamp=ceil(sampsize*tails/100);rightsamp=sampsize-floor(sampsize*tails/100);
            CI(ipeak,1)=bootdata(leftsamp);
            CI(ipeak,2)=bootdata(rightsamp);
            plot([CI(ipeak,1) CI(ipeak,2)],[0 0],'g','linew',5)
            if ipeak==2 && isigtau==2
                legend({'bootstrap dist.','mean boostrap value','best params from data',[num2str(CIp) '% confidence interval']},'fontsize',10)
            end
        end
    end
    saveas(gcf,[dataFolder filesep 'bootstrap_dist'],'pdf')       
%% Section 9: model each spread separately
%only Exp 3, model each frequency spread separately.
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [ dataFolder, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );

whichExp=3;
modelName = 'lme';
numt = 5;
saveFolder = [Folder filesep 'spreads' filesep];
mkdir(saveFolder)
ll=nan(length(sigmas),length(taus),3,2);
%    infos=cell(length(sigmas),length(taus));
spread2block = {1,[2,3],[4,5]};
for ispread=1:3
    uRA=bRA{whichExp}(bwhichSubjects{whichExp},spread2block{ispread},:,:,:);
    disp(ispread)
    for ipeak = 1:2
        disp(ipeak)
        sel=bT.Pot==ipeak & bT.ExpN==whichExp & ismember(bT.Block,spread2block{ispread});%fix
        %clear bRA
        vv=bT.val(sel);
        art=bT.Art(sel);
        seqInd=bT.seqInd(sel);
        subj=bT.Subject(sel);
        vv=vv(~art & seqInd>numt);
        subj=subj(~art & seqInd>numt);
        csubj=categorical(subj);

        tictot = tic;
        llst=nan(length(sigmas),length(taus));
        parfor it=1:length(taus)
            [ lls ] = calcLLsigmas( it, taus, sigmas, whichExp, uRA, art, seqInd, numt, vv, subj, csubj, modelName);
            llst(:,it)=lls;
        end
        ll(:,:,ispread,ipeak)=llst;
        disp(['Done ' whichpeaks{ipeak} ', spread' num2str(ispread) ' in ' num2str(toc(tictot))])

    end
end
save([saveFolder 'llspreads'],'ll','spread2block','whichExp')
    %% plot
    Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
%[ dataFolder, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );

    saveFolder = [Folder filesep 'spreads' filesep];
    load([saveFolder 'llspreads'],'ll','spread2block','whichExp')

    bestParams=nan(3,2,2);
    hf=ERPfigure;
    set(hf,'Position',[100 100 1000 600])
    isp=0;
    for ipeak=1:2
        for ispread=1:3
            isp=isp+1;
            sph{isp}=subplot(2,3,isp);
            llt=squeeze(ll(:,:,ispread,ipeak));
            imagesc(taus,sigmas,-exp(1)*(llt-max(max(llt))))
            [row col]=find(llt==max(max(llt)));    
            hold on 
            plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',1)
            plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',1)
            text(2,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',18);
            bestParams(ispread,ipeak,1)=sigmas(row);%sigma
            bestParams(ispread,ipeak,2)=taus(col);%tau
 
            colorbar
            %caxis([0 9.2])
            if ispread==1
                ylabel('\sigma')
            end
            if ipeak==2
                xlabel('\tau')
            end
            title([whichpeaks{ipeak} ' , condition ' num2str(ispread)])
            set(gca,'fontsize',16)
        end
    end
    suptitle(['Spreads'])
    saveas(gcf,[saveFolder 'llspreadsN1P2'],'png')
    saveas(gcf,[saveFolder 'llspreadsN1P2'],'pdf')
    
    for isp=1:length(sph)
        caxis(sph{isp},[0 20])
    end
    saveas(gcf,[saveFolder 'llspreadsN1P2_caxis'],'png')
    saveas(gcf,[saveFolder 'llspreadsN1P2_caxis'],'pdf')
   
    
%     for isp=1:length(sph)
%         caxis(sph{isp},[0 6])
%     end
%     saveas(gcf,[saveFolder 'llspreadsN1P2_CI95chi'],'png')
%     saveas(gcf,[saveFolder 'llspreadsN1P2_CI95chi'],'pdf')
%    
    save([saveFolder 'bestParams'],'bestParams');
    
    %plot just best sigmas
    freqspreads = [[10,10,9,10];[7 7 7 7];[4,4,3,3]];
    meanfreqspreads = mean(freqspreads');
    ERPfigure;
    hold on
    plot(meanfreqspreads,bestParams(:,1,1),'b','linewidth',2)
    plot(meanfreqspreads,bestParams(:,2,1),'r','linewidth',2)
    set(gca,'fontsize',20)
    legend({'N1','P2'},'fontsize',20,'Location','nw')
    plot(meanfreqspreads,bestParams(:,1,1),'.b','MarkerSize',50)
    plot(meanfreqspreads,bestParams(:,2,1),'.r','MarkerSize',50)
    title('\sigma as a function of frequency spread','fontsize',20)
    ylabel('Best fit \sigma (semitones)','fontsize',20)
    xlabel('Mean spread in the sequence (semitones)','fontsize',20)
    saveas(gcf,[saveFolder 'bestsigmas'],'png')
%% Section 10: bootstrap each spread separately
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [ dataFolder, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );

whichExp=3;
modelName = 'lme';
infoFlag=false;
numt = 5;
saveFolder = [Folder filesep 'spreads' filesep];
tag = '2';
mkdir(saveFolder)

np=100;
whichSubjects = bwhichSubjects{whichExp};
ns=length(whichSubjects);
wsp = cell(np+1,1);
for ip=1:np+1
    if ip==1
        wsp{ip}=whichSubjects;
    else
        wsp{ip}=sort(whichSubjects(randi(ns,size(whichSubjects))));
    end
end

ll=nan(np+1,length(sigmas),length(taus),3,2);
%    infos=cell(np,length(sigmas),length(taus));
spread2block = {1,[2,3],[4,5]};
bRAe=bRA{whichExp};
bTe=bT(bT.ExpN==whichExp,:);

parfor ip=1:np+1
    
    [ llp ] = bootstrapSpreads(ip, whichExp, bRAe, bTe, wsp, spread2block, whichpeaks, sigmas, taus, numt, modelName, infoFlag );
    ll(ip,:,:,:,:) = llp;
end
save([saveFolder 'llspreads_bootstrap_' tag],'ll','spread2block','whichExp','wsp')
    %% plot
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
saveFolder = [Folder filesep 'spreads' filesep];

load([saveFolder 'llspreads_bootstrap'])
%calc best params
bestausboot=nan(size(ll,1),3,2);
bestsigmasboot=nan(size(ll,1),3,2);
for ip=1:size(ll,1)
    for ispread=1:3
        for ipeak=1:2
            llp=squeeze(ll(ip,:,:,ispread,ipeak));
            [row col]=find(llp==max(max(llp)));
            bestausboot(ip,ispread,ipeak)=taus(col);
            bestsigmasboot(ip,ispread,ipeak)=sigmas(row);
        end
    end
end
save([saveFolder 'bestSigmasBoot'],'bestsigmasboot')
%plot histograms to check normality
load([saveFolder 'bestParams']);
freqspreads = [[10,10,9,10];[7 7 7 7];[4,4,3,3]];
meanfreqspreads = mean(freqspreads');

CIp=90;%confidence interval percent
tails=(100-CIp)/2;
ERPfigure;
isp=0;
for ispread=1:3
    for ipeak=1:2
        isp=isp+1;
        subplot(3,2,isp)
        %h1=histogram(bestausboot(:,ispread,ipeak),10);
        %hold on
        bootdata=sort(bestsigmasboot(:,ispread,ipeak));
        h2=histogram(bootdata,15);
        xlim([1 18])
        ylim([0 60])
        if ispread==1
            title(whichpeaks{ipeak})
        end
        hold on
        plot([mean(bootdata) mean(bootdata)],[0 60],'linew',2)
        plot([bestParams(ispread,ipeak,1) bestParams(ispread,ipeak,1)], [0 60],'linew',2)
        %CI
        sampsize=size(ll,1);
        leftsamp=ceil(sampsize*tails/100);rightsamp=sampsize-floor(sampsize*tails/100);
        CI(ispread,ipeak,1)=bootdata(leftsamp);
        CI(ispread,ipeak,2)=bootdata(rightsamp);
        plot([CI(ispread,ipeak,1) CI(ispread,ipeak,2)],[0 0],'g','linew',5)

        if ipeak==1
            ylabel([num2str(meanfreqspreads(ispread)) ' (semitones)'])
        end
        if ispread==3
            xlabel('best \sigma (semitones)')
            if ipeak==2
                legend({'bootstrap \sigma dist.','mean boostrap \sigma','best \sigma from data',[num2str(CIp) '% confidence interval']},'fontsize',10)
            end
        end
        set(gca,'fontsize',10)
    end
end
suptitle('bootstrap spreads')
saveas(gca,[saveFolder 'bootstrapSpreads'],'png')

%plot pairwise comparisons:
h=ERPfigure;
set(h,'Position',[100 100 500 800])
isp=0;
labels = {'spread 1-2 (2.75 semitone)','2-3 (3.5 semitone)','1-3 (6.25 semitone)'};
comparisons =[1,2;2,3;1,3];
for ic=1:3
    for ipeak=1:2
        isp=isp+1;
        subplot(3,2,isp)
        %h1=histogram(bestausboot(:,ispread,ipeak),10);
        %hold on
        diffdata=sort(bestsigmasboot(:,comparisons(ic,1),ipeak)-bestsigmasboot(:,comparisons(ic,2),ipeak));
        boxplot(diffdata)
        hold on
        plot(ones(size(diffdata))+0.02*randn(size(diffdata)),diffdata,'.','markersize',10)
        if ipeak==1
            ylabel(labels{ic})
        end
        if ic==1
            title(whichpeaks{ipeak})
        end
        set(gca,'fontsize',10)
        ylim([-10 10])
        hold on
        plot([0.5 1.5],[0 0],'k:')
    end
end
suptitle('bootstrap spreads pairwise comparisons')
saveas(gca,[saveFolder 'bootSpreads_pairComp'],'png')

%plot best fit sigma with error bars from bootstrap:
errorbars = nan(length(meanfreqspreads),2);
for ipeak = 1:2
    for ispread = 1:length(meanfreqspreads) 
        %CIs = Confidence(bestParams(:,ispread,ipeak,1));
        %errorbars(ispread,ipeak) = mean(bestParams(:,ispread,ipeak,1)) - CIs(1);
        x = bestsigmasboot(:,ispread,ipeak);
        errorbars(ispread,ipeak) = nanstd(x);
    end
end
ERPfigure;
hold on
plot(meanfreqspreads,bestParams(:,1,1),'b','linewidth',2)
plot(meanfreqspreads,bestParams(:,2,1),'r','linewidth',2)
set(gca,'fontsize',20)
legend({'N1','P2'},'fontsize',20,'Location','nw')
plot(meanfreqspreads,bestParams(:,1,1),'.b','MarkerSize',50,'HandleVisibility','off')
plot(meanfreqspreads,bestParams(:,2,1),'.r','MarkerSize',50,'HandleVisibility','off')
errorbar(meanfreqspreads,bestParams(:,1,1),errorbars(:,1),'b','HandleVisibility','off')
errorbar(meanfreqspreads,bestParams(:,2,1),errorbars(:,1),'r','HandleVisibility','off')
ylim([2 14])   
title('\sigma as a function of frequency spread','fontsize',20)
ylabel('Best fit \sigma (semitones)','fontsize',20)
xlabel('Mean spread in the sequence (semitones)','fontsize',20)
saveas(gcf,[saveFolder 'bestsigmas_booterrors'],'png')
%% anova of lme for the bootstrap
%create table
T=table;
for ip=1:size(bestsigmasboot,1)
    for ispread=1:size(bestsigmasboot,2)
        for ipeak=1:size(bestsigmasboot,3)
            tmp=struct;
            tmp.pot=ipeak;
            tmp.Cpot=categorical(ipeak);
            tmp.bootN=ip;
            tmp.CbootN=categorical(ip);
            tmp.bestsig=bestsigmasboot(ip,ispread,ipeak);
            tmp.cond=categorical(ispread);
            tmp.spreads=meanfreqspreads(ispread);
            T=[T; struct2table(tmp)];
        end
    end
end
formula1='bestsig~spreads+(spreads|bootN)';
lme1=fitlme(T,formula1);
anova(lme1)
formula2='bestsig~spreads*Cpot+(spreads|bootN)+(Cpot|bootN)';
lme2=fitlme(T,formula2);
anova(lme2)
compare(lme1,lme2)
T1=T(T.pot==1,:);
T2=T(T.pot==2,:);
lmeN1=fitlme(T1,formula1);
lmeP2=fitlme(T2,formula1);
%% scramble labels within spreads
Folder = ['S:\Lab-Shared\Experiments\N1P2\Analysis\Model\singleTrials'];
 [ dataFolder, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder );

whichExp=3;
modelName = 'lme';
infoFlag=false;
numt = 5;
saveFolder = [Folder filesep 'spreads' filesep];
np=100;
tag = num2str(np);
mkdir(saveFolder)

whichSubjects = bwhichSubjects{whichExp};
ns=length(whichSubjects);

ll=nan(np+1,length(sigmas),length(taus),3,2);
%    infos=cell(np,length(sigmas),length(taus));
spread2block = {1,[2,3],[4,5]};
bRAe=bRA{whichExp};
bTe=bT(bT.ExpN==whichExp,:);
%create scrambles
scrambles=cell(np+1,3);
for ispread=1:3
    ht=height(bTe(bTe.Pot==1 & ismember(bTe.Block,spread2block{ispread}),:));
    scrambles{1,ispread}=1:ht;
    for ip=2:np+1
        scrambles{ip,ispread}=randperm(ht);
    end
end

parfor ip=1:np+1
    disp(ip)
    [ llp ] = scrambleSpreads(ip, scrambles, whichExp, bRAe, bTe, whichSubjects, spread2block, whichpeaks, sigmas, taus, numt, modelName, infoFlag );
    ll(ip,:,:,:,:) = llp;
end
save([saveFolder 'llspreads_scramble_' tag],'ll','spread2block','whichExp','scrambles')
    %% plot
    saveFolder = [Folder filesep 'spreads' filesep];
    np=100;
    tag = num2str(np);
    load([saveFolder 'llspreads_scramble_' tag])
    bestParams=nan(size(ll,1)-1,3,2,2);
    ERPfigure;
   
    for ip=2:size(ll,1)
        isp=0;
        for ipeak=1:2
            for ispread=1:3
                isp=isp+1;
                subplot(2,3,isp)
                llt=squeeze(ll(ip,:,:,ispread,ipeak));
                imagesc(taus,sigmas,-2*(llt-max(max(llt))))
                [row col]=find(llt==max(max(llt)));    
                hold on 
                plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',1)
                plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',1)
                text(2,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',18);
                bestParams(ip-1,ispread,ipeak,1)=sigmas(row);%sigma
                bestParams(ip-1,ispread,ipeak,2)=taus(col);%tau
 
                %colorbar
                caxis([0 6])
                if ispread==1
                    ylabel('\sigma')
                end
                if ipeak==2
                    xlabel('\tau')
                end
                title([whichpeaks{ipeak} ' , condition ' num2str(ispread)])
            end
        end
        suptitle(['Spreads scramble #' num2str(ip)])
        drawnow
    end
    saveas(gcf,[saveFolder 'llspreads_scrambled'],'png')
    bestParamsScrambled=bestParams;
    save([saveFolder 'bestParamsScrambled'],'bestParamsScrambled');
    clear bestParams
    %% plot best sigmas real+scrambled
    load([saveFolder 'bestParams']);
    load([saveFolder 'bestSigmasBoot'])
    load([saveFolder 'bestParamsScrambled']);

    %First, the scrambled:
    errorbars = nan(length(meanfreqspreads),2);
    for ipeak = 1:2
        for ispread = 1:length(meanfreqspreads) 
            %CIs = Confidence(bestParams(:,ispread,ipeak,1));
            %errorbars(ispread,ipeak) = mean(bestParams(:,ispread,ipeak,1)) - CIs(1);
            x = bestParamsScrambled(:,ispread,ipeak,1);
            errorbars(ispread,ipeak) = nanstd(x)/sqrt(length(x(~isnan(x))));
        end
    end
    ERPfigure;
    hold on
    h1=plot(meanfreqspreads,mean(bestParamsScrambled(:,:,1,1)),'g-.','linewidth',2);
    h2=plot(meanfreqspreads,mean(bestParamsScrambled(:,:,2,1)),'m-.','linewidth',2);
    plot(meanfreqspreads,mean(bestParamsScrambled(:,:,1,1)),'.g','MarkerSize',30,'HandleVisibility','off')
    plot(meanfreqspreads,mean(bestParamsScrambled(:,:,2,1)),'.m','MarkerSize',30,'HandleVisibility','off')
    errorbar(meanfreqspreads,mean(bestParamsScrambled(:,:,1,1)),errorbars(:,1),'g:','HandleVisibility','off')
    errorbar(meanfreqspreads,mean(bestParamsScrambled(:,:,2,1)),errorbars(:,1),'m:','HandleVisibility','off')
    ylim([2 14])
    xlim([2 14])
    plot([2 14],[2 14],'k:')
    title('\sigma as a function of frequency spread')
    ylabel('Best fit \sigma (semitones)')
    xlabel('Mean spread in the sequence (semitones)')
    %now add the real data:
    %calc error bars from bootstrap
    errorbars = nan(length(meanfreqspreads),2);
    for ipeak = 1:2
        for ispread = 1:length(meanfreqspreads) 
            %CIs = Confidence(bestParams(:,ispread,ipeak,1));
            %errorbars(ispread,ipeak) = mean(bestParams(:,ispread,ipeak,1)) - CIs(1);
            x = bestsigmasboot(:,ispread,ipeak);
            errorbars(ispread,ipeak) = nanstd(x);
        end
    end
   
    h3=plot(meanfreqspreads,bestParams(:,1,1),'b','linewidth',2);
    h4=plot(meanfreqspreads,bestParams(:,2,1),'r','linewidth',2);
    plot(meanfreqspreads,bestParams(:,1,1),'.b','MarkerSize',50,'HandleVisibility','off')
    plot(meanfreqspreads,bestParams(:,2,1),'.r','MarkerSize',50,'HandleVisibility','off')
    errorbar(meanfreqspreads,bestParams(:,1,1),errorbars(:,1),'b','HandleVisibility','off')
    errorbar(meanfreqspreads,bestParams(:,2,1),errorbars(:,1),'r','HandleVisibility','off')
    title('\sigma as a function of frequency spread')
    ylabel('Best fit \sigma (semitones)')
    xlabel('Mean spread in the sequence (semitones)')
    hl=legend([h3,h4,h1,h2],{'N1','P2','N1 scrambled','P2 scrambled'},'Location','se');
    set(gca,'fontsize',20)
    set(hl,'fontsize',16)
    saveas(gcf,[saveFolder 'bestsigmasSpreads_boot+scrambled'],'png')
    saveas(gcf,[saveFolder 'bestsigmasSpreads_boot+scrambled'],'pdf')
    