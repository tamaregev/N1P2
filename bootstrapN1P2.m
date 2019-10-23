function [ llp, bestausp, bestsigmasp, rowsp, colsp ] = bootstrapN1P2(ip, whichExp, bT, wsp, bRA, numt, taus, sigmas, plotFlag, modelName  )
%BOOTSTRAPN1P2 estimates CI of best tau and sigma 

llp=nan(length(sigmas),length(taus),2);
bestausp=nan(1,2);
bestsigmasp=nan(1,2);
colsp=nan(1,2);
rowsp=nan(1,2);
% prep uRA as normal (first perm) and calc allwsp (which subjects in this permutations, coded for table of all exp)
switch whichExp
    case {1,2,3}
       % disp(['Exp' num2str(whichExp)])
        uRA=bRA{whichExp}(wsp{1,whichExp},:,:,:,:);
       allwsp=wsp{ip,whichExp}+whichExp*100;
    case 'all'
       % disp(whichExp)
        clear uRA
        allwsp=[];
        for iexp=1:3
            uRA{iexp}=bRA{iexp}(wsp{1,iexp},:,:,:,:);
            allwsp=[allwsp,wsp{ip,iexp}+iexp*100];
        end
end


for ipeak = 1:2
    %pT is a table for only this peak and experiment
    switch whichExp
        case {1,2,3}
            pT=bT(bT.Pot==ipeak & bT.ExpN==whichExp,:);
        case 'all'
            pT=bT(bT.Pot==ipeak,:);
    end
    %nT - new Table with the bootstrap participants:
    nT=table;
    % first keep only rows that exist in current group
    rowin=ismember(pT.Subject,allwsp);
    nT(rowin,:)=pT(rowin,:);
    nT(nT.Subject==0,:)=[];
    %Then add duplications
    for s = unique(allwsp)
        nsub=sum(allwsp==s);
        if nsub>1
            for in=1:nsub-1
                nT=[nT; pT(pT.Subject==s,:)];
            end
        end
    end

    vv=nT.val;
    art=nT.Art;
    seqInd=nT.seqInd;
    subj=nT.Subject;
    
    vv=vv(~art & seqInd>numt);
    subj=subj(~art & seqInd>numt);
    csubj=categorical(subj);
    tictot = tic;
    for it=1:length(taus)
        tictau=tic;
        for isigma=1:length(sigmas)
            ticsigtau=tic;
            
            switch whichExp
                case {1,2,3}
                    suRA=squeeze(uRA(:,:,:,isigma,it));
                    suRA=suRA(:); % this assumes that the table is arranged as above
                    suRA(isnan(suRA))=[];
                    if 0
                        %another way to calculate seqInd: (verifies that suRA(:) is correct
                        load([Folder filesep 'Exp' num2str(whichExp) filesep 'seqIndcat.mat'])   
                        useqIndcat=squeeze(seqIndcat(wsp{ip,whichExp},:,:));
                        useqIndcat(isnan(useqIndcat))=[];
                        %isequal(isnan(useqIndcatall),isnan(seqInd))
                    end
                case 'all'
                    suRA = [];
                    for iexp=1:3
                        tmpRA=squeeze(uRA{iexp}(:,:,:,isigma,it));
                        tmpRA=tmpRA(:);
                        suRA=[suRA; tmpRA]; 
                    end
                    suRA(isnan(suRA))=[];
                    if 0
                        %another way to calculate seqInd: (verifies that suRA(:) is correct
                        useqIndcatall=[];
                        for iexp = 1:3
                            load([Folder filesep 'Exp' num2str(iexp) filesep 'seqIndcat.mat'])   
                             useqIndcat=squeeze(seqIndcat(wsp{ip,iexp},:,:));
                             useqIndcat=useqIndcat(:);
                             useqIndcatall=[useqIndcatall;useqIndcat];
                        end
                        useqIndcatall(isnan(useqIndcatall))=[];
                        figure;plot(useqIndcatall);hold all;plot(seqInd)
                        isequal(seqInd,useqIndcatall)
                    end
            end
            % transform suRA as table due to bootstrap :
            nsuRA=nan(size(suRA));
            nsuRA(rowin,:)=suRA(rowin,:);
            nsuRA(isnan(nsuRA))=[];
            %Then add duplications
            for s = unique(allwsp)
                nsub=sum(allwsp==s);
                if nsub>1
                    for in=1:nsub-1
                        nsuRA=[nsuRA; suRA(pT.Subject==s,:)];
                    end
                end
            end

            
            nsuRA=nsuRA(~art & seqInd > numt);
            tt=table(subj,csubj,nsuRA,vv);
            switch modelName
                case 'lm'
                    formula = 'vv~nsuRA*csubj';
                    tmplme = fitlm(tt,formula);
                case 'lme'
                    formula = 'vv~nsuRA+(nsuRA|subj)';
                    tmplme = fitlme(tt,formula);
            end
            llp(isigma,it,ipeak)=tmplme.LogLikelihood;
            %infos{ip,isigma,it,ipeak} = extractInfolme(tmplme,0);
            %disp(['Done perm, itau, isigma = ' num2str([ip it isigma]) ' in ' num2str(toc(ticsigtau))])
        end
        disp(['Done perm, itau, = ' num2str([ip it]) ' in ' num2str(toc(tictau))])
    end
    llt=squeeze(llp(:,:,ipeak));
    [row col]=find(llt==max(max(llt)));
    bestausp(ipeak)=taus(col);
    bestsigmasp(ipeak)=sigmas(row);
    rowsp(ipeak)=row; colsp(ipeak)=col;

    if plotFlag
        subplot(1,2,ipeak)
        imagesc(taus,sigmas,llt-max(max(llt)));
        title([whichpeaks{ipeak} ])
        colorbar
        xlabel('taus')
        ylabel('sigmas')
        hold on
        plot([taus(col) taus(col)],[sigmas(1) sigmas(end)],'Color',[1 0 0],'linewidth',1)
        plot([taus(1) taus(end)],[sigmas(row) sigmas(row)],'Color',[1 0 0],'linewidth',1)
        text(2,15,['(' num2str(taus(col)) ',' num2str(sigmas(row)) ')'],'Color',[1,0,0],'fontsize',18);
    end
end
if plotFlag
    suptitle(['LLR, Exp' num2str(whichExp) '. Subgroup #' num2str(ip) ' ' num2str(percent) ' %'])
    drawnow
end    

end

