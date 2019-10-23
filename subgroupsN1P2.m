function [ llp, bestausp, bestsigmasp, rowsp, colsp ] = subgroupsN1P2(ip, whichExp, bT, wsp, bRA, numt, taus, sigmas, plotFlag, modelName  )
%SUBGROUPSN1P2 estimates tau and sigma for a random subgroup of all participants 

llp=nan(length(sigmas),length(taus),2);
bestausp=nan(1,2);
bestsigmasp=nan(1,2);
colsp=nan(1,2);
rowsp=nan(1,2);
% prep uRA and calc allwsp (which subjects in this permutations, coded for table of all exp)
switch whichExp
    case {1,2,3}
       % disp(['Exp' num2str(whichExp)])
       allwsp=wsp{ip,whichExp}+whichExp*100;
        uRA=bRA{whichExp}(wsp{ip,whichExp},:,:,:,:);
    case 'all'
       % disp(whichExp)
        clear uRA
        allwsp=[];
        for iexp=1:3
            uRA{iexp}=bRA{iexp}(wsp{ip,iexp},:,:,:,:);
            allwsp=[allwsp,wsp{ip,iexp}+iexp*100];
        end
end
for ipeak = 1:2
    switch whichExp
        case {1,2,3}
            sel=bT.Pot==ipeak & bT.ExpN==whichExp & ismember(bT.Subject,allwsp);
        case 'all'
            sel=bT.Pot==ipeak & ismember(bT.Subject,allwsp);
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
                        %TODO - check what happend with the squeeze!!!
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
            suRA=suRA(~art & seqInd > numt);
            tt=table(subj,csubj,suRA,vv);
            switch modelName
                case 'lm'
                    formula = 'vv~suRA*csubj';
                    tmplme = fitlm(tt,formula);
                case 'lme'
                    formula = 'vv~suRA+(suRA|subj)';
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

