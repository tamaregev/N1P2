function [ llp ] = bootstrapSpreads(ip, whichExp, bRA, bT, wsp, spread2block, whichpeaks, sigmas, taus, numt, modelName, infoFlag  )
%BOOTSTRAPSPREADS facilitates parfor in Section 10 of lmeModels_allExp
% spreads were manipulated only in exp 3 so all input variables are already
% specific to Exp 3 and whichExp = 3
    
llp=nan(length(sigmas),length(taus),3,2);
allwsp=wsp{ip}+whichExp*100;

for ispread=1:3
    
    ticspread = tic;
    % prep uRA as normal (first perm)
    uRA=bRA(wsp{1},spread2block{ispread},:,:,:);
    for ipeak = 1:2
        %pT is a table for only this peak, experiment, and blocks (spreads)
        pT=bT(bT.Pot==ipeak & ismember(bT.Block,spread2block{ispread}),:);
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

        llst=nan(length(sigmas),length(taus));
        for it=1:length(taus)
%             info=cell(length(sigmas),1);

            for isigma=1:length(sigmas)
                %disp([isigma,it])
                suRA=(uRA(:,:,:,isigma,it));
                suRA=suRA(:); % this assumes that the table is arranged as above
                suRA(isnan(suRA))=[];%the same way as when created bT
        
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
                llst(isigma,it)=tmplme.LogLikelihood;
                if infoFlag
                    info{isigma} = extractInfolme(tmplme,0);
                end
            end
        end
        llp(:,:,ispread,ipeak)=llst;
    end
    disp(['Done perm = ' num2str(ip) ' , spread ' num2str(ispread) ' in ' num2str(toc(ticspread))])
end

end
