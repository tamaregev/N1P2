function [ llst, infost ] = calcLL_scrambled( ip, scrambles, taus, sigmas, whichExp, bwhichSubjects, bRA, bT, numt, modelName, infoFlag)
%CALCLLSIGMAS facilitates converting the taus-for-loop to parfor  
llst=nan(length(sigmas),length(taus),2);
if infoFlag
    infost = cell(length(sigmas),length(taus),2);
else
    infost=nan;
end
scramble=scrambles{ip};

for ipeak = 1:2
    switch whichExp
        case {1,2,3}
            sel=bT.Pot==ipeak & bT.ExpN==whichExp;
            uRA=bRA{whichExp}(bwhichSubjects{whichExp},:,:,:,:);
        case 'all'
            sel=bT.Pot==ipeak;
            clear uRA
            for iexp=1:3
                uRA{iexp}=bRA{iexp}(bwhichSubjects{iexp},:,:,:,:);
            end
    end
    vv=bT.val(sel);
    art=bT.Art(sel);
    subj=bT.Subject(sel);
    seqInd=bT.seqInd(sel);
    subj=bT.Subject(sel);
    
    vv=vv(~art & seqInd > numt);
    subj=subj(~art & seqInd > numt);
    csubj=categorical(subj);
    
    for isigma = 1:length(sigmas)
        ticsig = tic;
        for itau = 1:length(taus)
            %disp([isigma itau]);
            switch whichExp
                case {1,2,3}
                    suRA=(uRA(:,:,:,isigma,itau));
                    suRA=suRA(:); % this assumes that the table is arranged as above
                    suRA(isnan(suRA))=[];%the same way as when created bT
                case 'all'
                    suRA = [];
                    for iexp=1:3
                        tmpRA=squeeze(uRA{iexp}(:,:,:,isigma,itau));
                        tmpRA=tmpRA(:);
                        suRA=[suRA; tmpRA]; 
                    end
                    suRA(isnan(suRA))=[];%the same way as when created bT
            end
            suRA=suRA(scramble);
            suRA=suRA(~art & seqInd > numt);
            tt=table(subj,csubj,suRA,vv);
            
            switch modelName
                case 'lme'
                    tmplme=fitlme(tt,'vv~suRA+(suRA|subj)');
                case 'lm'
                    tmplme=fitlm(tt,'vv~suRA*csubj');
            end
            llst(isigma,itau,ipeak)=tmplme.LogLikelihood;
            if infoFlag
                info = extractInfolme(tmplme);
                infost{isigma,itau,ipeak} = info;
            end
        end
        disp(['Done isigma = ' num2str(isigma) ' in ' num2str(toc(ticsig))])
    end
end

end

