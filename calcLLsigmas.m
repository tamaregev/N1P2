function [ lls, info ] = calcLLsigmas( it, taus, sigmas, whichExp, uRA, art, seqInd, numt, vv, subj, csubj, modelName, infoFlag)
%CALCLLSIGMAS facilitates converting the taus-for-loop to parfor  
lls=nan(length(sigmas),1);
info=cell(length(sigmas),1);

for isigma=1:length(sigmas)
    %disp([isigma,it])

    switch whichExp
        case {1,2,3}
            suRA=(uRA(:,:,:,isigma,it));
            suRA=suRA(:); % this assumes that the table is arranged as above
            suRA(isnan(suRA))=[];%the same way as when created bT
            if 0
                %another way to calculate seqInd: (verifies that suRA(:) is correct
                load([Folder filesep 'Exp' num2str(whichExp) filesep 'seqIndcat.mat'])   
                useqIndcat=squeeze(seqIndcat(bwhichSubjects{whichExp},:,:));
                useqIndcat(isnan(useqIndcat))=[];
                figure;plot(seqInd);hold all;plot(useqIndcat)
            end
        case 'all'
            suRA = [];
            for iexp=1:3
                tmpRA=squeeze(uRA{iexp}(:,:,:,isigma,it));
                tmpRA=tmpRA(:);
                suRA=[suRA; tmpRA]; 
            end
            suRA(isnan(suRA))=[];%the same way as when created bT
            if 0
                %another way to calculate seqInd: (verifies that suRA(:) is correct
                useqIndcatall=[];
                for iexp = 1:3
                    load([Folder filesep 'Exp' num2str(iexp) filesep 'seqIndcat.mat'])   
                     useqIndcat=squeeze(seqIndcat(bwhichSubjects{iexp},:,:));
                     useqIndcat=useqIndcat(:);
                     useqIndcatall=[useqIndcatall;useqIndcat];
                end
                useqIndcatall(isnan(useqIndcatall))=[];
                figure;plot(seqInd);hold all;plot(useqIndcatall)
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
    lls(isigma)=tmplme.LogLikelihood;
    if infoFlag
        info{isigma} = extractInfolme(tmplme,0);
    end
end
end

