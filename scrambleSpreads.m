function [ llp ] = scrambleSpreads(ip, scrambles, whichExp, bRA, bT, whichSubjects, spread2block, whichpeaks, sigmas, taus, numt, modelName, infoFlag  )
%SCRAMBLESPREADS facilitates parfor in Section 11 of lmeModels_allExp
% spreads were manipulated only in exp 3 so all input variables are already
% specific to Exp 3 and whichExp = 3
% scrambling lables per spread will show whether and to what extend we get
% sigmas modulated by the spread even when the model and data don't fit.
% Hopefully we could show that the modulation of best sigmas for the real
% data is above and beyond the modulation we would get here.
disp('I''m in')
llp=nan(length(sigmas),length(taus),3,2);

for ispread=1:3
    
    ticspread = tic;
    % prep uRA as normal (first perm)
    uRA=bRA(whichSubjects,spread2block{ispread},:,:,:);
    for ipeak = 1:2
        disp(['(ispread, ipeak) = ' num2str([ispread ipeak])]) 
        %pT is a table for only this peak, experiment, and blocks (spreads)
        pT=bT(bT.Pot==ipeak & ismember(bT.Block,spread2block{ispread}),:);
        %nT is the new scrambled table
        nT=pT(scrambles{ip,ispread},:);
        
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
