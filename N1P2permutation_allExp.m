function [ll, infos] = N1P2permutation_allExp(whichExp,perm,uRA,sigmas,taus,tt,art,seqInd,numt,otherpeak,modelName,infoFlag)
%N1P2PERMUTATION_ALLEXP performs one iteration of the permutation test
% to test for significance of difference between parameters of N1 and P2
 
%July 25 Tamar: Fixed the bug, added suRA(isnan(suRA))=[];

ll=nan(length(sigmas),length(taus),2);
if infoFlag
    infos = cell(length(sigmas),length(taus),2);
else
    infos=nan;
end
    for isig = 1:length(sigmas)
        ticsig = tic;
        for itau = 1:length(taus)
            %disp(itau)
            ttt=tt;
            op=otherpeak;
           
            switch whichExp
                case {1,2,3}
                    suRA=squeeze(uRA(:,:,:,isig,itau));
                    suRA=suRA(:); % this assumes that the table is arranged as above
                    suRA(isnan(suRA))=[];
                case 'all'
                    suRA = [];
                    for iexp=1:3
                        tmpRA=squeeze(uRA{iexp}(:,:,:,isig,itau));
                        tmpRA=tmpRA(:);
                        suRA=[suRA; tmpRA]; 
                    end
                    suRA(isnan(suRA))=[];
            end
            
            suRA=suRA(~art & seqInd > numt);
            for ipeak = 1:2
                btt=ttt{ipeak};
                btt.suRA=suRA;
                btt.vv(perm==1) = ttt{op(ipeak)}.vv(perm==1);
                
                switch modelName
                    case 'lme'
                        tmplme=fitlme(btt,'vv~suRA+(suRA|subj)');
                    case 'lm'
                        tmplme=fitlm(btt,'vv~suRA*csubj');
                end
                if infoFlag
                    info = extractInfolme(tmplme);
                    infos{isig,itau,ipeak} = info;
                end
                %nplme{ip,isig,itau,ipeak} = tmplme;
                ll(isig,itau,ipeak)=tmplme.LogLikelihood;
                
            end
        end
        disp(['Done isig = ' num2str(isig) ' in ' num2str(toc(ticsig))])
    end
    
end

