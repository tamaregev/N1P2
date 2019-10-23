function [ll, infos] = N1P2permutation(perm,uRA,sigmas,taus,tt,art,seqInd,numt,otherpeak,modelName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 nuRA=uRA;
 ll=nan(length(sigmas),length(taus),2);
 infos = cell(length(sigmas),length(taus),2);
    for isig = 1:length(sigmas)
        ticsig = tic;
        for itau = 1:length(taus)
            %disp(itau)
            ttt=tt;
            op=otherpeak;
            
            suRA=squeeze(nuRA(:,:,:,isig,itau));
            suRA=suRA(:);  
            
           
            
            
            
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
                info = extractInfolme(tmplme);
                %nplme{ip,isig,itau,ipeak} = tmplme;
                ll(isig,itau,ipeak)=tmplme.LogLikelihood;
                infos{isig,itau,ipeak} = info;
            end
        end
        disp(['Done isig = ' num2str(isig) ' in ' num2str(toc(ticsig))])
    end
    
end

