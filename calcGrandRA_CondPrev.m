function [ grandRA ] = calcGrandRA_CondPrev( bl, whichSubjects, RAs_SigTau, EventCodes, prevEventCodes, stimCodess, smplss, seqIndss, artIndss)
%calcGrandRA_CondPrev averages RA for all (current, previous) note pairs
%   based on calcGrandMatrix_condPrev_N1P2.m
%   
% Tamar, Dec 9 2018 

%TODO
% %assign EventCodes to relevantTrigs
% relevantTrigs = EventCodes;%historically built for correcting in some subjects
% 
% % previous trigs
% previousTrigs = prevEventCodes;

%is5first = squeeze(seqInd(it,is,:))<=5;
            
grandRA = nan(whichSubjects(end),1,size(RAs_SigTau,3),length(EventCodes),length(prevEventCodes),size(RAs_SigTau,6),size(RAs_SigTau,7));

for seq = 1:size(RAs_SigTau,3)
    for s=whichSubjects
        stims_prev = stimCodess(s,bl,seq,1:end-1);
        stims = stimCodess(s,bl,seq,2:end);    
        for con=1:length(EventCodes)
            for prevcon = 1:length(prevEventCodes)
                stimsInd_prev = stims_prev == prevEventCodes{prevcon};
                stimInd = stims == EventCodes{con};
                if any(stimsInd_prev & stimInd)
                    grandRA(s,1,seq,con,prevcon,:,:) = mean(RAs_SigTau(s,bl,seq,con,stimInd & stimsInd_prev,:,:),5);
                end
            end
        end
    end
end

%average across the 3 sequences:
grandRA_cat = [];
for is=1:size(RAs_SigTau,3)
    grandRA_cat = cat(8,grandRA_cat,grandRA(:,:,is,:,:,:,:));
end

grandRA = (mean(grandRA_cat,8));

