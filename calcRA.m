function [ RA, smpls, stimCode, seqInd ] = calcRA(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold, Mis)
%CALCRA calculates the expected RA (response adaptation) as in Herrmann
%et al. 2014
% 
% RA;  %types x sequences x channels x timepoints    
% is calculated for all blocks, per trial for a single subjects, 
% based on the specific train of stimuli presented during the experiment, 
% as documented in expdata.
% 
% Inputs:
%  phaseName = 'Passive';
%  sigma = 10;%MIDI
%  tau = 1;%seconds
%  SOA_threshold = 0.6;%seconds
% 
%   Written by Tamar Regev, Nov 2018
% 
block_list = expdata.(phaseName).block_list;
nBlocks = length(block_list);
blocks = expdata.(phaseName).blocks;
nBlockTypes = length(blocks);
trials = expdata.trials(strcmp({expdata.trials(:).phase_type},phaseName));
stimCoInd = eventCoInd(~ismember(eventCoInd(:,1),[210,200,100,110,254,255]),:);
newbli = find(eventCoInd(:,1)==100);%new block index in eventCoInd
if(~length(trials)==length(stimCoInd))
    warning('nTrials is not equal in edat and events indexes')
end
if(~length(newbli)==nBlocks)
    warning('nBlocks is not equal in edat and events indexes')
end

nTrialsPerBlock = expdata.(phaseName).blocks(1).num_trials;
nStim = length(expdata.(phaseName).blocks(1).CODE_NOTES);
whichBlocks = find(block_list==1);
if ~isnan(Mis)
    nChans = length(Mis);
    isnanMis = false;
else%Mis is nan
    isnanMis = true;
    nChans = nStim;
end
RA = nan(nBlockTypes,length(whichBlocks),nChans,nTrialsPerBlock);%types x sequences x channels x timepoints    
smpls = nan(nBlockTypes,length(whichBlocks),nTrialsPerBlock);%types x sequences x timepoints
stimCode = nan(nBlockTypes,length(whichBlocks),nTrialsPerBlock);%types x sequences x timepoints
seqInd = nan(nBlockTypes,length(whichBlocks),nTrialsPerBlock);%types x sequences x timepoints

for ibt = 1:nBlockTypes
    MIDIs = expdata.(phaseName).blocks(ibt).MIDIs;%MIDIs of all stimuli
    if isnanMis
        Mis = nan(length(MIDIs),1);
        for im = 1:length(MIDIs)
            Mis(im) = MIDIs{im};
        end
    end
    whichBlocks = find(block_list==ibt);        
    ib = 0;
    for bl = whichBlocks
        ib=ib+1;
        rj = 0;% 'real' j index - will reset to 1 when there is a large SOA and RA resent to R0
        for j = 1:nTrialsPerBlock
            %find current trial data:
            itrial = (bl-1)*nTrialsPerBlock + j;
            smpls(ibt,ib,j) = stimCoInd(trials(itrial).trial_number,2);
            stimCode(ibt,ib,j) = stimCoInd(trials(itrial).trial_number,1);
            if j==1
                rj = 1;
                RA(ibt,ib,:,j) = repmat(R0,nChans,1);
            else
                SOA = trials(itrial).time - trials(itrial-1).time;
                if SOA > SOA_threshold
                    rj = 1;
                    RA(ibt,ib,:,j) = repmat(R0,nChans,1);
                else
                    rj = rj + 1;
                    prevRA = squeeze(RA(ibt,ib,:,j-1));
                    note_code = trials(itrial-1).note_code;%current stimulus code
                    Mprev = expdata.(phaseName).blocks(ibt).MIDIs{expdata.(phaseName).blocks(ibt).CODE_NOTES==note_code};%current stimulus MIDI
                    RA(ibt,ib,:,j) = (prevRA + (1-prevRA) .* tuningCurve(Mis,Mprev,sigma)) .* exp(-SOA/tau);   
                end
            end
            seqInd(ibt,ib,j) = rj;
        end
    end
end

function tuning = tuningCurve(Mi,Mj,sigma)
    tuning = exp(-0.5*(((Mi-Mj)/sigma).^2));

