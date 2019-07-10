function [ RA, smpls, stimCode, seqInd ] = calcRA_MMNchroma(  R0, sigma, tau, expdata, eventCoInd, phaseName, SOA_threshold, Mis, bls, order)
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
%only relevant blocks:
if length(bls)>1
    %MMNchroma, Exp 1
    relBlockTypes = bls;
    phases = [1,3];
    block_list = [];
    for ip=1:length(phases)
        block_list = [block_list, expdata.phases(phases(ip)).block_list];%all experiment
    end
    nBlocks = length(block_list);
else
    %MMNchromaF, Exp 2
    relBlockTypes = 4;
    phases = 2;
    block_list = expdata.(phaseName).block_list;%all experiment
    nBlocks = length(block_list);
end
blocks = expdata.(phaseName);
nBlockTypes = length(relBlockTypes);
trials = expdata.trials(ismember([expdata.trials(:).phase_number],phases));
stimCoInd = eventCoInd(~ismember(eventCoInd(:,1),[210,200,100,110,254,255]),:);

newbli = find(eventCoInd(:,1)==100);%new block index in eventCoInd
if(~(length(trials)==length(stimCoInd(ismember(floor(stimCoInd(:,1)/1000),phases),1))))
    warning('nTrials is not equal in edat and events indexes')
end
if(~(12==nBlocks))
    warning('nBlocks is not equal to 12')
end

nTrialsPerBlock = blocks(1).num_trials;
%this is for MMNchroma, fix code notes from standards and deviant
CODE_STANDARDS = blocks(1).CODE_STANDARDS;
CODE_DEVIANT = blocks(1).CODE_DEVIANT;

CODE_NOTES = [CODE_STANDARDS, CODE_DEVIANT];
CODE_NOTES = CODE_NOTES(order);

nStim = length(CODE_NOTES);
whichBlocks = find(block_list==relBlockTypes(1));
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
    bt = relBlockTypes(ibt);
    MIDIs = [freq2MIDI([blocks(bt).toneSynth.Fstandards{:}]) , freq2MIDI([blocks(bt).toneSynth.Fdeviants{:}])];
    MIDIs = MIDIs(order);%MIDIs of all stimuli
    if isnanMis
        Mis = nan(length(MIDIs),1);
        for im = 1:length(MIDIs)
            Mis(im) = MIDIs(im);
        end
    end
    whichBlocks = find(block_list==bt);        
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
                    note_code = trials(itrial-1).trial_code;%current stimulus code
                    Mprev = MIDIs(CODE_NOTES==note_code);%current stimulus MIDI
                    RA(ibt,ib,:,j) = (prevRA + (1-prevRA) .* tuningCurve(Mis,Mprev,sigma)) .* exp(-SOA/tau);   
                end
            end
            seqInd(ibt,ib,j) = rj;
        end
    end
end

function tuning = tuningCurve(Mi,Mj,sigma)
    tuning = exp(-0.5*(((Mi-Mj)/sigma).^2));

