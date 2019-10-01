function [ singleAmps, isArtefact ] = calcAmpSingleTrials( ExpN )
%load continuous data, calculate voltage of N1 and P2 due to peak times of grand averages.
% The latter were computed in the past at: 
% Exp 3: Amplitudes_N1P2_GH

%params::
%general:
electrodeName = 'Cz';
load('L:\Experiments\MMNchroma\Analysis\Properties')
[ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );

addtag = '';
conditions = {'note 1','note 2','note 3','note 4','note 5'};
dt = 5;%ms before and after peak - window for calculating N1-P2 amp
srate = 512;
ds = round(srate*dt/1000);
ERPwin = -100:400; %for artefacts
ERPsamps = round(srate*ERPwin/1000);

%exp specific:
switch ExpN
    case 1
        addpath(['S:\Lab-Shared\Experiments\MMNchroma\Analysis'])
        Definitions_MMNchroma
        peakdate = '18-Nov-2018';
        
        nodeName = 'ImportMarkers_ReCoded';
        %load example expdata to calculate nTrialsPerBlock
        s=whichSubjects(1);
        load([EDATfolder ExpName '_' Subjects{s} '_' sessions{s}(1) '_expdata.mat' ])
        phaseName = 'stimParamsBlocksMMN';
        nTrialsPerBlock = expdata.(phaseName)(1).num_trials;
        block_list = [expdata.phases(1).block_list expdata.phases(3).block_list];
        nBlocksPerType = sum(block_list==block_list(1));
        STIMcodes = [20,30,11,40,50];
        order = cs;

    case 2
        addpath(['S:\Lab-Shared\Experiments\MMNchromaF\Analysis'])
        Definitions_MMNchromaF
        peakdate = '16-Nov-2018';
        nodeName = 'ImportMarkers_ReCoded';
        %load example expdata to calculate nTrialsPerBlock
        s=whichSubjects(1);
        load([EDATfolder ExpName '_' Subjects{s} '_' sessions{s}(1) '_expdata.mat' ])
        phaseName = 'MMN';
        nTrialsPerBlock = expdata.(phaseName).blocks(1).num_trials;
        block_list = expdata.(phaseName).block_list;
        nBlocksPerType = sum(block_list==block_list(1));
        
        STIMcodes = [20,30,11,40,50];
        order = cs;
        
        %codes and names table:
        names = {'Training_pure','Training_shep','MMN_F_pure','MMN_Fctrl_pure','MMN_classic_pure','MMN_F_shep','MMN_Fctrl_shep','MMN_classic_shep','Perception_pure','Perception_shep','MMN_classic_original_pure','MMN_classic_original_shep'}';
        PBcodes = [1200, 1400, 2300, 2400, 2500, 2800, 2900, 21000, 3200, 3400, 21100, 21200]';
        codes_names = table(names,PBcodes);

    case 3
        Definitions_N1P2
        
        %load data
        peakdate = '30-Oct-2018';
        nodeName = 'RDI_imported';
        STIMcodes = [1,2,3,4,5];
        phaseName = 'Passive';
        %load example expdata to calculate nTrialsPerBlock
        s=whichSubjects(1);
        load([EDATfolder ExpName '_' Subjects{s} '_' sessions{s}(1) '_expdata.mat' ])
        nTrialsPerBlock = expdata.(phaseName).blocks(1).num_trials;
        block_list = expdata.(phaseName).block_list;
        nBlocksPerType = sum(block_list==block_list(1));
end
Bp = [1 20];
low = Bp(1); high = Bp(2);
%load peaks 
load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate],'allGrandcon_times')

tictot = tic;
%for each subject, load raw data, and calc amp due to time of grand:
for s=whichSubjects
    %load raw data
    ticsubj = tic;
    Subj = Subjects{s};
    disp(Subj);

    FileName = [ExpName '_' Subj '_' sessions{s}(1) '_dt_' nodeName];
    %load and read markers
    if exist([ExportFolder FileName '.vmrk'], 'file') == 2
        VMRKfile = [ExportFolder FileName '.vmrk'];
    else
        VMRKfile = [AnalysisFolder Subj filesep FileName '.vmrk'];
    end
    [eventsCodesIndexes, artifactsIndexes]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed
    eventsCodesIndexes = eventsCodesIndexes(2:end,:);

    %load data
    %disp('loading dat file...');tic;
    if exist([ExportFolder FileName '.mat'], 'file') == 2
        load([ExportFolder FileName '.mat']);%disp(['Done loading in ' num2str(toc) ' sec.'])
    else
        load([AnalysisFolder Subj filesep FileName '.mat']);%disp(['Done loading in ' num2str(toc) ' sec.'])
    end
    
    %Band Pass
    if ~isnan(low) && ~isnan(high)
         allData = bandPassFilter(low,high,allData,srate);
    elseif ~isnan(low)
        allData = HPF(allData,srate,low);
    elseif ~isnan(high)
        allData = LPF(allData,srate,high);
    end

    %Relevant electrode
    data = allData(:,electrode);
    
    %init singleAmps:
    if s==whichSubjects(1)
        singleAmps = nan(whichSubjects(end),length(bls),nTrialsPerBlock*nBlocksPerType,2);%subjs x blockTypes x Trials x peaks
        isArtefact = nan(whichSubjects(end),length(bls),nTrialsPerBlock*nBlocksPerType);%subjs x blockTypes x Trials
        seqInds = nan(whichSubjects(end),length(bls),nTrialsPerBlock*nBlocksPerType);%subjs x blockTypes x Trials
        stimCodes = nan(whichSubjects(end),length(bls),nTrialsPerBlock*nBlocksPerType);%subjs x blockTypes x Trials
    end
    ib=0;
    for bl = bls
        ib=ib+1;
        blockName = blocks{bl};
        disp(blockName)
        %calc codes:
        switch ExpN
            case 1
               switch blockName 
                    case 'chroma'
                        EventCodes = {3120,3130,3140,3150,3111};
                    case 'chromactrl'
                        EventCodes = {3220,3230,3240,3250,3211};
                    case 'classic'
                        EventCodes = {[1320, 3320],[1330, 3330],[1340, 3340],[1350, 3350] ,[1311, 3311]};
                    case 'classicctrl'
                        EventCodes = {[1420, 3420],[1430, 3430],[1440, 3440],[1450, 3450],[1411, 3411]};
               end
                EventCodes = EventCodes(order);
            case 2
        
                PBcode = codes_names.PBcodes(strcmp(codes_names.names,blockName));
                EventCodes=cell(1,length(conditions));
                for i=1:length(conditions)
                    EventCodes{i} = PBcode + STIMcodes(i);
                end

            case 3
                Bcode = Bcodes_names.Bcodes(strcmp(Bcodes_names.Bnames,blockName));
                EventCodes=cell(1,length(conditions));
                for i=1:length(conditions)
                    EventCodes{i} = Bcode + STIMcodes(i);
                end

        end

        
        %
        %assign EventCodes to relevantTrigs
        relevantTrigs = EventCodes;%
        
        allBCodes = [];%only of current block type... 
        for ic=1:length(relevantTrigs)
            allBCodes = [allBCodes relevantTrigs{ic}];
        end
        blockCodesIndexes = eventsCodesIndexes(ismember(eventsCodesIndexes(:,1),allBCodes),:);
        
        for ipeak = 1:2
            Ptimes = allGrandcon_times{bl}(s,:,ipeak);
            
            for con=1:length(conditions)
                cond = conditions{con};
                Ptime = Ptimes(con);
                Psamp = round(srate*Ptime/1000);
                Psamps = (Psamp - ds) : (Psamp + ds);
                
                trigs = blockCodesIndexes(ismember(blockCodesIndexes(:,1),relevantTrigs{con}) ,2); 
                idx = find(ismember(blockCodesIndexes(:,1),relevantTrigs{con}));%serial index in sequence
       
                Plocs = trigs' + Psamps';%locations of peaks in experiment samples 
                ERPlocs = trigs' + ERPsamps';%ERPs in experiment samples
                
                wlen = length(Psamps);
                if length(idx)*5 ==  nTrialsPerBlock*nBlocksPerType
                else
                    warning([Subj ': number of found triggers is wrong'])
                end
                if trigs
                    %[AVG, SEGS, REJ] = SegAndAvg(allData(:,:),trigs,win,'reject',artifactsIndexes);
                    wdata = reshape(data(Plocs(:)),[wlen,length(idx)]);
                    singleAmps(s,ib,idx,ipeak) = mean(wdata);%subjs x blockTypes x Trials x peaks
                    artFlag = ismember(ERPlocs,artifactsIndexes);
                    isArtefact(s,ib,idx) = any(artFlag,1);
                end
            end
        end    
    end
    disp(['Done in ' num2str(toc(ticsubj))])
end
disp(['Done all in ' num2str(toc(tictot))])

end


