function [ grandMatrix, ntrials ] = calcGrandMatrix_condPrev_N1P2(Expinfo, whichSubjects, conditions, prevConditions, EventCodes, prevEventCodes, nelectrodes, win, bsl, Bp  )
%12/3/2018 - Tamar, adapted from calcGrandMatrixAll_condPrevAllComb for
%N1P2 experiment
%added Bp as in calcGrandMatrix_N1P2

tictot = tic;

ExpName = Expinfo.ExpName;
ExportFolder = Expinfo.ExportFolder;
AnalysisFolder = Expinfo.AnalysisFolder;
Subjects = Expinfo.Subjects;
sessions = Expinfo.sessions;
srate=Expinfo.srate;
low = Bp(1); high = Bp(2);
grandMatrix = nan(size(win,2), max(whichSubjects), length(conditions), length(prevConditions), nelectrodes);
%time x subjects x (current)condition x prevConditions x electrodes
ntrials = nan(max(whichSubjects),length(conditions),length(prevConditions),2);
nodeName = 'RDI_imported';

    for s = whichSubjects
        ticsubj = tic;
        Subj = Subjects{s};
        disp(Subj);
        %disp(['Processing ' Subj '...'])

        FileName = [ExpName '_' Subj '_' sessions{s}(1) '_dt_' nodeName];
        %load and read markers
        if exist([ExportFolder FileName '.vmrk'], 'file') == 2
            VMRKfile = [ExportFolder FileName '.vmrk'];
        else
            VMRKfile = [AnalysisFolder Subj filesep FileName '.vmrk'];
        end
        [eventsCodesIndexes, artifactsIndexes]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed
        eventsCodesIndexes_prev = eventsCodesIndexes(1:end-1,:);
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
%             allData = HPF(allData,srate,low);
%             allData = LPF(allData,srate,high);
        elseif ~isnan(low)
            allData = HPF(allData,srate,low);
        elseif ~isnan(high)
            allData = LPF(allData,srate,high);
        end
        
        %assign EventCodes to relevantTrigs
        relevantTrigs = EventCodes;%historically built for correcting in some subjects
        
        % previous trigs
        previousTrigs = prevEventCodes;
         
        for con=1:length(conditions)
            cond = conditions{con};
            for prevcon = 1:length(prevConditions)
                prevcond = conditions{prevcon};
                
                trigs = eventsCodesIndexes(ismember(eventsCodesIndexes(:,1),relevantTrigs{con}) & ismember(eventsCodesIndexes_prev(:,1),previousTrigs{prevcon}) ,2); 

                if trigs
                    [AVG, SEGS, REJ] = SegAndAvg(allData(:,:),trigs,win,'reject',artifactsIndexes);
                    if bsl
                        AVG_bc = bsxfun(@minus,AVG,mean(AVG(bsl,:),1));
                    else
                        AVG_bc = AVG;
                    end
                    clear AVG
                    %assign into grandMatrix:
                    grandMatrix(:,s,con,prevcon,:) = AVG_bc;
                    ntrials(s,con,prevcon,1) = size(SEGS,3);
                    ntrials(s,con,prevcon,2) = numel(REJ);
                else
                    grandMatrix(:,s,con,prevcon,:) = nan;
                end
            end
        end
        %disp(['Done subject ' Subj ' in ' num2str(toc(ticsubj)) ' sec.'])
    end
    
    disp(['Done all in ' num2str(toc(tictot)) ' sec.' ])

end
% matrixdate - 19-Mar-2018 Bp1-20 Bl1 - Bl3a


