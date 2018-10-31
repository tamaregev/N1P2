function [ grandMatrix, ntrials ] = calcGrandMatrix_N1P2(Expinfo, whichSubjects, conditions, EventCodes, nelectrodes, win, bsl, Bp  )
%
% Sep. 27 2016 - Tamar - adapted from calcGrandMatrixAll of MMNchroma
tictot = tic;

ExpName = Expinfo.ExpName;
ExportFolder = Expinfo.ExportFolder;
AnalysisFolder = Expinfo.AnalysisFolder;
Subjects = Expinfo.Subjects;
sessions = Expinfo.sessions;
srate=Expinfo.srate;
low = Bp(1); high = Bp(2);
grandMatrix = nan(size(win,2), max(whichSubjects), length(conditions),nelectrodes);
ntrials = nan(max(whichSubjects),length(conditions),2);

    for s = whichSubjects
        ticsubj = tic;
        Subj = Subjects{s}
        disp(['Processing ' Subj '...'])

        FileName = [ExpName '_' Subj '_' sessions{s}(1) '_dt_ImportMarkers_ReCoded'];
        %load and read markers
        if exist([ExportFolder FileName '.vmrk'], 'file') == 2
            VMRKfile = [ExportFolder FileName '.vmrk'];
        else
            VMRKfile = [AnalysisFolder Subj filesep FileName '.vmrk'];
        end
        [eventsCodesIndexes, artifactsIndexes]=read_markers_artifacts(VMRKfile,15);%check that this is the row of Mk1 indeed
        %load data
        disp('loading dat file...');tic;
        if exist([ExportFolder FileName '.mat'], 'file') == 2
            load([ExportFolder FileName '.mat']);disp(['Done loading in ' num2str(toc) ' sec.'])
        else
            load([AnalysisFolder Subj filesep FileName '.mat']);disp(['Done loading in ' num2str(toc) ' sec.'])
        end
        
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

        for con=1:length(conditions)
            cond = conditions{con};
            %disp(['Condition # ' num2str(con) ': ' cond])
            trigs = eventsCodesIndexes(ismember(eventsCodesIndexes(:,1),relevantTrigs{con}),2); 
            if trigs
                [AVG, SEGS, REJ] = SegAndAvg(allData(:,:),trigs,win,'reject',artifactsIndexes);
                if bsl
                    AVG_bc = bsxfun(@minus,AVG,mean(AVG(bsl,:),1));
                else
                    AVG_bc = AVG;
                end
                clear AVG
                %assign into grandMatrix:
                grandMatrix(:,s,con,:) = AVG_bc;
                ntrials(s,con,1) = size(SEGS,3);
                ntrials(s,con,2) = numel(REJ);
            else
                grandMatrix(:,s,con,:) = nan;
            end
        end
        disp(['Done subject ' Subj ' in ' num2str(toc(ticsubj)) ' sec.'])
    end
    
    disp(['Done all in ' num2str(toc(tictot)) ' sec.' ])

end

