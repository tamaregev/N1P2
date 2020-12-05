%MasterScript_N1P2
%
%14/1/2018 - Tamar. Adapted from MasterScript_Ushape
% Preprocessing of data from experiment N1P2

%% Definitions:
addpath('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
%% BDF Detrend

for i=31:37
    tic
    Subj = Subjects{i}
    session = sessions{i};
    BDFfolder = ['L' ResultsFolder(2:end) 'EEG\BDF_before_dtrend\'];
    
    multipleSes = false;
    if strcmp(session, '')
        FileName = [ExpName '_' Subj ];
    elseif length(session) == 1
        FileName = [ExpName '_' Subj '_' session];
    elseif length(session) > 1
        multipleSes = true;
    end
    
    disp(['Detrending ' Subj ' ...'])
    if multipleSes
        for s=1:length(session)
            ses=session(s);
            FileName = [ExpName '_' Subj '_' session(s)]; 
            detrendBDFblocks([BDFfolder FileName '.bdf']);
            movefile([cd filesep FileName '_dt.dat'],[RawDataFolder FileName  '_dt.dat']);
            movefile([cd filesep FileName '_dt.vhdr'],[RawDataFolder FileName '_dt.vhdr']);
            movefile([cd filesep FileName '_dt.vmrk'],[RawDataFolder FileName '_dt.vmrk']);
        end
    else
        detrendBDFblocks([BDFfolder FileName '.bdf']);
        movefile([cd filesep FileName '_dt.dat'],['L' RawDataFolder(2:end) FileName  '_dt.dat']);
        movefile([cd filesep FileName '_dt.vhdr'],['L' RawDataFolder(2:end) FileName '_dt.vhdr']);
        movefile([cd filesep FileName '_dt.vmrk'],['L' RawDataFolder(2:end) FileName '_dt.vmrk']);
    end
    disp(['Done dtrending ' Subj 'in ' num2str(toc) ' sec.'])
end
% done until 130
% 
%% Go to - Analyzer
% until generic data export

%% move exported files to relevant folder
% after importing recoded markers to analyzer and generic data export
ExportFolderNew = ['L' AnalyzerExportFolder(2:end) 'beforeSegmentations\'];
lastnode = 'RDI_imported';
for s = 7
    %[7,18:19,21:28,30:33,35:37]
    tic
    Subj = Subjects{s}
    disp([ Subj '...'])
    if strcmp(sessions{s},'')
        FileName = [ExpName '_' Subj '_dt_' lastnode];
    else
        FileName = [ExpName '_' Subj '_' sessions{s}(1) '_dt_' lastnode];
    end
    datFile = [ FileName '.dat'];
    vmrkFile = [ FileName '.vmrk'];
    vhdrFile = [ FileName '.vhdr'];
    movefile(['L' AnalyzerExportFolder(2:end) datFile],[ExportFolderNew datFile]);
    movefile(['L' AnalyzerExportFolder(2:end) vmrkFile],[ExportFolderNew vmrkFile]);
    movefile(['L' AnalyzerExportFolder(2:end) vhdrFile],[ExportFolderNew vhdrFile]);
    disp('Done')
end
%Done until 105, 8-17

%% copy files to D:
% copy data files
for s = 7
    %[7,18:19,21:28,30:33,35:37]
    Subj = Subjects{s};
    disp(['Copying ' Subj '...'])
    sourceFolder = ['L:\Experiments\' ExpName '\Experiment\Results\EEG\export\beforeSegmentations\'];
    destFolder = ['D:\Experiments\' ExpName '\Experiment\Results\EEG\export\beforeSegmentations\'];
    mkdir(destFolder);
    fileName = [ExpName '_' Subj '_' sessions{s}(1) '_dt_' lastnode];
    source ={[sourceFolder fileName '.dat'],[sourceFolder fileName '.vhdr'],[sourceFolder fileName '.vmrk']};
    destination = {[destFolder fileName '.dat'],[destFolder fileName '.vhdr'],[destFolder fileName '.vmrk']};
    for i=1:length(source)
        copyfile(source{i},destination{i})
    end
    % EDAT
    sourceFolder = ['L:\Experiments\' ExpName '\Experiment\Results\EDAT\'];
    destFolder = ['D:\Experiments\' ExpName '\Experiment\Results\EDAT\'];
    mkdir(destFolder);
    fileName = [ExpName '_' Subj '_' sessions{s}(1) '_expdata'];
    source =[sourceFolder fileName '.mat'];
    destination = [destFolder fileName '.mat'];
    copyfile(source,destination)
    disp(['Done Subj ' Subj])
end
%Done until 105, 8-17
%% read data
%do this in D and copy to L
for s = [7,18:19,21:28,30:33,35:37]
    tic
    Subj = Subjects{s};
    mkdir(['D' AnalysisFolder(2:end) Subj])
    disp(['Loading ' Subj '...'])
    if strcmp(sessions{s},'')
        FileName = [ExpName '_' Subj '_dt_' lastnode];
    else
        FileName = [ExpName '_' Subj '_' sessions{s}(1) '_dt_' lastnode];
    end
    datFile = ['D' ExportFolder(2:end) FileName '.dat'];
    [allData] = read_analyzer_UnsegmentedData(datFile);
    save(['D' AnalysisFolder(2:end) Subj filesep FileName],'allData')
    source = ['D' AnalysisFolder(2:end) Subj filesep FileName '.mat'];
    destFolder = ['L' AnalysisFolder(2:end) Subj filesep];
    mkdir(destFolder);
    destination = [destFolder FileName '.mat'];
    copyfile(source,destination)
    disp(['Done ' Subj ' in ' num2str(toc) 'sec.'])
end
%Done until 106

%% for grand matrix Go to - 
% grand_N1P2.m
% 
