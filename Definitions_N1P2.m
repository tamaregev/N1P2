 %Definitions_N1P2
ExpName = 'N1P2';
Subjects =  {'s101','s102','s103','s104','s105','s106','s107','s108','s109','s110','s111','s112','s113','s114','s115','s116','s117','s118','s119','s120','s121','s122','s123','s124','s125','s126','s127','s128','s129','s130','s131','s132','s133','s134','s135','s136','s137'};
sessions = {'1','3','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1'};
nEDATS = ones(size(Subjects));
versions = {'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',''};

% 1 - Tamar
% 20 - Prozac
% 34 - didn't record first 3 blocks
% 29 - too many artifacts
badSubjects = [1, 20, 29, 34];
longSubjects = [11, 12, 17, 35];%1/2 infinite artifacts - bug in read_markers_artifacts

mode = [];
%drive = 'D:';%local
drive = 'L:';%network

ResultsFolder = [drive '\Experiments\' ExpName '\Experiment\Results\'];
AnalyzerExportFolder = ['L' ResultsFolder(2:end) 'EEG\export\'];
RawDataFolder = [ResultsFolder 'EEG\raw\'];
AnalysisFolder = [drive '\Experiments\' ExpName '\Analysis\'];
EDATfolder = [ResultsFolder 'EDAT\'];

ExportFolder = [AnalyzerExportFolder 'beforeSegmentations\'];

Expinfo.ExpName = ExpName;
Expinfo.Subjects = Subjects;Expinfo.sessions = sessions;
Expinfo.versions = versions;Expinfo.AnalysisFolder = AnalysisFolder;
Expinfo.EDATfolder = EDATfolder;Expinfo.ExportFolder = ExportFolder;
srate = 512;
%%
