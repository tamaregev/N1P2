%% functional Plots for model prediction
%% definitions
Definitions_N1P2

grandFolder = [AnalysisFolder 'grandAverage'];
mixedFolder = [AnalysisFolder 'MixedModel\'];
modelFolder = [AnalysisFolder 'Model\'];

srate = 512;
Expinfo.srate = srate;

%codes and names table:
Bnames = {'1','2a','2b','3a','3b'}';
Bcodes = [10,20,30,40,50]';
Bcodes_names = table(Bnames,Bcodes);
blocks = Bnames;
bls=1:5;
%% plot all cond, prevcond

sigma = 10;
tau = 0.4;

if 1 % as in functionalPlots_GH
    if 1 %params
        electrodeName = 'Cz';%'Fz'/'Cz'/'central cluster'/'GFP'
        addtag = '';
        include = false;%include subjects for which there was a manual inspection
    end
    %load tables
    peakdate = '30-Oct-2018';
    load([mixedFolder 'N1P2_' electrodeName '_' addtag peakdate])%loaded for allGrandcon_as_median
    nCond = size(allGrandcon_times{1},2);
    savedate = '31-Oct-2018';
    load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
    blocks = {'b1','b2a','b2b','b3a','b3b'};
    ibs = 1:5;
    relCols = 2:7;
    %exclude subjects with no clear N1 or P2:
    if include
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
    else
        mss = [];
        for bl=1:length(allGrandcon_as_median)
            for ipeak = 1:2%here I exclude from all analysis any subject that either had a problem in the N1 or in the P2s
                [r, c]=find(allGrandcon_as_median{bl}(:,:,ipeak));
                mss = [mss, r'];
            end
        end
        mss = unique(mss);
        allSubjects = 1:length(Subjects);
        includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
    end
    for ib=ibs
        allTables.(blocks{ib}) = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
    end
end

% prepare allTables for the model
%load Model
RAdate = '19-Dec-2018';
loadFolder = [modelFolder RAdate filesep];
load([loadFolder 'RApeaks'])
load([loadFolder 'Params'])
sigmas = floor(sigmas*1000)/1000;taus=floor(taus*1000)/1000;
load([loadFolder 'scalingFactors'])  
RATables = struct;
nlines = height(allTables.b1);
N1 = nan(nlines,1);P2 = nan(nlines,1);
for ib=ibs
    Table = [allTables.(blocks{ib})(:,1:7), table(N1), table(P2)];
    for iline = 1:nlines
        subj = Table.subject(iline);
        curr = str2double(Table.currNote{iline}(6));
        prev = str2double(Table.prevNote{iline}(6));
        RA = RApeaks(subj,ib,curr,prev,sigmas==sigma,taus==tau);
        Predict = 1-RA;
        Table.N1(iline) = Predict.*SFs(sigma==sigmas,tau==taus,1);
        Table.P2(iline) = Predict.*SFs(sigma==sigmas,tau==taus,2);
    end
    RATables.(blocks{ib}) = Table;
end

%plot here:
hf = plotFunctional_N1P2( RATables, relCols, blocks, ibs );

suptitle(['Exp: ' ExpName '. Model prediction sigma = ' num2str(sigma) ', tau = ' num2str(tau)])
