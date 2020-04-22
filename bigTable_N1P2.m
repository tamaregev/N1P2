%%create biTable for general lme models
%% All experiments together - prep big table
% One BIG model!!
include = false;%include subjects for which time of grandcon peak was determined due to median
%Experiment 1
N = 1;
cd('S:\Lab-Shared\Experiments\MMNchroma\Analysis')
Definitions_MMNchroma;%new idea for the first time! move definitions into a separate script!
savedate = '13-May-2019';%peaks calculated at the Amplitudes_MMNchroma script at:
load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
if include
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
    mss = [];%known
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end
ib=1;
Table = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
%add block variable
blockN = (categorical(ones(size(Table,1),1)));
spread = ((max(Table.currMIDI) -  min(Table.currMIDI))*ones(size(Table,1),1))/4;
bigTable.(['E' num2str(N)]) = [Table,table(blockN),table(spread)];
if length(bls)>1 
    for ib=2:length(bls)
        Table = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
        blockN = (categorical(ib*ones(size(Table,1),1)));
        spread = ((max(Table.currMIDI) -  min(Table.currMIDI))*ones(size(Table,1),1))/4;
        Table = [Table,table(blockN),table(spread)];
        bigTable.(['E' num2str(N)]) = [bigTable.(['E' num2str(N)]);Table];
    end
end
ExpN = categorical(N*ones(size(bigTable.(['E' num2str(N)]),1),1));
bigTable.(['E' num2str(N)]) = [bigTable.(['E' num2str(N)]), table(ExpN)];

%Experiment 2
N = 2;
cd('S:\Lab-Shared\Experiments\MMNchromaF\Analysis')
Definitions_MMNchromaF;
savedate = '13-May-2019';%peaks calculated at the Amplitudes_MMNchroma script at:
load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
if include
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
    mss = 31;%known
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end
ib=1;
Table = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
%add block variable
blockN = (categorical(ones(size(Table,1),1)));
spread = ((max(Table.currMIDI) -  min(Table.currMIDI))*ones(size(Table,1),1))/4;
bigTable.(['E' num2str(N)]) = [Table,table(blockN),table(spread)];
if length(bls)>1 
    for ib=2:length(bls)
        Table = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
        blockN = (categorical(ib*ones(size(Table,1),1)));
        spread = ((max(Table.currMIDI) -  min(Table.currMIDI))*ones(size(Table,1),1))/4;
        Table = [Table,table(blockN),table(spread)];
        bigTable.(['E' num2str(N)]) = [bigTable.(['E' num2str(N)]);Table];
    end
end
ExpN = categorical(N*ones(size(bigTable.(['E' num2str(N)]),1),1));
bigTable.(['E' num2str(N)]) = [bigTable.(['E' num2str(N)]), table(ExpN)];

% Experiment 3
N = 3;
cd('S:\Lab-Shared\Experiments\N1P2\Analysis\N1P2_GH')
Definitions_N1P2
cd(GHfolder)
blocks = {'b1','b2a','b2b','b3a','b3b'};
savedate = '31-Oct-2018';
load([mixedFolder 'LMEtables_combined_' savedate ], 'allTables')
if include
    includeSubjects = allSubjects(~ismember(allSubjects,badSubjects));
else
    mss = [5,14];%known
    includeSubjects = allSubjects(~ismember(allSubjects,[badSubjects,mss]));
end
ib=1;
Table = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);
%add block variable
blockN = (categorical(ones(size(Table,1),1)));
spread = ((max(Table.currMIDI) -  min(Table.currMIDI))*ones(size(Table,1),1))/4;
bigTable.(['E' num2str(N)]) = [Table,table(blockN),table(spread)];
if length(bls)>1 
    for ib=2:length(bls)
        Table = allTables.(blocks{ib})(ismember(allTables.(blocks{ib}).subject,includeSubjects),:);;
        blockN = (categorical(ib*ones(size(Table,1),1)));
        spread = ((max(Table.currMIDI) -  min(Table.currMIDI))*ones(size(Table,1),1))/4;
        Table = [Table,table(blockN),table(spread)];
        bigTable.(['E' num2str(N)]) = [bigTable.(['E' num2str(N)]);Table];
    end
end
ExpN = categorical(N*ones(size(bigTable.(['E' num2str(N)]),1),1));
bigTable.(['E' num2str(N)]) = [bigTable.(['E' num2str(N)]), table(ExpN)];

%all Experiments together 
ie = 1;
theTable = bigTable.(['E' num2str(ie)]);
for ie=2:3
   theTable = [theTable; bigTable.(['E' num2str(ie)])];
end
theTableNoP = theTable(:,[1:7,11:13]);
theTableP = theTable(:,8:9);
Voltage = zscore(-1*theTableP.N1);N1V = table(Voltage);
Voltage = zscore(theTableP.P2);P2V = table(Voltage);
Ns = cell(size(theTable,1),1);for i=1:length(Ns),Ns{i,1} = 'N1';end
Ps = cell(size(theTable,1),1);for i=1:length(Ps),Ps{i,1} = 'P2';end
Potential = categorical([Ns;Ps]);
theTableN = [theTableNoP,N1V];
theTableP = [theTableNoP,P2V];
theTable = [theTableN;theTableP];
theTable = [theTable,table(Potential)];
%fix subjects with different numbers. Oct 22
theTable.subject = theTable.subject + grp2idx(theTable.ExpN)*100;
theTable.ExpN = grp2idx(theTable.ExpN);
isN1 = zeros(height(theTable),1);isP2 = zeros(height(theTable),1);
isN1(theTable.Potential=='N1')=1;isP2(theTable.Potential=='P2')=1;
theTable=[theTable,table(isN1),table(isP2)];
%save into mixed folder of Exp3 N1P2
save([mixedFolder filesep 'theTable'],'theTable')