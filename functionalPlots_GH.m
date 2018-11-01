% functional plots Nov 1 2018
% plotting N1P2 amplitudes as a function of the explaining variables
%% definitions
Definitions_N1P2

grandFolder = [AnalysisFolder 'grandAverage'];
mixedFolder = [AnalysisFolder 'MixedModel\'];

srate = 512;
Expinfo.srate = srate;

%codes and names table:
Bnames = {'1','2a','2b','3a','3b'}';
Bcodes = [10,20,30,40,50]';
Bcodes_names = table(Bnames,Bcodes);
blocks = Bnames;
%% plot

%load peak amplitudes




