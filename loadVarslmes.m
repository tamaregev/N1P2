function [ dataFolder, bT, bRA, bwhichSubjects, taus, sigmas, whichpeaks ] = loadVarslmes( Folder )
%LOADVARSLMES is instead of running Section 2 in the script
%lmeModels_allExp before every Section..
dataFolder = [Folder filesep 'allExp'];
cd(dataFolder)

load dataTable bT
tic;load bRA;toc%about 15 sec on Sassi
load bwhichSubjects.mat
load taus
load sigmas
whichpeaks = {'N1','P2'};

end

