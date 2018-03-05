
%this is meant to be master code to run analysis code!
%Establishes the home folder as the base
masterFolder = pwd;
masterDir = dir;
masterDir = {masterDir.name};
%this code here is to remove folders reported in dir like "." or ".."
masterIndex = strfind(masterDir,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
%masterFolders is a listing of all folders that I want to target for batch
%analysis
masterFolders = masterDir(masterIndex);
numFolders = size(masterFolders,2);

for masterCount = 1:numFolders
%generates text that signifies path to the target folder
targetPath = strcat(masterFolder,'\',masterFolders{masterCount});
cd(targetPath)
analysisWrapperPhotPavlov(targetPath);
close all

cd(masterFolder)
end






