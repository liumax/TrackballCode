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
masterStruct = struct;
for masterCount = 1:numFolders
    %generates text that signifies path to the target folder
    targetPath = strcat(masterFolder,'\',masterFolders{masterCount});
    cd(targetPath)
    %what lists matlab files in the folder. can extract based on file type.
    fileName = what;
    %pulls the first .mat file. this should be the actual sound data file. 
    fileName = fileName.mat{1};
    %finds first period, makes a reduced file name based on that.
    periodFinder = strfind(fileName,'.');
    fileName = fileName(1:periodFinder-1);
    saveName = genvarname(fileName);
    %performs analysis
    [s] = analysisBasicTuning(fileName);
    %determine number of units
    unitNum = length(s.DesignationName);
    %in a for loop, generate array space for each unit
    for u = 1:unitNum
        minorName = strcat(saveName,'_',s.DesignationName{u});
        masterStruct.(minorName) = s.(s.DesignationName{u});
        masterStruct.(minorName).SoundData = s.SoundData;
    end
    cd(masterFolder)
end

saveName = strcat('BatchTuningAnalysis',date,'.mat');
save (fullfile(pwd,saveName),'masterStruct');