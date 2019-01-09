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

[findString] = functionCellStringFind(masterFolders,'.pdf');
masterFolders(findString) = [];

numFolders = length(masterFolders);
for masterCount = 1:numFolders
    %generates text that signifies path to the target folder
    targetPath = strcat(masterFolder,'\',masterFolders{masterCount});
    cd(targetPath)

    fileName = masterFolders{masterCount};
%     fileName = fileName(1:end-4);
    
%     periodFinder = strfind(fileName,'.');
%     fileName = fileName(1:periodFinder-1);
    try
        analysisTuningAltLaser(fileName);
    catch
        disp('DIO FUCKUP')
    end
    cd(masterFolder)
end