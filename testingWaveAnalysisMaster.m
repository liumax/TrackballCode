


%This is meant to control wave analysis

waveStruct = struct;

%finds info about hte current directory
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

for folderCount = 1:numFolders
    %generates text that signifies path to the target folder
    targetPath = strcat(masterFolder,'\',masterFolders{folderCount});
    %moves to target path
    cd(targetPath)
    %pulls all matlab files
    fileNames = what;
    %extracts names of .mat files specifically
    fileNames = fileNames.mat;
    %pulls first file, assuming this will be the soundfile. Uses this to
    %create substructure in structured array. 
    nameHolder = fileNames{1};
    nameHolder = nameHolder(1:strfind(nameHolder,'.')-1);
    nameHolder = genvarname(nameHolder);
    waveStruct.(nameHolder) = [];
    %searches for files containing "analysis" in the name, pulls index.
    fileIndex = strfind(fileNames,'Analysis');
    fileIndex = find(not(cellfun('isempty', fileIndex)));
    %extracts specific file
    files = (fileNames{fileIndex});
    matclustHolder = open((files));
    %extracts name of structured array inside. This is to deal with issues
    %that may relate to matclustStruct vs masterStruct.
    temp = fieldnames(matclustHolder);
    matclustHolder = matclustHolder.matclustStruct;

    structFields = fieldnames(matclustHolder);
    fieldIndex = strfind(structFields,'matclust');
    fieldIndex = find(not(cellfun('isempty', fieldIndex)));
    for trodesCount = 1:size(fieldIndex,1)
       averageWaveData = matclustHolder.(structFields{trodesCount}).AverageWaveForms(:,:,2);
       averageRate = matclustHolder.(structFields{trodesCount}).AverageFiringRate;
       waveStruct.(nameHolder).(genvarname(num2str(trodesCount))).AverageFiringRate = averageRate;
       for clusterCount = 1:matclustHolder.(structFields{trodesCount}).Clusters
           
       end
    end
    
end