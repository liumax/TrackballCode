


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

counterIndex = 1;
fullWaveInfo = zeros(4,1000);
allInfoHolder = cell(1000,1);

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
       averageWaves = matclustHolder.(structFields{trodesCount}).AverageWaveForms(:,:,2);
       averageRate = matclustHolder.(structFields{trodesCount}).AverageFiringRate;
       rpvs = matclustHolder.(structFields{trodesCount}).RPVs;
       waveStruct.(nameHolder).(genvarname(num2str(trodesCount))).AverageFiringRate = averageRate;
       [output,allPeaks] = functionWavePropertyExtraction(averageWaves,averageRate);
       outputSizer = size(output,2);
       fullWaveInfo(1:3,counterIndex:counterIndex+outputSizer-1) = output;
       fullWaveInfo(4,counterIndex:counterIndex+outputSizer-1) = rpvs;
       counterIndex = counterIndex + outputSizer;
    end
    
    cd(masterFolder)
end

%clears out zero values!
fullWaveInfo = fullWaveInfo(:,fullWaveInfo(1,:)>0);