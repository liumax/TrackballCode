

%This code is meant to batch process hardware files from noldus. 

masterFolder = pwd;
masterDir = dir;
masterDir = {masterDir.name};
%this code here is to remove folders reported in dir like "." or ".."
masterIndex = strfind(masterDir,'Hardware');
secondIndex = strfind(masterDir,'.txt');
masterIndex = find(not(cellfun('isempty', masterIndex)));
secondIndex = find(not(cellfun('isempty', secondIndex)));
trueIndex = intersect(masterIndex,secondIndex);
%masterFolders is a listing of all folders that I want to target for batch
%analysis
masterFolders = masterDir(trueIndex);
numFiles = length(masterFolders);
for i = 1:numFiles
    %find true name
    dashFind = strfind(masterFolders{i},'-');
    spaceFind = strfind(masterFolders{i},' ');
    fname = masterFolders{i};
    trialNum = str2num(fname(spaceFind(end-1)+1:dashFind(end)-1));
    arenaNum = str2num(fname(spaceFind(end)+1));
    %pull data
    [timeStore,stateStore] = functionNoldusHardwareExtractCommand(fname);
    %save data
    newName = strcat('HardwareTrial',num2str(trialNum),'Arena',num2str(arenaNum),'.mat');
    save(newName,'timeStore','stateStore')
end