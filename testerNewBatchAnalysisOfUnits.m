%This code is meant to go through, run analysis scripts, and then show
%figures to allow for approval/disapproval. 

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

%establishes master array for storage of data

masterStruct = struct;
masterStruct.masterNames = cell(0,0);
masterStruct.TuningDesignation = [];
masterStruct.LaserDesignation = [];

for masterCount = 1:numFolders
    %generates text that signifies path to the target folder
    targetPath = strcat(masterFolder,'\',masterFolders{masterCount});
    cd(targetPath)
    %what lists matlab files in the folder. can extract based on file type.
    fileName = what;
    fileName = fileName.mat{1};
    periodFinder = strfind(fileName,'.');
    fileName = fileName(1:periodFinder-1);
    [s] = functionBasicTuning(fileName);
    %opens tool which will allow for selection of tuned units, stores these
    %units in decisionTuning (1 for tuned, 0 for untuned).
    [decisionTuning] = functionTuningSelectionTool(s,fileName);
    
    %adjust designationNames to include filename
    desigNames = genvarname(strcat(fileName,s.DesignationName));
    %store these data! 
    masterStruct.masterNames = [masterStruct.masterNames,desigNames];
    masterStruct.TuningDesignation = [masterStruct.TuningDesignation;decisionTuning];
    %store unit analysis data.
    for i = 1:length(s.DesignationName);
        masterStruct.(desigNames{i}) = s.(s.DesignationName{i});
    end
    %store sound data.
    masterStruct.(strcat('SoundData',fileName)) = s.SoundData;
    
    %calculate laser data. 
    [s] = functionBasicLaserResponse(fileName);
    %opens tool which will allow for selection of laser units, stores these
    %units in decisionTuning (1 for laserID, 0 for nonID).
    [decisionLaser] = functionLaserSelectionTool(s,fileName);
    masterStruct.LaserDesignation = [masterStruct.TuningDesignation;decisionLaser];
    


    cd(masterFolder)
end