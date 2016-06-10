%This function is meant to calculate information of tuning curves based on
%selected subset of spikes. 

%NOTE: spikeName and ttlName must be in quotations.
function [masterStruct] = functionPairingToneAnalysis(masterStruct,truncatedNames,...
    spikeName,ttlName,soundName,rasterWindow,histBin,clusterWindow,...
    clims1,rpvTime,trodesDesignation,fileName); 
%pulls out things into variables to make them easier to call.

soundData = masterStruct.SoundData.(soundName);
%pulls out freq and db data for easy access.
targetSet = soundData.TargetSet;
controlSet = soundData.ControlSet;
%pulls frequencies.
targetFreq = targetSet(1);
controlFreq = controlSet(1);
%pulls corresponding trials of each.
targetTrials = find(soundData.Frequencies == targetFreq);
controlTrials = find(soundData.Frequencies == controlFreq);
%set up array of names for naming things!
toneNames = {'Target','Control'};
%pullls master array.
master = masterStruct.SoundData.(soundName).MasterArray;
%finds all TTL times and feeds them into the master array.
TTLTimes = masterStruct.TTLs.(ttlName);
master(:,1) = TTLTimes;
%divide master array into two arrays: one for control, one for target.
masterTarget = master(targetTrials,:);
masterControl = master(controlTrials,:);

%pulls unique frequency and DB information
uniqueFreqs = unique(master(:,2));
uniqueDBs = unique(master(:,3));
numDBs = size(unique(master(:,3)),1);
numFreqs = size(unique(master(:,2)),1);
%converts raster window from ratio to actual time in seconds.
rasterWindow = rasterWindow*soundData.ToneDur;
%generates an axis for raster at ms resolution
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
%sets LFP window to the same as teh raster window.
lfpWindow = rasterWindow;
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2];
numTrodes = size(truncatedNames,2);
%for loop going through individual matclust files.
for i = 1:numTrodes
    %preps individual spike sets and cluster.
    spikeTimes = masterStruct.(truncatedNames{i}).(spikeName);
    clusters = masterStruct.(truncatedNames{i}).Clusters;
    for j = 1:size(toneNames,2) %for control vs target
        structureTarget = struct; %generates target structure for holding info.
        nameHolder = strcat('master',toneNames{j}); %generates name automatically! avoids hardcoding
        [structureTarget] = functionPairingRasterHistExtraction(i,clusters,...
    (nameHolder),structureTarget,spikeTimes,rasterWindow,histBin,histBinVector);
        %store the data in the appropriate structured area.
        masterStruct.(truncatedNames{i}).(soundName).(toneNames{j}) = structureTarget;
        structureTarget = [];
    end  
end
end