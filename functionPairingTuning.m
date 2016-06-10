%This function is meant to calculate information of tuning curves based on
%selected subset of spikes. 

%NOTE: spikeName and ttlName must be in quotations.
function [masterStruct] = functionPairingTuning(masterStruct,truncatedNames,...
    spikeName,ttlName,soundName,rasterWindow,histBin,clusterWindow,...
    clims1,rpvTime,trodesDesignation,fileName); 
%pulls out things into variables to make them easier to call.

soundData = masterStruct.SoundData.(soundName);
master = masterStruct.SoundData.(soundName).MasterArray;
TTLTimes = masterStruct.TTLs.(ttlName);
master(:,1) = TTLTimes;
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
    structureTarget = struct;
    %extracts overall rasters and histograms
    [structureTarget] = functionPairingRasterHistExtraction(i,clusters,...
    master,numDBs,numFreqs,structureTarget,truncatedNames,spikeTimes,...
    rasterWindow,histBin,histBinVector);
    %Extracts frequency specific traces.
    [structureTarget] = functionPairingFreqAmpRasterHist(i, clusters,...
    structureTarget,histBinNum,histBinVector,histBin,truncatedNames,...
    master,soundData,masterStruct);
    %stores into master structure
    masterStruct.(truncatedNames{i}).(soundName) = structureTarget;
    structureTarget = [];
    %plots in standard format
    [masterStruct] = functionPairingTuningPlot(i,...
    truncatedNames,rpvTime,clusterWindow,rasterWindow,histBinVector,...
    clims1,numFreqs,numDBs,uniqueFreqs,uniqueDBs,...
    trodesDesignation,masterStruct,fileName,clusters,soundName);
    

%     
end


end