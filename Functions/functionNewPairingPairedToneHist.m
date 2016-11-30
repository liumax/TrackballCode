%This function will generate information about pairing (histograms and
%rasters)

%Inputs:
%s: large structured array for storage
%desigNames: names of single units (ntXclusterY)
%spikeName: single name of the target spikes. name should be spikesPairing,
%or something like this
%soundName: single name of sound information to extract, ex: Pairing
%params: structured array with appropriate parameters

% Outputs:
% s: large structured array. Stores data in s.designNames.fieldName (ex: Tuning1Analysis)

%NOTE: spikeName and ttlName must be in quotations.
function [s] = functionNewPairingPairedToneHist(s,desigNames,...
    spikeName,soundName,params); 

%generates name for new field in structured array
fieldName = strcat(soundName,'Analysis');

%pulls out things into variables to make them easier to call.
rasterWindow = params.pairingWindow;
histBin = params.histBin;
baselineBin = params.baselineBin;
smoothingBins = params.smoothingBins;
defaultBins = params.defaultBins;
calcWindow = params.calcWindow;
zLimit = params.zLimit;
firstSpikeWindow = params.firstSpikeWindow;

soundData = s.SoundData.(soundName);
alignTimes = s.TTLs.(soundName); %times for aligning rasters and histograms.
toneReps = soundData.ToneRepetitions;
%converts raster window from ratio to actual time in seconds.
rasterWindow = rasterWindow*soundData.ToneDuration;
baselineBin = baselineBin*soundData.ToneDuration;
firstSpikeWindow = firstSpikeWindow*soundData.ToneDuration;
%generates an axis for raster at ms resolution
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
%sets LFP window to the same as teh raster window.
lfpWindow = rasterWindow;
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2];

numUnits = size(s.DesignationName,2);

%for loop to go through each individual unit
for i = 1:numUnits
    %prep individual sets of spikes:
    spikeTimes = s.(desigNames{i}).(spikeName);
    [rasters] = functionBasicRaster(spikeTimes,alignTimes,rasterWindow);
    fullRasterData = rasters;
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    for k = 1:length(alignTimes)
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > baselineBin(1) & rasters(:,1) < baselineBin(2)),1);
    end
    
    averageRate = mean(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
    averageSTD = std(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
    
    [histCounts histCenters] = hist(rasters(:,1),histBinVector);
    fullHistData = histCounts'/length(alignTimes)/histBin;
    
    [firstSpikeTimes,firstSpikeStats,binSpikeTimes,binSpikeStats] = ...
    functionBasicFirstSpikeTiming(firstSpikeWindow,fullRasterData,toneReps,2,1:1:toneReps);
    
    %calculate significant response for combined histogram of all responses
    %to all sounds.
    inputRaster = rasters(rasters(:,1) > calcWindow(1) & rasters(:,1) < calcWindow(2),1);
    [respStore] = ...
    functionBasicResponseSignificance(smoothingBins,defaultBins,calcWindow,...
    averageRate,averageSTD,inputRaster,zLimit,length(alignTimes));
    fullResp = respStore;
    %if there is a significant response, stores important values: start,
    %duration, and peak size.
    if ~isempty(respStore{1})
        fullRespGraph(1) = respStore{1}(1,1);
        fullRespGraph(2) = respStore{1}(1,2)-respStore{1}(1,1);
        fullRespGraph(3) = respStore{1}(1,5);
    else
        fullRespGraph = zeros(3,1);
    end

    s.(desigNames{i}).(fieldName).AllRasters = fullRasterData;
    s.(desigNames{i}).(fieldName).AllHistograms = fullHistData;
    s.(desigNames{i}).(fieldName).AverageRate = averageRate;
    s.(desigNames{i}).(fieldName).AverageSTD = averageSTD;
    s.(desigNames{i}).(fieldName).FullResponseStats = fullResp;
    s.(desigNames{i}).(fieldName).FullResponseGraphs = fullRespGraph;
    s.(desigNames{i}).(fieldName).FirstSpikeTimes = firstSpikeTimes;
    s.(desigNames{i}).(fieldName).FirstSpikeStats = firstSpikeStats;
    s.(desigNames{i}).(fieldName).BinSpikeTimes = binSpikeTimes;
    s.(desigNames{i}).(fieldName).BinSpikeStats = binSpikeStats;
end



end