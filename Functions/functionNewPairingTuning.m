%This function will generate information about tuning curves for the
%purposes of calculations of pairing. 

%Inputs:
%s: large structured array for storage
%desigNames: names of single units (ntXclusterY)
%spikeName: single name of the target spikes. name should be spikesTuning1,
%or something like this
%soundName: single name of sound information to extract, ex: Tuning1
%params: structured array with appropriate parameters

% Outputs:
% s: large structured array. Stores data in s.designNames.fieldName (ex: Tuning1Analysis)

%NOTE: spikeName and ttlName must be in quotations.
function [s] = functionNewPairingTuning(s,desigNames,...
    spikeName,soundName,params); 

%generates name for new field in structured array
fieldName = strcat(soundName,'Analysis');

%pulls out things into variables to make them easier to call.
rasterWindow = params.rasterWindow;
histBin = params.histBin;
baselineBin = params.baselineBin;
smoothingBins = params.smoothingBins;
defaultBins = params.defaultBins;
calcWindow = params.calcWindow;
zLimit = params.zLimit;
firstSpikeWindow = params.firstSpikeWindow;

soundData = s.SoundData.(soundName);
master = s.SoundData.(soundName).MasterArray;
alignTimes = s.TTLs.(soundName); %times for aligning rasters and histograms.
master(:,1) = alignTimes;
s.SoundData.(soundName).MasterArray(:,1) = alignTimes;
uniqueFreqs = unique(master(:,2));
uniqueDBs = unique(master(:,3));
numDBs = size(unique(master(:,3)),1);
numFreqs = size(unique(master(:,2)),1);
toneReps = soundData.ToneReps;
%converts raster window from ratio to actual time in seconds.
rasterWindow = rasterWindow*soundData.ToneDur;
baselineBin = baselineBin*soundData.ToneDur;
firstSpikeWindow = firstSpikeWindow*soundData.ToneDur;
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
    rasters(:,3) = master(rasters(:,2),5); %adds information about frequency/amplitude
    fullRasterData = rasters;
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    for k = 1:length(master)
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > baselineBin(1) & rasters(:,1) < baselineBin(2)),1);
    end
    
    averageRate = mean(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
    averageSTD = std(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
    
    [histCounts histCenters] = hist(rasters(:,1),histBinVector);
    fullHistData = histCounts'/length(master)/histBin;
    
    %calculate significant response for combined histogram of all responses
    %to all sounds.
    inputRaster = rasters(rasters(:,1) > calcWindow(1) & rasters(:,1) < calcWindow(2),1);
    [respStore] = ...
    functionBasicResponseSignificance(smoothingBins,defaultBins,calcWindow,...
    averageRate,averageSTD,inputRaster,zLimit,length(master));
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
    
    %allocates empty array.
    sigResp = cell(numFreqs,numDBs);
    
    for k = 1:numFreqs
        for l = 1:numDBs
            targetTrials = master(master(:,2) == uniqueFreqs(k) & master(:,3) == uniqueDBs(l),4); %finds the trial number for all trials of given frequency and amplitude
            findMatches = find(ismember(fullRasterData(:,2),targetTrials)); %uses these trial numbers to pull raster index
            targetRasters = fullRasterData(findMatches,:); %extracts target rasters from all rasters using the above index.
            organizedRasters{k,l} = targetRasters; %saves to organized rasters
            [histCounts,histCenters] = hist(targetRasters(:,1),histBinVector); %calculates histogram with defined bin size
            organizedHist(k,l,:) = histCounts/toneReps/histBin; %saves histogram
            [firstSpikeTimes,firstSpikeStats,binSpikeTimes,binSpikeStats] = ...
                functionBasicFirstSpikeTiming(firstSpikeWindow,targetRasters,toneReps,2,targetTrials); %calculates information about first spike timing
            firstSpikeTimeHolder{k,l} = firstSpikeTimes; %saves first spike times
            firstSpikeStatsHolder(k,l,:,:) = firstSpikeStats; %saves statistics about first spikes
            binSpikeHolder{k,l} = binSpikeTimes; %binned spikes from the defined window.
            binSpikeStatsHolder(k,l,:,:) = binSpikeStats; %stats about those spikes
            inputRaster = targetRasters(targetRasters(:,1) > calcWindow(1) & targetRasters(:,1) < calcWindow(2),1);
            [respStore] = ...
            functionBasicResponseSignificance(smoothingBins,defaultBins,calcWindow,...
            averageRate,averageSTD,inputRaster,zLimit,toneReps); %calculates significance of responses to specific frequencies and dBs
            sigResp{k,l} = respStore;
            if ~isempty(respStore{1}) %if there is significant response, stores this for later display.
                sigRespGraph(k,l,1) = respStore{1}(1,1);
                sigRespGraph(k,l,2) = respStore{1}(1,2)-respStore{1}(1,1);
                sigRespGraph(k,l,3) = respStore{1}(1,5);
            else
                sigRespGraph(k,1,:) = 0;
            end
        end
        freqSpecHist(k,:) = mean(squeeze(organizedHist(k,:,:)));
    end
    s.(desigNames{i}).(fieldName).AllRasters = fullRasterData;
    s.(desigNames{i}).(fieldName).AllHistograms = fullHistData;
    s.(desigNames{i}).(fieldName).FreqDBRasters = organizedRasters;
    s.(desigNames{i}).(fieldName).FreqDBHistograms = organizedHist;
    s.(desigNames{i}).(fieldName).FirstSpikeTimes = firstSpikeTimeHolder;
    s.(desigNames{i}).(fieldName).FirstSpikeStats = firstSpikeStatsHolder;
    s.(desigNames{i}).(fieldName).BinSpikes = binSpikeHolder;
    s.(desigNames{i}).(fieldName).BinSpikeStats = binSpikeStatsHolder;
    s.(desigNames{i}).(fieldName).FrequencyHistograms = freqSpecHist;
    s.(desigNames{i}).(fieldName).AverageRate = averageRate;
    s.(desigNames{i}).(fieldName).AverageSTD = averageSTD;
    s.(desigNames{i}).(fieldName).ResponseStats = sigResp;
    s.(desigNames{i}).(fieldName).ResponseStatsGraph = sigRespGraph;
    s.(desigNames{i}).(fieldName).FullResponseStats = fullResp;
    s.(desigNames{i}).(fieldName).FullResponseGraphs = fullRespGraph;
end



end