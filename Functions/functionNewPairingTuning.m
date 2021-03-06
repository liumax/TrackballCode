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
toneWindow = params.toneWindow;
genWindow = params.genWindow;
latBin = params.latBin;
percentCutoff = params.percentCutoff;
baselineCutoff = params.baseCutoff;
histBin = params.histBin;
baselineBin = params.baselineBin;
calcWindow = params.calcWindow;
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
toneWindow = toneWindow * soundData.ToneDur;
genWindow = genWindow * soundData.ToneDur;
calcWindow = calcWindow*soundData.ToneDur;
baselineBin = baselineBin*soundData.ToneDur;
firstSpikeWindow = firstSpikeWindow*soundData.ToneDur;
%generates an axis for raster at ms resolution
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
%sets LFP window to the same as teh raster window.
lfpWindow = rasterWindow;
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2];

numUnits = size(s.DesignationName,2);

s.BaselineSpikes.(spikeName) = zeros(numUnits,1);
s.SignificantSpikes.(spikeName) = zeros(numUnits,1);

%for loop to go through each individual unit
for i = 1:numUnits
    %prep individual sets of spikes:
    spikeTimes = s.(desigNames{i}).(spikeName);
    %rasterize spike data!
    [rasters] = functionBasicRaster(spikeTimes,alignTimes,rasterWindow);
    rasters(:,3) = master(rasters(:,2),5); %adds information about frequency/amplitude
    fullRasterData = rasters;
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    %Also perform histogram on each trial.
    indivHist = zeros(length(histBinVector),length(master));
    for k = 1:length(master)
        indivHist(:,k) = hist(rasters(rasters(:,2) == k,1),histBinVector)/histBin;
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > baselineBin(1) & rasters(:,1) < baselineBin(2)),1);
    end
    
    averageRate = mean(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
    averageSTD = std(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
    
    %make compensatory binning relative to average rate.
    compBinTone = averageRate*(toneWindow(2)-toneWindow(1));
    compBinGen = averageRate*(genWindow(2)-genWindow(1));
    
    [histCounts histCenters] = hist(rasters(:,1),histBinVector);
    fullHistData = histCounts'/length(master)/histBin;
    
    %calculate significant response for combined histogram of all responses
    %to all sounds.
    
    [generalResponseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,length(master));
    
    disp(strcat('Baseline Spikes:',num2str(generalResponseHist.SpikeNumber),' Unit:',(desigNames{i})))
    s.BaselineSpikes.(spikeName)(i) = generalResponseHist.SpikeNumber;
    if generalResponseHist.Warning == 0 & generalResponseHist.SigSpike == 1
        s.SignificantSpikes.(spikeName)(i) = 1;
    end
    [generalLatPeakBinOut] = functionLatPeakBinCalculation(toneWindow,genWindow,rasterWindow,rasters,...
                length(master),2,[1:1:length(master)],latBin,histBin,percentCutoff,baselineCutoff);
    
    %allocates empty array.
    responseHistHolder = cell(numFreqs,numDBs);
    organizedHist = zeros(numFreqs,numDBs,length(histBinVector));
    organizedRasters = cell(numFreqs,numDBs);
    histErr = zeros(numFreqs,numDBs,length(histBinVector));
    latStore = zeros(numFreqs,numDBs);
    peakStoreTone = zeros(numFreqs,numDBs);
    peakStoreGen = zeros(numFreqs,numDBs);
    binStoreTone = zeros(numFreqs,numDBs);
    binStoreGen = zeros(numFreqs,numDBs);
    probStoreTone = zeros(numFreqs,numDBs);
    probStoreGen = zeros(numFreqs,numDBs);
    
    for k = 1:numFreqs
        for l = 1:numDBs
            targetTrials = master(master(:,2) == uniqueFreqs(k) & master(:,3) == uniqueDBs(l),4); %finds the trial number for all trials of given frequency and amplitude
            findMatches = find(ismember(fullRasterData(:,2),targetTrials)); %uses these trial numbers to pull raster index
            targetRasters = fullRasterData(findMatches,:); %extracts target rasters from all rasters using the above index.
            organizedRasters{k,l} = targetRasters; %saves to organized rasters
            [histCounts,histCenters] = hist(targetRasters(:,1),histBinVector); %calculates histogram with defined bin size
            organizedHist(k,l,:) = histCounts/toneReps/histBin; %saves histogram
            %Calculate standard error based on individual histograms
            specHist = indivHist(:,targetTrials);
            histErr(k,l,:) = std(specHist,0,2)/sqrt(length(targetTrials));
            [latPeakBinOut] = functionLatPeakBinCalculation(toneWindow,genWindow,rasterWindow,targetRasters,...
                length(targetTrials),2,targetTrials,latBin,histBin,percentCutoff,baselineCutoff);
            latPeakBinOut.BinnedSpikesToneComp = latPeakBinOut.BinnedSpikesTone - compBinTone;
            latPeakBinOut.BinnedSpikesGenComp = latPeakBinOut.BinnedSpikesGen - compBinGen;
            latePeakBinStore{k,l} = latPeakBinOut;
            latStore(k,l) = latPeakBinOut.ResponseLatency;
            peakStoreTone(k,l) = latPeakBinOut.PeakRespTone;
            peakStoreGen(k,l) = latPeakBinOut.PeakRespGen;
            binStoreTone(k,l) = mean(latPeakBinOut.BinnedSpikesTone);
            binStoreGen(k,l) = mean(latPeakBinOut.BinnedSpikesGen);
            probStoreTone(k,l) = latPeakBinOut.ProbSpikeTone;
            probStoreGen(k,l) = latPeakBinOut.ProbSpikeGen;
            latPeakBinOut = [];
            inputRaster = targetRasters(:,1);
            [responseHist] = functionBasicResponseSignificance(s,calcWindow,...
            spikeTimes,alignTimes,length(master));
            responseHistHolder{k,l} = responseHist;
            responseHist = [];
        end
        freqSpecHist(k,:) = mean(squeeze(organizedHist(k,:,:))');
    end
    s.(desigNames{i}).(fieldName).AllRasters = fullRasterData;
    s.(desigNames{i}).(fieldName).AllHistograms = fullHistData;
    s.(desigNames{i}).(fieldName).IndivHistograms = indivHist;
    s.(desigNames{i}).(fieldName).FreqDBRasters = organizedRasters;
    s.(desigNames{i}).(fieldName).FreqDBHistograms = organizedHist;
    s.(desigNames{i}).(fieldName).FreqDBHistogramErrors = histErr;

    s.(desigNames{i}).(fieldName).OverallLatBinPeakCalcs = generalLatPeakBinOut;
    s.(desigNames{i}).(fieldName).IndividualLatbinPeakCalcs = latePeakBinStore;
    s.(desigNames{i}).(fieldName).LatencyStore = latStore;
    s.(desigNames{i}).(fieldName).PeakStoreTone = peakStoreTone;
    s.(desigNames{i}).(fieldName).PeakStoreGen = peakStoreGen;
    s.(desigNames{i}).(fieldName).BinStoreTone = binStoreTone;
    s.(desigNames{i}).(fieldName).BinStoreGen = binStoreGen;
    s.(desigNames{i}).(fieldName).ProbStoreTone = probStoreTone;
    s.(desigNames{i}).(fieldName).ProbStoreGen = probStoreGen;
    
    s.(desigNames{i}).(fieldName).FrequencyHistograms = freqSpecHist;
    s.(desigNames{i}).(fieldName).AverageRate = averageRate;
    s.(desigNames{i}).(fieldName).AverageSTD = averageSTD;
    s.(desigNames{i}).(fieldName).AllHistogramSig = generalResponseHist;
    s.(desigNames{i}).(fieldName).SpecHistogramSig = responseHistHolder;
    s.(desigNames{i}).(fieldName).HistBinVector = histBinVector;
    
    s.(desigNames{i}).(fieldName).CompBinTone = compBinTone;
    s.(desigNames{i}).(fieldName).CompBinGen = compBinGen;
    %fuck it, lets just calculate all the shit here before the plotting
    %code
    
    s.(desigNames{i}).(fieldName).BinStoreToneComp = binStoreTone - compBinTone;
    s.(desigNames{i}).(fieldName).BinStoreGenComp = binStoreGen - compBinGen;
    s.(desigNames{i}).(fieldName).PeakStoreGenComp = peakStoreGen - averageRate;
    
end



end