%This code will perform analysis of tuning, laser response, and try and
%remove duplicate cells from overlap across tetrodes. It will also try and
%remove artifact "units". 




function [s] = analysisIDTuningCombined(fileName);

%% Constants and things you might want to tweak
%Tuning Parameters
s.Parameters.RasterWindow = [-1 3]; %ratio for raster window. will be multiplied by toneDur
s.Parameters.RPVTime = 0.001; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.SampleRate = 30000;%trodes sampling rate
s.Parameters.DefaultBins = 0.001;% bin size for calculating significant responses
s.Parameters.SmoothingBins = [0.01 0.001];%bins for smoothing
s.Parameters.CalcWindow = [0 2]; %window for calculating significant responses
s.Parameters.zLimit = 3; %zlimit for calculating significant responses
% s.Parameters.FirstSpikeWindow = [0 0.5 1 1.5]; %ratios! need to be multiplied by tone duration.
s.Parameters.FirstSpikeWindow = [0 1];
s.Parameters.ChosenSpikeBin = 1; %delineates which spike window I will graph.
s.Parameters.BaselineBin = [-1 0]; %ratio for bin from which baseline firing rate will be calculated
s.Parameters.LFPWindow = [-1 3];

%stuff for significance
s.Parameters.calcWindow = [0 2]; %defines period for looking for responses, based on toneDur
s.Parameters.zLimit = [0.05 0.01 0.001];
s.Parameters.numShuffle = 1000;
s.Parameters.firstSpikeWindow = [0 1];%defines period for looking for first spike, based on toneDur
s.Parameters.chosenSpikeBin = 1; %spike bin selected in binSpike (in the event of multiple spike bins)
s.Parameters.minSpikes = 100; %minimum number of spikes to do spike shuffling
s.Parameters.minSigSpikes = 2; %minimum number of significant points to record a significant response.
s.Parameters.BaselineWindow = [-0.4 0]; %window for counting baseline spikes, in SECONDS. NOTE THIS IS DIFFERENT FROM RASTER WINDOW
s.Parameters.BaselineCalcBins = 1; %bin size in seconds if there are insufficient baseline spikes to calculate a baseline rate.

%Parameters for laser ID
s.Parameters.LaserWindow = [-0.1 0.1]; %duration for laser rasters/histograms, in seconds
s.Parameters.laserBin = 0.001; %histogram bin size for lasers
s.Parameters.laserCutoffLate = 15; %time in ms after which responses are not counted
s.Parameters.laserCutoffEarly = 2; %time in ms before which responses are not counted


%% sets up file saving stuff
saveName = strcat(fileName,'FullTuningAndIDAnalysis','.mat');
fname = saveName;
pname = pwd;

%% Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%% Find and ID Matclust Files for Subsequent Analysis. Generates Structured Array for Data Storage
%pull matclust file names
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');
%generate placeholder structure
% s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Parameters.RPVTime,...
    matclustFiles,s,s.Parameters.ClusterWindow);

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%% Extracts Sound Data from soundFile, including freq, db, reps.
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;

%pull important sound and trial information
uniqueFreqs = unique(soundData.Frequencies);
uniqueDBs = unique(soundData.dBs);
numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);
%store these in structured array
s.SoundData.UniqueFrequencies = uniqueFreqs;
s.SoundData.UniqueDBs = uniqueDBs;
s.SoundData.NumFreqs = numFreqs;
s.SoundData.NumDBs = numDBs;

toneDur = soundData.ToneDuration;
toneReps = soundData.ToneRepetitions;
totalTrialNum = length(soundData.Frequencies);

%Recalculate raster window and such
s.Parameters.RasterWindow = s.Parameters.RasterWindow * toneDur;
s.Parameters.FirstSpikeWindow = s.Parameters.FirstSpikeWindow * toneDur;
s.Parameters.BaselineBin = s.Parameters.BaselineBin * toneDur;
calcWindow = s.Parameters.calcWindow*toneDur;
rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes
s.Parameters.CalcWindow = s.Parameters.CalcWindow * toneDur;
s.Parameters.LFPWindow = s.Parameters.LFPWindow * toneDur;

% Generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for i =1:numFreqs
    freqNameHolder{i} = num2str(uniqueFreqs(i));
end

% Finds the number of octaves, makes array of octave steps. This will be used for imagesc graphing applications
totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));
%This then makes an array of the full octave steps I've made
octaveRange = zeros(totalOctaves + 1,2);
octaveRange(1,1) = uniqueFreqs(1);
for i = 1:totalOctaves
    octaveRange (i+1,1) = octaveRange(i,1)*2;
end
%next, I find the positions from uniqueFreqs that match octaveRange
for i = 1:size(octaveRange,1);
    octaveRange(i,2) = find(round(uniqueFreqs) == octaveRange(i,1));
end

% Does the same for dBs. 
dbSteps = uniqueDBs(2) - uniqueDBs(1);
totalDBs = (uniqueDBs(end) - uniqueDBs(1))/dbSteps;
dbRange = zeros(totalDBs + 1,2);
dbRange(:,1) = uniqueDBs(1):dbSteps:uniqueDBs(end);
for i = 1:size(dbRange,1)
    dbRange(i,2) = find(uniqueDBs == dbRange(i,1));
end

%stores into s.
s.Parameters.OctaveRange = octaveRange;
s.Parameters.DBRange = dbRange;

%Sets up master array so that we have a size comparison for DIO
%information. 
master = zeros(size(soundFile.soundData.Frequencies,1),5);

%master(:,1) is reserved for actual input times
%master(:,2) is frequency
master(:,2) = soundData.Frequencies;
%master(:,3) is dB
master(:,3) = soundData.dBs;
%master(:,4) is trial number (chronological)
master(:,4) = 1:1:totalTrialNum;
%master(:,5) is trial num, arranging trials in order from small dB to large
%dB, and low freq to high freq. frequency is larger category.
sortingCounter = 1;
for i = 1:numFreqs
    for j = 1:numDBs
        sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == uniqueDBs(j));
        master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end

%% Extract DIO information. Tuning curve should rely on just one DIO output, DIO1. 
%Laser response should rely on DIO2

%find DIO folder and D1/D2 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.SampleRate);

%insert to master. check for errors
if length(dioTimes) ~= length(master)
    length(dioTimes)
    length(master)
    error('dioTimes and master mismatched')
elseif length(dioTimes) == length(master)
    master(:,1) = dioTimes;
end

%extracts DIO stuffs! this code is written to extract inputs for d2
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%pulls out DIO2 up state onsets.
[dio2Times,dio2TimeDiff] = functionBasicDIOCheck(DIO2Data,s.Parameters.SampleRate);



%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db
for i = 1:numUnits
    %allocate a bunch of empty arrays.
    fullHistData = zeros(histBinNum,1);
    
    organizedRasters = cell(numFreqs,numDBs,1); %allocates cell array for rasters
    organizedHist = zeros(numFreqs,numDBs,histBinNum,1);
    
    freqSpecHist = zeros(numFreqs,histBinNum,1);
    
    firstSpikeTimeHolder = cell(numFreqs,numDBs,1);
    firstSpikeStatsHolder = zeros(numFreqs,numDBs,4,size(s.Parameters.FirstSpikeWindow,2)-1);
    
    binSpikeHolder = cell(numFreqs,numDBs,1);
    binSpikeStatsHolder = zeros(numFreqs,numDBs,2,size(s.Parameters.FirstSpikeWindow,2)-1);
    
    averageSpikeHolder = zeros(totalTrialNum,1);
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    alignTimes = master(:,1);
    
    %calculates rasters based on spike information. 
    [rasters] = functionBasicRaster(spikeTimes,alignTimes,s.Parameters.RasterWindow);
    rasters(:,3) = master(rasters(:,2),5); %adds information about frequency/amplitude
    fullRasterData = rasters;
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    for k = 1:totalTrialNum
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > s.Parameters.BaselineBin(1) & rasters(:,1) < s.Parameters.BaselineBin(2)),1);
    end
    
    averageRate = mean(averageSpikeHolder/(s.Parameters.BaselineBin(2)-s.Parameters.BaselineBin(1)));
    averageSTD = std(averageSpikeHolder/(s.Parameters.BaselineBin(2)-s.Parameters.BaselineBin(1)));
    averageSTE = averageSTD/(sqrt(totalTrialNum-1));
    
    %calculates histograms of every single trial. 
    fullHistHolder = zeros(length(histBinVector),totalTrialNum);
    for k =1:totalTrialNum
        if size(rasters(rasters(:,2) == k,:),1) > 0
            [counts centers] = hist(rasters(rasters(:,2) == k,1),histBinVector);
            fullHistHolder(:,k) = counts;
        end
    end
    
    %calculate standard deviation for each bin
    histSTE = (std(fullHistHolder'))';
    
    [histCounts histCenters] = hist(rasters(:,1),histBinVector);
    fullHistData = histCounts'/totalTrialNum/s.Parameters.histBin;
    
    %calculate significant response for combined histogram of all responses
    %to all sounds.    
    [generalResponseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,length(master));

    disp(strcat('Baseline Spikes:',num2str(generalResponseHist.SpikeNumber),' Unit:',(desigNames{i})))
    s.BaselineSpikes = generalResponseHist.SpikeNumber;
    if generalResponseHist.Warning == 0 & generalResponseHist.SigSpike == 1
        s.SignificantSpikes(i) = 1;
    end
    
    %allocates empty array.
    organizedHist = zeros(numFreqs,numDBs,length(histBinVector));
    organizedRasters = cell(numFreqs,numDBs);
    responseHistHolder = cell(numFreqs,numDBs);
    histErr = zeros(numFreqs,numDBs,length(histBinVector));
    
    for k = 1:numFreqs
        for l = 1:numDBs
            targetTrials = master(master(:,2) == uniqueFreqs(k) & master(:,3) == uniqueDBs(l),4); %finds the trial number for all trials of given frequency and amplitude
            findMatches = find(ismember(fullRasterData(:,2),targetTrials)); %uses these trial numbers to pull raster index
            targetRasters = fullRasterData(findMatches,:); %extracts target rasters from all rasters using the above index.
            organizedRasters{k,l} = targetRasters; %saves to organized rasters
            [histCounts,histCenters] = hist(targetRasters(:,1),histBinVector); %calculates histogram with defined bin size
            organizedHist(k,l,:) = histCounts/toneReps/s.Parameters.histBin; %saves histogram
            specHist = fullHistHolder(:,targetTrials);
            histErr(k,l,:) = std(specHist,0,2)/sqrt(length(targetTrials));
            [firstSpikeTimes,firstSpikeStats,binSpikeTimes,binSpikeStats] = ...
                functionBasicFirstSpikeTiming(s.Parameters.FirstSpikeWindow,targetRasters,toneReps,2,targetTrials); %calculates information about first spike timing
            firstSpikeTimeHolder{k,l} = firstSpikeTimes; %saves first spike times
            firstSpikeStatsHolder(k,l,:,:) = firstSpikeStats; %saves statistics about first spikes
            binSpikeHolder{k,l} = binSpikeTimes; %binned spikes from the defined window.
            binSpikeStatsHolder(k,l,:,:) = binSpikeStats; %stats about those spikes
            [responseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,toneReps);
            responseHistHolder{k,l} = responseHist;
        end
        freqSpecHist(k,:) = mean(squeeze(organizedHist(k,:,:)));
    end
        disp(strcat('Raster and Histogram Extraction Complete: Unit ',num2str(i),' Out of ',num2str(numUnits)))
    s.(desigNames{i}).AllRasters = fullRasterData;
    s.(desigNames{i}).AllHistograms = fullHistData;
    s.(desigNames{i}).IndividualHistograms = fullHistHolder; 
    s.(desigNames{i}).HistogramStandardDeviation = histSTE;
    s.(desigNames{i}).FreqDBRasters = organizedRasters;
    s.(desigNames{i}).FreqDBHistograms = organizedHist;
    s.(desigNames{i}).FreqDBHistogramErrors = histErr;
    s.(desigNames{i}).FirstSpikeTimes = firstSpikeTimeHolder;
    s.(desigNames{i}).FirstSpikeStats = firstSpikeStatsHolder;
    s.(desigNames{i}).BinSpikes = binSpikeHolder;
    s.(desigNames{i}).BinSpikeStats = binSpikeStatsHolder;
    s.(desigNames{i}).FrequencyHistograms = freqSpecHist;
    s.(desigNames{i}).AverageRate = averageRate;
    s.(desigNames{i}).AverageSTD = averageSTD;
    s.(desigNames{i}).AverageSTE = averageSTE;
    s.(desigNames{i}).HistBinVector = histBinVector;
    s.(desigNames{i}).AllHistogramSig = generalResponseHist;
    s.(desigNames{i}).SpecHistogramSig = responseHistHolder;
end





