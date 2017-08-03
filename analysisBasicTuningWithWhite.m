%This code is meant for the basic analysis of tuning curve data, examining
%and plotting basic properties of the tuning curve and displaying them as a
%figure. 

%%Inputs
%fileName: name of the target file, without file extension. 

%%Outputs
%s: structured array that has all the data from the analysis. s will
%contain the following: 
%
%n number of units, which are independent of clusters/trodes. This is to
%say that each unit will have its own named field. For example, 2 clusters
%on channel 10 will have two separate fields.

%DesignationArray/DesignationName: The array indicates which probe sites
%and clusters were used, designation names should be the names of all
%units.

%SoundData: all the sound data from the tuning curve, with addition of
%unique frequencies/dbs and number of frequencies/dbs.

function [s] = analysisBasicTuningWithWhite(fileName);
%% Constants and things you might want to tweak
%TOGGLES FOR ENABLING/DISABLING FEATURES
s.Parameters.toggleRPV = 0; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 0; %1 means you want to eliminate duplicates.
toggleROC = 0; %toggle for tuning on/off ROC analysis

%PARAMETERS FOR BASIC ARRANGEMENT OF DATA
s.Parameters.RasterWindow = [-4 3]; %ratio for raster window. will be multiplied by toneDur
s.Parameters.ToneWindow = [0 1];
s.Parameters.GenWindow = [0 2];
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate
s.Parameters.LFPWindow = [-1 3];

%stuff for significance
s.Parameters.calcWindow = [0 2]; %defines period for looking for responses, based on toneDur
s.Parameters.zLimit = [0.05 0.01 0.001];
s.Parameters.minSpikes = 100; %minimum number of spikes in baseline, if lower, triggers a warning
s.Parameters.minSigSpikes = 5; %minimum number of significant points to record a significant response.
s.Parameters.PercentCutoff = 99.9; %for significance in latency calculations
s.Parameters.BaselineCutoff = 95; %for the onset in latency calculations
s.Parameters.latBin = 0.001; %histogram bins for latency and significance calculations
s.Parameters.SigSmoothWindow = 11; %window of smoothing for calculations of significance

%for duplicate elimination
s.Parameters.DownSampFactor = 10; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr

%for rotary encoder:
s.Parameters.InterpolationStepRotary = 0.01; %interpolation steps in seconds.

%for edr
s.Parameters.EDRdownsamp = 20; %number of samples to downsample by. Smoothing is likely unnecessary
s.Parameters.EDRTimeCol = 1;
s.Parameters.EDRTTLCol = 3;
s.Parameters.EDRPiezoCol = 2;

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%set other things
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
% format short
disp('Parameters Set')
%% sets up file saving stuff
saveName = strcat(fileName,'FullTuningAnalysis','.mat');
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
%pull lfp file names. this allows detection of the number of trodes. 
try
    [lfpFiles] = functionFileFinder(subFoldersCell,'LFP','LFP');
    s.NumberTrodes = length(lfpFiles);
catch
    [paramFiles] = functionFileFinder(subFoldersCell,'matclust','param');
    s.NumberTrodes = length(paramFiles) - length(matclustFiles);
end
disp('Matclust Files Extracted')
%generate placeholder structure
% s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Parameters.RPVTime,...
    matclustFiles,s,s.Parameters.ClusterWindow,subFoldersCell);

if toggleDuplicateElimination ==1
    if length(s.DesignationName) > 1
        disp('Now Selecting Based on xCORR')
        [s] = functionDuplicateElimination(s,s.Parameters.DownSampFactor,...
            s.Parameters.corrSlide,s.Parameters.ThresholdComparison,s.Parameters.trodesFS,s.Parameters.RPVTime,s.Parameters.ClusterWindow,...
            s.ShankDesignation,s.ShankMap,s.Shanks);
    end
else
    disp('NOT EXECUTING DUPLICATE ELIMINATION')
end
%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%% Extracts Sound Data from soundFile, including freq, db, reps.
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;

%check for presence of white noise
try
    whiteStatus = soundData.WhiteNoise;
catch
    whiteStatus = 0;
end

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
s.Parameters.ToneWindow = s.Parameters.ToneWindow * toneDur;
s.Parameters.GenWindow = s.Parameters.GenWindow * toneDur;
calcWindow = s.Parameters.calcWindow*toneDur;
rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes
% s.Parameters.CalcWindow = s.Parameters.CalcWindow * toneDur;
s.Parameters.LFPWindow = s.Parameters.LFPWindow * toneDur;

% Generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for i =1:numFreqs
    freqNameHolder{i} = num2str(uniqueFreqs(i));
end

if whiteStatus ==0
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
elseif whiteStatus == 1;
    % Finds the number of octaves, makes array of octave steps. This will be used for imagesc graphing applications
    totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(2)));
    %This then makes an array of the full octave steps I've made
    octaveRange = zeros(totalOctaves + 1,2);
    octaveRange(1,1) = uniqueFreqs(2);
    for i = 1:totalOctaves
        octaveRange (i+1,1) = octaveRange(i,1)*2;
    end
    %next, I find the positions from uniqueFreqs that match octaveRange
    for i = 1:size(octaveRange,1);
        octaveRange(i,2) = find(round(uniqueFreqs) == octaveRange(i,1));
    end
    octaveRange(2:end+1,:) = octaveRange(1:end,:);
    octaveRange(1,1) = 0;
    octaveRange(1,2) = 1;
end

% Does the same for dBs. 
if length(uniqueDBs) == 1
    dbSteps = 1;
    totalDBs = 1;
    dbRange = [100,1];
else
    dbSteps = uniqueDBs(2) - uniqueDBs(1);
    totalDBs = (uniqueDBs(end) - uniqueDBs(1))/dbSteps;
    dbRange = zeros(totalDBs + 1,2);
    dbRange(:,1) = uniqueDBs(1):dbSteps:uniqueDBs(end);
    for i = 1:size(dbRange,1)
        dbRange(i,2) = find(uniqueDBs == dbRange(i,1));
    end
end


%stores into s.
s.Parameters.OctaveRange = octaveRange;
s.Parameters.DBRange = dbRange;

disp('Sound Data Extracted')
%% Extract DIO information. Tuning curve should rely on just one DIO output, DIO1.
disp('Beginning DIO Extraction...')
%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
if length(D1FileName) == 0
    [D1FileName] = functionFileFinder(subFoldersCell,'DIO','Din1');
end
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);
s.SoundData.ToneTimes = dioTimes;
s.SoundData.ToneTimeDiff = dioTimeDiff;

disp('DIO Data Extracted, Checking for Errors')
%insert to master. check for errors
if length(dioTimes) == length(soundData.Frequencies)
    %Sets up master array so that we have a size comparison for DIO
    %information. 
    master = zeros(size(soundFile.soundData.Frequencies,1),5);

    master(:,1) = dioTimes;
    %master(:,2) is frequency
    master(:,2) = s.SoundData.Frequencies;
    %master(:,3) is dB
    master(:,3) = s.SoundData.dBs;
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
    
    repArray = ones(numFreqs,numDBs);
    repArray = repArray * toneReps;
elseif length(dioTimes) ~= length(s.SoundData.Frequencies) %error case
    disp('ERROR IN DIOTIMES vs PROJECTED TONES')
    %the below code is optimized for single errors, and assumes that the
    %first TTL is correct. This code may need to be rewritten in the
    %future to account for failures of the first TTL. 
    predTimes = s.SoundData.Delays;
    minPred = min(predTimes);
    
    repArray  = zeros(numFreqs,numDBs); %make array of zeros!
    alertArray = zeros(numFreqs,numDBs);
    %assume that the first TTL is fine
    try
        freqFind = find(uniqueFreqs == s.SoundData.TrialMatrix(1,2));
        dbFind = find(uniqueDBs == s.SoundData.TrialMatrix(1,3));
    catch
        freqFind = find(uniqueFreqs == s.SoundData.Frequencies(1));
        dbFind = find(uniqueDBs == s.SoundData.dBs(1));
    end
    repArray(freqFind,dbFind) = repArray(freqFind,dbFind) + 1;
    whileTrig = 1;
    whileCounter = 1;
    alertCounter = 1;
    while whileTrig == 1
        %first see if I've exceeded the bounds
        if whileCounter > length(predTimes) %if hit bound of predicted times]
            disp('Went through all predicted times')
            break
            
        elseif whileCounter > length(dioTimeDiff) %if hit bound of received times
            disp('Went through all received times')
            totalTrialNum = whileCounter;
            break
            
        else
            %check the difference between the current values. If exceeds
            %threshold, flag. Otherwise, need to record the tuning
            %information so I get accurate representation of repetitions
            diffCheck = dioTimeDiff(whileCounter) - predTimes(whileCounter);
            if diffCheck < minPred & diffCheck > 0
                %now lets fill in rep numbers
                try
                    freqFind = find(uniqueFreqs == s.SoundData.TrialMatrix(whileCounter + 1,2));
                    dbFind = find(uniqueDBs == s.SoundData.TrialMatrix(whileCounter + 1,3));
                catch
                    freqFind = find(uniqueFreqs == s.SoundData.Frequencies(1));
                    dbFind = find(uniqueDBs == s.SoundData.dBs(1));
                end
                repArray(freqFind,dbFind) = repArray(freqFind,dbFind) + 1;
                whileCounter = whileCounter + 1;
            elseif diffCheck >= minPred
                disp('Error Found! Insufficient DIO')
                try
                    freqFind = find(uniqueFreqs == s.SoundData.TrialMatrix(whileCounter + 1,2));
                    dbFind = find(uniqueDBs == s.SoundData.TrialMatrix(whileCounter + 1,3));
                catch
                    freqFind = find(uniqueFreqs == s.SoundData.Frequencies(1));
                    dbFind = find(uniqueDBs == s.SoundData.dBs(1));
                end
                %store information about the failed DIO
                alertArray(freqFind,dbFind) = alertArray(freqFind,dbFind) + 1;
                alertDesig(alertCounter) = whileCounter;
                %first fix predicted times
                predTimes(whileCounter+1) = predTimes(whileCounter+1) + predTimes(whileCounter);
                predTimes(whileCounter) = [];
                try
                    s.SoundData.TrialMatrix(whileCounter,:) = [];
                    s.SoundData.Frequencies(whileCounter) = [];
                    s.SoundData.dBs(whileCounter) = [];
                catch
                    s.SoundData.Frequencies(whileCounter) = [];
                    s.SoundData.dBs(whileCounter) = [];
                end
                %dont need to update whileCounter, since this will cycle
                %through again.
            elseif diffCheck < 0 %this seems to work more or less!
                disp('Error Found! Excess DIO')
                try
                    freqFind = find(uniqueFreqs == s.SoundData.TrialMatrix(whileCounter + 1,2));
                    dbFind = find(uniqueDBs == s.SoundData.TrialMatrix(whileCounter + 1,3));
                catch
                    freqFind = find(uniqueFreqs == s.SoundData.Frequencies(1));
                    dbFind = find(uniqueDBs == s.SoundData.dBs(1));
                end
                %delete the excess DIO
                dioTimes(whileCounter + 1) = [];
                dioTimeDiff = diff(dioTimes);
            else
                error('Code Failure In DIO Repair')
            end
        end
    end
    master = zeros(totalTrialNum,5);

    master(:,1) = dioTimes;
    %master(:,2) is frequency
    master(:,2) = s.SoundData.Frequencies(1:totalTrialNum);
    %master(:,3) is dB
    master(:,3) = s.SoundData.dBs(1:totalTrialNum);
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
    
end
disp('DIO data successfully extracted and stored.')

%% Extract data from rotary encoder.
[funcOut] = functionRotaryExtraction(s.Parameters.trodesFS,s.Parameters.InterpolationStepRotary,subFoldersCell);
s.RotaryData = funcOut;
disp('Rotary Encoder Data Extracted')
%rasterize this data
jumpsBack = round(s.Parameters.RasterWindow(1)/s.Parameters.InterpolationStepRotary);
jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);
velRaster = zeros(jumpsForward-jumpsBack+1,totalTrialNum);
length(s.RotaryData.Velocity);
for i = 1:totalTrialNum
    %find the time from the velocity trace closest to the actual stim time
    targetInd = find(s.RotaryData.Velocity(:,1)-dioTimes(i) > 0,1,'first');
    %pull appropriate velocity data
    velRaster(:,i) = s.RotaryData.Velocity([targetInd+jumpsBack:targetInd+jumpsForward],2);
end
%make average trace:
averageVel = mean(velRaster,2);
velVector = [s.Parameters.RasterWindow(1):s.Parameters.InterpolationStepRotary:s.Parameters.RasterWindow(2)];
velZero = find(velVector >= 0,1,'first');

%Change toggle for ROC analysis if insufficient running is found
if s.RotaryData.RawDistance(end,2) < 10;
    toggleROC = 0;
    disp('Resetting ROC Toggle Due to Lack of Movement')
end


figure
subplot(2,1,1)
plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2))
xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
title('Velocity Over Session (cm/s)')
subplot(2,1,2)
hold on
plot(velVector,averageVel,'r','LineWidth',2)
plot([0 0],[ylim],'b');
xlim([velVector(1) velVector(end)])
title('Average Velocity Traces')

%see if EDR data exists
fileNames = dir(homeFolder);
fileNames = {fileNames.name};
targetFileFinder = strfind(fileNames,'.EDR'); %examines names for D1
targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
if ~isempty(targetFileFinder)
    %flip toggle for graphing
    edrToggle = 1;
    %pull filename
    targetFile = fileNames{targetFileFinder};%pulls out actual file name
    %extract data
    [s] = functionEDRPull(s,targetFile,s.Parameters.EDRTimeCol,s.Parameters.EDRTTLCol,s.Parameters.EDRPiezoCol);
    %downsample for sake of space and time
    s.EDR.NewData = downsample(s.EDR.FullData,s.Parameters.EDRdownsamp);
    %confirm ttls line up
    if length(s.EDR.TTLTimes) ~= totalTrialNum
        error('Incorrect match of EDR TTLs and trial number')
    end
    %figure out raster window
    edrTimeStep = mean(diff(s.EDR.FullData(:,s.Parameters.EDRTimeCol)));
    jumpsBack = round(s.Parameters.RasterWindow(1)/(edrTimeStep*s.Parameters.EDRdownsamp));
    jumpsForward = round(s.Parameters.RasterWindow(2)/(edrTimeStep*s.Parameters.EDRdownsamp));
    edrRaster = zeros(jumpsForward-jumpsBack+1,totalTrialNum);
    for i = 1:totalTrialNum
        %find time closest to actual stim time
        targetInd = find(s.EDR.NewData(:,s.Parameters.EDRTimeCol) - s.EDR.TTLTimes(i) > 0,1,'first');
        edrRaster(:,i) = s.EDR.NewData([targetInd + jumpsBack:targetInd + jumpsForward],s.Parameters.EDRPiezoCol);
    end
    
    %calculate overall mean
    edrMean = mean(edrRaster');
    edrAbsMean = mean(abs(edrRaster'));
    edrVector = [s.Parameters.RasterWindow(1):edrTimeStep*s.Parameters.EDRdownsamp:s.Parameters.RasterWindow(2)];
    edrZero = find(edrVector >= 0,1,'first');
else
    disp('NO EDR FILE FOUND')
    edrToggle = 0;
end

if edrToggle == 1
    figure
    hold on
    plot(edrVector,edrMean,'b')
    plot(edrVector,edrAbsMean,'r')
    plot([0 0],[ylim],'b')
    xlim([edrVector(1) edrVector(end)])
    title('Mean and AbsMean Piezo')
end


%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db
for i = 1:numUnits
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    alignTimes = master(:,1);
    
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    
    %calculates rasters based on spike information. 
    [rasters] = functionBasicRaster(spikeTimes,alignTimes,s.Parameters.RasterWindow);
    rasters(:,3) = master(rasters(:,2),5); %adds information about frequency/amplitude
    fullRasterData = rasters;
    
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    averageSpikeHolder = zeros(totalTrialNum,1);
    for k = 1:totalTrialNum
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > s.Parameters.RasterWindow(1) & rasters(:,1) < 0),1);
    end
    averageRate = mean(averageSpikeHolder/(-s.Parameters.RasterWindow(1)));
    averageSTD = std(averageSpikeHolder/(-s.Parameters.RasterWindow(1)));
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
%     [generalResponseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,length(master));
    [generalResponseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,totalTrialNum,...
        s.Parameters.minSpikes,s.Parameters.latBin,[s.Parameters.RasterWindow(1),0],s.Parameters.zLimit,s.Parameters.minSigSpikes,s.Parameters.SigSmoothWindow);
    
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
    latPeakBinStore = cell(numFreqs,numDBs);
    latStore = zeros(numFreqs,numDBs);
    peakStoreGen = zeros(numFreqs,numDBs);
    binStoreTone = zeros(numFreqs,numDBs);
    binStoreGen = zeros(numFreqs,numDBs);
    probStoreTone = zeros(numFreqs,numDBs);
    probStoreGen = zeros(numFreqs,numDBs);
    freqSpecHist = zeros(numFreqs,histBinNum,1);
    
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
            [latPeakBinOut] = functionLatPeakBinCalculation(s.Parameters.ToneWindow,s.Parameters.GenWindow,s.Parameters.RasterWindow,...
                targetRasters,length(targetTrials),2,targetTrials,s.Parameters.latBin,s.Parameters.histBin,s.Parameters.PercentCutoff,s.Parameters.BaselineCutoff);
            latePeakBinStore{k,l} = latPeakBinOut;
            latStore(k,l) = latPeakBinOut.ResponseLatency;
            peakStoreGen(k,l) = latPeakBinOut.PeakRespGen;
            binStoreTone(k,l) = mean(latPeakBinOut.BinnedSpikesTone);
            binStoreGen(k,l) = mean(latPeakBinOut.BinnedSpikesGen);
            probStoreTone(k,l) = latPeakBinOut.ProbSpikeTone;
            probStoreGen(k,l) = latPeakBinOut.ProbSpikeGen;
            latPeakBinOut = [];
%             [responseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,toneReps);
            [responseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes(targetTrials),length(targetTrials),...
        s.Parameters.minSpikes,s.Parameters.latBin,[s.Parameters.RasterWindow(1),0],s.Parameters.zLimit,s.Parameters.minSigSpikes,s.Parameters.SigSmoothWindow);
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
    s.(desigNames{i}).FrequencyHistograms = freqSpecHist;
    s.(desigNames{i}).AverageRate = averageRate;
    s.(desigNames{i}).SessionFiring = sessionFiring;
    s.(desigNames{i}).AverageSTD = averageSTD;
    s.(desigNames{i}).AverageSTE = averageSTE;
    s.(desigNames{i}).HistBinVector = histBinVector;
    s.(desigNames{i}).AllHistogramSig = generalResponseHist;
    s.(desigNames{i}).SpecHistogramSig = responseHistHolder;
    s.(desigNames{i}).LatPeakBin = latePeakBinStore;
    s.(desigNames{i}).LatencyMap = latStore;
    s.(desigNames{i}).PeakMap = peakStoreGen;
    s.(desigNames{i}).BinTone = binStoreTone;
    s.(desigNames{i}).BinGen = binStoreGen;
    s.(desigNames{i}).ProbTone = probStoreTone;
    s.(desigNames{i}).ProbGen = probStoreGen;
    
    if toggleROC == 1
        targetName = desigNames{i};
        [s] = functionLocomotionROC(s,targetName);
    end
    
end

%calculate and plot LFP information
% [lfpStruct] = functionLFPaverage(master, s.Parameters.LFPWindow, s,homeFolder,fileName, uniqueFreqs, uniqueDBs, numFreqs, numDBs);
% s.LFP = lfpStruct;

%% Plotting

if toggleTuneSelect == 1 %if you want tuning selection...
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    decisionTuning = zeros(numUnits,1);
    for i = 1:numUnits
        %plots average waveform
        subplot(4,6,1)
        hold on
        plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
        title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
        %plots ISI
        subplot(4,6,2)
        hist(s.(desigNames{i}).ISIGraph,1000)
        histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
        line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
        xlim(s.Parameters.ClusterWindow)
        title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
            strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
        %Plot binned response during tone period
        subplot(4,3,4)
        imagesc(s.(desigNames{i}).BinTone')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Mean Binned Response (tone)')
        %Plot binned response during general period
        subplot(4,3,7)
        imagesc(s.(desigNames{i}).BinGen')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Mean Binned Response (general)')
        %Plot peak response during general period
        subplot(4,3,10)
        imagesc(s.(desigNames{i}).PeakMap')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Peak Response (general)')
        %plot velocity data
        subplot(4,3,6)
        hold on
        plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'b')
        plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
        xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
        ylim([-0.1,1])
        if toggleROC == 1
            title(strcat('Vel & Firing Rate. AUC:',num2str(s.(desigNames{i}).TrueAUC),'99%Range',num2str(prctile(s.(desigNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(desigNames{i}).ShuffleAUC,1))))
        else
            title('Vel & Firing Rate')
        end
        %plot probability of response (tone)
        subplot(4,3,9)
        imagesc(s.(desigNames{i}).ProbTone')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Probability of Response (tone)')
        %plot probability of response (gen)
        subplot(4,3,12)
        imagesc(s.(desigNames{i}).ProbGen')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Probability of Response (general)')

        %plots rasters (chronological)
        subplot(3,3,2)
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,2),'k.','markersize',4)
        hold on
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        title({fileName;desigNames{i}},'fontweight','bold')
        set(0, 'DefaulttextInterpreter', 'none')
        %plots rasters (frequency and amplitude organized)
        subplot(3,3,5)
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,3),'k.','markersize',4)
        hold on
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        rasterFreqLines = zeros(numFreqs,2);
        rasterFreqLines(:,1) = toneReps*size(uniqueDBs,1)/2:toneReps*size(uniqueDBs,1):totalTrialNum;
        rasterFreqLines(:,2) = uniqueFreqs;
        %this generates green lines separating by Frequency
        for k = 1:size(uniqueFreqs,1)
            plot(s.Parameters.RasterWindow,[toneReps*numDBs*k toneReps*numDBs*k],'g','LineWidth',1)
        end
        set(gca,'YTick',rasterFreqLines(:,1));
        set(gca,'YTickLabel',rasterFreqLines(:,2));
        set(gca,'Ydir','reverse')
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Descending = increase in amplitude and freq')
        %plot heatmap organized by frequency
        subplot(3,3,8)
        imagesc(s.(desigNames{i}).FrequencyHistograms(:,:))
        colormap(parula)
        colorbar
        set(gca,'YTick',octaveRange(:,2));
        set(gca,'YTickLabel',octaveRange(:,1));
        set(gca,'XTick',[1:10:size(histBinVector,2)]);
        set(gca,'XTickLabel',histBinVector(1:20:end));
        histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
        histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),toneDur);
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
    %         title('Heatmap by Frequency and Time Max')
        title('Frequency Arranged Heatmap')
        % plot histogram.
        subplot(4,3,3)
        plot(histBinVector,s.(desigNames{i}).AllHistograms,'k','LineWidth',2)
        hold on
        plot(histBinVector,s.(desigNames{i}).AllHistograms - s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
        plot(histBinVector,s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
        %plot significant values
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1,1),...
            'b*')
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1,1),...
            'bo')
        %plot negative values for first tuning
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
            'k*')
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Histogram')

        hold off
        spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')

        %ask for input! 
        promptCounter = 1; %This is used to run the while loop.
        whileCounter = 0; %this is the counter that gets updated to exit the loop

        while whileCounter ~= promptCounter
            try
                prompt = 'How is this unit tuned? (excite(e)/inhib(i)/both(b)/none(n))';
                str = input(prompt,'s');
                if str~='n' & str~='e' & str~='i' & str~='b'
                    error
                else
                    whileCounter = 1;
                end
            catch
            end
        end
        if strfind(str,'e') | strfind(str,'i') | strfind(str,'b')
            decisionTuning(i) = 1;
        elseif strfind(str,'n')
            decisionTuning(i) = 0;
        end
        
        tuningType{i} = str;
        %clear figure.
        clf

    end
    s.TuningDecision = decisionTuning;
    s.TuningType = tuningType;
    close
else %in the case you dont want to do tuning selection, default to normal system
    for i = 1:numUnits
        hFig = figure;
        set(hFig, 'Position', [10 80 1240 850])
        %plots average waveform
        subplot(4,6,1)
        hold on
        plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
        title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
        %plots ISI
        subplot(4,6,2)
        hist(s.(desigNames{i}).ISIGraph,1000)
        histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
        line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
        xlim(s.Parameters.ClusterWindow)
        title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
            strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
        %Plot binned response during tone period
        subplot(4,3,4)
        imagesc(s.(desigNames{i}).BinTone')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Mean Binned Response (tone)')
        %Plot binned response during general period
        subplot(4,3,7)
        imagesc(s.(desigNames{i}).BinGen')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Mean Binned Response (general)')
        %Plot peak response during general period
        subplot(4,3,10)
        imagesc(s.(desigNames{i}).PeakMap')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Peak Response (general)')
        %plot velocity data
        subplot(4,3,6)
        hold on
        plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'b')
        plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
        xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
        ylim([-0.1,1])
        if toggleROC == 1
            title(strcat('Vel & Firing Rate. AUC:',num2str(s.(desigNames{i}).TrueAUC),'99%Range',num2str(prctile(s.(desigNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(desigNames{i}).ShuffleAUC,1))))
        else
            title('Vel & Firing Rate')
        end
        %plot probability of response (tone)
        subplot(4,3,9)
        imagesc(s.(desigNames{i}).ProbTone')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Probability of Response (tone)')
        %plot probability of response (gen)
        subplot(4,3,12)
        imagesc(s.(desigNames{i}).ProbGen')
        colormap(parula)
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Probability of Response (general)')

        %plots rasters (chronological)
        subplot(3,3,2)
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,2),'k.','markersize',4)
        hold on
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        title({fileName;desigNames{i}},'fontweight','bold')
        set(0, 'DefaulttextInterpreter', 'none')
        %plots rasters (frequency and amplitude organized)
        subplot(3,3,5)
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,3),'k.','markersize',4)
        hold on
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        rasterFreqLines = zeros(numFreqs,2);
        rasterFreqLines(:,1) = cumsum(sum(repArray'));
        rasterFreqLines(:,2) = uniqueFreqs;
        %this generates green lines separating by Frequency
        for k = 1:size(uniqueFreqs,1)
            plot(s.Parameters.RasterWindow,[toneReps*numDBs*k toneReps*numDBs*k],'g','LineWidth',1)
        end
        set(gca,'YTick',rasterFreqLines(:,1));
        set(gca,'YTickLabel',rasterFreqLines(:,2));
        set(gca,'Ydir','reverse')
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Descending = increase in amplitude and freq')
        %plot heatmap organized by frequency
        subplot(3,3,8)
        imagesc(s.(desigNames{i}).FrequencyHistograms(:,:))
        colormap(parula)
        colorbar
        set(gca,'YTick',octaveRange(:,2));
        set(gca,'YTickLabel',octaveRange(:,1));
        set(gca,'XTick',[1:10:size(histBinVector,2)]);
        set(gca,'XTickLabel',histBinVector(1:20:end));
        histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
        histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),toneDur);
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
    %         title('Heatmap by Frequency and Time Max')
        title('Frequency Arranged Heatmap')
        % plot histogram.
        subplot(4,3,3)
        plot(histBinVector,s.(desigNames{i}).AllHistograms,'k','LineWidth',2)
        hold on
        plot(histBinVector,s.(desigNames{i}).AllHistograms - s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
        plot(histBinVector,s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
        %plot significant values
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1,1),...
            'b*')
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1,1),...
            'bo')
        %plot negative values for first tuning
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
            'k*')
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Histogram')

        hold off
        spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')

    end
end

%% Saving
save(fullfile(pname,fname),'s');

end