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

function [s] = analysisTuningAltLaser(fileName);
%% Constants and things you might want to tweak

%TOGGLES FOR ENABLING/DISABLING FEATURES
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.
toggleROC = 0; %toggle for tuning on/off ROC analysis

%PARAMETERS FOR BASIC ARRANGEMENT OF DATA
s.Parameters.RasterWindow = [-6 4]; %ratio for raster window. will be multiplied by toneDur
s.Parameters.BaseWindow = [-4 0];
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

%for cell type
s.Parameters.PVLim = [0.0004 0.0005];
s.Parameters.ChatLim = 1.1;

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%for plotting the laser related stuff
s.Parameters.LaserWindow = [-0.3 0.4]; %plotting window for rasters
s.Parameters.LaserBin = 0.01; %histogram bin size
s.Parameters.LaserAnalysis = [-0.2,0;0.1,0.3];
% s.Parameters.LaserLim = 0.015; %maximum lag value for calculation.

%for looking at target window
toneTarget = 2; %this selects for which part of response I care about. 1 is fast, 2 is tone, 3 is general
sigCutoff = 0.05;
%set other things
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
% format short
disp('Parameters Set')
%% sets up file saving stuff
saveName = strcat(fileName,'FullTuningAnalysis','.mat');
fname = saveName;
pname = pwd;
%check to see if analysis file already exists.
searchNames = what;
searchNames = searchNames.mat;
searchResult = strcmp(searchNames,saveName);
if find(searchResult == 1) 
    prevFile = 1;
    disp('Previous Analysis Found!!')
else
    prevFile = 0;
end

%% Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
if ismac
    subFoldersCell = strsplit(subFolders,':')';
elseif ispc
    subFoldersCell = strsplit(subFolders,';')';
end

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

%if a previous analysis file exists, I want to default to that selection of
%duplicates, to avoid repeatedly choosing which units to save. 
if prevFile == 0
disp('No Previous Analysis, Performing Duplicate Elimination')
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
elseif prevFile == 1
    disp('Using Existing Analysis for Duplicate Elimination')
    n = load(saveName);
    %load deleted units
    if isfield(n.s,'DeletedUnits')
        s.DeletedUnits = n.s.DeletedUnits;
        %now delete bad units
        delUnits = fields(s.DeletedUnits);
        for indCount = 1:length(delUnits)
            s = rmfield(s,delUnits{indCount});
        end
        %replace designation array and names. 
        s.DesignationArray = n.s.DesignationArray;
        s.DesignationName = n.s.DesignationName;
    else
        %replace designation array and names. 
        s.DesignationArray = n.s.DesignationArray;
        s.DesignationName = n.s.DesignationName;
    end
    
    disp('Finished! Closing original file')
    clear n
end
    
%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;


%generate master array for 2-d storage of important values
masterData = zeros(numUnits,10);
masterHeader = cell(10,1);
masterInd = 1;

%% Produce Indices for sorting by position on shank
for reInd = 1:length(s.DesignationName) %go through every unit
    %extract information about tetrode
    testTet = s.DesignationArray(reInd,1);
    %extract biggest waveform
    [bigVal bigInd] = max(max(s.(s.DesignationName{reInd}).AverageWaveForms));
    %calculate overall index
    overInd(reInd) = (testTet - 1)*4 + bigInd;
end

s.PeakWaveIndex = overInd;
[s.SortedPeakWaveIndex s.SortedPeakWaveOrder] = sort(overInd);
s.ShankLength = max(s.ShankMap(:,1)/s.Shanks);
%now lets try and map this onto real space. We will have to do something to
%allow things that are on the same electrode to map out separately. 

%Go through each electrode and find matching sites. if exist, then store.
%if extras, then figure out spacing. 
cellCount = 1;
posArray = zeros(numUnits,2);
for i = 1:length(s.ShankMap)
    %see if there are cells centered on given electrode
    cellFind = find(s.SortedPeakWaveIndex == s.ShankMap(i,1));
    if length(cellFind) == 0
        disp(strcat('No Cells For Electrode',num2str(s.ShankMap(i,1))))
    elseif length(cellFind) == 1
        disp(strcat('Single Cell For Electrode',num2str(s.ShankMap(i,1))))
        posArray(cellCount,1) = s.ShankMap(i,1);
        posArray(cellCount,2) = s.ShankMap(i,3);
        cellCount = cellCount + 1;
    elseif length(cellFind) > 1
        disp(strcat('Multiple Cells For Electrode',num2str(s.ShankMap(i,1))))
        incAmount = 1/length(cellFind); %amount to increment different cells by
        for incInd = 1:length(cellFind)
            posArray(cellCount,1) = s.ShankMap(i,1) + (incInd-1)*incAmount;
            posArray(cellCount,2) = s.ShankMap(i,3);
            cellCount = cellCount + 1;
        end
    end
end
%make things negative so they plot out correctly
posArray(:,1) = -1*posArray(:,1);
posArray(posArray(:,2) == 2,1) = posArray(posArray(:,2) == 2,1)+s.ShankLength;

s.ElectrodePositionDisplayArray = posArray;

%store this into master
masterData(:,masterInd) = posArray(s.SortedPeakWaveOrder,1); masterHeader{masterInd} = 'Distance from Top of Shank'; masterInd = masterInd + 1;
masterData(:,masterInd) = posArray(s.SortedPeakWaveOrder,2); masterHeader{masterInd} = 'Shank Designation'; masterInd = masterInd + 1;
masterData(:,masterInd) = s.SortedPeakWaveOrder; masterHeader{masterInd} = 'SimpleShankOrder'; masterInd = masterInd + 1;

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
uniqueDBs = unique(soundData.TrialMatrix(soundData.TrialMatrix(:,2)==uniqueFreqs(1),3));
uniqueDBSteps = unique(soundData.TrialMatrix(soundData.TrialMatrix(:,2)==uniqueFreqs(1),3));
numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);
trialMatrix = soundData.TrialMatrix;
trialMatrixLaser = soundData.TrialMatrixPaired;
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
s.Parameters.BaseWindow = s.Parameters.BaseWindow * toneDur;
calcWindow = s.Parameters.calcWindow*toneDur;
rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin
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
    octaveRange(:,2) = interp1(round(uniqueFreqs),[1:length(uniqueFreqs)],octaveRange(:,1));
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
    dbSteps = uniqueDBSteps(2) - uniqueDBSteps(1);
    totalDBs = (uniqueDBSteps(end) - uniqueDBSteps(1))/dbSteps;
    dbRange = zeros(totalDBs + 1,2);
    dbRange(:,1) = uniqueDBSteps(1):dbSteps:uniqueDBSteps(end);
    for i = 1:size(dbRange,1)
        dbRange(i,2) = find(uniqueDBSteps == dbRange(i,1));
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

%now we need to finish processing the DIO times to separate out laser vs
%non-laser. First, find all laser points
diffFind = find(dioTimeDiff < 0.01);
diffFindInd = diffFind + 1;
diffFindNew = sort([diffFind;diffFindInd]);
toneTimeLaser = diffFindInd + 1;
toneTime = diffFindInd + 2;
allToneTime = [toneTime;toneTimeLaser];

%clean up, since sometimes you'll end up with excess pulses.
finalPulse = length(dioTimes);
toneTimeLaser(toneTimeLaser > finalPulse) = [];
toneTime(toneTime>finalPulse) = [];

%replace with actual times.
toneTimeLaser = dioTimes(toneTimeLaser);
toneTime = dioTimes(toneTime);

figure
plot(toneTime - toneTimeLaser)
hold on
plot(soundData.DelaysPaired,'r')
title('Delay between non-laser and laser tones, expected in red')
% pause

%pull D2 info (laser ID)
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
if length(D2FileName) == 0
    [D2FileName] = functionFileFinder(subFoldersCell,'DIO','Din2');
end
D2FileName = D2FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%pulls out DIO up state onsets.
[dio2Times,dio2TimeDiff] = functionBasicDIOCheck(DIO2Data,s.Parameters.trodesFS);

laserOnsetTimes = dio2Times;

timeFinder = find(toneTimeLaser > s.TimeFilterRange(2));
toneTimeLaser(timeFinder) = [];

timeFinder = find(laserOnsetTimes > s.TimeFilterRange(2));
laserOnsetTimes(timeFinder) = [];

timeFinder = find(toneTime > s.TimeFilterRange(2));
toneTime(timeFinder) = [];

trialNum = length(toneTime);
trialNumLaser = length(toneTimeLaser);
% totalTrialNum = length(master);
disp('Time Filter Applied')

%now we want to determine how many trials were present for each tone. This
%should be a fairly simple nested for loop.

matrixTrialNum = zeros(numFreqs,numDBs);
for freqInd = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(freqInd),3));
    for dbInd = 1:numDBs
        matrixTrialNum(freqInd,dbInd) = length(find(trialMatrix(:,2) == uniqueFreqs(freqInd) & trialMatrix(:,3) == subUniqueDB(dbInd)));
    end
end

matrixTrialNumLaser = zeros(numFreqs,numDBs);
for freqInd = 1:numFreqs
    subUniqueDB = unique(trialMatrixLaser(trialMatrixLaser(:,2) == uniqueFreqs(freqInd),3));
    for dbInd = 1:numDBs
        matrixTrialNumLaser(freqInd,dbInd) = length(find(trialMatrixLaser(:,2) == uniqueFreqs(freqInd) & trialMatrixLaser(:,3) == subUniqueDB(dbInd)));
    end
end

% s.TrialMatrix = master;

s.TrialNumbers = matrixTrialNum;
s.TrialNumbersLaser = matrixTrialNumLaser;
% s.RepArray = repArray;

%% Generate Masters

%form master matrix
master = zeros(trialNum,5);

master(:,1) = toneTime;
%master(:,2) is frequency
master(:,2) = s.SoundData.Frequencies(1:trialNum);
%master(:,3) is dB
master(:,3) = s.SoundData.dBs(1:trialNum);
%master(:,4) is trial number (chronological)
master(:,4) = 1:1:trialNum;
%master(:,5) is trial num, arranging trials in order from small dB to large
%dB, and low freq to high freq. frequency is larger category.
sortingCounter = 1;
for i = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(i),3));
    for j = 1:numDBs
        sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == subUniqueDB(j));
        master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end


%form masterLaser matrix
masterLaser = zeros(trialNumLaser,5);

masterLaser(:,1) = toneTimeLaser;
%masterLaser(:,2) is frequency
masterLaser(:,2) = s.SoundData.FrequenciesPaired(1:trialNumLaser);
%masterLaser(:,3) is dB
masterLaser(:,3) = s.SoundData.dBsPaired(1:trialNumLaser);
%masterLaser(:,4) is trial number (chronological)
masterLaser(:,4) = 1:1:trialNumLaser;
%masterLaser(:,5) is trial num, arranging trials in order from small dB to large
%dB, and low freq to high freq. frequency is larger category.
sortingCounter = 1;
for i = 1:numFreqs
    subUniqueDB = unique(trialMatrixLaser(trialMatrixLaser(:,2) == uniqueFreqs(i),3));
    for j = 1:numDBs
        sortingFinder = find(masterLaser(:,2) == uniqueFreqs(i) & masterLaser(:,3) == subUniqueDB(j));
        masterLaser(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end

%% Extract data from rotary encoder.
[funcOut] = functionRotaryExtraction(s.Parameters.trodesFS,s.Parameters.InterpolationStepRotary,subFoldersCell);
s.RotaryData = funcOut;
disp('Rotary Encoder Data Extracted')
%rasterize this data
jumpsBack = round(s.Parameters.RasterWindow(1)/s.Parameters.InterpolationStepRotary);
jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);
velRaster = zeros(jumpsForward-jumpsBack+1,totalTrialNum);
length(s.RotaryData.Velocity);
for i = 1:length(allToneTime)
    %find the time from the velocity trace closest to the actual stim time
    targetInd = find(s.RotaryData.Velocity(:,1) - dioTimes(allToneTime(i)) > 0,1,'first');
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

%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db

[overOut,indivOut,s,masterData,masterHeader,masterInd] = functionTuningDataExtraction(numUnits,numDBs,numFreqs,uniqueFreqs,s,masterData,masterHeader,masterInd,histBinVector,trialNum,master,desigNames,calcWindow,histBinNum,whiteStatus);
s.NonLaserOverall = overOut;
for i = 1:numUnits
    fn = fieldnames(indivOut.(desigNames{i}));
    for j = 1:length(fn)
        s.(desigNames{i}).(fn{j}) = indivOut.(desigNames{i}).(fn{j});
    end
end

[overOut,indivOut,s,masterData,masterHeader,masterInd] = functionTuningDataExtraction(numUnits,numDBs,numFreqs,uniqueFreqs,s,masterData,masterHeader,masterInd,histBinVector,trialNumLaser,masterLaser,desigNames,calcWindow,histBinNum,whiteStatus);
s.LaserOverall = overOut;
for i = 1:numUnits
    fn = fieldnames(indivOut.(desigNames{i}));
    for j = 1:length(fn)
        s.(desigNames{i}).(strcat(fn{j},'Laser')) = indivOut.(desigNames{i}).(fn{j});
    end
end

%now lets try and plot out slope of binned responses. Use tone period only

%pull out significant values only, just for tones.  
% sigCutoff = 0.05;
minResp = 5;
analysisWindow = 2;
disp('Analyzing Binned Responses for Linear Regression')
for i = 1:numUnits
    disp(strcat('Analyzing Unit-',desigNames{i}))
    valStore = [];
    if s.SoundData.WhiteNoise == 0
        tester = s.(desigNames{i}).BinTone;
        valStore(:,1) = reshape(tester,1,[]);
        tester = s.(desigNames{i}).BinToneLaser;
        valStore(:,2) = reshape(tester,1,[]);

        [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigCutoff); %b1 is for y intercept, b2 is slope.
    elseif s.SoundData.WhiteNoise == 1
        tester = s.(desigNames{i}).BinTone(2:end,:);
        valStore(:,1) = reshape(tester,1,[]);
        tester = s.(desigNames{i}).BinToneLaser(2:end,:);
        valStore(:,2) = reshape(tester,1,[]);

        [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigCutoff); %b1 is for y intercept, b2 is slope.
    end
    %for bintr, first row is range for intercept, second row is for range of slope
%         b
%         bintr
    s.(desigNames{i}).RegRawVals = valStore;
    s.(desigNames{i}).RegVals = b;
    s.(desigNames{i}).RegValSig = bintr;
    bigRegStore(:,i) = b;
    s.(desigNames{i}).RegValSigAlt = bintjm;
    %store if significant changes.
    if bintr(1,1)*bintr(1,2) > 0
        disp('Significant Y Intercept Change!')
        s.RegressionValueSig(i,1) = 1; %store 1 for significant y intercept change
        s.RegressionValueSig(i,2) = sign(bintr(1,1)); %store sign of change
    else
        disp('Insignificant Y Intercept')
        s.RegressionValueSig(i,1) = 0; %store 1 for significant y intercept change
        s.RegressionValueSig(i,2) = 0;
    end
    %now look at slope. 
    if bintr(2,1) > 1 %this indicates range is above 1
        disp('Significant Positive Slope Change')
        s.RegressionValueSig(i,3) = 1;
        s.RegressionValueSig(i,4) = 1;
    elseif bintr(2,2) < 1 %this indicates range is below 1
        disp('Significant Negative Slope Change')
        s.RegressionValueSig(i,3) = 1;
        s.RegressionValueSig(i,4) = -1;
    else
        disp('No Slope Change')
        s.RegressionValueSig(i,3) = 0;
        s.RegressionValueSig(i,4) = 0;
    end
%     findSig = find(s.(desigNames{i}).BinSigVals(2:end,:,analysisWindow)<sigVal);
%     findSigLaser = find(s.(desigNames{i}).BinSigValsLaser(2:end,:,analysisWindow)<sigVal);
%     %find intersection
%     sigInter = intersect(findSig,findSigLaser);
%     sigInter
    %now pull those values from binned amounts. 
%     if length(sigInter) < minResp
%         disp('Insufficient Points for Linear Regression')
%         s.(desigNames{i}).RegRawVals = [];
%         s.(desigNames{i}).RegVals = [];
%         s.(desigNames{i}).RegValSig = [];
%         s.(desigNames{i}).RegValSigAlt = [];
%         s.RegressionValueSig(i,1) = 0;
%         s.RegressionValueSig(i,2) = 0;
%         s.RegressionValueSig(i,3) = 0;
%         s.RegressionValueSig(i,4) = 0;
%     else
%         disp('Sufficient Points for Linear Regression')
%         valStore = [];
%         tester = s.(desigNames{i}).BinTone(2:end,:);
%         valStore(:,1) = tester(sigInter);
%         tester = s.(desigNames{i}).BinToneLaser(2:end,:);
%         valStore(:,2) = tester(sigInter);
%         
%         [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal); %b1 is for y intercept, b2 is slope. 
%         %for bintr, first row is range for intercept, second row is for range of slope
% %         b
% %         bintr
%         s.(desigNames{i}).RegRawVals = valStore;
%         s.(desigNames{i}).RegVals = b;
%         s.(desigNames{i}).RegValSig = bintr;
%         s.(desigNames{i}).RegValSigAlt = bintjm;
%         %store if significant changes.
%         if bintr(1,1)*bintr(1,2) > 0
%             disp('Significant Y Intercept Change!')
%             s.RegressionValueSig(i,1) = 1; %store 1 for significant y intercept change
%             s.RegressionValueSig(i,2) = sign(bintr(1,1)); %store sign of change
%         else
%             disp('Insignificant Y Intercept')
%             s.RegressionValueSig(i,1) = 0; %store 1 for significant y intercept change
%             s.RegressionValueSig(i,2) = 0;
%         end
%         %now look at slope. 
%         if bintr(2,1) > 1 %this indicates range is above 1
%             disp('Significant Positive Slope Change')
%             s.RegressionValueSig(i,3) = 1;
%             s.RegressionValueSig(i,4) = 1;
%         elseif bintr(2,2) < 1 %this indicates range is below 1
%             disp('Significant Negative Slope Change')
%             s.RegressionValueSig(i,3) = 1;
%             s.RegressionValueSig(i,4) = -1;
%         else
%             disp('No Slope Change')
%             s.RegressionValueSig(i,3) = 0;
%             s.RegressionValueSig(i,4) = 0;
%         end
%             
%     end
end
masterData(:,masterInd:masterInd+3) = zeros(numUnits,4);
masterData(:,masterInd:masterInd+3) = s.RegressionValueSig; 
masterData(:,masterInd+4:masterInd+5) = zeros(numUnits,2);
masterData(:,masterInd + 4:masterInd + 5) =  bigRegStore';
masterHeader{masterInd} = 'RegYIntSig'; 
masterHeader{masterInd+1} = 'RegYIntDir';
masterHeader{masterInd+2} = 'RegSlopeSig';
masterHeader{masterInd+3} = 'RegSlopeDir';
masterHeader{masterInd + 4} = 'YintVal';
masterHeader{masterInd + 5} = 'SlopeVal';
masterInd = masterInd + 6;
%% Laser Analysis

laserBinVect = [s.Parameters.LaserWindow(1):s.Parameters.LaserBin:s.Parameters.LaserWindow(2)-s.Parameters.LaserBin];

for i = 1:numUnits
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    alignTimes = dio2Times;
    %now align spikes trial by trial
    [laserRasters] = functionBasicRaster(spikeTimes,alignTimes,s.Parameters.LaserWindow);
    %generate overall histogram and compensate
    [histCounts histCenters] = hist(laserRasters(:,1),laserBinVect);
    laserHistData = histCounts'/length(alignTimes)/s.Parameters.LaserBin;
    %now pull trial by trial to see if there are spikes during the
    %laser period.
%     laserResp = zeros(length(alignTimes),1);
%     for j = 1:length(alignTimes)
%         %find if there are any spikes
%         laseFind = find(laserRasters(:,2) == j & laserRasters(:,1) > 0 & laserRasters(:,1) < s.Parameters.LaserLim);
%         if laseFind
%             laserResp(j) = 1;
%         end
%     end
    s.(desigNames{i}).LaserRasters = laserRasters;
    s.(desigNames{i}).LaserHist = laserHistData;
%     s.(desigNames{i}).LaserResps = laserResp;
    %now get some simple response data from laser. 
    preBin = length(find(laserRasters(:,1) >= s.Parameters.LaserAnalysis(1,1) & laserRasters(:,1) <= s.Parameters.LaserAnalysis(1,2)));
    postBin = length(find(laserRasters(:,1) >= s.Parameters.LaserAnalysis(2,1) & laserRasters(:,1) <= s.Parameters.LaserAnalysis(2,2)));
    modInd = (postBin - preBin)/(postBin + preBin);
    masterData(i,masterInd) = modInd; masterHeader{masterInd} = 'ModulationIndex'; 
end
masterInd = masterInd + 1;

%% Saving
save(fullfile(pname,fname),'s','masterData','masterHeader');

%% Plotting!!

%% first plot general figure. 
[indCellType] = functionCellStringFind(masterHeader,'CellType');
indCellType = indCellType(1);
[indModInd] = functionCellStringFind(masterHeader,'ModulationIndex');
indModInd = indModInd(1);
%determine if there are units in each category
findPVs = find(masterData(:,indCellType) == 1);
if findPVs
    for i = 1:length(findPVs)
        pvStores(:,i) = s.(s.DesignationName{findPVs(i)}).SessionFiring;
    end
    avPV = mean(pvStores');
end

findMSNs = find(masterData(:,indCellType) == 0);
if findMSNs
    for i = 1:length(findMSNs)
        msnStores(:,i) = s.(s.DesignationName{findMSNs(i)}).SessionFiring;
    end
    avMSN = mean(msnStores');
end

[indPkTr] = functionCellStringFind(masterHeader,'PeakTrough');
[indISI] = functionCellStringFind(masterHeader,'isiCov');
[indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
[indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');


hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
%Column 1
%plot spike width vs coefficient of variation
subplot(6,4,1)
hold on

plot(masterData(:,indPkTr),masterData(:,indISI),'k.')
plot(masterData(masterData(:,indCellType) == 1,indPkTr),masterData(masterData(:,indCellType) == 1,indISI),'r.')
plot(masterData(masterData(:,indCellType) == 2,indPkTr),masterData(masterData(:,indCellType) == 2,indISI),'g.')
xlabel('Peak Trough (ms)')
ylabel('ISI Coefficient of Variation')
title(fileName,'fontweight','bold', 'Interpreter', 'none');

%plot loco speed vs firing rates
subplot(6,4,5)
hold on
plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'g')
if findMSNs
    plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],avMSN/max(avMSN),'k')
end
if findPVs
    plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],avPV/max(avPV),'r')
end
xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
ylim([-0.1,1])

%plot out proportion of each cell type

subplot(3,4,5)
cellDist = [length(findMSNs),length(findPVs)];
pie(cellDist)
labels = {'MSNs','PVs'};
detZero = find(cellDist == 0);
labels(detZero) = [];
legend(labels,'Location','southoutside','Orientation','horizontal')


%Column 2
if findMSNs
    subplot(3,4,2)

    hold on
    holder = masterData(findMSNs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    posResp = find(holder(:,1) == 1);


    plot(s.NonLaserOverall.PosWidths(:,findMSNs(posResp),1),'k.')
    plot(mean(s.NonLaserOverall.PosWidths(:,findMSNs(posResp),1)'),'k-')
    plot(s.NonLaserOverall.PosWidths(:,findMSNs(posResp),2),'b.')
    plot(mean(s.NonLaserOverall.PosWidths(:,findMSNs(posResp),2)'),'b-')
    plot(s.NonLaserOverall.PosWidths(:,findMSNs(posResp),3),'m.')
    plot(mean(s.NonLaserOverall.PosWidths(:,findMSNs(posResp),3)'),'m-')
    xlim([0 size(s.NonLaserOverall.PosWidths,1) + 1])
    title(strcat(num2str(length(findMSNs)),'-MSN Tuning Width Responses fast(k) tone(b) gen(m)'))


    %plot distribution of positive, negative, both, and untuned units
    subplot(3,4,6)

    holder = masterData(findMSNs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    det = holder(:,1) + holder(:,2);
    det = hist(det,[-2:1:1]);

    pie(det)
    labels = {'Neg','Mix','None','Pos'};
    detZero = find(det == 0);
    labels(detZero) = [];
    legend(labels,'Location','southoutside','Orientation','horizontal')



    subplot(3,4,10)
    hold on
    plot(s.NonLaserOverall.PosWidths(:,findMSNs(posResp),2),s.LaserOverall.PosWidths(:,findMSNs(posResp),2),'r.')
    plot([0 max(max(max(s.LaserOverall.PosWidths)))],[0 max(max(max(s.LaserOverall.PosWidths)))],'k')
    title('MSN Width Change with Laser x normal y laser')

    %plot out modulation index. 
    subplot(3,4,3)
    hold on
    hist(masterData(findMSNs,indModInd),[-1:0.1:1])
    title('Laser Mod Index MSNs')

    %plot out average change to firing rate of responses
    subplot(3,4,9)
    hold on
    for i = 1:length(findMSNs)
        %first find significant responses in baseline
        findBaseSig = find(s.(desigNames{findMSNs(i)}).BinSigVals(:,:,3) <= 0.05); 
        sigFilter = NaN(size(s.(desigNames{findMSNs(i)}).BinSigVals(:,:,3)));
        sigFilter(findBaseSig) = 1;
        preResp = s.(desigNames{findMSNs(i)}).BinDiff(:,:,3).*sigFilter;
        postResp = s.(desigNames{findMSNs(i)}).BinDiffLaser(:,:,3).*sigFilter;
        ratioVal(i) = nanmean(nanmean((postResp-preResp)./(postResp+preResp)));
        sigValNum(i) = nansum(nansum(sigFilter));
    end
    plot(ratioVal,sigValNum,'ko')
    xlim([-1.2 1.2])
    title('MSN ModIndexSigResponse Gen Window (x) vs NumSigResp (y)')
end
%Column 3 PV CELLS
if findPVs
    subplot(3,4,4)
    hold on
    holder = masterData(findPVs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    posResp = find(holder(:,1) == 1);

    plot(s.NonLaserOverall.PosWidths(:,findPVs(posResp),1),'k.')
    plot(mean(s.NonLaserOverall.PosWidths(:,findPVs(posResp),1)'),'k-')
    plot(s.NonLaserOverall.PosWidths(:,findPVs(posResp),2),'b.')
    plot(mean(s.NonLaserOverall.PosWidths(:,findPVs(posResp),2)'),'b-')
    plot(s.NonLaserOverall.PosWidths(:,findPVs(posResp),3),'m.')
    plot(mean(s.NonLaserOverall.PosWidths(:,findPVs(posResp),3)'),'m-')
    xlim([0 size(s.NonLaserOverall.PosWidths,1) + 1])
    title(strcat(num2str(length(findPVs)),'-PV Tuning Width Responses fast(k) tone(b) gen(m)'))

    %plot distribution of positive, negative, both, and untuned units
    subplot(3,4,8)
    holder = masterData(findPVs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    det = holder(:,1) + holder(:,2);
    det = hist(det,[-2:1:1]);
    pie(det)
    labels = {'Neg','Mix','None','Pos'};
    detZero = find(det == 0);
    labels(detZero) = [];
    legend(labels,'Location','southoutside','Orientation','horizontal')
    
    subplot(3,4,12)
    hold on
    plot(s.NonLaserOverall.PosWidths(:,findPVs(posResp),2),s.LaserOverall.PosWidths(:,findPVs(posResp),2),'r.')
    plot([0 max(max(max(s.LaserOverall.PosWidths)))],[0 max(max(max(s.LaserOverall.PosWidths)))],'k')
    title('PV Width Change with Laser x normal y laser')
    
    subplot(3,4,7)
    hold on
    hist(masterData(findPVs,indModInd),[-1:0.1:1])
    title('Laser Mod Index PVs')
    
    subplot(3,4,11)
    hold on
    ratioVal = [];
    sigValNum = [];
    for i = 1:length(findPVs)
        %first find significant responses in baseline
        findBaseSig = find(s.(desigNames{findPVs(i)}).BinSigVals(:,:,3) <= 0.05); 
        sigFilter = NaN(size(s.(desigNames{findPVs(i)}).BinSigVals(:,:,3)));
        sigFilter(findBaseSig) = 1;
        preResp = s.(desigNames{findPVs(i)}).BinDiff(:,:,3).*sigFilter;
        postResp = s.(desigNames{findPVs(i)}).BinDiffLaser(:,:,3).*sigFilter;
        ratioVal(i) = nanmean(nanmean((postResp-preResp)./(postResp+preResp)));
        sigValNum(i) = nansum(nansum(sigFilter));
    end
    plot(ratioVal,sigValNum,'k.')
    xlim([-1.2 1.2])
    title('PV ModIndexSigResponse Gen Window (x) vs NumSigResp (y)')
    
    
end


spikeGraphName = strcat(fileName,'GeneralFigureCellTypes');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%% Plot differentiating by shank
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])

%plot response to white noise as heatmaps for both shanks
whiteStore = [];
genStore = [];
zeroPoint = find(histBinVector < 0,1,'last');
%z score these responses so that things will display a bit more nicely. Or
%alternatively use normalized rates. switch to normalized rates. 
for i = 1:numUnits
%     if std(s.(desigNames{i}).FrequencyHistograms(1,1:zeroPoint)) > 0
    if mean(s.(desigNames{i}).FrequencyHistograms(1,1:zeroPoint)) > 0
%         whiteStore(:,i) = (s.(desigNames{i}).FrequencyHistograms(1,:)-mean(s.(desigNames{i}).FrequencyHistograms(1,1:zeroPoint)))/std(s.(desigNames{i}).FrequencyHistograms(1,1:zeroPoint));
        whiteStore(:,i) = (s.(desigNames{i}).FrequencyHistograms(1,:)-mean(s.(desigNames{i}).FrequencyHistograms(1,1:zeroPoint)))/(max((s.(desigNames{i}).FrequencyHistograms(1,:))-mean(s.(desigNames{i}).FrequencyHistograms(1,1:zeroPoint))));
    else
        whiteStore(:,i) = NaN(length(s.(desigNames{i}).FrequencyHistograms(1,:)),1);
    end
%     if std(mean(s.(desigNames{i}).FrequencyHistograms(2:end,1:zeroPoint)))
    if mean(mean(s.(desigNames{i}).FrequencyHistograms(2:end,1:zeroPoint)))
%         genStore(:,i) = (mean(s.(desigNames{i}).FrequencyHistograms(2:end,:))-mean(mean(s.(desigNames{i}).FrequencyHistograms(2:end,1:zeroPoint))))/std(mean(s.(desigNames{i}).FrequencyHistograms(2:end,1:zeroPoint)));
        genStore(:,i) = (mean(s.(desigNames{i}).FrequencyHistograms(2:end,:))-mean(mean(s.(desigNames{i}).FrequencyHistograms(2:end,1:zeroPoint))))/(max(mean(s.(desigNames{i}).FrequencyHistograms(2:end,:)))-mean(mean(s.(desigNames{i}).FrequencyHistograms(2:end,1:zeroPoint))));
    else
        genStore(:,i) = NaN(length(s.(desigNames{i}).FrequencyHistograms(1,:)),1);
    end
    %pull overall maps
    fullDiffs(:,:,i) = squeeze(mean(s.(desigNames{i}).BinDiff,2));
end

[indShankDes] = functionCellStringFind(masterHeader,'Shank Designation');
findShank1 = find(masterData(:,indShankDes) == 1);
map = [0 0 1];
for i = 2:21
    if i < 11
        map(i,:) = map(i-1,:);
        map(i,1:2) = map(i,1:2) + 0.1;
    elseif i == 11
        map(i,:) = [1 1 1];
    elseif i > 11
        map(i,:) = map(i-1,:);
        map(i,2:3) = map(i,2:3) - 0.1;
    end
end

%find PV cells on shank 1
findPVshank1 = intersect(findPVs,findShank1);
% findCHATshank1 = intersect(findCHATs,findShank1);

subplot(2,5,1)
imagesc((whiteStore(:,s.SortedPeakWaveOrder(findShank1))'),[-1 1])
hold on
for i = 1:length(findPVshank1)
    plot([0 length(whiteStore)],[s.SortedPeakWaveOrder(findPVshank1(i)) s.SortedPeakWaveOrder(findPVshank1(i))],'r','LineWidth',2)
end
% for i = 1:length(findCHATshank1)
%     plot([0 length(whiteStore)],[s.SortedPeakWaveOrder(findCHATshank1(i)) s.SortedPeakWaveOrder(findCHATshank1(i))],'g','LineWidth',2)
% end
colormap(map)
colorbar
title('Shank 1 Resp to WN, BaseNorm FR LogScale')

subplot(2,5,2)
imagesc((genStore(:,s.SortedPeakWaveOrder(findShank1))'),[-1 1])
colormap(map)
colorbar
title('Resp to Tones, BaseNorm FR LogScale')

subplot(2,5,3)
imagesc(squeeze(fullDiffs(:,1,s.SortedPeakWaveOrder(findShank1)))')
colormap parula
colorbar
title(strcat(fileName,'Bin Subtracted Fast Period'))
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

subplot(2,5,4)
imagesc(squeeze(fullDiffs(:,2,s.SortedPeakWaveOrder(findShank1)))')
colormap parula
colorbar
title('Bin Subtracted Tone Period')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

subplot(2,5,5)
imagesc(squeeze(fullDiffs(:,3,s.SortedPeakWaveOrder(findShank1)))')
colormap parula
colorbar
title('Bin Subtracted Gen Period')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));



findShank2 = find(masterData(:,indShankDes) == 2);

%find PV cells on shank 2
findPVshank2 = intersect(findPVs,findShank2);
% findCHATshank2 = intersect(findCHATs,findShank2);

subplot(2,5,6)
imagesc((whiteStore(:,s.SortedPeakWaveOrder(findShank2))'),[-1 1])
hold on
for i = 1:length(findPVshank2)
    plot([0 length(whiteStore)],[s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1) s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1)],'r','LineWidth',2)
end
% for i = 1:length(findCHATshank2)
%     plot([0 length(whiteStore)],[s.SortedPeakWaveOrder(findCHATshank2(i))-length(findShank1) s.SortedPeakWaveOrder(findCHATshank2(i))-length(findShank1)],'g','LineWidth',2)
% end
colormap(map)
colorbar
title('Shank 2 Resp to WN, BaseNorm FR LogScale')

subplot(2,5,7)
imagesc((genStore(:,s.SortedPeakWaveOrder(findShank2))'),[-1 1])
colormap(map)
colorbar
title('Resp to Tones, BaseNorm FR LogScale')

subplot(2,5,8)
imagesc(squeeze(fullDiffs(:,1,s.SortedPeakWaveOrder(findShank2)))')
colormap parula
colorbar
title('Bin Subtracted Fast Period')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

subplot(2,5,9)
imagesc(squeeze(fullDiffs(:,2,s.SortedPeakWaveOrder(findShank2)))')
colormap parula
colorbar
title('Bin Subtracted Tone Period')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

subplot(2,5,10)
imagesc(squeeze(fullDiffs(:,3,s.SortedPeakWaveOrder(findShank2)))')
colormap parula
colorbar
title('Bin Subtracted Gen Period')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

spikeGraphName = strcat(fileName,'ProbeTopology');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% now plot individual cells. 
for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])

    % Column 1
    %plots average waveform
    subplot(4,8,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    if ismember(i,findPVs)
        title(strcat('PV AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    elseif ismember(i,findMSNs)
        title(strcat('MSN AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
%     elseif ismember(i,findCHATs)
%         title(strcat('CHAT AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    else
        title(strcat('UNK AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    end

    %plots ISI
    subplot(4,8,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})

    % plot histogram.
    subplot(4,4,5)
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
    histLims = max([max(s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation),max(s.(desigNames{i}).AllHistogramsLaser + s.(desigNames{i}).HistogramStandardDeviationLaser)]);
    ylim([0 histLims])
    title('Histogram')

    %plot out rasters, organized!
    subplot(2,4,5)
    plot(s.(desigNames{i}).AllRasters(:,1),...
        s.(desigNames{i}).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    rasterFreqLines = zeros(numFreqs,2);
    if numDBs ==1
        rasterFreqLines(:,1) = cumsum((matrixTrialNum'));
    elseif numDBs > 1
        rasterFreqLines(:,1) = cumsum(sum(matrixTrialNum'));
    end

    rasterFreqLines(:,2) = uniqueFreqs;
    %this generates green lines separating by Frequency
    tempHold = 1;
    for k = 1:size(uniqueFreqs,1)
        plot(s.Parameters.RasterWindow,[tempHold+sum(matrixTrialNum(k,:)) tempHold+sum(matrixTrialNum(k,:))],'g','LineWidth',1)
        tempHold = tempHold + sum(matrixTrialNum(k,:));
    end
    set(gca,'YTick',rasterFreqLines(:,1));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 totalTrialNum])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    title('Descending = increase in amplitude and freq')

    % Column 2
    %plot FR and velocity
    subplot(4,4,2)
    hold on
    plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'b')
    plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
    xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
    ylim([-0.1,1])
    title({fileName;desigNames{i}},'fontweight','bold', 'Interpreter', 'none');
    
    % plot histogram for LASER trials
    subplot(4,4,6)
    plot(histBinVector,s.(desigNames{i}).AllHistogramsLaser,'k','LineWidth',2)
    hold on
    plot(histBinVector,s.(desigNames{i}).AllHistogramsLaser - s.(desigNames{i}).HistogramStandardDeviationLaser,'b','LineWidth',1)
    plot(histBinVector,s.(desigNames{i}).AllHistogramsLaser + s.(desigNames{i}).HistogramStandardDeviationLaser,'b','LineWidth',1)
    %plot significant values
    plot(s.(desigNames{i}).AllHistogramSigLaser.Centers(...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(:,3) == 1),...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(:,3) == 1,1),...
        'b*')
    plot(s.(desigNames{i}).AllHistogramSigLaser.Centers(...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(:,4) == 1),...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(:,4) == 1,1),...
        'bo')
    %plot negative values for first tuning
    plot(s.(desigNames{i}).AllHistogramSigLaser.Centers(...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(:,6) == 1),...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(...
        s.(desigNames{i}).AllHistogramSigLaser.Histogram(:,6) == 1,1),...
        'k*')
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    ylim([0 histLims])
    title('Histogram')

    %plot out rasters, organized! FOR LASER TRIALS
    subplot(2,4,6)
    plot(s.(desigNames{i}).AllRastersLaser(:,1),...
        s.(desigNames{i}).AllRastersLaser(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    rasterFreqLines = zeros(numFreqs,2);
    if numDBs ==1
        rasterFreqLines(:,1) = cumsum((matrixTrialNum'));
    elseif numDBs > 1
        rasterFreqLines(:,1) = cumsum(sum(matrixTrialNum'));
    end

    rasterFreqLines(:,2) = uniqueFreqs;
    %this generates green lines separating by Frequency
    tempHold = 1;
    for k = 1:size(uniqueFreqs,1)
        plot(s.Parameters.RasterWindow,[tempHold+sum(matrixTrialNum(k,:)) tempHold+sum(matrixTrialNum(k,:))],'g','LineWidth',1)
        tempHold = tempHold + sum(matrixTrialNum(k,:));
    end
    set(gca,'YTick',rasterFreqLines(:,1));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 totalTrialNum])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    title('Descending = increase in amplitude and freq')
    
    % Column 3: binned responses and tuning curves
%     clims = [min([min(min(s.(desigNames{i}).BinTone)),min(min(s.(desigNames{i}).BinToneLaser))]),...
%         max([max(max(s.(desigNames{i}).BinTone)),max(max(s.(desigNames{i}).BinToneLaser))])];
%     %Plot binned response during tone period %%%%HERE WE WANT OUTLINE
%     subplot(4,4,3)
%     imagesc(s.(desigNames{i}).BinDiff(:,:,2)',clims)
%     colormap(parula)
%     colorbar
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     set(gca,'YTick',dbRange(:,2));
%     set(gca,'YTickLabel',dbRange(:,1));
%     title('Mean Binned Response (tone), nL');
%     
%     %Plot binned response during tone period LASER
%     subplot(4,4,7)
%     imagesc(s.(desigNames{i}).BinDiffLaser(:,:,2)',clims)
%     colormap(parula)
%     colorbar
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     set(gca,'YTick',dbRange(:,2));
%     set(gca,'YTickLabel',dbRange(:,1));
%     title('Mean Binned Response (tone), LASER');
%     
%     clims = [min([min(min(s.(desigNames{i}).BinGen)),min(min(s.(desigNames{i}).BinGenLaser))]),...
%         max([max(max(s.(desigNames{i}).BinGen)),max(max(s.(desigNames{i}).BinGenLaser))])];
%     %Plot binned response during general period
%     subplot(4,4,11)
%     imagesc(s.(desigNames{i}).BinGen',clims)
%     colormap(parula)
%     colorbar
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     set(gca,'YTick',dbRange(:,2));
%     set(gca,'YTickLabel',dbRange(:,1));
%     title('Mean Binned Response (general), nL')
%     
%     %Plot binned response during general period LASER
%     subplot(4,4,15)
%     imagesc(s.(desigNames{i}).BinGenLaser',clims)
%     colormap(parula)
%     colorbar
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     set(gca,'YTick',dbRange(:,2));
%     set(gca,'YTickLabel',dbRange(:,1));
%     title('Mean Binned Response (general), LASER')
%plot out responses by frequency, with x axis being amplitude
    subplot(2,4,3)
    %plot things linked by frequency
    hold on
    plotThresh = 0.01; %significance threshold for plotting stars. 
    maxval = max([max(max(s.(desigNames{i}).BinDiff(:,:,toneTarget))),max(max(s.(desigNames{i}).BinDiffLaser(:,:,toneTarget)))]);
    plot([0 maxval],[0 maxval],'k','LineWidth',2)
    for j = 1:numFreqs
        %first determine if by linear fit there is a change to the slope or
        %yint. Only use significant values and positive values
        sigVals = find(s.(desigNames{i}).BinSigVals(j,:,toneTarget) < sigCutoff);
        sigValsLaser = find(s.(desigNames{i}).BinSigValsLaser(j,:,toneTarget) < sigCutoff);
        posVals = find(s.(desigNames{i}).BinDiff(j,:,toneTarget) > 0);
        posValsLaser = find(s.(desigNames{i}).BinDiffLaser(j,:,toneTarget) > 0);
        %now find intersect of these!
        intersectNorm = intersect(sigVals,posVals);
        intersectLaser = intersect(sigValsLaser,posValsLaser);
        fullIntersect = intersect(intersectNorm,intersectLaser);
        if length(fullIntersect) > 3
            %now pull values
            disp('At Least 3 Significant Positive Responses')
            normVals = s.(desigNames{i}).BinDiff(j,:,toneTarget);
            laserVals = s.(desigNames{i}).BinDiffLaser(j,:,toneTarget);
            normVals = normVals(fullIntersect);
            laserVals = laserVals(fullIntersect);
            [b,bintr,bintjm] = gmregress(normVals,laserVals,sigCutoff); %b1 is for y intercept, b2 is slope. 
            %for bintr, first row is range for intercept, second row is for range of slope

            %if significant slope, then plot the line thicker
            if 1 >= bintr(2,1) && 1 <= bintr(2,2) %this is the case in which 1 falls within the confidence bounds, indicating insignificant change in slope
                 plot(normVals,laserVals,'g.-','Color',[0 j/numFreqs 0])
            else
                plot(normVals,laserVals,'g.-','LineWidth',2,'Color',[0 j/numFreqs 0])
            end

            %if intercept is significant change, mark this
            if 0 >= bintr(1,1) && 0 <= bintr(1,2) %this is the case in which 1 falls within the confidence bounds, indicating insignificant change in slope
                 plot(s.(desigNames{i}).BinDiff(j,1,toneTarget),s.(desigNames{i}).BinDiffLaser(j,1,toneTarget),'o','Color',[0 j/numFreqs 0])
            else
                plot(s.(desigNames{i}).BinDiff(j,1,toneTarget),s.(desigNames{i}).BinDiffLaser(j,1,toneTarget),'o','Color',[0 j/numFreqs 0])
            end
            %now see if differences significant
            for k = 1:numDBs
                nlVal = s.(desigNames{i}).LatPeakBin{j, k}.BinnedSpikesTone - s.(desigNames{i}).LatPeakBin{j, k}.BinnedSpikesToneBase;
                laserVal = s.(desigNames{i}).LatPeakBinLaser{j, k}.BinnedSpikesTone - s.(desigNames{i}).LatPeakBinLaser{j, k}.BinnedSpikesToneBase;
                testVal = ranksum(nlVal,laserVal);
                if testVal < plotThresh
                    plot(s.(desigNames{i}).BinDiff(j,k,toneTarget),s.(desigNames{i}).BinDiffLaser(j,k,toneTarget),'r*')
                end
            end
        end
    end
    title('Plotting Responses linked by Frequency Value')
    
    %plot things linked by amplitude
    subplot(2,4,7)
    hold on
    plotThresh = 0.01; %significance threshold for plotting stars. 
    maxval = max([max(max(s.(desigNames{i}).BinDiff(:,:,toneTarget))),max(max(s.(desigNames{i}).BinDiffLaser(:,:,toneTarget)))]);
    plot([0 maxval],[0 maxval],'k','LineWidth',2)
    for j = 1:numDBs
        %first determine if by linear fit there is a change to the slope or
        %yint. Only use significant values and positive values
        sigVals = find(s.(desigNames{i}).BinSigVals(:,j,toneTarget) < sigCutoff);
        sigValsLaser = find(s.(desigNames{i}).BinSigValsLaser(:,j,toneTarget) < sigCutoff);
        posVals = find(s.(desigNames{i}).BinDiff(:,j,toneTarget) > 0);
        posValsLaser = find(s.(desigNames{i}).BinDiffLaser(:,j,toneTarget) > 0);
        %now find intersect of these!
        intersectNorm = intersect(sigVals,posVals);
        intersectLaser = intersect(sigValsLaser,posValsLaser);
        fullIntersect = intersect(intersectNorm,intersectLaser);
        if length(fullIntersect) > 3
            disp('At Least 3 Positive Significant Responses!')
            %now pull values
            normVals = s.(desigNames{i}).BinDiff(:,j,toneTarget);
            laserVals = s.(desigNames{i}).BinDiffLaser(:,j,toneTarget);
            normVals = normVals(fullIntersect);
            laserVals = laserVals(fullIntersect);
            [b,bintr,bintjm] = gmregress(normVals,laserVals,sigCutoff); %b1 is for y intercept, b2 is slope. 
            %for bintr, first row is range for intercept, second row is for range of slope

            %if significant slope, then plot the line thicker
            if 1 >= bintr(2,1) && 1 <= bintr(2,2) %this is the case in which 1 falls within the confidence bounds, indicating insignificant change in slope
                 plot(normVals,laserVals,'g.-','Color',[0 j/numDBs 0])
            else
                plot(normVals,laserVals,'g.-','LineWidth',2,'Color',[0 j/numDBs 0])
            end

            %if intercept is significant change, mark this
            if 0 >= bintr(1,1) && 0 <= bintr(1,2) %this is the case in which 1 falls within the confidence bounds, indicating insignificant change in slope
                 plot(normVals,laserVals,'o','Color',[0 j/numDBs 0])
            else
                plot(normVals,laserVals,'o','Color',[0 j/numDBs 0])
            end
            %now see if differences significant
            for k = 1:numFreqs
                nlVal = s.(desigNames{i}).LatPeakBin{k,j}.BinnedSpikesTone - s.(desigNames{i}).LatPeakBin{k,j}.BinnedSpikesToneBase;
                laserVal = s.(desigNames{i}).LatPeakBinLaser{k,j}.BinnedSpikesTone - s.(desigNames{i}).LatPeakBinLaser{k,j}.BinnedSpikesToneBase;
                testVal = signrank(nlVal,laserVal);
                if testVal < plotThresh
                    plot(s.(desigNames{i}).BinDiff(k,j,toneTarget),s.(desigNames{i}).BinDiffLaser(k,j,toneTarget),'r*')
                end
            end
        end
    end
    title('Plotting Responses linked by DB Level')

    % Column 4

%     %plot out binned spikes (fast)
%     subplot(4,4,4)
%     hold on
%     for cInd = 1:numDBs
%         plot(s.(desigNames{i}).BinDiff(:,cInd,1),'Color',[cInd/numDBs 0 0])
%         %find significant points, by p < 0.05
%         findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,1)<0.05);
%         plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,1),'r*')
%         %find significant points, by p < 0.01
%         findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,1)<0.01);
%         plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,1),'ro')
%         %now plot laser related. 
%         plot(s.(desigNames{i}).BinDiffLaser(:,cInd,1),'Color',[0 cInd/numDBs 0])
%         %find significant points, by p < 0.05
%         findSigs = find(s.(desigNames{i}).BinSigValsLaser(:,cInd,1)<0.05);
%         plot(findSigs,s.(desigNames{i}).BinDiffLaser(findSigs,cInd,1),'g*')
%         %find significant points, by p < 0.01
%         findSigs = find(s.(desigNames{i}).BinSigValsLaser(:,cInd,1)<0.01);
%         plot(findSigs,s.(desigNames{i}).BinDiffLaser(findSigs,cInd,1),'go')
%     end
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     ylabel('Binned Spikes/Fast Period')
%     title(strcat('Curve Fast Period, nl width:',num2str(s.NonLaserOverall.PosWidths(3,i,1)),',laser width:',num2str(s.LaserOverall.PosWidths(3,i,1))))
    subplot(4,4,4)
    hold on
    for j = 1:numFreqs
        plot((s.(desigNames{i}).BinDiffLaser(j,:,2)-s.(desigNames{i}).BinDiff(j,:,2)),'LineWidth',2,'Color',[j/numFreqs 0 0])
    end
    %plot out binned spikes (tone)
    subplot(4,4,8)
    hold on
    for cInd = 1:numDBs
        plot(s.(desigNames{i}).BinDiff(:,cInd,2),'Color',[cInd/numDBs 0 0])
        %find significant points, by p < 0.05
        findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,2)<0.05);
        plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,2),'r*')
        %find significant points, by p < 0.01
        findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,2)<0.01);
        plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,2),'ro')
        %now plot laser related. 
        plot(s.(desigNames{i}).BinDiffLaser(:,cInd,2),'Color',[0 cInd/numDBs 0])
        %find significant points, by p < 0.05
        findSigs = find(s.(desigNames{i}).BinSigValsLaser(:,cInd,2)<0.05);
        plot(findSigs,s.(desigNames{i}).BinDiffLaser(findSigs,cInd,2),'g*')
        %find significant points, by p < 0.01
        findSigs = find(s.(desigNames{i}).BinSigValsLaser(:,cInd,2)<0.01);
        plot(findSigs,s.(desigNames{i}).BinDiffLaser(findSigs,cInd,2),'go')
    end
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    ylabel('Binned Spikes/Fast Period')
    title(strcat('Curve Tone Period, nl width:',num2str(s.NonLaserOverall.PosWidths(end,i,2)),',laser width:',num2str(s.LaserOverall.PosWidths(end,i,2))))
    
    %plot out binned responses to general period, non laser vs laser
    subplot(4,4,12)
    hold on
    plot(reshape(s.(desigNames{i}).BinDiff(:,:,3),[],1),reshape(s.(desigNames{i}).BinDiffLaser(:,:,3),[],1),'r.')
    plot([0 max(max(s.(desigNames{i}).BinDiffLaser(:,:,3)))],[0 max(max(s.(desigNames{i}).BinDiffLaser(:,:,3)))],'k')
    %now plot regression line?
    if ~isnan(s.(desigNames{i}).RegVals(1)) & ~isnan(s.(desigNames{i}).RegVals(2))
        plot([0 max(max(s.(desigNames{i}).BinDiffLaser(:,:,3)))],[s.(desigNames{i}).RegVals(1) max(max(s.(desigNames{i}).BinDiffLaser(:,:,3)))*s.(desigNames{i}).RegVals(2)],'c')
        if s.(desigNames{i}).RegValSig(1,1) > 0
            plot(0,s.(desigNames{i}).RegVals(1),'c*')
        elseif s.(desigNames{i}).RegValSig(1,2) < 0
            plot(0,s.(desigNames{i}).RegVals(1),'c*')
        end
    end
    
    if s.(desigNames{i}).RegValSig(2,1) > 1
        title('Bin Tone nl(x) vs lsr(y) Sig Pos Slope')
    elseif s.(desigNames{i}).RegValSig(2,2) <1
        title('Bin Tone nl(x) vs lsr(y) Sig Neg Slope')
    end
%     title('normal (x) vs Laser (y)')
    axis equal
%     %plot rasters to laser
%     subplot(4,4,12)
%     hold on
%     plot(s.(desigNames{i}).LaserRasters(:,1),...
%         s.(desigNames{i}).LaserRasters(:,2),'k.','markersize',4)
%     xlim(s.Parameters.LaserWindow)
%     ylim([0 trialNumLaser])
%     title('Raster Aligned to Laser Onset')
    
%     %plot histogram to laser onset.
%     subplot(4,4,16)
%     plot(laserBinVect,s.(desigNames{i}).LaserHist)
%     xlim(s.Parameters.LaserWindow)
%     title('Histogram Aligned to Laser Onset')
    
    subplot(4,4,16)
    hold on
    newWin = [-0.02 0.04];
    newVector = [newWin(1):0.0005:newWin(2)];
    testRast = s.(desigNames{i}).AllRasters(s.(desigNames{i}).AllRasters >= newWin(1) & s.(desigNames{i}).AllRasters <= newWin(2) ,1);
    plot(newVector,smooth(hist(testRast,newVector),3),'r')
    testRast = s.(desigNames{i}).AllRastersLaser(s.(desigNames{i}).AllRastersLaser >= newWin(1) & s.(desigNames{i}).AllRastersLaser <= newWin(2) ,1);
    plot(newVector,smooth(hist(testRast,newVector),3),'b')
    xlim(newWin)
    title('Plot of Latency')
    
    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
    close

end



end