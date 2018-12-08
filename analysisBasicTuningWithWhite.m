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
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.
toggleROC = 0; %toggle for tuning on/off ROC analysis

%PARAMETERS FOR BASIC ARRANGEMENT OF DATA
s.Parameters.RasterWindow = [-5 6]; %ratio for raster window. will be multiplied by toneDur
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
s.Parameters.LaserWindow = [-0.1 0.1]; %plotting window for rasters
s.Parameters.LaserBin = 0.001; %histogram bin size
s.Parameters.LaserLim = 0.015; %maximum lag value for calculation.

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
    s.DeletedUnits = n.s.DeletedUnits;
    %now delete bad units
    delUnits = fields(s.DeletedUnits);
    for indCount = 1:length(delUnits)
        s = rmfield(s,delUnits{indCount});
    end
    %replace designation array and names. 
    s.DesignationArray = n.s.DesignationArray;
    s.DesignationName = n.s.DesignationName;
    
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
%     for i = 1:size(octaveRange,1);
%         octaveRange(i,2) = find(round(uniqueFreqs) == octaveRange(i,1));
%     end
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

%check for tons of DIO2. if alot this should be an ID tuning.
if length(dio2Times) > 10
    idToggle = 1;
    disp('Laser Pulses Detected: ID Session')
else
    idToggle = 0;
    disp('No Laser Pulses Detected')
end

%now lets see how many changes there are to DIO2

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
        subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(i),3));
        for j = 1:numDBs
            sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == subUniqueDB(j));
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
        subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(i),3));
        for j = 1:numDBs
            sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == subUniqueDB(j));
            master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
            sortingCounter = sortingCounter + size(sortingFinder,1);
        end
    end
    
end
disp('DIO data successfully extracted and stored.')

%170807 adding edit to allow for windowing by the time filter
%We want to remove points beyond the range of the time filter. Do this
%here.

timeFinder = find(master(:,1) > s.TimeFilterRange(2));
master(timeFinder,:) = [];
totalTrialNum = length(master);
disp('Time Filter Applied')

%now we want to determine how many trials were present for each tone. This
%should be a fairly simple nested for loop.

matrixTrialNum = zeros(numFreqs,numDBs);
for freqInd = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(freqInd),3));
    for dbInd = 1:numDBs
        matrixTrialNum(freqInd,dbInd) = length(find(master(:,2) == uniqueFreqs(freqInd) & master(:,3) == subUniqueDB(dbInd)));
    end
end

s.TrialMatrix = master;

s.TrialNumbers = matrixTrialNum;
s.RepArray = repArray;

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
    targetInd = find(s.RotaryData.Velocity(:,1)-master(i,1) > 0,1,'first');
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

%% see if EDR data exists
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

s.PosWidths = zeros(numDBs,numUnits,3); %store width of response for positive responses here. 
s.NegWidths = zeros(numDBs,numUnits,3); %store width of response for negative responses here. 
s.PosTots = zeros(numDBs,numUnits,3); %store total number of significant responses here
s.NegTots = zeros(numDBs,numUnits,3);  %store total number of significant responses here
s.PosCont = zeros(numDBs,numUnits,3); %store whether response is contiguous or not here
s.NegCont = zeros(numDBs,numUnits,3); %store whether response is contiguous or not here
s.PosEdgeWarn = zeros(numDBs,numUnits,3); %store whether significant responses abut an edge
s.NegEdgeWarn = zeros(numDBs,numUnits,3); %store whether significant responses abut an edge
s.PosGaussWidth = zeros(numDBs,numUnits,3); % will perform gaussian fit on binned spikes if there is sufficient significant values. Store half-peak width
s.RespWidthPos = zeros(numFreqs,numDBs,numUnits,2); %width of firing rate, positive
widthStore = zeros(numFreqs,numDBs,numUnits);
bigWidth = [];
for i = 1:numUnits
    masterHolder = masterInd;
    
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
    
    %store average rate into master
    masterData(i,masterHolder) = averageRate;
    masterHeader{masterHolder} = 'Average Rate';
    masterHolder = masterHolder + 1;
    
    %now pull waveform data and isi data to categorize cells
    waveForms = s.(desigNames{i}).AverageWaveForms;
    
    %find size of waveforms
    waveSize = size(waveForms);

    %find maximum waveform by max peak size
    maxWave = max(waveForms);
    [maxVal maxInd] = max(maxWave);

    %chose the big wave, interpolate to fine degree
    chosenWave = waveForms(:,maxInd);
    interpVect = [1:0.1:40];
    interpWave = interp1(1:40,chosenWave,interpVect,'spline');
    
    %now we need to find the peak. Find this starting at point 10. 

    [pkVal pkInd] = max(interpWave(100:end));
    pkInd = pkInd + 100 - 1;
    %now we need to find the minimum following the peak

    [troughVal troughInd] = min(interpWave(pkInd:end));
    troughInd = troughInd + pkInd - 1;

    peakTrough = (troughInd - pkInd)/300000;
    
    %find ISIs
    isiTimes = diff(spikeTimes);
    isiCov = std(isiTimes)/mean(isiTimes);
    
    masterData(i,masterHolder) = peakTrough;
    masterHeader{masterHolder} = 'PeakTrough(ms)';
    masterHolder = masterHolder + 1;
    
    masterData(i,masterHolder) = isiCov;
    masterHeader{masterHolder} = 'isiCov';
    masterHolder = masterHolder + 1;
    
    
    if peakTrough < s.Parameters.PVLim(1) & isiCov > s.Parameters.ChatLim
        masterData(i,masterHolder) = 1; %pv cell
    elseif peakTrough > s.Parameters.PVLim(2) & isiCov < s.Parameters.ChatLim & averageRate > 2
        masterData(i,masterHolder) = 2; %ChAT Cell
    elseif peakTrough > s.Parameters.PVLim(2) & isiCov > s.Parameters.ChatLim
        masterData(i,masterHolder) = 0; %MSN
    else
        masterData(i,masterHolder) = NaN; %label as unknown
    end
    masterHeader{masterHolder} = 'CellType';
    masterHolder = masterHolder + 1;
    
    
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
    if generalResponseHist.Warning == 0 & generalResponseHist.SigSpikePos == 1 
        s.SignificantSpikes(i) = 1;
    end
    
    masterData(i,masterHolder) = generalResponseHist.MeanBaseline;
    masterHeader{masterHolder} = 'BaselineRate';
    masterHolder = masterHolder + 1;
    
    masterData(i,masterHolder) = generalResponseHist.SigSpikePos;
    masterHeader{masterHolder} = 'PosSigGenHist';
    masterHolder = masterHolder + 1;
    
    masterData(i,masterHolder) = generalResponseHist.SigSpikeNeg;
    masterHeader{masterHolder} = 'NegSigGenHist';
    masterHolder = masterHolder + 1;
    
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
        subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(k),3));
        for l = 1:numDBs
            targetTrials = master(master(:,2) == uniqueFreqs(k) & master(:,3) == subUniqueDB(l),4); %finds the trial number for all trials of given frequency and amplitude
            findMatches = find(ismember(fullRasterData(:,2),targetTrials)); %uses these trial numbers to pull raster index
            targetRasters = fullRasterData(findMatches,:); %extracts target rasters from all rasters using the above index.
            organizedRasters{k,l} = targetRasters; %saves to organized rasters
            [histCounts,histCenters] = hist(targetRasters(:,1),histBinVector); %calculates histogram with defined bin size
            organizedHist(k,l,:) = histCounts/toneReps/s.Parameters.histBin; %saves histogram
            specHist = fullHistHolder(:,targetTrials);
            histErr(k,l,:) = std(specHist,0,2)/sqrt(length(targetTrials));
            [latPeakBinOut] = functionLatPeakBinCalculation(s.Parameters.ToneWindow,s.Parameters.GenWindow,s.Parameters.RasterWindow,...
                targetRasters,length(targetTrials),2,targetTrials,s.Parameters.latBin,s.Parameters.PercentCutoff,s.Parameters.BaselineCutoff);
            latePeakBinStore{k,l} = latPeakBinOut;
            latStore(k,l) = latPeakBinOut.ResponseLatency;
            peakStoreFast(k,l) = latPeakBinOut.PeakRespFast;
            peakStoreTone(k,l) = latPeakBinOut.PeakRespTone;
            peakStoreGen(k,l) = latPeakBinOut.PeakRespGen;
            binStoreFast(k,l) = mean(latPeakBinOut.BinnedSpikesFast);
            binStoreTone(k,l) = mean(latPeakBinOut.BinnedSpikesTone);
            binStoreGen(k,l) = mean(latPeakBinOut.BinnedSpikesGen);
            probStoreTone(k,l) = latPeakBinOut.ProbSpikeTone;
            probStoreGen(k,l) = latPeakBinOut.ProbSpikeGen;
            binSigVals(k,l,:) = latPeakBinOut.BinSigVals;
            binDiff(k,l,:) = latPeakBinOut.BinDiff;
            latPeakBinOut = [];
%             [responseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,toneReps);
            [responseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes(targetTrials),length(targetTrials),...
        s.Parameters.minSpikes,s.Parameters.latBin,[s.Parameters.RasterWindow(1),0],s.Parameters.zLimit,s.Parameters.minSigSpikes,s.Parameters.SigSmoothWindow);
            responseHistHolder{k,l} = responseHist;
            if responseHist.SigSpikePos == 1%if have a significant response, check for width
                [widthOut] = functionResponseWidth(responseHist);
                if length(widthOut.Widths) > 0
                    widthStore(k,l,i) = widthOut.Widths(1);
                    widthLat(k,l,i) = widthOut.StartsEnds(1);
                    bigWidth{k,l,i} = widthOut;
                end
            else
                widthOut = [];
            end
        end
        if numDBs == 1
            freqSpecHist(k,:) = (squeeze(organizedHist(k,:,:)));
        elseif numDBs > 1
            freqSpecHist(k,:) = mean(squeeze(organizedHist(k,:,:)));
        end
        
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
    s.(desigNames{i}).PeakMapGen = peakStoreGen;
    s.(desigNames{i}).PeakMapTone = peakStoreTone;
    s.(desigNames{i}).BinFast = binStoreFast;
    s.(desigNames{i}).BinTone = binStoreTone;
    s.(desigNames{i}).BinGen = binStoreGen;
    s.(desigNames{i}).BinDiff = binDiff;
    s.(desigNames{i}).ProbTone = probStoreTone;
    s.(desigNames{i}).ProbGen = probStoreGen;
    s.(desigNames{i}).BinSigVals = binSigVals;
    s.(desigNames{i}).WidthData = widthOut;
    
    %store some information about width
    CFTrig = 0;
    BFTrig = 0;
    for j = 1:numDBs
        sigThresh = 0.01;
        %pull out sign of change, and significance value
        if whiteStatus == 1
            targetSig = squeeze(binSigVals(2:end,j,:));
            targetSign = squeeze(binDiff(2:end,j,:));
            targetBins = squeeze(binStoreFast(2:end,j,:));
            
        else
            targetSig = squeeze(binSigVals(:,j,:));
            targetSign = squeeze(binDiff(:,j,:));
            targetBins = squeeze(binStoreFast(:,j,:));
        end

        for kk = 1:size(targetSig,2)
            %first, deal with positives
            findPos =  find(targetSig(:,kk) < sigThresh & targetSign(:,kk) > 0);
            if findPos%in the event of positive and significant events
                
                if CFTrig == 0
                    sigValLength = length(find(targetSig(:,kk) < sigThresh));
                    if whiteStatus == 1
                        if sigValLength == 1
                            masterData(i,masterHolder) = uniqueFreqs(find(targetSig(:,kk)<sigThresh)+1);
                        else
                            masterData(i,masterHolder) = mean(uniqueFreqs(find(targetSig(:,kk)<sigThresh)+1));
                        end
                    else
                        if sigValLength == 1
                            masterData(i,masterHolder) = uniqueFreqs(find(targetSig(:,kk)<sigThresh));
                        else
                            masterData(i,masterHolder) = mean(uniqueFreqs(find(targetSig(:,kk)<sigThresh)));
                        end
                    end
                    masterHeader{masterHolder} = 'CF';
                    masterHolder = masterHolder + 1;
                    disp('CF Detected')
                    CFTrig = 1;
                end
                %determine number of significant values
                if find(findPos == 1) | find(findPos == length(targetSig)) %determine if anything on the edge of the tuning range
                    s.PosEdgeWarn(j,i,kk) = 1;
                end
                if length(findPos) > 1
                    try
                        x = [1:length(targetSig)];
                        y = targetBins;
                        f = fit(x',y,'gauss1'); %perform single gaussian fit
                        newVect = [1:0.1:length(targetSig)];
                        fVals = f(newVect);
                        %find peak
                        [peakVal peakInd] = max(fVals);

                        %find half peak width
                        firstBound = find(fVals(1:peakInd) - peakVal/2 > 0,1,'first');
                        lastBound = find(fVals(peakInd:end) - peakVal/2 > 0,1,'last');
                        s.PosGaussWidth(j,i,kk) = (lastBound + peakInd - firstBound)/10;
                    catch
                        s.PosGaussWidth(j,i,kk) = 0;
                    end
                end
                s.PosTots(j,i,kk) = length(findPos);
                diffFind = diff(findPos);
                if diffFind == 1
                    disp('All Consecutive!')
                    s.PosCont(j,i,kk) = 1;
                    %since all consecutive, width equals length of findPos
                    s.PosWidths(j,i,kk) = length(findPos);
                else
                    disp('Not Consecutive, Finding Widest')
                    s.PosCont(j,i,kk) = 0;
%                     findSpaces = length(find(diffFind ~= 1)); %determine how many different segments there are
                    findGaps = find(diffFind ~= 1);
                    findGaps(2:end+1) = findGaps;
                    findGaps(1) = 0;
                    findGaps(end+1) = length(findPos);
                    diffLengths = diff(findGaps);
                    disp('Widest Point Found')
                    s.PosWidths(j,i,kk) = max(diffLengths);
                end
                %find maximum value. This is BF. 
                if j == numDBs & BFTrig == 0
                    [M bfFind] = max(findPos);
                    if whiteStatus == 1
                        masterData(i,masterHolder) = uniqueFreqs(findPos(bfFind)+1);
                    else
                        masterData(i,masterHolder) = uniqueFreqs(findPos(bfFind));
                    end
                    masterHeader{masterHolder} = 'BF';
                    masterHolder = masterHolder + 1;
                    disp('BF Detected')
                    BFTrig = 1;
                end
            else
                disp('No Positives Found')
            end
            
            %now do negative
            sigThresh = 0.01;
            findNeg =  find(targetSig(:,kk) < sigThresh & targetSign(:,kk) < 0);
            if findNeg%in the event of positive and significant events
                if find(findNeg == 1) | find(findNeg == length(targetSig))
                    s.NegEdgeWarn(j,i,kk) = 1;
                end
                s.NegTots(j,i,kk) = length(findNeg);
                diffFind = diff(findNeg);
                if diffFind == 1
                    disp('All Consecutive!')
                    s.NegCont(j,i,kk) = 1;
                    %since all consecutive, width equals length of findNeg
                    s.NegWidths(j,i,kk) = length(findNeg);
                else
                    disp('Not Consecutive, Finding Widest')
                    s.NegCont(j,i,kk) = 0;
%                     findSpaces = length(find(diffFind ~= 1)); %determine how many different segments there are
                    findGaps = find(diffFind ~= 1);
                    findGaps(2:end+1) = findGaps;
                    findGaps(1) = 0;
                    findGaps(end+1) = length(findNeg);
                    diffLengths = diff(findGaps);
                    disp('Widest Point Found')
                    s.NegWidths(j,i,kk) = max(diffLengths);
                end
            else
                disp('No Negatives Found')
            end
        end
        if j == numDBs & CFTrig == 0
            masterData(i,masterHolder) = NaN;
            masterHeader{masterHolder} = 'CF';
            masterHolder = masterHolder + 1;
            disp('NO CF')
        end
        
        if j == numDBs & BFTrig == 0
            masterData(i,masterHolder) = NaN;
            masterHeader{masterHolder} = 'BF';
            masterHolder = masterHolder + 1;
            disp('NO BF')
        end
    end
    
    if toggleROC == 1
        [velOut] = functionLocomotionROC(spikeTimes,s.RotaryData.Velocity);
        s.(desigNames{i}).TrueAUC = velOut.TrueAUC;
        s.(desigNames{i}).ShuffleAUC = velOut.ShuffleAUC;
        masterData(i,masterHolder) = velOut.TrueAUC;%find best frequency, ignores white noise
        masterHeader{masterHolder} = 'LocoAUC';
        masterHolder = masterHolder + 1;
    else
        s.(desigNames{i}).TrueAUC = 0;
        s.(desigNames{i}).ShuffleAUC = zeros(1000,1);
    end
    
    AUCLims(1) = prctile(s.(desigNames{i}).ShuffleAUC,0.05);
    AUCLims(2) = prctile(s.(desigNames{i}).ShuffleAUC,99.95);
    if s.(desigNames{i}).TrueAUC > AUCLims(2) | s.(desigNames{i}).TrueAUC < AUCLims(1)
        s.(desigNames{i}).AUCSig = 1;
    else
        s.(desigNames{i}).AUCSig = 0;
    end
    masterData(i,masterHolder) = s.(desigNames{i}).AUCSig;%find best frequency, ignores white noise
    masterHeader{masterHolder} = 'LocoAUCSig';
    masterHolder = masterHolder + 1;
    
end
s.WidthData = widthStore;
s.WidthLatData = widthLat;
s.FullWidth = bigWidth;

masterInd = masterHolder;

%% Laser Analysis

laserBinVect = [s.Parameters.LaserWindow(1):s.Parameters.LaserBin:s.Parameters.LaserWindow(2)-s.Parameters.LaserBin];
if idToggle == 1
    
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
        laserResp = zeros(length(alignTimes),1);
        for j = 1:length(alignTimes)
            %find if there are any spikes
            laseFind = find(laserRasters(:,2) == j & laserRasters(:,1) > 0 & laserRasters(:,1) < s.Parameters.LaserLim);
            if laseFind
                laserResp(j) = 1;
            end
        end
        s.(desigNames{i}).LaserRasters = laserRasters;
        s.(desigNames{i}).LaserHist = laserHistData;
        s.(desigNames{i}).LaserResps = laserResp;
        
    end
else
    for i = 1:numUnits
        s.(desigNames{i}).LaserRasters = [];
        s.(desigNames{i}).LaserHist = zeros(length(laserBinVect),1);
        s.(desigNames{i}).LaserResps = zeros(length(dio2Times),1);
    end
    
end

%calculate and plot LFP information
% [lfpStruct] = functionLFPaverage(master, s.Parameters.LFPWindow, s,homeFolder,fileName, uniqueFreqs, uniqueDBs, numFreqs, numDBs);
% s.LFP = lfpStruct;

%% Saving
save(fullfile(pname,fname),'s','masterData','masterHeader');

%% Plotting!!

%% first plot general figure. 
[indCellType] = functionCellStringFind(masterHeader,'CellType');


%determine if there are units in each category
findPVs = find(masterData(:,indCellType) == 1);
if findPVs
    for i = 1:length(findPVs)
        pvStores(:,i) = s.(s.DesignationName{findPVs(i)}).SessionFiring;
    end
    avPV = mean(pvStores');
end

findMSNs = find(masterData(:,indCellType) == 0);
for i = 1:length(findMSNs)
    msnStores(:,i) = s.(s.DesignationName{findMSNs(i)}).SessionFiring;
end
avMSN = mean(msnStores');

findCHATs = find(masterData(:,indCellType) == 2);
if findCHATs
    for i = 1:length(findCHATs)
        chatStores(:,i) = s.(s.DesignationName{findCHATs(i)}).SessionFiring;
    end
    avCHAT = mean(chatStores');
end

[indPkTr] = functionCellStringFind(masterHeader,'PeakTrough');
[indISI] = functionCellStringFind(masterHeader,'isiCov');
[indPosSig] = functionCellStringFind(masterHeader,'PosSigGenHist');
[indNegSig] = functionCellStringFind(masterHeader,'NegSigGenHist');




hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
% Column 1
%plot spike width vs coefficient of variation
subplot(3,4,1)
hold on

plot(masterData(:,indPkTr),masterData(:,indISI),'k.')
plot(masterData(masterData(:,indCellType) == 1,indPkTr),masterData(masterData(:,indCellType) == 1,indISI),'r.')
plot(masterData(masterData(:,indCellType) == 2,indPkTr),masterData(masterData(:,indCellType) == 2,indISI),'g.')
xlabel('Peak Trough (ms)')
ylabel('ISI Coefficient of Variation')
title(fileName,'fontweight','bold', 'Interpreter', 'none');

%plot loco speed vs firing rates
subplot(3,4,5)
hold on
plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'g')
plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],avMSN/max(avMSN),'k')
if findPVs
    plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],avPV/max(avPV),'r')
end
xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
ylim([-0.1,1])

%plot out proportion of each cell type

subplot(3,4,9)
cellDist = [length(findMSNs),length(findPVs),length(findCHATs)];
pie(cellDist)
labels = {'MSNs','PVs','ChATs'};
detZero = find(cellDist == 0);
labels(detZero) = [];
legend(labels,'Location','southoutside','Orientation','horizontal')


% Column 2

subplot(3,4,2)
hold on

hold on
holder = masterData(findMSNs,[indPosSig,indNegSig]);
holder(:,2) = holder(:,2) * -2;
posResp = find(holder(:,1) == 1);


plot(s.PosWidths(:,findMSNs(posResp),1),'k.')
plot(mean(s.PosWidths(:,findMSNs(posResp),1)'),'k-')
plot(s.PosWidths(:,findMSNs(posResp),2),'b.')
plot(mean(s.PosWidths(:,findMSNs(posResp),2)'),'b-')
plot(s.PosWidths(:,findMSNs(posResp),3),'m.')
plot(mean(s.PosWidths(:,findMSNs(posResp),3)'),'m-')
xlim([0 size(s.PosWidths,1) + 1])
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

% Column 3 PV CELLS
if findPVs
    subplot(3,4,4)
    hold on
    holder = masterData(findPVs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    posResp = find(holder(:,1) == 1);

    plot(s.PosWidths(:,findPVs(posResp),1),'k.')
    plot(mean(s.PosWidths(:,findPVs(posResp),1)'),'k-')
    plot(s.PosWidths(:,findPVs(posResp),2),'b.')
    plot(mean(s.PosWidths(:,findPVs(posResp),2)'),'b-')
    plot(s.PosWidths(:,findPVs(posResp),3),'m.')
    plot(mean(s.PosWidths(:,findPVs(posResp),3)'),'m-')
    xlim([0 size(s.PosWidths,1) + 1])
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
end

% Column 4 CHATs
if findCHATs
    subplot(3,4,3)
    hold on
    holder = masterData(findCHATs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    posResp = find(holder(:,1) == 1);

    plot(s.PosWidths(:,findCHATs(posResp),1),'k.')
    plot(mean(s.PosWidths(:,findCHATs(posResp),1)'),'k-')
    plot(s.PosWidths(:,findCHATs(posResp),2),'b.')
    plot(mean(s.PosWidths(:,findCHATs(posResp),2)'),'b-')
    plot(s.PosWidths(:,findCHATs(posResp),3),'m.')
    plot(mean(s.PosWidths(:,findCHATs(posResp),3)'),'m-')
    xlim([0 size(s.PosWidths,1) + 1])
    title(strcat(num2str(length(findCHATs)),'-CHAT Tuning Width Responses fast(k) tone(b) gen(m)'))

    %plot distribution of positive, negative, both, and untuned units
    subplot(3,4,7)
    holder = masterData(findCHATs,[indPosSig,indNegSig]);
    holder(:,2) = holder(:,2) * -2;
    det = holder(:,1) + holder(:,2);
    det = hist(det,[-2:1:1]);
    pie(det)
    labels = {'Neg','Mix','None','Pos'};
    detZero = find(det == 0);
    labels(detZero) = [];
    legend(labels,'Location','southoutside','Orientation','horizontal')
end

spikeGraphName = strcat(fileName,'GeneralFigureCellTypes');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

close

% %% Now plot figure of all cross correlograms. 
% %first, calculate all xcorrs. 
% crossWindow = [-0.01 0.01];
% xCorrVect = [crossWindow(1):0.001:crossWindow(2)];
% %divide by shank
% for m = 1:2
%     disp(strcat('xCorr-Shank-',num2str(m)))
%     targUnits = find(masterData(:,2) == m);
% %     targPVs = find(masterData(targUnits,7) == 1);
%     targNames = desigNames(targUnits);
%     numTargs = length(targUnits);
%     for i = 1:numTargs
%         disp(i)
%         for j = 1:numTargs
%             %compute xcorr
%             spikeTimes1 = s.(targNames{i}).SpikeTimes;
%             spikeTimes2 = s.(targNames{j}).SpikeTimes;
%             spikeStore = ones(100000,1);
%             sCounter = 1;
%             for spikeInd = 1:length(spikeTimes1)
%                 %subtract spike time from first train out of all of second train
%                 subSpikes = spikeTimes2 - spikeTimes1(spikeInd);
%                 %remove things outside the window of interest
%                 subSpikes(subSpikes>crossWindow(2) | subSpikes<crossWindow(1)) = [];
%                 spikeStore(sCounter:(sCounter + length(subSpikes)-1)) = subSpikes;
%                 sCounter = sCounter + length(subSpikes);
%             end
%             spikeStore(sCounter:end) = [];
%             %now we need to standardize so we can store large array
%             tempStore = [];
%             tempStore = hist(spikeStore,xCorrVect);
%             bigXCorrStore(m,i,j,:) = tempStore;
%         end
% 
%     end
% end
% 
% subplot = @(m,n,p) subtightplot (m, n, p, [0.01], [0.01 0.01], [0.03 0.01]);
% 
% % for m = 1:2
% %     targUnits = find(masterData(:,2) == m);
% %     targNames = desigNames(targUnits);
% %     numTargs = length(targUnits);
% %     hFig = figure;
% %     set(hFig, 'Position', [10 80 1900 1000])
% %     for i = 1:numTargs
% %         for j = 1:numTargs
% %             subplot(numTargs,numTargs,(i-1)*numTargs+j)
% %             bar(xCorrVect,squeeze(bigXCorrStore(m,i,j,:)))
% %             hold on
% %             plot([0 0],[0 max(squeeze(bigXCorrStore(m,i,j,:)))],'r')
% %             plot(crossWindow,[0 0],'b')
% %             axis off
% %         end
% %     end
% % end
% 
% 
% for m = 1:2
%     targUnits = find(masterData(:,2) == m);
%     targPVs = find(masterData(targUnits,7) == 1);
%     targNames = desigNames(targUnits);
%     numTargs = length(targUnits);
%     if targPVs
%         hFig = figure;
%         set(hFig, 'Position', [10 80 1600 800])
%         for i = 1:length(targPVs)
%             for j = 1:numTargs
%                 subplot(length(targPVs),numTargs,(i-1)*numTargs+j)
%                 bar(xCorrVect,squeeze(bigXCorrStore(m,targPVs(i),j,:)))
%                 hold on
%                 plot([0 0],[0 max(squeeze(bigXCorrStore(m,targPVs(i),j,:)))],'r')
%                 plot(crossWindow,[0 0],'b')
%                 axis off
%                 xlim(crossWindow)
%             end
%         end
%         
%         spikeGraphName = strcat('CrossCorrPVvsOthersShank',num2str(m));
%         savefig(hFig,spikeGraphName);
% 
%         %save as PDF with correct name
%         set(hFig,'Units','Inches');
%         pos = get(hFig,'Position');
%         set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%         print(hFig,spikeGraphName,'-dpdf','-r0')
%     end
% end

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
findCHATshank1 = intersect(findCHATs,findShank1);

subplot(2,5,1)
imagesc((whiteStore(:,s.SortedPeakWaveOrder(findShank1))'),[-1 1])
hold on
for i = 1:length(findPVshank1)
    plot([0 length(whiteStore)],[s.SortedPeakWaveOrder(findPVshank1(i)) s.SortedPeakWaveOrder(findPVshank1(i))],'r','LineWidth',2)
end
for i = 1:length(findCHATshank1)
    plot([0 length(whiteStore)],[s.SortedPeakWaveOrder(findCHATshank1(i)) s.SortedPeakWaveOrder(findCHATshank1(i))],'g','LineWidth',2)
end
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
findCHATshank2 = intersect(findCHATs,findShank2);

subplot(2,5,6)
imagesc((whiteStore(:,s.SortedPeakWaveOrder(findShank2))'),[-1 1])
hold on
for i = 1:length(findPVshank2)
    plot([0 length(whiteStore)],[s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1) s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1)],'r','LineWidth',2)
end
for i = 1:length(findCHATshank2)
    plot([0 length(whiteStore)],[s.SortedPeakWaveOrder(findCHATshank2(i))-length(findShank1) s.SortedPeakWaveOrder(findCHATshank2(i))-length(findShank1)],'g','LineWidth',2)
end
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

close


%% now plot individual cells. 
for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])

    %% Column 1
    %plots average waveform
    subplot(4,8,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    if ismember(i,findPVs)
        title(strcat('PV AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    elseif ismember(i,findMSNs)
        title(strcat('MSN AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    elseif ismember(i,findCHATs)
        title(strcat('CHAT AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
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

    %% Column 2
    %Plot binned response during tone period
    subplot(4,4,2)
    imagesc(s.(desigNames{i}).BinTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title({fileName;desigNames{i}},'fontweight','bold', 'Interpreter', 'none');

    %Plot binned response during general period
    subplot(4,4,6)
    imagesc(s.(desigNames{i}).BinGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Response (general)')
    %Plot peak response during tone period
    subplot(4,4,10)
    imagesc(s.(desigNames{i}).PeakMapTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Peak Response (tone)')
    %Plot peak response during general period
    subplot(4,4,14)
    imagesc(s.(desigNames{i}).PeakMapGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Peak Response (general)')

    %% Column 3
%         %plot probability of response (tone)
%         subplot(4,4,3)
%         imagesc(s.(desigNames{i}).ProbTone')
%         colormap(parula)
%         colorbar
%         set(gca,'XTick',octaveRange(:,2));
%         set(gca,'XTickLabel',octaveRange(:,1));
%         set(gca,'YTick',dbRange(:,2));
%         set(gca,'YTickLabel',dbRange(:,1));
%         title('Probability of Response (tone)')
%         %plot probability of response (gen)
%         subplot(4,4,7)
%         imagesc(s.(desigNames{i}).ProbGen')
%         colormap(parula)
%         colorbar
%         set(gca,'XTick',octaveRange(:,2));
%         set(gca,'XTickLabel',octaveRange(:,1));
%         set(gca,'YTick',dbRange(:,2));
%         set(gca,'YTickLabel',dbRange(:,1));
%         title('Probability of Response (general)')

    %plot out binned spikes (tone)
    subplot(4,4,7)
    hold on
    for cInd = 1:numDBs
        plot(s.(desigNames{i}).BinDiff(:,cInd,1),'Color',[cInd/numDBs 0 0])
        %find significant points, by p < 0.05
        findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,1)<0.05);
        plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,1),'g*')
        %find significant points, by p < 0.01
        findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,1)<0.01);
        plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,1),'go')
    end
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    ylabel('Binned Spikes/Fast Period')
    title('Tuning Curves Across Fast Period')

    %plot out binned spikes (tone)

    subplot(4,4,11)
    hold on
    for cInd = 1:numDBs
        plot(s.(desigNames{i}).BinDiff(:,cInd,2),'Color',[cInd/numDBs 0 0])
        %find significant points, by p < 0.05
        findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,2)<0.05);
        plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,2),'g*')
        %find significant points, by p < 0.01
        findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,2)<0.01);
        plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,2),'go')
    end
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    ylabel('Binned Spikes/Tone Period')
    title('Tuning Curves Across Tone Period')
    %plot out binned spikes (general)
    subplot(4,4,15)
    hold on
    for cInd = 1:numDBs
        plot(s.(desigNames{i}).BinDiff(:,cInd,3),'Color',[cInd/numDBs 0 0])
        %find significant points, by p < 0.05
        findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,3)<0.05);
        plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,3),'g*')
        %find significant points, by p < 0.01
        findSigs = find(s.(desigNames{i}).BinSigVals(:,cInd,3)<0.01);
        plot(findSigs,s.(desigNames{i}).BinDiff(findSigs,cInd,3),'go')
    end
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    ylabel('Binned Spikes/Gen Period')
    title('Tuning Curves Across General Period')

    %% Column 4
    %plot velocity data
    subplot(4,4,4)
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


    %plot heatmap organized by frequency
    subplot(4,4,8)
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

    if idToggle == 0

    elseif idToggle == 1
        %plot laser raster
        subplot(4,4,12)
        plot(s.(desigNames{i}).LaserRasters(:,1),s.(desigNames{i}).LaserRasters(:,2),'k.')
        xlim([s.Parameters.LaserWindow(1),s.Parameters.LaserWindow(2)])
        title(strcat('Laser Raster, Response %:',num2str(sum(s.(desigNames{i}).LaserResps)/length(dio2Times))))
        subplot(4,4,16)
        plot(laserBinVect,s.(desigNames{i}).LaserHist)
        xlim([s.Parameters.LaserWindow(1),s.Parameters.LaserWindow(2)])
        title('Laser Histogram')
    end



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