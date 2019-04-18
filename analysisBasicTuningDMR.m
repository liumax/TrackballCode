
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

%for dmr!
s.Parameters.STRFWin = [-0.2 0.01];
s.Parameters.dmrTiming = 6;
s.Parameters.expectedDur = 10; %expected time duration. 
numLags = 100;

%we have multiple datasets, so lets figure out which one is which. 
dateVal = str2num(fileName(1:6));
if dateVal < 190211
%     sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR10mindsT5dsF5.mat';
%     sprFile =
%     'Z:\Max\dmrOutputFiles\extractedDMR400-6400-approx1kHz.mat'; %This is
%     the 1 kHz example. 
    sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR400-6400-approx1kHz.mat';
    % dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
    ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimesCORRECTFreqSampling.mat';

%     ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimes.mat';
else 
    sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR4kHz-16kHz_10mindsT30dsF5.mat';
    % dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
    ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-4000flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimes.mat';
end
s.Parameters.MinFreq = 4000;
s.Parameters.MaxFreq = 64000;
% s.Parameters.DMRtimeDilation = 0.9995; %dilation of time. Multiply Trodes time by this to get DMR time.  
% s.Parameters.DMRtimeBin = 1/6000;
% newWin = round(s.Parameters.STRFWin/s.Parameters.DMRtimeBin);
preDelay = 0.05;
postDelay = 1/6;

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
        disp('Duplicates Found')
        s.DeletedUnits = n.s.DeletedUnits;
        %now delete bad units
        delUnits = fields(s.DeletedUnits);
        for indCount = 1:length(delUnits)
            s = rmfield(s,delUnits{indCount});
        end
        %now we need to convert names into numbers for finding?
        for i = 1:length(s.DesignationName)
            findt = strfind(s.DesignationName{i},'t');
            findt = findt(1);
            findc = strfind(s.DesignationName{i},'c');
            findr = strfind(s.DesignationName{i},'r');
            numVals(i) = str2num([s.DesignationName{i}(findt+1:findc-1),'9',s.DesignationName{i}(findr+1:end)]);
        end
        
        for i = 1:length(n.s.DesignationName)
            findt = strfind(n.s.DesignationName{i},'t');
            findt = findt(1);
            findc = strfind(n.s.DesignationName{i},'c');
            findr = strfind(n.s.DesignationName{i},'r');
            numValsPrev(i) = str2num([n.s.DesignationName{i}(findt+1:findc-1),'9',n.s.DesignationName{i}(findr+1:end)]);
        end
        %find intersect
        [C,ia,ib] = intersect(numVals,numValsPrev);
        %ia is the set of names that I want. 
        
        
        %replace designation array and names. 
        s.DesignationArray = s.DesignationArray(ia,:);
        s.DesignationName = s.DesignationName(ia);
    else
        disp('No Duplicates Found')
        %replace designation array and names. 
        s.DesignationArray = s.DesignationArray;
        s.DesignationName = s.DesignationName;
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

%now find long pause!
findPause = find(dioTimeDiff > 10);
if length(findPause) > 1
    error('MULTIPLE LONG PAUSES')
else
    disp('Found Long Pause in DIO Times')
end

%now lets divide the data by the long pause
tuningDIO = dioTimes(1:findPause);
tuningDIODiff = diff(tuningDIO);
dmrDIO = dioTimes(findPause+1:end);

%% Clean up DIO for tuning. 
disp('DIO Data Extracted, Checking for Errors')
%insert to master. check for errors
if length(tuningDIO) == length(soundData.Frequencies)
    %Sets up master array so that we have a size comparison for DIO
    %information. 
    master = zeros(size(soundFile.soundData.Frequencies,1),5);

    master(:,1) = tuningDIO;
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
elseif length(tuningDIO) ~= length(s.SoundData.Frequencies) %error case
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
            
        elseif whileCounter > length(tuningDIODiff) %if hit bound of received times
            disp('Went through all received times')
            totalTrialNum = whileCounter;
            break
            
        else
            %check the difference between the current values. If exceeds
            %threshold, flag. Otherwise, need to record the tuning
            %information so I get accurate representation of repetitions
            diffCheck = tuningDIODiff(whileCounter) - predTimes(whileCounter);
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
                tuningDIO(whileCounter + 1) = [];
                tuningDIODiff = diff(tuningDIO);
            else
                error('Code Failure In DIO Repair')
            end
        end
    end
    master = zeros(totalTrialNum,5);

    master(:,1) = tuningDIO;
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

%% now clean up DIO for DMR
%now lets go through more carefully in the DMR block and look for failures
%load the data file.
load(ttlInfo)
lengthTrue = length(ttlOnsetTime);
if length(dmrDIO) == lengthTrue
    disp('DMR has correct number of TTLs')
else
    error('DMR has incorrect number of TTLs')
end
%now go through an examine for deviations of timing. 
devCounter = 1; %this counts number of deviating TTLs
devStore = [];
%now pull up appropriate ITIs
blockITIs = diff(dmrDIO);
%compare against true values with subtraction
itiMod = ttlDiffTime - blockITIs; %ttlDiffTime is in the ttl file. 
%find all remainders that are too big
bigMod = find(abs(itiMod) > 0.001);
if bigMod
    disp(strcat(num2str(length(bigMod)),'-deviations found, marking...'))
    for j = 1:length(bigMod)
        devStore(devCounter,1) = j; %store the chunk
        devStore(devCounter,2) = bigMod(j) - 1; %stores the TTLs before and after. 
        devStore(devCounter,3) = bigMod(j) +1;
%         devStore(devCounter,4) = ttlChunks(i,1) - 1 + bigMod(j) - 1; %stores the TTLs before and after. 
%         devStore(devCounter,5) = ttlChunks(i,1) - 1 + bigMod(j) + 1;

        devCounter = devCounter + 1;
    end
else
    disp('No Deviations Found! :D')
end

disp('Applying TTL QC to Spiking Data by Removing Spikes')
for i = 1:numUnits
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    for j = 1:size(devStore,1)
        findSpikes = find(spikeTimes >= dioTimes(devStore(j,2)) & spikeTimes <= dioTimes(devStore(j,3)));
        disp(strcat(num2str(length(findSpikes)),'-spikes removed'))
        spikeTimes(findSpikes) = [];
        s.(desigNames{i}).SpikeTimes = spikeTimes;
    end
end

s.SoundData.DMRPulses = dmrDIO;

%% Generate master
%form master matrix
trialNum = length(tuningDIO);
master = zeros(trialNum,5);

master(:,1) = tuningDIO;
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

%% Extract data from rotary encoder.


[D3FileName] = functionFileFinder(subFoldersCell,'DIO','D3');
if length(D3FileName) == 0
    [D3FileName] = functionFileFinder(subFoldersCell,'DIO','Din3');
end
D3FileName = D3FileName{1};

[D4FileName] = functionFileFinder(subFoldersCell,'DIO','D4');
if length(D4FileName) == 0
    [D4FileName] = functionFileFinder(subFoldersCell,'DIO','Din4');
end
D4FileName = D4FileName{1};

[funcOut] = functionNewRotaryExtraction(D3FileName,D4FileName);
s.RotaryData = funcOut;
funcOut = [];
disp('Rotary Encoder Data Extracted')
%181217 data output is now in trodes samples, at 1ms intervals. Need to fix this! 
% newTimes = [(round(timeMin*(1/interpStep)))*interpStep:interpStep:(round(timeMax*(1/interpStep)))*interpStep];
newTimeVector = [s.RotaryData.Distance(1,1)/s.Parameters.trodesFS:s.Parameters.InterpolationStepRotary:s.RotaryData.Distance(end,1)/s.Parameters.trodesFS];
newDistVector = interp1(s.RotaryData.Distance(:,1),s.RotaryData.Distance(:,2),newTimeVector);
newVelVector = interp1(s.RotaryData.Distance(1:end-1,1)/s.Parameters.trodesFS,s.RotaryData.Velocity,newTimeVector);

%rasterize this data
jumpsBack = round(s.Parameters.RasterWindow(1)/s.Parameters.InterpolationStepRotary);
jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);
velRaster = zeros(jumpsForward-jumpsBack+1,totalTrialNum);
length(s.RotaryData.Velocity);
for i = 1:length(tuningDIO)
    %find the time from the velocity trace closest to the actual stim time
    targetInd = find(newTimeVector - tuningDIO(i) > 0,1,'first');
    %pull appropriate velocity data
    if targetInd < length(newTimeVector) - jumpsForward
        velRaster(:,i) = newVelVector([targetInd+jumpsBack:targetInd+jumpsForward]);
    else
        velRaster(:,i) = zeros(jumpsForward-jumpsBack+1,1);
    end
end
%make average trace:
averageVel = mean(velRaster,2);
velVector = [s.Parameters.RasterWindow(1):s.Parameters.InterpolationStepRotary:s.Parameters.RasterWindow(2)];
velZero = find(velVector >= 0,1,'first');

%Change toggle for ROC analysis if insufficient running is found
if s.RotaryData.Distance(end,2) < 10;
    toggleROC = 0;
    disp('Resetting ROC Toggle Due to Lack of Movement')
end


figure
subplot(2,1,1)
plot(newTimeVector,newVelVector)
xlim([newTimeVector(1),newTimeVector(end)])
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


%% Now process DMR

load(sprFile)

%downsample to 1 ms. 
if dateVal < 190211
else
%     stimulus = stimulus(:,1:6:end);
end

% Figure out how DMR file maps onto trodes times
tempVect = linspace(dmrTimes(1),dmrTimes(2),length(stimulus));

trueDMRtimes = interp1(ttlOnsetTime,dmrDIO,tempVect);
dmrStore = [];
dmrStore = trueDMRtimes;
trueDMRtimes = [];

dmrStep = mean(mode(diff(dmrStore)));
newWin = round(s.Parameters.STRFWin/dmrStep);

spikeHistVect = dmrStore+dmrStep/2;
for i = 1:numUnits
    spikeStore = s.(desigNames{i}).SpikeTimes(s.(desigNames{i}).SpikeTimes > spikeHistVect(1) + s.Parameters.STRFWin(1) & s.(desigNames{i}).SpikeTimes < spikeHistVect(end));
    spikeArray(i,:) = hist(spikeStore,spikeHistVect);
end

%store number of spikes used for DMR.
masterData(:,masterInd) = sum(spikeArray'); masterHeader{masterInd} = 'SpikesPerDMR'; masterInd = masterInd + 1;

[sta, stabigmat, spkcountvec] = quick_calc_sta(stimulus, spikeArray, numLags);
[sta_sig, ptd, siglevel] = ne_sig_sta_from_stim_obs_resp(sta, spikeArray, stimulus, 10, numLags, 95);

s.STAs = sta;
s.STASig = sta_sig;
%lets now generate some STA based data by reshaping and pulling time and
%frequency information. We will also extract the RTF.

%first, reshape the STAs!
for i = 1:numUnits
    tmpStore = reshape(sta_sig(i,:),numLags,[]);
    newSTA(:,:,i) = tmpStore;
    newTMP = tmpStore;
    newTMP(tmpStore >=0) = 0;
    negSTA(:,:,i) = newTMP;
    tmpStore(tmpStore <= 0) = 0;
    posSTA(:,:,i) = tmpStore;
end

%to determine peak and width, lets compress along time axis or freq axis.. 
smoothWind = 1;
for i = 1:numUnits
    compTimingPos(:,i) = sum(posSTA(:,:,i));
    compTimingNeg(:,i) =  sum(negSTA(:,:,i));
    compWidthPos(:,i) = smooth(sum(posSTA(:,:,i)'),smoothWind);
    compWidthNeg(:,i) = smooth(sum(negSTA(:,:,i)'),smoothWind);
end

%now lets test significance of the responses. Since this requires many
%loops to get things done, I think the best idea is to actually downsample
%down to a 5ms time step, so that I can still get through things, show
%significance 

%the way in which I will do that is by smoothing by 5, then downsampling by
%5. 
dsLags = 20;
dsStimulus = [];
for i = 1:length(faxis)
    dsStimulus(i,:) = downsample(smooth(stimulus(i,:),5),5)*5;
end


for i = 1:numUnits
    dsSpikes(i,:) = downsample(smooth(spikeArray(i,:),5),5)*5;
end


tempVect = linspace(dmrTimes(1),dmrTimes(2),length(dsStimulus));

dsDMRtimes = interp1(ttlOnsetTime,dmrDIO,tempVect);

dsDmrStep = mean(mode(diff(dsDMRtimes)));
% newWin = round(s.Parameters.STRFWin/dsDmrStep);

dsSpikeHistVect = dsDMRtimes+dsDmrStep/2;

%Now lets generate split STAs to calculate correlation. 
randVals = randperm(length(dsSpikes),floor(length(dsSpikes)/2));
splitSpikes1 = zeros(size(dsSpikes));
splitSpikes1(:,randVals) = dsSpikes(:,randVals);
randVals2 = [1:1:length(dsSpikes)];
randVals2(randVals) = [];
splitSpikes2 = zeros(size(dsSpikes));
splitSpikes2(:,randVals2) = dsSpikes(:,randVals2);
disp('Calculating Correlation Coefficient of Split STA')
[sta1, stabigmat, spkcountvec] = quick_calc_sta(dsStimulus, splitSpikes1, dsLags);
[sta2, stabigmat, spkcountvec] = quick_calc_sta(dsStimulus, splitSpikes2, dsLags);

%now lets generate correlation coefficients. 
for i = 1:numUnits
    tester = corrcoef(sta1(i,:),sta2(i,:));
    realCorrStore(i) = tester(2);
end

disp('Extracting Nonlinearities in Response')
%now lets extract nonlinearities
nonLinearBin = 15;
for i = 1:numUnits
    [xprior, xposterior] = ne_sta_stimulus_projection(reshape(s.STAs(i,:),length(faxis),[]), spikeArray(i,:),stimulus);
    [px{i},pspk{i},pxspk{i},xbincenter] = calc_px_pspk_pxspk(xprior,xposterior, nonLinearBin);
end
disp('Finished extracting nonlinearities...')

s.RealCorrStore = realCorrStore;

s.DMRTimes = trueDMRtimes;
s.DMRStep = dmrStep;
s.DMRfaxis = faxis;

s.SpikeArray = spikeArray;
s.NonLinear.px = px;
s.NonLinear.pspk = pspk;
s.NonLinear.pxspk = pxspk;
s.NonLinear.xbins = xbincenter;

%% Pull waveforms from different time periods, so we can plot them out.

for i = 1:numUnits
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    findTuningSpikes = [1:find(spikeTimes < dmrStore(1),1,'last')];
    findDMRSpikes = [find(spikeTimes < dmrStore(1),1,'last'):find(spikeTimes < dmrStore(end) - s.Parameters.STRFWin(2) ,1,'last')];
    numTuningSpikes(i) = length(findTuningSpikes);
    numDMRSpikes(i) = length(findDMRSpikes);
    %now lets pull average waves. 
    waveForms = s.(desigNames{i}).AllWaves;
    tuningWaves = waveForms(:,:,findTuningSpikes);
    avWaveTuning = reshape(mean(tuningWaves,3),[],1);
    dmrWaves = waveForms(:,:,findDMRSpikes);
    avWaveDMR = reshape(mean(dmrWaves,3),[],1);
    s.(desigNames{i}).MeanTuningWave = avWaveTuning;
    s.(desigNames{i}).MeanDMRWave = avWaveDMR;
end

%% Saving
save(fullfile(pname,fname),'s','masterData','masterHeader');
% save(fullfile(pname,strcat(fileName,'DMRData')),'sta','sta_sig','spikeArray','stimulus','faxis')
save(fullfile(pname,strcat(fileName,'DMRData')),'sta','spikeArray','stimulus','faxis','spikeHistVect')

%% Plotting shanks!
hFig = figure;
[indCellType] = functionCellStringFind(masterHeader,'CellType');
indCellType = indCellType(1);
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
set(hFig, 'Position', [10 80 1900 1000])

%plot response to all tones
respStore = [];
zeroPoint = find(histBinVector < 0,1,'last');
%z score these responses so that things will display a bit more nicely.
tarWindows = [0.05,s.Parameters.ToneWindow(2),s.Parameters.GenWindow(2)];
zlims = [-20 20];
targetX = [0:40:length(histBinVector)];
targetX(1) = 1;
for i = 1:numUnits
    respStore(:,i) = (s.(desigNames{i}).AllHistograms-mean(s.(desigNames{i}).AllHistograms(1:zeroPoint)))/std(s.(desigNames{i}).AllHistograms(1:zeroPoint));
    respStore(:,i) = respStore(:,i)/max(respStore(:,i));
    %pull overall maps
    tarBinDiff = squeeze(mean(s.(desigNames{i}).BinDiff,2));
    for j = 1:size(s.(desigNames{i}).BinDiff,3)
        %fast window is hardcoded for 50 ms. Others are target defined. 
        fullDiffs(:,j,i) = (tarBinDiff(:,j)/tarWindows(j) - s.(desigNames{i}).AverageRate)/s.(desigNames{i}).AverageSTD;
    end
%     fullDiffs(:,:,i) = bsxfun(@rdivide, squeeze(mean(s.(desigNames{i}).BinDiff,2)), squeeze(mean(s.(desigNames{i}).BinDiff,2)));
%     fullDiffs(:,:,i) = squeeze(mean(s.(desigNames{i}).BinDiff,2))./max(abs(squeeze(mean(s.(desigNames{i}).BinDiff,2))))';
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

subplot(2,4,1)
imagesc((respStore(:,s.SortedPeakWaveOrder(findShank1))'),[-1 1])
colormap(map)
colorbar
hold on
for i = 1:length(findPVshank1)
    plot([0 length(respStore)],[s.SortedPeakWaveOrder(findPVshank1(i)) s.SortedPeakWaveOrder(findPVshank1(i))],'g','LineWidth',2)
end
ylim([0.5 length(findShank1)+0.5])
xlim([0.5 length(respStore)+0.5])
set(gca, 'Ydir', 'reverse')
set(gca,'XTick',targetX);
set(gca,'XTickLabel',histBinVector(targetX));
title('Shank 1 PSTH, ZScore')

subplot(2,4,2)
hold on
imagesc(squeeze(fullDiffs(:,1,s.SortedPeakWaveOrder(findShank1)))',zlims)
for i = 1:length(findPVshank1)
    plot([0 size(fullDiffs,1)],[s.SortedPeakWaveOrder(findPVshank1(i)) s.SortedPeakWaveOrder(findPVshank1(i))],'g','LineWidth',2)
end
% colormap parula
colorbar
ylim([0.5 length(findShank1)+0.5])
xlim([0.5 size(fullDiffs,1)+0.5])
set(gca, 'Ydir', 'reverse')
title(strcat(fileName,'BinSubFast (Z)'), 'Interpreter', 'none')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

subplot(2,4,3)
hold on
imagesc(squeeze(fullDiffs(:,2,s.SortedPeakWaveOrder(findShank1)))',zlims)
for i = 1:length(findPVshank1)
    plot([0 size(fullDiffs,1)],[s.SortedPeakWaveOrder(findPVshank1(i)) s.SortedPeakWaveOrder(findPVshank1(i))],'g','LineWidth',2)
end
% colormap parula
colorbar
ylim([0.5 length(findShank1)+0.5])
xlim([0.5 size(fullDiffs,1)+0.5])
title('BinSubTone (Z)')
set(gca, 'Ydir', 'reverse')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

subplot(2,4,4)
hold on
imagesc(squeeze(fullDiffs(:,3,s.SortedPeakWaveOrder(findShank1)))',zlims)
for i = 1:length(findPVshank1)
    plot([0 size(fullDiffs,1)],[s.SortedPeakWaveOrder(findPVshank1(i)) s.SortedPeakWaveOrder(findPVshank1(i))],'g','LineWidth',2)
end
% colormap parula
colorbar
title('BinSubGeneral (Z)')
ylim([0.5 length(findShank1)+0.5])
xlim([0.5 size(fullDiffs,1)+0.5])
set(gca, 'Ydir', 'reverse')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));


findShank2 = find(masterData(:,indShankDes) == 2);

%find PV cells on shank 2
findPVshank2 = intersect(findPVs,findShank2);
% findCHATshank2 = intersect(findCHATs,findShank2);

subplot(2,4,5)
hold on
imagesc((respStore(:,s.SortedPeakWaveOrder(findShank2))'),[-1 1])
for i = 1:length(findPVshank2)
    plot([0 length(respStore)],[s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1) s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1)],'g','LineWidth',2)
end
% colormap(map)
colorbar
ylim([0.5 length(findShank2)+0.5])
xlim([0.5 length(respStore)+0.5])
set(gca,'XTick',targetX);
set(gca,'XTickLabel',histBinVector(targetX));
set(gca, 'Ydir', 'reverse')
title('Shank 2 PSTH, ZScore')

subplot(2,4,6)
hold on
imagesc(squeeze(fullDiffs(:,1,s.SortedPeakWaveOrder(findShank2)))',zlims)
for i = 1:length(findPVshank2)
    plot([0 size(fullDiffs,1)],[s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1) s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1)],'g','LineWidth',2)
end
% colormap parula
colorbar
title('BinSubFast (Z)')
ylim([0.5 length(findShank2)+0.5])
xlim([0.5 size(fullDiffs,1)+0.5])
set(gca, 'Ydir', 'reverse')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

subplot(2,4,7)
hold on
imagesc(squeeze(fullDiffs(:,2,s.SortedPeakWaveOrder(findShank2)))',zlims)
for i = 1:length(findPVshank2)
    plot([0 size(fullDiffs,1)],[s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1) s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1)],'g','LineWidth',2)
end
% colormap parula
colorbar
title('BinSubTone (Z)')
ylim([0.5 length(findShank2)+0.5])
xlim([0.5 size(fullDiffs,1)+0.5])
set(gca, 'Ydir', 'reverse')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

subplot(2,4,8)
hold on
imagesc(squeeze(fullDiffs(:,3,s.SortedPeakWaveOrder(findShank2)))',zlims)
for i = 1:length(findPVshank2)
    plot([0 size(fullDiffs,1)],[s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1) s.SortedPeakWaveOrder(findPVshank2(i))-length(findShank1)],'g','LineWidth',2)
end
% colormap parula
colorbar
title('BinSubGeneral (Z)')
ylim([0.5 length(findShank2)+0.5])
xlim([0.5 size(fullDiffs,1)+0.5])
set(gca, 'Ydir', 'reverse')
set(gca,'XTick',octaveRange(:,2));
set(gca,'XTickLabel',octaveRange(:,1));

spikeGraphName = strcat(fileName,'ProbeTopology');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%% Plotting individual units

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])
    % Column 1
    %plots average waveform
    subplot(4,3,1)
    hold on
    plot(s.(desigNames{i}).MeanTuningWave,'k','LineWidth',2)
    plot(s.(desigNames{i}).MeanDMRWave,'r','LineWidth',2)
    if ismember(i,findPVs)
        title(strcat('PV AvRate:',num2str(s.(desigNames{i}).AverageRate),'Tuning',num2str(numTuningSpikes(i)),'DMR',num2str(numDMRSpikes(i))))
    elseif ismember(i,findMSNs)
        title(strcat('MSN AvRate:',num2str(s.(desigNames{i}).AverageRate),'Tuning',num2str(numTuningSpikes(i)),'DMR',num2str(numDMRSpikes(i))))
%     elseif ismember(i,findCHATs)
%         title(strcat('CHAT AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    else
        title(strcat('UNK AvRate:',num2str(s.(desigNames{i}).AverageRate),'Tuning',num2str(numTuningSpikes(i)),'DMR',num2str(numDMRSpikes(i))))
    end

    
    % plot histogram.
    subplot(4,3,4)
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
    histLims = max(s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation);
    ylim([0 histLims])
    title('Histogram')

    %plot out rasters, organized!
    subplot(2,3,4)
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
    %plots ISI
    subplot(4,3,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})

    %plot FR and velocity
    subplot(4,3,5)
    hold on
    plot(newTimeVector,newVelVector/max(newVelVector),'b')
    plot([newTimeVector(1):s.Parameters.SpeedFiringBins:newTimeVector(end)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
    xlim([newTimeVector(1),newTimeVector(end)])
    ylim([-0.1,1])
    title({fileName;desigNames{i}},'fontweight','bold', 'Interpreter', 'none');
    
    %plot heatmap tuning curves
    subplot(4,3,8)
    imagesc(s.(desigNames{i}).BinTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Resp (tone)')
    
    subplot(4,3,11)
    imagesc(s.(desigNames{i}).BinGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Resp (general)')
    
    
    % Column 3
    
    %plot DMR response
    subplot(3,3,3)
    imagesc(reshape(sta(i,:),length(faxis),[]))
    colormap('parula')
    colorbar
    format short
    set(gca,'XTick',[0:20:numLags]);
    set(gca,'XTickLabel',[dmrStep*100:-20*dmrStep:dmrStep]);
    set(gca,'YTick',[1:10:length(faxis)]);
    set(gca,'YTickLabel',[faxis([1:10:end])]);
    title(strcat('DMR STA Based on',num2str(sum(spikeArray(i,:)))))
    
    %plot out thresholded DMR
    subplot(3,3,6)
    imagesc(reshape(sta_sig(i,:),length(faxis),[]))
    colormap('parula')
    colorbar
    format short
    set(gca,'XTick',[0:20:numLags]);
    set(gca,'XTickLabel',[dmrStep*100:-20*dmrStep:dmrStep]);
    set(gca,'YTick',[1:10:length(faxis)]);
    set(gca,'YTickLabel',[faxis([1:10:end])]);
    title(strcat('Thresholded DMR a=0.95 SelfCorr = ',num2str(s.RealCorrStore(i))))
    
    %plot nonlinearity
    subplot(3,3,9)
    hold on
    nl = pspk{i} .* pxspk{i} ./ px{i};

    maxmax = max(nl);
    plot(xbincenter, nl, 'ko-', 'markerfacecolor', 'k');
    minx = min(xbincenter);
    maxx = max(xbincenter);
    plot([minx-1 maxx+1], [pspk{i} pspk{i}], 'k--');
    xlim([minx-2 maxx+2]);
    if ~isnan(maxmax)
        ylim([0 maxmax]);
    else
        ylim([0 1])
    end
    xlabel('Projection (SD)');
    ylabel('P(spk|x)');
%     freqs
    
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


