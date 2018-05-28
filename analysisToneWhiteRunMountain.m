





function [s] = analysisToneWhiteRunMountain(fileName);

%set other things
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%TOGGLES FOR ENABLING/DISABLING FEATURES
s.Parameters.toggleRPV = 0; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 0; %1 means you want to eliminate duplicates.
toggleROC = 1; %toggle for tuning on/off ROC analysis

% %Toggles for different aspects of analysis
% toggleAuditory = 0;
% toggleLight = 0;

%basic parameters
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.trodesFS = 30000;%trodes sampling rate

%for duplicate elimination
s.Parameters.DownSampFactor = 10; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr

%for rotary encoder:
s.Parameters.InterpolationStepRotary = 0.01; %interpolation steps in seconds.

%for edr
% s.Parameters.EDRdownsamp = 20; %number of samples to downsample by. Smoothing is likely unnecessary
s.Parameters.EDRTimeCol = 1;
s.Parameters.EDRTTLCol = [3:5];
s.Parameters.EDRPiezoCol = 2;

%for cell type
s.Parameters.PVLim = [0.0004 0.0005];
s.Parameters.ChatLim = 1.1;

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 


%Parameters for auditory stimuli, visual stimuli, locomotion
s.Parameters.RasterWindow = [-0.5 1];
s.Parameters.histBin = 0.005;

%parameters specifically for auditory
calcWindow = [0 0.2]; %defines period for looking for responses, based on toneDur
s.Parameters.zLimit = [0.05 0.01 0.001];
s.Parameters.minSpikes = 100; %minimum number of spikes in baseline, if lower, triggers a warning
s.Parameters.minSigSpikes = 5; %minimum number of significant points to record a significant response.
s.Parameters.PercentCutoff = 99.9; %for significance in latency calculations
s.Parameters.BaselineCutoff = 95; %for the onset in latency calculations
s.Parameters.latBin = 0.001; %histogram bins for latency and significance calculations
s.Parameters.SigSmoothWindow = 11; %window of smoothing for calculations of significance
rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
s.Parameters.ToneWindow = [0 0.1];
s.Parameters.GenWindow = [0 0.2];

%% sets up file saving stuff
saveName = strcat(fileName,'LightToneAnalysis','.mat');
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


%% Extract trodes timestamps
disp('Extracting Mountainsort Timestamps')
fileNames = dir('*timestamps*');
timeName = fileNames.name;
timeStamps = readmda(timeName);
%all we really care about is the first time point. extract this!
firstTime = timeStamps(1);

%% Find mountainsort output. Use Manual Input

counter = 1;
unitData = zeros(2,3);
clusterVect = [s.Parameters.ClusterWindow(1):0.0001:s.Parameters.ClusterWindow(2)];
for i = 1:2
    disp('Select Shank')
    tarFileName = uigetfile('*.mda',strcat('Select Shank-',num2str(i)));
    if tarFileName
        %extract file
        holder = readmda(tarFileName);
        %this produces a three row matrix, with row 1 being the electrode with
        %the biggest signal, row 2 being the sample # relative to the first sample
        %, and row 3 being the unit number. We now want to change this to
        %seconds
        holder(2,:) = (holder(2,:)+firstTime)/s.Parameters.trodesFS;
        %determine number of units
        uniqueUnits = unique(holder(3,:));
        numUnits = length(uniqueUnits);
        disp(strcat('Found ',num2str(numUnits),' Units'))
        %determine biggest channel for each unit
        channelDesig = [];
        for j = 1:numUnits
            chanHold = holder(1,holder(3,:) == uniqueUnits(j));
            channelDesig(j) = mode(chanHold);
        end
        %store shank, electrode, and nt name
        unitData(counter:counter + numUnits - 1,1) = i;
        unitData(counter:counter + numUnits - 1,2) = channelDesig; %store biggest channel
        unitData(counter:counter + numUnits - 1,3) = uniqueUnits; %store channel name
        %store this information in combined manner
        s.OrderData(counter:counter+length(unitData)-1,:) = unitData;
        clipName = uigetfile('*.mda',strcat('Select Clips for-',num2str(i)));
        recordClips = readmda(clipName);
        %finally, pull spikes
        for j = 1:numUnits
            targetName = strcat('nt',num2str(i),'_',num2str(uniqueUnits(j)));
            desigNames{counter+j-1} = targetName;
            nameHold{counter - 1 + j} = targetName;
            targetSpikes = holder(2,holder(3,:) == uniqueUnits(j));
            s.(targetName).SpikeTimes = targetSpikes;
%             s.(targetName).BigChan = 
            %make crude autocorrelogram
            pruned = diff(targetSpikes);
            pruned(pruned>s.Parameters.ClusterWindow(2)) = [];
            s.(targetName).ISIGraph = hist(pruned,clusterVect);
            %store rpv violation number
            s.(targetName).RPVNum = length(find(pruned < s.Parameters.RPVTime));
%             %now store all spikes belonging to the target unit
%             s.(targetName).AllSpikes = -recordClips(:,:,holder(3,:) == uniqueUnits(j));
            %now pull template
            s.(targetName).TemplateSpike = recordClips(:,:,j);
        end

        counter = counter + numUnits;
    end
end
numUnits = size(unitData,1);
holder = [];


s.DesignationName = desigNames;


%generate master array for 2-d storage of important values
masterData = zeros(numUnits,10);
masterHeader = cell(10,1);
masterInd = 1;

%store shank designation info into master
% masterData(:,masterInd) = posArray(s.SortedPeakWaveOrder,1); masterHeader{masterInd} = 'Distance from Top of Shank'; masterInd = masterInd + 1;
masterData(:,masterInd) = s.OrderData(:,1); masterHeader{masterInd} = 'Shank Designation'; masterInd = masterInd + 1;
% masterData(:,masterInd) = s.SortedPeakWaveOrder; masterHeader{masterInd} = 'SimpleShankOrder'; masterInd = masterInd + 1;


%% Now pull DIO data! This will determine what things need analysis

disp('Beginning DIO Extraction for Audio')
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
length(dioTimes)
disp('Auditory DIOs detected and Saved')


%pull D2 info (syncing pulse)
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
if length(D2FileName) == 0
    [D2FileName] = functionFileFinder(subFoldersCell,'DIO','Din2');
end
D2FileName = D2FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%pulls out DIO up state onsets.
[dio2Times,dio2TimeDiff] = functionBasicDIOCheck(DIO2Data,s.Parameters.trodesFS);
length(dio2Times)
disp('Synchronization DIOs detected')

%pull D5 info (light pulse)
[D5FileName] = functionFileFinder(subFoldersCell,'DIO','D5');
if length(D5FileName) == 0
    [D5FileName] = functionFileFinder(subFoldersCell,'DIO','Din5');
end
D5FileName = D5FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIO5Data] = readTrodesExtractedDataFile(D5FileName);

%pulls out DIO up state onsets.
[dio5Times,dio5TimeDiff] = functionBasicDIOCheck(DIO5Data,s.Parameters.trodesFS);
length(dio5Times)
disp('Light DIOs detected')

%% Now lets pull out auditory related information
soundName = strcat(fileName,'Sound.mat');
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

% %Recalculate raster window and such
% s.Parameters.RasterWindow = s.Parameters.RasterWindow * toneDur;
% s.Parameters.ToneWindow = s.Parameters.ToneWindow * toneDur;
% s.Parameters.GenWindow = s.Parameters.GenWindow * toneDur;
% calcWindow = s.Parameters.calcWindow*toneDur;
% rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
% histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
% %histBinVector is for the purposes of graphing. This provides a nice axis
% %for graphing purposes
% % s.Parameters.CalcWindow = s.Parameters.CalcWindow * toneDur;
% s.Parameters.LFPWindow = s.Parameters.LFPWindow * toneDur;

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

%% Now try and fix DIO data, see if there are problems
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


% %170807 adding edit to allow for windowing by the time filter
% %We want to remove points beyond the range of the time filter. Do this
% %here.

% timeFinder = find(master(:,1) > s.TimeFilterRange(2));
% master(timeFinder,:) = [];
% totalTrialNum = length(master);
% disp('Time Filter Applied')

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
%clean up positive binary
[locoOut] = locoBinFix(s.RotaryData.SimpleLocomotion,1,s.Parameters.InterpolationStepRotary,1,3);
s.RotaryData.PosStartsEnds = [locoOut.starts,locoOut.ends];
%clean up negative binary
[locoOut] = locoBinFix(locoOut.newSimp,1,s.Parameters.InterpolationStepRotary,-2,-3);
s.RotaryData.NegStartsEnds = [locoOut.starts,locoOut.ends];
s.RotaryData.NewSimp = locoOut.newSimp;
disp('Fixed Negative and Positive Starts/Ends')


%% Extract EDR Data
fileNames = dir(homeFolder);
fileNames = {fileNames.name};
targetFileFinder = strfind(fileNames,'.EDR'); %examines names for D1
targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file

%pull filename
targetFile = fileNames{targetFileFinder};%pulls out actual file name
%extract data
[EDROut] = functionEDRPull(targetFile,s.Parameters.EDRTimeCol,s.Parameters.EDRPiezoCol,s.Parameters.EDRTTLCol);

%now we need to sync edr with trodes. 

edrSyncDiff = diff(EDROut.TTLs{3});
%check for absurdly short things. if they are present remove and adjust edr
%data as well. Do this by iteratively going through every diff value and
%determining whats going on. 

whileCounter = 1;
whileTrigger = 0;
while whileTrigger == 0;
    if whileCounter == length(edrSyncDiff)
        break
    end
    if edrSyncDiff(whileCounter) < 0.5
        EDROut.TTLs{3}(whileCounter+1) = [];
        edrSyncDiff = diff(EDROut.TTLs{3});
    else
        whileCounter = whileCounter + 1;
    end
end

%now we align ttls with each other
[xcf,lags,bounds]  = crosscorr(edrSyncDiff,dio2TimeDiff,100);
[xcMax maxInd] = max(xcf);
xcLag = lags(maxInd);
disp('CorrLag')
disp(xcLag)
disp('MaxCorr')
disp(xcMax)

if xcLag == 0
    lengthEDR = length(EDROut.TTLs{3});
    lengthMCU = length(dio2Times);
    if lengthEDR == lengthMCU
        disp('MCU and EDR TTLs matched and aligned')
    elseif lengthEDR > lengthMCU
        disp('Sync TTLs Aligned, Excess EDR, adjusting...')
        EDROut.TTLs{3}(lengthMCU+1:end)=[];
        edrSyncDiff = diff(EDROut.TTLs{3});
        disp('Excess EDR Pulses Removed')
    elseif lengthEDR < lengthMCU
        disp('Sync TTLs Aligned, Excess MCU, adjusting...')
        dio2Times(lengthEDR+1:end)=[];
        dio2TimeDiff = diff(dio2Times);
        disp('Excess MCU Pulses Removed')
    end
else
    figure
    plot(edrSyncDiff,'r')
    hold on
    plot(dio2TimeDiff,'b')
    title('EDR Sync (r) vs MCU Sync (b)')
    error('FAILURE TO ALIGN SYNCING PULSES')
end

edrMagTimes = interp1(EDROut.TTLs{3},dio2Times,EDROut.PiezoPower(:,1));

%now lets try and going through and find onsets. 
magData = EDROut.PiezoPower(:,2);
magDiff = diff(magData);
threshHiPass = 1000;
threshReset = 200;
whileTrig = 0;
whileInd = 1;
whileCounter = 1;
onsetStore = [];
while whileTrig == 0
    %first, establish break for when exceeds length
    testFind = find(magData(whileCounter:end) > threshHiPass);
    if length(testFind) == 0
        break
    end
    %now lets do our iterative search. first, find the next point exceeding
    %the high threshold
    findNext = find(magData(whileCounter:end) > threshHiPass,1,'first');
    %next, find the earliest point before with a slope of less than 0
    %(bottom)
    findBottom = find(magDiff(1:whileCounter + findNext-2) < 0, 1, 'last');
    onsetStore(whileInd) = findBottom + 2;
    %store magnitude
    peakSizeStore(whileInd) = magData(findNext + whileCounter - 1) - magData(findBottom+2);
    %now that we've stored the value, we need to find the next ok point. 
    findLim = find(magData(whileCounter + findNext:end) < threshReset,1,'first');
    whileCounter = whileCounter + findNext + findLim;
    wholeStore(whileInd,1) = findNext;
    wholeStore(whileInd,2) = findBottom;
    wholeStore(whileInd,3) = findLim;
    wholeStore(whileInd,4) = whileCounter;
    
    whileInd = whileInd + 1;  
end

figure
plot(magData,'k')
hold on
plot(onsetStore,magData(onsetStore),'r*')
plot([1 length(magData)],[1000 1000],'r')
plot([1 length(magData)],[200 200],'c')
plot(magDiff,'g')
plot([1 length(magData)],[0 0],'k')
title('EDR Peak Detection')

edrOnsetTimes = edrMagTimes(onsetStore);
edrOnsetTimes(isnan(edrOnsetTimes)) = [];


%% Now process spiking data aligned to various points
%now that we have the desired alignment times, we can basically make simple
%rasters for every single thing. Tone related will be a bit more
%complicated...


for i = 1:numUnits
    %%Process general information
    masterHolder = masterInd;
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes';
    alignTimes = master(:,1);
    
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    
    %make a fine grained one as well!
    fineSessFire = hist(spikeTimes,[s.RotaryData.Velocity(1,1):0.01:s.RotaryData.Velocity(end,1)]);
    fineSessFire(end) = 0; %this is to compensate for problems with spikes coming after the period
    fineSessFire(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    
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
    waveForms = s.(desigNames{i}).TemplateSpike;
    
    %pull out biggest waveform
    maxWave = max(-waveForms');
    [maxVal findMax] = max(maxWave);
   
    %chose the big wave, interpolate to fine degree
    chosenWave = waveForms(findMax,:);
    interpVect = [1:0.1:100];
    interpWave = interp1(1:100,chosenWave,interpVect,'spline');
    
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
    bigHistStore(:,i) = fullHistData;
    
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
    %% Process Tone Frequency Related Info
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
                targetRasters,length(targetTrials),2,targetTrials,s.Parameters.latBin,s.Parameters.histBin,s.Parameters.PercentCutoff,s.Parameters.BaselineCutoff);
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
    s.(desigNames{i}).FineSession = fineSessFire;
    s.(desigNames{i}).FineSessionVector = [s.RotaryData.Velocity(1,1):0.01:s.RotaryData.Velocity(end,1)];
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
    
    %% Now do rasters aligned to light delivery
    [rastersLight] = functionBasicRaster(spikeTimes,dio5Times,s.Parameters.RasterWindow);
    [histCounts histCenters] = hist(rastersLight(:,1),histBinVector);
    histLight = histCounts'/length(dio5Times)/s.Parameters.histBin;
    s.(desigNames{i}).LightRaster = rastersLight;
    s.(desigNames{i}).LightHist = histLight;
    bigHistLight(:,i) = histLight; %for plotting purposes!
    
    %% Now do rasters aligned to significant deviations based on EDR
    
    [rastersEDR] = functionBasicRaster(spikeTimes,edrOnsetTimes,s.Parameters.RasterWindow);
    [histCounts histCenters] = hist(rastersEDR(:,1),histBinVector);
    histEDR = histCounts'/length(edrOnsetTimes)/s.Parameters.histBin;
    s.(desigNames{i}).EDRRaster = rastersEDR;
    s.(desigNames{i}).EDRHist = histEDR;
    bigHistEDR(:,i) = histEDR; %for plotting purposes!
    
    
    %% Now do locomotion overall correlation
    
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

masterInd = masterHolder;

%now generate normalized overall histograms for light and movement
findZero = find(histBinVector < 0,1,'last');
for i = 1:numUnits
    bigHistLightNorm(:,i) = bigHistLight(:,i) / mean(bigHistLight([1:findZero],i));
    bigHistEDRNorm(:,i) = bigHistEDR(:,i) / mean(bigHistEDR([1:findZero],i));
    bigHistStoreNorm(:,i) = bigHistStore(:,i) / mean(bigHistStore([1:findZero],i));
end

%first make EDR rasters for different events
lightEDRalign = [];
extractWindow = [round(s.Parameters.RasterWindow(1)/nanmean(diff(edrMagTimes))) round(s.Parameters.RasterWindow(2)/nanmean(diff(edrMagTimes)))-1];
edrVect = [s.Parameters.RasterWindow(1):nanmean(diff(edrMagTimes)):s.Parameters.RasterWindow(2)];
for j = 1:length(dio5Times)
    targetTime = find(edrMagTimes - dio5Times(j)>0,1,'first');
    lightEDRalign(:,j) = EDROut.PiezoPower(extractWindow(1)+targetTime:extractWindow(2)+targetTime,2);
end

soundEDRAlign = [];
extractWindow = [round(s.Parameters.RasterWindow(1)/nanmean(diff(edrMagTimes))) round(s.Parameters.RasterWindow(2)/nanmean(diff(edrMagTimes)))-1];
edrVect = [s.Parameters.RasterWindow(1):nanmean(diff(edrMagTimes)):s.Parameters.RasterWindow(2)];
for j = 1:length(master)
    targetTime = find(edrMagTimes - master(j,1)>0,1,'first');
    soundEDRAlign(:,j) = EDROut.PiezoPower(extractWindow(1)+targetTime:extractWindow(2)+targetTime,2);
end

%now lets also set this to different frequencies and generate averages. 
soundEDRMean = [];
for j = 1:numFreqs
    targetTrials = master(master(:,2) == uniqueFreqs(j),4);
    soundEDRMean(:,j) = mean(soundEDRAlign(:,targetTrials)');
end

%also perform FFT of locomotion
X = s.RotaryData.Velocity(:,2);
Fs = 100;
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1floco = P2(1:L/2+1);
P1floco(2:end-1) = 2*P1floco(2:end-1);
floco = Fs*(0:(L/2))/L;


%% Plotting!!

%% Want to plot by shanks, overall responses to light, sound, locomotion. Do heatmaps?

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
findFirst = find(masterData(:,1) == 1);
findSecond = find(masterData(:,1) == 2);
%plot out responses to blue light
subplot(2,4,1)
imagesc(bigHistLightNorm(:,findFirst)')
colormap('parula')
colorbar
title('Shank1 Resp Blue Light')

subplot(2,4,5)
imagesc(bigHistLightNorm(:,findSecond)')
colormap('parula')
colorbar

%plot out responses to auditory overall
subplot(2,4,2)
imagesc(bigHistStoreNorm(:,findFirst)')
title('Shank1 Resp Auditory')
colormap('parula')
colorbar

subplot(2,4,6)
imagesc(bigHistStoreNorm(:,findSecond)')
colormap('parula')
colorbar

%plot out responses to loco start
subplot(2,4,3)
imagesc(bigHistEDRNorm(:,findFirst)')
title('Shank1 Resp EDR')
colormap('parula')
colorbar

subplot(2,4,7)
imagesc(bigHistEDRNorm(:,findSecond)')
colormap('parula')
colorbar

%plot out responses to auditory overall, reserved for A2A stim in the
%future. 
subplot(2,4,4)
title('Shank1 Resp A2A')


subplot(2,4,8)



%save

spikeGraphName = strcat(fileName,'OverallResponses');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
%% Plot out positions of units on shank


hFig = figure;
subplot(2,1,1)
firstShank = find(s.OrderData(:,1) == 1);
hist(s.OrderData(firstShank,2),[1:32])
title('Distribution of Max Wave on Shank 1')
xlim([1 32])
subplot(2,1,2)
secondShank = find(s.OrderData(:,1) == 2);
hist(s.OrderData(secondShank,2),[1:32])
title('Distribution of Max Wave on Shank 2')
xlim([1 32])
spikeGraphName = strcat(fileName,'UnitDistribution');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%% Plot out individual cells more specifically?


for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])
    %First column: templates, autocorrelogram, EDR raster, and velocity vs
    %EDR vs EDR power vs smoothed firing rate
    
    %template
    subplot(2,8,1)
    hold on
    fullMin = 0;
    fullMax = 0;
    for j = 1:32
        plot(s.(desigNames{i}).TemplateSpike(j,:)-j)
        tempMin = min(s.(desigNames{i}).TemplateSpike(j,:)-j);
        tempMax = max(s.(desigNames{i}).TemplateSpike(j,:)-j);
        if tempMin < fullMin
            fullMin = tempMin;
        end
        if tempMax > fullMax
            fullMax = tempMax;
        end
    end
    ylim([fullMin,fullMax])
    title(strcat((desigNames{i}),'Shank',num2str(s.OrderData(i)),'Template'), 'Interpreter', 'none');
    
    
    %autocorrelogram
    
    subplot(4,8,2)
    hold on
    bar(clusterVect,s.(desigNames{i}).ISIGraph)
    histMax = max(s.(desigNames{i}).ISIGraph);
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVNum/length(s.(desigNames{i}).SpikeTimes)));...
        strcat(num2str(s.(desigNames{i}).RPVNum),'/',num2str(length(s.(desigNames{i}).SpikeTimes)))})
    
    %do fft of fine firing
    
    subplot(4,8,10)
    hold on
    X = s.(desigNames{i}).FineSession;
    Fs = 100;
    L = length(X);
    Y = fft(X);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,smooth(P1,41)/max(smooth(P1,41))) 
    plot(f,smooth(P1floco,41)/max(smooth(P1floco,41)),'r')
    xlim([0.5 20])
    title('norm FFT of 10ms FR (b) loco (r)')
    
    %EDR raster
    subplot(4,4,9)
    hold on
    rasterPlot(s.(desigNames{i}).EDRRaster(:,1),s.(desigNames{i}).EDRRaster(:,2))
    xlim(s.Parameters.RasterWindow)
    ylim([0 length(edrOnsetTimes)])
    title('EDR RASTER')
    
    %big velocity etc plot
    subplot(4,1,4)
    hold on
    %plot velocity
    plot(downsample(s.RotaryData.Velocity(:,1),10),downsample(s.RotaryData.Velocity(:,2)./max(s.RotaryData.Velocity(:,2)),10),'k')
    %plot EDR power
    plot(downsample(edrMagTimes,10),downsample(EDROut.PiezoPower(:,2)./max(EDROut.PiezoPower(:,2)),10),'b')
    %plot firing rate
    plot(downsample(s.(desigNames{i}).FineSessionVector,5),downsample(smooth(s.(desigNames{i}).FineSession,31)./max(smooth(s.(desigNames{i}).FineSession,31)),5),'r')
    %restrict timeframe to EDR times.
    edrStart = edrMagTimes(find(~isnan(edrMagTimes),1,'first'));
    edrEnd = edrMagTimes(find(~isnan(edrMagTimes),1,'last'));
    xlim([edrStart edrEnd])
    if toggleROC == 1
        title(strcat('Normalized Vel (k), EDR power (b) and Firing Rate (r) 10ms bin 300ms smooth AUC:',num2str(s.(desigNames{i}).TrueAUC),'99%Range',num2str(prctile(s.(desigNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(desigNames{i}).ShuffleAUC,1))))
    else
        title('Normalized Vel (k), EDR power (b) and Firing Rate (r) 10ms bin 300ms smooth')
    end
    
    
    %Column 2: rasters for lights and sound, edr histogram
    
    subplot(4,4,2)
    hold on
    rasterPlot(s.(desigNames{i}).LightRaster(:,1),s.(desigNames{i}).LightRaster(:,2))
    xlim(s.Parameters.RasterWindow)
    ylim([0 length(dio5Times)])
    title('Light Raster')
    
    subplot(4,4,6)
    hold on
    rasterPlot(s.(desigNames{i}).AllRasters(:,1),s.(desigNames{i}).AllRasters(:,3))
    xlim(s.Parameters.RasterWindow)
    ylim([0 length(master)])
    set(gca,'Ydir','reverse')
    title('Sound Raster Organized by Freq Ascending top to bottom')
    
    subplot(4,4,10)
    plot(histBinVector,s.(desigNames{i}).EDRHist)
    xlim(s.Parameters.RasterWindow)
    title('Histogram of EDR Response')
    
    %Column 3: Histogram for light, frequency averaged heatmap for sound
    
    subplot(4,4,3)
    plot(histBinVector,s.(desigNames{i}).LightHist)
    xlim(s.Parameters.RasterWindow)
    title('Histogram of Light Response')
    
    subplot(4,4,7)
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
    title('Frequency Responses')
    
    %column 4, plot out EDR aligned to light, edr aligned to different
    %frequencies?
    
    subplot(4,4,4)
    
%     plot(edrVect,lightEDRalign)
%     hold on
    plot(edrVect,mean(lightEDRalign'),'LineWidth',2)
    title('Mean EDR Response to Light')
    
    
    subplot(4,4,8)
    imagesc(soundEDRMean')
    colormap('parula')
    colorbar
%     plot(edrVect,soundEDRAlign)
%     hold on
%     plot(edrVect,mean(soundEDRAlign'),'LineWidth',2)
    
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    close
end


%% Saving
save(fullfile(pname,fname),'s','masterData','masterHeader');


end
