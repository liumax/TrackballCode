%This code is meant to apply tuning analysis to selected kilosorted
%recordings. The goal here is to determine which units have cross
%correlogram relationships, and to examine their sorting. 



%for testing purposes
fileName = '190205ML181105C_RAudStr3667pen2rec1tuningDMR'

%% Constants


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

%% Set up file saving stuff. 

saveName = strcat(fileName,'Analysis','.mat');
fname = saveName;
pname = pwd;

%% Establishes folders
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%% Extract files from kilosort. 

%load kilosort spikes
clustfolder ='kilosortoutput\';
spfolder = 'kilosortoutput\';
%get clusters
[cids, cgs] = readClusterGroupsCSV([clustfolder 'cluster_group.tsv']);
phyGoodClusters = cids(cgs==2);
clu=phyGoodClusters;
numClu = numel(clu);

%get the full sp trains from kilosortfinal directory
spAmp = readNPY([spfolder 'amplitudes.npy']);
spTime = readNPY([spfolder 'spike_times.npy']);
spTemp = readNPY([spfolder 'spike_templates.npy']);
spClu = readNPY([spfolder 'spike_clusters.npy']);

spTimeSec=double(spTime)/30000;

%get peak channel of each template and cluster
load([spfolder 'rez.mat']);
templatesOrig = double(readNPY([spfolder 'templates.npy']))*200;
templates = zeros(size(templatesOrig,1), size(templatesOrig,2),length(rez.connected));
templates(:,:,rez.connected) = templatesOrig;
tempPerClu = findTempForEachClu(spClu, spTemp);
[~, tempPkChan] = max(squeeze(max(templates,[],2))-squeeze(min(templates,[],2)),[],2);
cluPkChan=tempPerClu;
cluPkChan(:)=NaN;
for i=1:length(cluPkChan)
   if ~isnan(tempPerClu(i))
      cluPkChan(i)=tempPkChan(tempPerClu(i)+1); 
   end
end


%This outputs a bunch of spikes and their designations. we need to
%appropriately sort these now. First lets generate the right names
cluNames = [];
numUnits = length(clu);
for i = 1:numUnits
    cluNames{i} = strcat('clu',num2str(clu(i)));
end
s.DesignationName = cluNames;
for i = 1:numUnits
    tarSpikes = find(spClu == clu(i));
    s.(cluNames{i}).SpikeTimes = spTimeSec(tarSpikes);
end

%determine peak channel for given cluster. remember that cluster names are
%zero indexed!!

s.PeakChanVals = cluPkChan(clu+1);

%generate master array for 2-d storage of important values
masterData = zeros(numUnits,10);
masterHeader = cell(1,1);
masterInd = 1;

% datPath = strcat(fileName);
datPath = strcat(fileName,'.phy.dat');
% load(datPath)
filenamestruct = dir(datPath);
dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(64*dataTypeNBytes);  % Number of samples per channel
%get channel map
chanMap = getfield(rez.ops,'chanMap');
%extract only the clusters I want
tester = ismember(spClu,clu);
newClu = spClu(tester);
newST = spTimeSec(tester);

tarPath = strcat(pwd,'\',fileName,'.phy.dat');
medWFs = extractMedianWFs(newClu, newST, 30000,tarPath, 'int16', [64,nSamp], chanMap, 1);

%store the waveforms! also extract characteristics
for i = 1:numUnits
    tarWave = squeeze(medWFs(i,:,:));
    s.(cluNames{i}).medianWave = tarWave;
    %find minima for each wave. Note that waves are now downwards pointing
    for j = 1:64
        tarWave(j,:) = tarWave(j,:) - mean(tarWave(j,1:5));
    end
    [pks inds] = min(tarWave');
    %find the minimal minima
    [metpks metinds] = min(pks);
    %lets pull the specific channel!
    selWave = tarWave(metinds,:);
    %generate smoothed version
    selWaveInt = interp1([1:length(selWave)],selWave,[1:0.1:length(selWave)],'spline');
    %now find minimum (peak)
    [selpk selind] = min(selWaveInt(1:250));
    %now find trough, which determines peak trough time
    [trghpk trphind] = max(selWaveInt(selind:end));
    pktrough = trphind/30000/10;
    
    %now lets determine half max width.
    halfFront = find(selWaveInt(1:selind) < selpk/2,1,'first');
    halfBack = find(selWaveInt(halfFront:end) > selpk/2,1,'first');
    halfWidth = halfBack/30000/10;
    masterData(i,masterInd) = pktrough;
    masterHeader{masterInd} = 'PeakTroughms';
    masterData(i,masterInd+1) = halfWidth;
    masterHeader{masterInd+1} = 'HalfMaxms';
    
end

masterInd = masterInd + 2;

% %lets try plotting things
% figure
% hold on
% for i = 1:32
% plot([1:54],squeeze(medWFs(2,i,:))/max(max(medWFs(1,:,:)))-i)
% end
% for i = 33:64
% plot([55:108],squeeze(medWFs(2,i,:))/max(max(medWFs(1,:,:)))-i+32)
% end

%extract the waveform properties from median waveforms. Use waveform with
%largest amplitude. 
 
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
histBinNum = round((s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin);
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

targetFolderSearch = dir;%pulls dir from DIO folder
targetFileNames = {targetFolderSearch.name}';%pulls names section
targetFileFinder = strfind(targetFileNames,'D1'); %examines names for D1
targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name

if length(targetFiles) == 0
    targetFolderSearch = dir;%pulls dir from DIO folder
    targetFileNames = {targetFolderSearch.name}';%pulls names section
    targetFileFinder = strfind(targetFileNames,'Din1'); %examines names for D1
    targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
    targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name
end
D1FileName = targetFiles{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);
timeFirst = double(DIOData.fields(1).data(1))/s.Parameters.trodesFS;
timeLast = double(DIOData.fields(1).data(end))/s.Parameters.trodesFS;
%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);
dioTimes = dioTimes - timeFirst; %need to do this because kilosort times are from zero, not from trodes time. 

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
    spikeTimes = s.(cluNames{i}).SpikeTimes;
    for j = 1:size(devStore,1)
        findSpikes = find(spikeTimes >= dioTimes(devStore(j,2)) & spikeTimes <= dioTimes(devStore(j,3)));
        disp(strcat(num2str(length(findSpikes)),'-spikes removed'))
        spikeTimes(findSpikes) = [];
        s.(cluNames{i}).SpikeTimes = spikeTimes;
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

targetFolderSearch = dir;%pulls dir from DIO folder
targetFileNames = {targetFolderSearch.name}';%pulls names section
targetFileFinder = strfind(targetFileNames,'D3'); %examines names for D1
targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name

if length(targetFiles) == 0
    targetFolderSearch = dir;%pulls dir from DIO folder
    targetFileNames = {targetFolderSearch.name}';%pulls names section
    targetFileFinder = strfind(targetFileNames,'Din3'); %examines names for D1
    targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
    targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name
end
D3FileName = targetFiles{1};

%pull D2 info (laser ID)
targetFolderSearch = dir;%pulls dir from DIO folder
targetFileNames = {targetFolderSearch.name}';%pulls names section
targetFileFinder = strfind(targetFileNames,'D4'); %examines names for D1
targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name

if length(targetFiles) == 0
    targetFolderSearch = dir;%pulls dir from DIO folder
    targetFileNames = {targetFolderSearch.name}';%pulls names section
    targetFileFinder = strfind(targetFileNames,'Din4'); %examines names for D1
    targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
    targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name
end
D4FileName = targetFiles{1};

[funcOut] = functionNewRotaryExtraction(D3FileName,D4FileName);
%now we need to fix the timing.
funcOut.Distance(:,1) = funcOut.Distance(:,1) - timeFirst*s.Parameters.trodesFS;
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

[overOut,indivOut,s,masterData,masterHeader,masterInd] = functionTuningDataExtraction(numUnits,numDBs,numFreqs,uniqueFreqs,s,masterData,masterHeader,masterInd,histBinVector,trialNum,master,cluNames,calcWindow,histBinNum,whiteStatus);
s.NonLaserOverall = overOut;
for i = 1:numUnits
    fn = fieldnames(indivOut.(cluNames{i}));
    for j = 1:length(fn)
        s.(cluNames{i}).(fn{j}) = indivOut.(cluNames{i}).(fn{j});
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
    spikeStore = s.(cluNames{i}).SpikeTimes(s.(cluNames{i}).SpikeTimes > spikeHistVect(1) + s.Parameters.STRFWin(1) & s.(cluNames{i}).SpikeTimes < spikeHistVect(end));
    spikeArray(i,:) = hist(spikeStore,spikeHistVect);
end

[sta, stabigmat, spkcountvec] = quick_calc_sta(stimulus, spikeArray, numLags);
% [sta_sig, ptd, siglevel] = ne_sig_sta_from_stim_obs_resp(sta, spikeArray, stimulus, 10, numLags, 95);

s.STAs = sta;
% s.STASig = sta_sig;

% %now lets generate jittered spike trains and see what this produces?
% jitterNum = 100;
% jitterTime = 0.3; %+/- 150 ms. 
% for j = 1:jitterNum
%     disp(j)
%     for i = 1:numUnits
%         
%         spikeStore = s.(cluNames{i}).SpikeTimes(s.(cluNames{i}).SpikeTimes > spikeHistVect(1) + s.Parameters.STRFWin(1) & s.(cluNames{i}).SpikeTimes < spikeHistVect(end));
% 
%         jitterVals = rand(length(spikeStore),1)-0.5;
%         fakeSpikes = spikeStore + jitterVals;
%         tempSpikes(i,:) = hist(fakeSpikes,spikeHistVect);
%         
%     end
%     [sta, stabigmat, spkcountvec] = quick_calc_sta(stimulus, tempSpikes, numLags);
%     jitSTAStore{j} = sta;
% end
% 
% %also generate split spike STAs to have internal comparison. 
% for i = 1:numUnits
%     spikeStore = s.(cluNames{i}).SpikeTimes(s.(cluNames{i}).SpikeTimes > spikeHistVect(1) + s.Parameters.STRFWin(1) & s.(cluNames{i}).SpikeTimes < spikeHistVect(end));
%     sel1 = randperm(length(spikeStore),round(length(spikeStore)/2));
%     sel2 = [1:length(spikeStore)];
%     sel2(ismember(sel1,sel2)) = [];
%     spikes1(i,:) = hist(spikeStore(sel1),spikeHistVect);
%     spikes2(i,:) = hist(spikeStore(sel2),spikeHistVect);
% end
% [sta1, stabigmat, spkcountvec] = quick_calc_sta(stimulus, spikes1, numLags);
% [sta2, stabigmat, spkcountvec] = quick_calc_sta(stimulus, spikes2, numLags);
% 
% s.JitteredSTAs = jitSTAStore;
% s.SplitSTA1 = sta1;
% s.SplitSTA2 = sta2;


s.DMRTimes = trueDMRtimes;
s.DMRStep = dmrStep;
s.DMRfaxis = faxis;

s.SpikeArray = spikeArray;

% %going to go through each thing and perform STRFs
% for i = 1:numUnits
%     disp(strcat('Working on Unit ',num2str(i)))
%     %pulls spike times and times for alignment
%     spikeTimes = s.(cluNames{i}).SpikeTimes;
%     %lets pull the stuff based on all spikes
%     
%     indivStore = [];
%     avStore = [];
%     %pull out only the spikes during the desired window
%     firstSpike = find(spikeTimes > dmrStore(1) - s.Parameters.STRFWin(1) ,1,'first');
%     lastSpike = find(spikeTimes < dmrStore(end) - s.Parameters.STRFWin(2) ,1,'last');
%     targetSpikes = spikeTimes(firstSpike:lastSpike);
%     if targetSpikes
%         masterData(i,masterInd) = length(targetSpikes);
%         
%         for k = 1:length(targetSpikes)
%             findTarget = find(dmrStore - targetSpikes(k) > 0,1,'first');
%             takeChunk = stimulus(:,findTarget+newWin(1):findTarget+newWin(2));
%             indivStore(:,:,k) = takeChunk;
%         end
% %         avStore(:,:,i) = mean(indivStore,3);
% %         s.(cluNames{i}).IndivStim = indivStore;
%         s.(cluNames{i}).DMRAverage = mean(indivStore,3);
% %         bigAvStore(:,:,i) = mean(squeeze(avStore(:,:,i,:)),3);
% %         s.(cluNames{i}).FullAverage = mean(squeeze(avStore(:,:,i,:)),3);
%     else
% %         avStore(:,:,i) = ;
%         s.(cluNames{i}).DMRAverage = zeros(size(stimulus,1),newWin(2)-newWin(1)+1);
% %         bigAvStore(:,:,i) = mean(squeeze(avStore(:,:,i)),3);
% %         s.(cluNames{i}).FullAverage = mean(squeeze(avStore(:,:,i)),3);
%     end
% end
% masterHeader{masterInd} = 'NumSpikesDMR';
% masterInd = masterInd + 1;




%% Saving
save(fullfile(pname,fname),'s','masterData','masterHeader');
% save(fullfile(pname,strcat(fileName,'DMRData')),'sta','sta_sig','spikeArray','stimulus','faxis')
save(fullfile(pname,strcat(fileName,'DMRData')),'sta','spikeArray','stimulus','faxis')

%now lets remove other stuff
rez = [];
spikeArray = [];
stimulus = [];
indivOut = [];

%% Plotting

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

%First, lets plot waveforms. Split by shank.
plotVect = [0:1/30000:53/30000];
[indPkTr] = functionCellStringFind(masterHeader,'PeakTrough');
indPkTr = indPkTr(1);
for i = 1:2
    findTars = find(s.PeakChanVals <= 32*i & s.PeakChanVals > 32*(i-1)+1);
    if findTars
        hFig = figure;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.07], [0.01 0.01]);
        set(hFig, 'Position', [10 10 1900 1000])
        for j = 1:length(findTars)
            subplot(1,length(findTars),j)
            hold on
            for k = 1:32
                if masterData(findTars(j),indPkTr) < s.Parameters.PVLim(1) %is FSI
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/max(max(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[1 0 0]);
                elseif masterData(findTars(j),indPkTr) > s.Parameters.PVLim(2) %is MSN
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/max(max(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[0 0 0]);
                else %unknown
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/max(max(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[.7 .7 .7]);
                end
            end
            title({(cluNames{findTars(j)});strcat('Num Spikes:',num2str(length(s.(cluNames{findTars(j)}).SpikeTimes)));strcat('Max Amp:',num2str(max(max(s.(cluNames{findTars(j)}).medianWave))))})
            xlim([plotVect(1) plotVect(end)])
            ylim([-33 1])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
        end
        spikeGraphName = strcat(fileName,'Shank',num2str(i),'Waves');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end


%first do cross correlograms. We want to first generate some kind of
%shuffled comparison to calculate significance by which we can extract
%putatively monosynaptic linked units. 
shuffSpikes = cell(numUnits,1);
shuffNum = 1000;
shuffVal = 0.01; %since multiply by 0.5, need to double here. https://www.nature.com/articles/nn.2134#methods this paper uses +/- 5 ms
for i = 1:numUnits
    %determine number of spikes
    tarSpikes = s.(cluNames{i}).SpikeTimes;
    %generate random array
    randArray = rand(length(tarSpikes),shuffNum) - 0.5;
    shuffSpikes{i} = randArray*shuffVal + repmat(tarSpikes,1,shuffNum);
end

s.JitteredSpikes = shuffSpikes;

%now lets actually do the targeted cross correlograms! Try to only target FSIs and
%align MSN spikes to them. 
findFSIs =find(masterData(:,indPkTr) < s.Parameters.PVLim(1));
findMSNs = find(masterData(:,indPkTr) > s.Parameters.PVLim(2));

synchWind = 0.001; %time difference within which spikes are considered synchronous.
rasterWindow = [-0.01025 0.01025];
rasterVector = [-0.01:0.0005:0.01];


for i = 1:numUnits
    for j = 1:numUnits
        disp(strcat('Running Cross Corr for Unit',num2str(cluNames{j}),'-aligned to',num2str(cluNames{i})))
        [rasters] = crosscorrelogram(s.(cluNames{i}).SpikeTimes,s.(cluNames{j}).SpikeTimes,rasterWindow);
        histStore = hist(rasters,rasterVector);
        rasters = [];
        trueStore{(i),(j)} = histStore;
    end
end

s.CrossCorrs = trueStore;

for i = 1:length(findFSIs)
    for j = 1:length(findMSNs)
%         disp(strcat('Running Cross Corr for Unit',num2str(cluNames{findMSNs(j)}),'-aligned to',num2str(cluNames{findFSIs(i)})))
%         [rasters] = crosscorrelogram(s.(cluNames{findFSIs(i)}).SpikeTimes,s.(cluNames{findMSNs(j)}).SpikeTimes,rasterWindow);
%         histStore = hist(rasters,rasterVector);
%         rasters = [];
%         trueStore{findFSIs(i),findMSNs(j)} = histStore;
        %now lets try and remove synchronous FSIs spikes. This way, we can
        %eliminate any peak in the FSI aligned MSN spikes. 
        diffArray = [];
        newFSIspike = s.(cluNames{findFSIs(i)}).SpikeTimes;
        for k = length(s.(cluNames{findMSNs(j)}).SpikeTimes):-1:1
            diffArray = s.(cluNames{findFSIs(i)}).SpikeTimes - s.(cluNames{findMSNs(j)}).SpikeTimes(k);
            spikeFind = find(abs(diffArray) < synchWind);
            newFSIspike(spikeFind) = [];
        end
        [rasters] = crosscorrelogram(newFSIspike,s.(cluNames{findMSNs(j)}).SpikeTimes,rasterWindow);
        histStore = hist(rasters,rasterVector);
        rasters = [];
        remStore{findFSIs(i),findMSNs(j)} = histStore;
        
        
        %now lets run the jittered spikes
        tempHist = [];
        shuffTar = shuffSpikes{findMSNs(j)};
        disp(strcat('Running Shuffle Cross Corr for Unit',num2str(cluNames{findMSNs(j)}),'-aligned to',num2str(cluNames{findFSIs(i)})))
        for k = 1:shuffNum
            [rasters] = crosscorrelogram(s.(cluNames{findFSIs(i)}).SpikeTimes,shuffTar(:,k),rasterWindow);
            rasters(rasters < rasterWindow(1) | rasters > rasterWindow(2)) = [];
            tempHist(:,k) = hist(rasters,rasterVector);
        end
        shuffStore{findFSIs(i),findMSNs(j)} = tempHist;
    end
end

s.CrossCorrJitter = shuffStore;
s.CrossCorrRemSynch = remStore;

hFig = figure;
[B,I] = sort(s.PeakChanVals);
rasterWindow = [-0.01025 0.01025];
rasterVector = [-0.01:0.0005:0.01];
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.02], [0.01 0.01]);
set(hFig, 'Position', [10 10 1800 1000])
for i = 1:numUnits
    disp(strcat('unit',num2str(i)))
    %determine units greater than or equal to current unit
%     tarUnits = i:numUnits;
    for j = 1:numUnits
        disp('Next Unit!')
        histStore = s.CrossCorrs{i,j};
        if j == i
            histStore(21) = 0;
        end
        if ismember(I(i),findFSIs) && ismember(I(j),findMSNs);
            subplot(numUnits,numUnits,numUnits*(i-1)+j)
            hold on
            plot(rasterVector,mean(s.CrossCorrJitter{I(i),I(j)}'),'Color',[0.7 0.7 0.7],'LineWidth',2)
            Y = prctile(s.CrossCorrJitter{I(i),I(j)},[.5 99.5],2);
            plot(rasterVector,Y(:,1),'Color',[0.7 0.7 0.7])
            plot(rasterVector,Y(:,2),'Color',[0.7 0.7 0.7])
            plot(rasterVector,smooth(s.CrossCorrRemSynch{I(i),I(j)},3),'c','LineWidth',2)
            plot(rasterVector,smooth(histStore,3),'r','LineWidth',2)
            hold on
            plot([0 0],[0 max(histStore)],'r')
            if max(histStore) < 1
                ylim([0 1])
            else
                ylim([0 max(histStore)])
            end
            xlim([-0.01 0.01])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'TickLength',[0 0])
            rasters = [];
        elseif ismember(I(i),findFSIs) && ~ismember(I(j),findMSNs);
            subplot(numUnits,numUnits,numUnits*(i-1)+j)
            hold on
%             plot(rasterVector,mean(shuffStore{I(i),I(j)}'),'Color',[0.7 0.7 0.7],'LineWidth',2)
%             plot(rasterVector,mean(shuffStore{I(i),I(j)}')+std(shuffStore{I(i),I(j)}'/sqrt(shuffNum)),'Color',[0.7 0.7 0.7])
%             plot(rasterVector,mean(shuffStore{I(i),I(j)}')-std(shuffStore{I(i),I(j)}'/sqrt(shuffNum)),'Color',[0.7 0.7 0.7])
            plot(rasterVector,smooth(histStore,3),'r','LineWidth',2)
            hold on
            plot([0 0],[0 max(histStore)],'r')
            if max(histStore) < 1
                ylim([0 1])
            else
                ylim([0 max(histStore)])
            end
            xlim([-0.01 0.01])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'TickLength',[0 0])
            rasters = [];
        else
            subplot(numUnits,numUnits,numUnits*(i-1)+j)
            plot(rasterVector,smooth(histStore,3),'k','LineWidth',2)
            hold on
            plot([0 0],[0 max(histStore)],'k')
            if max(histStore) < 1
                ylim([0 1])
            else
                ylim([0 max(histStore)])
            end
            xlim([-0.01 0.01])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'TickLength',[0 0])
            rasters = [];
        end
        
        %grey out non-identified units
        if ~ismember(I(j),[findFSIs;findMSNs])
            plot([0 0],[0 max(histStore)],'LineWidth',2,'Color',[0.7 0.7 0.7])
            plot(rasterVector,smooth(histStore,3),'LineWidth',2,'Color',[0.7 0.7 0.7])
        end
        
        %now lets try and plot out titles only for top row.
        if i == 1
            title(cluNames{I(j)})
        end
        
    end
end
spikeGraphName = strcat(fileName,'AllCrossCorr');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and separate out tone and non-tone period spikes
spikesNonTone = [];
spikesTone = [];
cluNonTone = [];
cluTone = [];
allSpikes = spTimeSec;
allClu = spClu;
counterTone = 1;
counterNonTone = 1; 
for i = 1:length(dioTimes)
    if rem(i,500) == 0 
        disp('500 trials')
    end
    %subtract tone time
    toneSub = allSpikes - dioTimes(i);
    %extract all spikes up to the end of the tone
    subSpikes = toneSub(toneSub < 0.1);
    %find negative (non tone) and positive (tone) values
    findNeg = find(subSpikes <= 0);
    findPos = find(subSpikes > 0);
    %store both spike times and cluster designations
    spikesNonTone(counterNonTone:counterNonTone - 1 + length(findNeg)) = allSpikes(findNeg);
    spikesTone(counterTone:counterTone - 1 + length(findPos)) = allSpikes(findPos);
    cluNonTone(counterNonTone:counterNonTone - 1 + length(findNeg)) = allClu(findNeg);
    cluTone(counterTone:counterTone - 1 + length(findPos)) = allClu(findPos);
    %now we need to subtract out these spikes
    allSpikes(1:length(subSpikes)) = [];
    allClu(1:length(subSpikes)) = [];
    counterTone = counterTone + length(findPos);
    counterNonTone = counterNonTone + length(findNeg);
end
spikesNonTone(counterNonTone:counterNonTone - 1 + length(allSpikes)) = allSpikes;
cluNonTone(counterNonTone:counterNonTone - 1 + length(allSpikes)) = allClu;

%now we can separate out the clusters.

% try doing this with either tone window spikes and non-tone spikes
hFig = figure;
[B,I] = sort(s.PeakChanVals);
rasterWindow = [-0.01 0.01];
rasterVector = [-0.01:0.0005:0.01];
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.02], [0.01 0.01]);
set(hFig, 'Position', [10 80 1000 1000])
for i = 1:numUnits
    %determine units greater than or equal to current unit
%     tarUnits = i:numUnits;
    tarSpikesNonTone = spikesNonTone(cluNonTone == clu(I(i)));
    tarSpikesTone = spikesTone(cluTone == clu(I(i)));
    for j = 1:numUnits
        disp('Next Unit!')
        %plot out non-tone crossCorr
        [rasters] = crosscorrelogram(tarSpikesNonTone,spikesNonTone(cluNonTone == clu(I(j))),rasterWindow);
        histStore = hist(rasters,rasterVector);
        if j == i
            histStore(21) = 0;
        end
        nonToneStore{i,j} = histStore;
        subplot(numUnits,numUnits,numUnits*(i-1)+j)
        plot(rasterVector,histStore,'k','LineWidth',2)
        max1 =max(histStore);
        hold on
        %plot out tone cross corr
        [rasters] = crosscorrelogram(tarSpikesTone,spikesTone(cluTone == clu(I(j))),rasterWindow);
        histStore = hist(rasters,rasterVector);
        ToneStore{i,j} = histStore;
        if j == i
            histStore(21) = 0;
        end
        
        plot(rasterVector,histStore,'m','LineWidth',2)
        max2 = max(histStore);
        plot([0 0],[0 max([max1,max2])],'r')
        if max([max1,max2]) < 1
            ylim([0 1])
        else
            ylim([0 max([max1,max2])])
        end
        xlim([-0.01 0.01])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'TickLength',[0 0])
        rasters = [];
    end
end
spikeGraphName = strcat(fileName,'TonevsBaseCrossCorr');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

s.CrossCorrTone = ToneStore;
s.CrossCorrNonTone = nonToneStore;

%now try and pull spikes using a shuffle in which I shift by one trial (in
%terms of the same stimulus). We need to first check and find these trials,
%and figure out how to index it. 

%first thing to do is basically streamline the dBs.
newDBs = s.SoundData.dBs;
for i = 1:s.SoundData.NumFreqs
    findTars = find(s.SoundData.Frequencies == s.SoundData.UniqueFrequencies(i));
    tarDBs = s.SoundData.dBs(findTars);
    uniqueDBVals = unique(tarDBs);
    for j = 1:length(uniqueDBVals)
        findFreq = find(s.SoundData.Frequencies == s.SoundData.UniqueFrequencies(i));
        findDBs = find(s.SoundData.dBs == uniqueDBVals(j));
        specFind = intersect(findFreq,findDBs);
        newDBs(specFind) = s.SoundData.UniqueDBs(j);
    end
end

%now lets try and find the presentations of all the unique stimuli. First,
%make a lookup table. 
oldInd = [1:length(dioTimes)];
newInd = [];
counter = 1;
for i= 1:s.SoundData.NumFreqs
    for j = 1:s.SoundData.NumDBs
        %find trials of particular frequency and db. 
        findFreq = find(s.SoundData.Frequencies == s.SoundData.UniqueFrequencies(i));
        findDBs = find(newDBs == s.SoundData.UniqueDBs(j));
        specFind = intersect(findFreq,findDBs);
        newInd(specFind,1) = counter;
        newInd(specFind,2) = [1:length(specFind)];
        counter = counter + 1;
    end
end

subInd = [];
for i = 1:counter - 1
    testFind = find(newInd(:,1) == i);
    subInd(testFind(1:end-1)) = testFind(2:end);
    subInd(testFind(end)) = testFind(1);
end
subInd = subInd';
if length(subInd) > length(dioTimes)
    subInd(length(dioTimes)+1:end) = [];
end

hFig = figure;
[B,I] = sort(s.PeakChanVals);
tarWindow = [-0.4 0.4];
rasterWindow = [-0.01 0.01];
rasterVector = [-0.01:0.0005:0.01];
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.02], [0.01 0.01]);
set(hFig, 'Position', [10 80 1000 1000])
for i = 1:numUnits
    for j = 1:numUnits
        histStore = zeros(1,length(rasterVector));
        disp('Next Unit!')
        %plot out non-tone crossCorr
        for k = 1:length(dioTimes)
            tarSpikes = s.(cluNames{I(i)}).SpikeTimes - dioTimes(k);
            tarSpikes(tarSpikes < tarWindow(1) | tarSpikes > tarWindow(2)) = [];
            tarSpikesMatch = s.(cluNames{I(j)}).SpikeTimes - dioTimes(subInd(k));
            tarSpikesMatch(tarSpikesMatch < tarWindow(1) | tarSpikesMatch > tarWindow(2)) = [];
            if length(tarSpikes)>0 & length(tarSpikesMatch) > 0
                [rasters] = crosscorrelogram(tarSpikes,tarSpikesMatch,rasterWindow);
                tempHist = hist(rasters,rasterVector);
                histStore = histStore + tempHist;
                trialScrambleStore{i,j} = histStore;
            end
            
        end

%         if j == i
%             histStore(21) = 0;
%         end
        subplot(numUnits,numUnits,numUnits*(i-1)+j)
        plot(rasterVector,histStore,'k','LineWidth',2)
        max1 =max(histStore);
        hold on
        %plot out tone cross corr
        
        plot([0 0],[0 max1],'r')
        if max1 < 1
            ylim([0 1])
        else
            ylim([0 max1])
        end
        xlim([-0.01 0.01])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'TickLength',[0 0])
        rasters = [];
    end
end
spikeGraphName = strcat(fileName,'TrialShuffCrossCorr');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

s.CrossCorrTrialShuffle = trialScrambleStore;

 
%% Plotting individual units


for i = 1:numUnits
    spikeTimes = s.(cluNames{i}).SpikeTimes;
    findTuningSpikes = [1:find(spikeTimes < dmrStore(1),1,'last')];
    findDMRSpikes = [find(spikeTimes < dmrStore(1),1,'last'):find(spikeTimes < dmrStore(end) - s.Parameters.STRFWin(2) ,1,'last')];
    numTuningSpikes(i) = length(findTuningSpikes);
    numDMRSpikes(i) = length(findDMRSpikes);
end

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])
    % Column 1
    %plots average waveform
    subplot(4,4,1)
    hold on
    plot(s.(cluNames{i}).medianWave,'k','LineWidth',2) 
%     plot(s.(cluNames{i}).MeanDMRWave,'r','LineWidth',2)
    if ismember(i,findFSIs)
        title(strcat('PV AvRate:',num2str(s.(cluNames{i}).AverageRate),'Tuning',num2str(numTuningSpikes(i)),'DMR',num2str(numDMRSpikes(i))))
    elseif ismember(i,findMSNs)
        title(strcat('MSN AvRate:',num2str(s.(cluNames{i}).AverageRate),'Tuning',num2str(numTuningSpikes(i)),'DMR',num2str(numDMRSpikes(i))))
%     elseif ismember(i,findCHATs)
%         title(strcat('CHAT AverageFiringRate:',num2str(s.(cluNames{i}).AverageRate)))
    else
        title(strcat('UNK AvRate:',num2str(s.(cluNames{i}).AverageRate),'Tuning',num2str(numTuningSpikes(i)),'DMR',num2str(numDMRSpikes(i))))
    end

    
    % plot histogram.
    subplot(4,4,5)
    plot(histBinVector,s.(cluNames{i}).AllHistograms,'k','LineWidth',2)
    hold on
    plot(histBinVector,s.(cluNames{i}).AllHistograms - s.(cluNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
    plot(histBinVector,s.(cluNames{i}).AllHistograms + s.(cluNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
    %plot significant values
    plot(s.(cluNames{i}).AllHistogramSig.Centers(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,3) == 1),...
        s.(cluNames{i}).AllHistogramSig.Histogram(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,3) == 1,1),...
        'b*')
    plot(s.(cluNames{i}).AllHistogramSig.Centers(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,4) == 1),...
        s.(cluNames{i}).AllHistogramSig.Histogram(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,4) == 1,1),...
        'bo')
    %plot negative values for first tuning
    plot(s.(cluNames{i}).AllHistogramSig.Centers(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
        s.(cluNames{i}).AllHistogramSig.Histogram(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
        'k*')
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    histLims = max(s.(cluNames{i}).AllHistograms + s.(cluNames{i}).HistogramStandardDeviation);
    ylim([0 histLims])
    title('Histogram')

    %plot out rasters, organized!
    subplot(2,4,5)
    plot(s.(cluNames{i}).AllRasters(:,1),...
        s.(cluNames{i}).AllRasters(:,3),'k.','markersize',4)
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
%     %plots ISI
%     subplot(4,4,2)
%     hist(s.(cluNames{i}).ISIGraph,1000)
%     histMax = max(hist(s.(cluNames{i}).ISIGraph,1000));
%     line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
%     xlim(s.Parameters.ClusterWindow)
%     title({strcat('ISI RPV %: ',num2str(s.(cluNames{i}).RPVPercent));...
%         strcat(num2str(s.(cluNames{i}).RPVNumber),'/',num2str(s.(cluNames{i}).TotalSpikeNumber))})

    %plot FR and velocity
    subplot(4,4,6)
    hold on
    plot(newTimeVector,newVelVector/max(newVelVector),'b')
    plot([newTimeVector(1):s.Parameters.SpeedFiringBins:newTimeVector(end)],s.(cluNames{i}).SessionFiring/max(s.(cluNames{i}).SessionFiring),'r')
    xlim([newTimeVector(1),newTimeVector(end)])
    ylim([-0.1,1])
    title({fileName;cluNames{i}},'fontweight','bold', 'Interpreter', 'none');
    
    %plot heatmap tuning curves
    subplot(4,4,10)
    imagesc(s.(cluNames{i}).BinTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Resp (tone)')
    
    subplot(4,4,14)
    imagesc(s.(cluNames{i}).BinGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Resp (general)')
    
    
    % Column 3
    
    %plot DMR response
    subplot(2,2,2)
    imagesc(reshape(sta(i,:),length(faxis),[]))
    colormap('parula')
    colorbar
    format short
    set(gca,'XTick',[0:20:numLags]);
    set(gca,'XTickLabel',[dmrStep*100:-20*dmrStep:dmrStep]);
    set(gca,'YTick',[1:10:length(faxis)]);
    set(gca,'YTickLabel',[faxis([1:10:end])]);
    title(strcat('DMR STA Based on',num2str(numDMRSpikes(i))))
    
%     %plot out thresholded DMR
%     subplot(2,2,4)
%     imagesc(reshape(sta_sig(i,:),length(faxis),[]))
%     colormap('parula')
%     colorbar
%     format short
%     set(gca,'XTick',[0:20:numLags]);
%     set(gca,'XTickLabel',[dmrStep*100:-20*dmrStep:dmrStep]);
%     set(gca,'YTick',[1:10:length(faxis)]);
%     set(gca,'YTickLabel',[faxis([1:10:end])]);
%     title('Thresholded DMR a=0.95')
    
%     freqs
    
    hold off
    spikeGraphName = strcat(fileName,cluNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
    close

end
























































%% Pull response data!
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
    spikeTimes = s.(cluNames{i}).SpikeTimes;
    alignTimes = master(:,1);
    
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
    
    disp(strcat('Baseline Spikes:',num2str(generalResponseHist.SpikeNumber),' Unit:',(cluNames{i})))
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
    s.(cluNames{i}).AllRasters = fullRasterData;
    s.(cluNames{i}).AllHistograms = fullHistData;
    s.(cluNames{i}).IndividualHistograms = fullHistHolder; 
    s.(cluNames{i}).HistogramStandardDeviation = histSTE;
    s.(cluNames{i}).FreqDBRasters = organizedRasters;
    s.(cluNames{i}).FreqDBHistograms = organizedHist;
    s.(cluNames{i}).FreqDBHistogramErrors = histErr;
    s.(cluNames{i}).FrequencyHistograms = freqSpecHist;
    s.(cluNames{i}).AverageRate = averageRate;
%     s.(cluNames{i}).SessionFiring = sessionFiring;
    s.(cluNames{i}).AverageSTD = averageSTD;
    s.(cluNames{i}).AverageSTE = averageSTE;
    s.(cluNames{i}).HistBinVector = histBinVector;
    s.(cluNames{i}).AllHistogramSig = generalResponseHist;
    s.(cluNames{i}).SpecHistogramSig = responseHistHolder;
    s.(cluNames{i}).LatPeakBin = latePeakBinStore;
    s.(cluNames{i}).LatencyMap = latStore;
    s.(cluNames{i}).PeakMapGen = peakStoreGen;
    s.(cluNames{i}).PeakMapTone = peakStoreTone;
    s.(cluNames{i}).BinFast = binStoreFast;
    s.(cluNames{i}).BinTone = binStoreTone;
    s.(cluNames{i}).BinGen = binStoreGen;
    s.(cluNames{i}).BinDiff = binDiff;
    s.(cluNames{i}).ProbTone = probStoreTone;
    s.(cluNames{i}).ProbGen = probStoreGen;
    s.(cluNames{i}).BinSigVals = binSigVals;
    s.(cluNames{i}).WidthData = widthOut;
    
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
    
    
end
s.WidthData = widthStore;
s.WidthLatData = widthLat;
s.FullWidth = bigWidth;

masterInd = masterHolder;



%% Plotting

%First, lets plot waveforms. Split by shank.
plotVect = [0:1/30000:53/30000];
[indPkTr] = functionCellStringFind(masterHeader,'PeakTrough');
indPkTr = indPkTr(1);
for i = 1:2
    findTars = find(s.PeakChanVals <= 32*i & s.PeakChanVals > 32*(i-1)+1);
    if findTars
        hFig = figure;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.07], [0.01 0.01]);
        set(hFig, 'Position', [10 10 1900 1000])
        for j = 1:length(findTars)
            subplot(1,length(findTars),j)
            hold on
            for k = 1:32
                if masterData(findTars(j),indPkTr) < s.Parameters.PVLim(1) %is FSI
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/max(max(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[1 0 0]);
                elseif masterData(findTars(j),indPkTr) > s.Parameters.PVLim(2) %is MSN
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/max(max(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[0 0 0]);
                else %unknown
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/max(max(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[.7 .7 .7]);
                end
            end
            title({(cluNames{findTars(j)});strcat('Num Spikes:',num2str(length(s.(cluNames{findTars(j)}).SpikeTimes)));strcat('Max Amp:',num2str(max(max(s.(cluNames{findTars(j)}).medianWave))))})
            xlim([plotVect(1) plotVect(end)])
            ylim([-33 1])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
        end
        spikeGraphName = strcat(fileName,'Shank',num2str(i),'Waves');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end


%first do cross correlograms. We want to first generate some kind of
%shuffled comparison to calculate significance by which we can extract
%putatively monosynaptic linked units. 
shuffSpikes = cell(numUnits,1);
shuffNum = 1000;
shuffVal = 0.01; %since multiply by 0.5, need to double here. https://www.nature.com/articles/nn.2134#methods this paper uses +/- 5 ms
for i = 1:numUnits
    %determine number of spikes
    tarSpikes = s.(cluNames{i}).SpikeTimes;
    %generate random array
    randArray = rand(length(tarSpikes),shuffNum) - 0.5;
    shuffSpikes{i} = randArray*shuffVal + repmat(tarSpikes,1,shuffNum);
end

s.JitteredSpikes = shuffSpikes;

%now lets actually do the targeted cross correlograms! Try to only target FSIs and
%align MSN spikes to them. 
findFSIs =find(masterData(:,indPkTr) < s.Parameters.PVLim(1));
findMSNs = find(masterData(:,indPkTr) > s.Parameters.PVLim(2));

synchWind = 0.001; %time difference within which spikes are considered synchronous.
rasterWindow = [-0.01025 0.01025];
rasterVector = [-0.01:0.0005:0.01];


for i = 1:numUnits
    for j = 1:numUnits
        disp(strcat('Running Cross Corr for Unit',num2str(cluNames{j}),'-aligned to',num2str(cluNames{i})))
        [rasters] = crosscorrelogram(s.(cluNames{i}).SpikeTimes,s.(cluNames{j}).SpikeTimes,rasterWindow);
        histStore = hist(rasters,rasterVector);
        rasters = [];
        trueStore{(i),(j)} = histStore;
    end
end

s.CrossCorrs = trueStore;
trueStore = [];

for i = 1:length(findFSIs)
    for j = 1:length(findMSNs)
%         disp(strcat('Running Cross Corr for Unit',num2str(cluNames{findMSNs(j)}),'-aligned to',num2str(cluNames{findFSIs(i)})))
%         [rasters] = crosscorrelogram(s.(cluNames{findFSIs(i)}).SpikeTimes,s.(cluNames{findMSNs(j)}).SpikeTimes,rasterWindow);
%         histStore = hist(rasters,rasterVector);
%         rasters = [];
%         trueStore{findFSIs(i),findMSNs(j)} = histStore;
        %now lets try and remove synchronous FSIs spikes. This way, we can
        %eliminate any peak in the FSI aligned MSN spikes. 
        diffArray = [];
        newFSIspike = s.(cluNames{findFSIs(i)}).SpikeTimes;
        for k = length(s.(cluNames{findMSNs(j)}).SpikeTimes):-1:1
            diffArray = s.(cluNames{findFSIs(i)}).SpikeTimes - s.(cluNames{findMSNs(j)}).SpikeTimes(k);
            spikeFind = find(abs(diffArray) < synchWind);
            newFSIspike(spikeFind) = [];
        end
        [rasters] = crosscorrelogram(newFSIspike,s.(cluNames{findMSNs(j)}).SpikeTimes,rasterWindow);
        histStore = hist(rasters,rasterVector);
        rasters = [];
        remStore{findFSIs(i),findMSNs(j)} = histStore;
        
        
        %now lets run the jittered spikes
        tempHist = [];
        shuffTar = shuffSpikes{findMSNs(j)};
        disp(strcat('Running Shuffle Cross Corr for Unit',num2str(cluNames{findMSNs(j)}),'-aligned to',num2str(cluNames{findFSIs(i)})))
        for k = 1:shuffNum
            [rasters] = crosscorrelogram(s.(cluNames{findFSIs(i)}).SpikeTimes,shuffTar(:,k),rasterWindow);
            rasters(rasters < rasterWindow(1) | rasters > rasterWindow(2)) = [];
            tempHist(:,k) = hist(rasters,rasterVector);
        end
        shuffStore{findFSIs(i),findMSNs(j)} = tempHist;
    end
end

s.CrossCorrJitter = shuffStore;
s.CrossCorrRemSynch = remStore;

hFig = figure;
[B,I] = sort(s.PeakChanVals);
rasterWindow = [-0.01025 0.01025];
rasterVector = [-0.01:0.0005:0.01];
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.02], [0.01 0.01]);
set(hFig, 'Position', [10 10 1800 1000])
for i = 1:numUnits
    %determine units greater than or equal to current unit
%     tarUnits = i:numUnits;
    for j = 1:numUnits
        disp('Next Unit!')
        
        [rasters] = crosscorrelogram(s.(cluNames{I(i)}).SpikeTimes,s.(cluNames{I(j)}).SpikeTimes,rasterWindow);
        histStore = hist(rasters,rasterVector);
        if j == i
            histStore(21) = 0;
        end
        if ismember(I(i),findFSIs) && ismember(I(j),findMSNs);
            subplot(numUnits,numUnits,numUnits*(i-1)+j)
            hold on
            plot(rasterVector,mean(s.CrossCorrJitter{I(i),I(j)}'),'Color',[0.7 0.7 0.7],'LineWidth',2)
            Y = prctile(s.CrossCorrJitter{I(i),I(j)},[.5 99.5],2);
            plot(rasterVector,Y(:,1),'Color',[0.7 0.7 0.7])
            plot(rasterVector,Y(:,2),'Color',[0.7 0.7 0.7])
            plot(rasterVector,smooth(s.CrossCorrRemSynch{I(i),I(j)},3),'c','LineWidth',2)
            plot(rasterVector,smooth(histStore,3),'r','LineWidth',2)
            hold on
            plot([0 0],[0 max(histStore)],'r')
            if max(histStore) < 1
                ylim([0 1])
            else
                ylim([0 max(histStore)])
            end
            xlim([-0.01 0.01])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'TickLength',[0 0])
            rasters = [];
        elseif ismember(I(i),findFSIs) && ~ismember(I(j),findMSNs);
            subplot(numUnits,numUnits,numUnits*(i-1)+j)
            hold on
%             plot(rasterVector,mean(shuffStore{I(i),I(j)}'),'Color',[0.7 0.7 0.7],'LineWidth',2)
%             plot(rasterVector,mean(shuffStore{I(i),I(j)}')+std(shuffStore{I(i),I(j)}'/sqrt(shuffNum)),'Color',[0.7 0.7 0.7])
%             plot(rasterVector,mean(shuffStore{I(i),I(j)}')-std(shuffStore{I(i),I(j)}'/sqrt(shuffNum)),'Color',[0.7 0.7 0.7])
            plot(rasterVector,smooth(histStore,3),'r','LineWidth',2)
            hold on
            plot([0 0],[0 max(histStore)],'r')
            if max(histStore) < 1
                ylim([0 1])
            else
                ylim([0 max(histStore)])
            end
            xlim([-0.01 0.01])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'TickLength',[0 0])
            rasters = [];
        else
            subplot(numUnits,numUnits,numUnits*(i-1)+j)
            plot(rasterVector,smooth(histStore,3),'k','LineWidth',2)
            hold on
            plot([0 0],[0 max(histStore)],'k')
            if max(histStore) < 1
                ylim([0 1])
            else
                ylim([0 max(histStore)])
            end
            xlim([-0.01 0.01])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'TickLength',[0 0])
            rasters = [];
        end
        
        %grey out non-identified units
        if ~ismember(I(j),[findFSIs;findMSNs])
            plot([0 0],[0 max(histStore)],'LineWidth',2,'Color',[0.7 0.7 0.7])
            plot(rasterVector,smooth(histStore,3),'LineWidth',2,'Color',[0.7 0.7 0.7])
        end
        
        %now lets try and plot out titles only for top row.
        if i == 1
            title(cluNames{I(j)})
        end
        
    end
end
spikeGraphName = strcat(fileName,'AllCrossCorr');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and separate out tone and non-tone period spikes
spikesNonTone = [];
spikesTone = [];
cluNonTone = [];
cluTone = [];
allSpikes = spTimeSec;
allClu = spClu;
counterTone = 1;
counterNonTone = 1; 
for i = 1:length(dioTimes)
    if rem(i,500) == 0 
        disp('500 trials')
    end
    %subtract tone time
    toneSub = allSpikes - dioTimes(i);
    %extract all spikes up to the end of the tone
    subSpikes = toneSub(toneSub < 0.1);
    %find negative (non tone) and positive (tone) values
    findNeg = find(subSpikes <= 0);
    findPos = find(subSpikes > 0);
    %store both spike times and cluster designations
    spikesNonTone(counterNonTone:counterNonTone - 1 + length(findNeg)) = allSpikes(findNeg);
    spikesTone(counterTone:counterTone - 1 + length(findPos)) = allSpikes(findPos);
    cluNonTone(counterNonTone:counterNonTone - 1 + length(findNeg)) = allClu(findNeg);
    cluTone(counterTone:counterTone - 1 + length(findPos)) = allClu(findPos);
    %now we need to subtract out these spikes
    allSpikes(1:length(subSpikes)) = [];
    allClu(1:length(subSpikes)) = [];
    counterTone = counterTone + length(findPos);
    counterNonTone = counterNonTone + length(findNeg);
end
spikesNonTone(counterNonTone:counterNonTone - 1 + length(allSpikes)) = allSpikes;
cluNonTone(counterNonTone:counterNonTone - 1 + length(allSpikes)) = allClu;

%now we can separate out the clusters.

% try doing this with either tone window spikes and non-tone spikes
hFig = figure;
[B,I] = sort(s.PeakChanVals);
rasterWindow = [-0.01 0.01];
rasterVector = [-0.01:0.0005:0.01];
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.02], [0.01 0.01]);
set(hFig, 'Position', [10 80 1000 1000])
for i = 1:numUnits
    %determine units greater than or equal to current unit
%     tarUnits = i:numUnits;
    tarSpikesNonTone = spikesNonTone(cluNonTone == clu(I(i)));
    tarSpikesTone = spikesTone(cluTone == clu(I(i)));
    for j = 1:numUnits
        disp('Next Unit!')
        %plot out non-tone crossCorr
        [rasters] = crosscorrelogram(tarSpikesNonTone,spikesNonTone(cluNonTone == clu(I(j))),rasterWindow);
        histStore = hist(rasters,rasterVector);
        if j == i
            histStore(21) = 0;
        end
        nonToneStore{i,j} = histStore;
        subplot(numUnits,numUnits,numUnits*(i-1)+j)
        plot(rasterVector,histStore,'k','LineWidth',2)
        max1 =max(histStore);
        hold on
        %plot out tone cross corr
        [rasters] = crosscorrelogram(tarSpikesTone,spikesTone(cluTone == clu(I(j))),rasterWindow);
        histStore = hist(rasters,rasterVector);
        ToneStore{i,j} = histStore;
        if j == i
            histStore(21) = 0;
        end
        
        plot(rasterVector,histStore,'m','LineWidth',2)
        max2 = max(histStore);
        plot([0 0],[0 max([max1,max2])],'r')
        if max([max1,max2]) < 1
            ylim([0 1])
        else
            ylim([0 max([max1,max2])])
        end
        xlim([-0.01 0.01])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'TickLength',[0 0])
        rasters = [];
    end
end
spikeGraphName = strcat(fileName,'TonevsBaseCrossCorr');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

s.CrossCorrTone = ToneStore;
s.CrossCorrNonTone = nonToneStore;

%now try and pull spikes using a shuffle in which I shift by one trial (in
%terms of the same stimulus). We need to first check and find these trials,
%and figure out how to index it. 

%first thing to do is basically streamline the dBs.
newDBs = s.SoundData.dBs;
for i = 1:s.SoundData.NumFreqs
    findTars = find(s.SoundData.Frequencies == s.SoundData.UniqueFrequencies(i));
    tarDBs = s.SoundData.dBs(findTars);
    uniqueDBVals = unique(tarDBs);
    for j = 1:length(uniqueDBVals)
        findFreq = find(s.SoundData.Frequencies == s.SoundData.UniqueFrequencies(i));
        findDBs = find(s.SoundData.dBs == uniqueDBVals(j));
        specFind = intersect(findFreq,findDBs);
        newDBs(specFind) = s.SoundData.UniqueDBs(j);
    end
end

%now lets try and find the presentations of all the unique stimuli. First,
%make a lookup table. 
oldInd = [1:length(dioTimes)];
newInd = [];
counter = 1;
for i= 1:s.SoundData.NumFreqs
    for j = 1:s.SoundData.NumDBs
        %find trials of particular frequency and db. 
        findFreq = find(s.SoundData.Frequencies == s.SoundData.UniqueFrequencies(i));
        findDBs = find(newDBs == s.SoundData.UniqueDBs(j));
        specFind = intersect(findFreq,findDBs);
        newInd(specFind,1) = counter;
        newInd(specFind,2) = [1:length(specFind)];
        counter = counter + 1;
    end
end

subInd = [];
for i = 1:counter - 1
    testFind = find(newInd(:,1) == i);
    subInd(testFind(1:end-1)) = testFind(2:end);
    subInd(testFind(end)) = testFind(1);
end
subInd = subInd';
if length(subInd) > length(dioTimes)
    subInd(length(dioTimes)+1:end) = [];
end

hFig = figure;
[B,I] = sort(s.PeakChanVals);
tarWindow = [-0.4 0.4];
rasterWindow = [-0.01 0.01];
rasterVector = [-0.01:0.0005:0.01];
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.02], [0.01 0.01]);
set(hFig, 'Position', [10 80 1000 1000])
for i = 1:numUnits
    for j = 1:numUnits
        histStore = zeros(1,length(rasterVector));
        disp('Next Unit!')
        %plot out non-tone crossCorr
        for k = 1:length(dioTimes)
            tarSpikes = s.(cluNames{I(i)}).SpikeTimes - dioTimes(k);
            tarSpikes(tarSpikes < tarWindow(1) | tarSpikes > tarWindow(2)) = [];
            tarSpikesMatch = s.(cluNames{I(j)}).SpikeTimes - dioTimes(subInd(k));
            tarSpikesMatch(tarSpikesMatch < tarWindow(1) | tarSpikesMatch > tarWindow(2)) = [];
            if length(tarSpikes)>0 & length(tarSpikesMatch) > 0
                [rasters] = crosscorrelogram(tarSpikes,tarSpikesMatch,rasterWindow);
                tempHist = hist(rasters,rasterVector);
                histStore = histStore + tempHist;
                trialScrambleStore{i,j} = histStore;
            end
            
        end

%         if j == i
%             histStore(21) = 0;
%         end
        subplot(numUnits,numUnits,numUnits*(i-1)+j)
        plot(rasterVector,histStore,'k','LineWidth',2)
        max1 =max(histStore);
        hold on
        %plot out tone cross corr
        
        plot([0 0],[0 max1],'r')
        if max1 < 1
            ylim([0 1])
        else
            ylim([0 max1])
        end
        xlim([-0.01 0.01])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'TickLength',[0 0])
        rasters = [];
    end
end
spikeGraphName = strcat(fileName,'TrialShuffCrossCorr');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

s.CrossCorrTrialShuffle = trialScrambleStore;

%% SAVE DATA


save(fullfile(pname,fname),'s','masterData','masterHeader');


%% NOW PLOT INDIVIDUAL UNITS


for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])

    %% Column 1
    %plot out median wave form
    subplot(1,4,1)
    hold on
    if s.PeakChanVals(i) <=32
        for j = 1:32
            if masterData(i,indPkTr) < s.Parameters.PVLim(1) %is FSI
                plot(plotVect,s.(cluNames{i}).medianWave(j,:)/max(max(s.(cluNames{i}).medianWave)) - j + 1,'LineWidth',2,'Color',[1 0 0]);
            elseif masterData(i,indPkTr) > s.Parameters.PVLim(2) %is MSN
                plot(plotVect,s.(cluNames{i}).medianWave(j,:)/max(max(s.(cluNames{i}).medianWave)) - j + 1,'LineWidth',2,'Color',[0 0 0]);
            else %unknown
                plot(plotVect,s.(cluNames{i}).medianWave(j,:)/max(max(s.(cluNames{i}).medianWave)) - j + 1,'LineWidth',2,'Color',[.7 .7 .7]);
            end
        end
        ylim([-33 1])
        xlim([plotVect(1) plotVect(end)])
        title(strcat('Num Spikes:',num2str(length(s.(cluNames{i}).SpikeTimes)),'Max Amp:',num2str(max(max(s.(cluNames{i}).medianWave)))))
    else
        for j = 1:32
            if masterData(i,indPkTr) < s.Parameters.PVLim(1) %is FSI
                plot(plotVect,s.(cluNames{i}).medianWave(32+j,:)/max(max(s.(cluNames{i}).medianWave)) - j + 1,'LineWidth',2,'Color',[1 0 0]);
            elseif masterData(i,indPkTr) > s.Parameters.PVLim(2) %is MSN
                plot(plotVect,s.(cluNames{i}).medianWave(32+j,:)/max(max(s.(cluNames{i}).medianWave)) - j + 1,'LineWidth',2,'Color',[0 0 0]);
            else %unknown
                plot(plotVect,s.(cluNames{i}).medianWave(32+j,:)/max(max(s.(cluNames{i}).medianWave)) - j + 1,'LineWidth',2,'Color',[.7 .7 .7]);
            end
        end
        ylim([-33 1])
        xlim([plotVect(1) plotVect(end)])
        title(strcat('Num Spikes:',num2str(length(s.(cluNames{i}).SpikeTimes)),'Max Amp:',num2str(max(max(s.(cluNames{i}).medianWave)))))        
    end
    % plot histogram.
    subplot(4,4,6)
    plot(histBinVector,s.(cluNames{i}).AllHistograms,'k','LineWidth',2)
    hold on
    plot(histBinVector,s.(cluNames{i}).AllHistograms - s.(cluNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
    plot(histBinVector,s.(cluNames{i}).AllHistograms + s.(cluNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
    %plot significant values
    plot(s.(cluNames{i}).AllHistogramSig.Centers(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,3) == 1),...
        s.(cluNames{i}).AllHistogramSig.Histogram(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,3) == 1,1),...
        'b*')
    plot(s.(cluNames{i}).AllHistogramSig.Centers(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,4) == 1),...
        s.(cluNames{i}).AllHistogramSig.Histogram(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,4) == 1,1),...
        'bo')
    %plot negative values for first tuning
    plot(s.(cluNames{i}).AllHistogramSig.Centers(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
        s.(cluNames{i}).AllHistogramSig.Histogram(...
        s.(cluNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
        'k*')
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    title('Histogram')

    %plot out rasters, organized!
    subplot(2,4,6)
    plot(s.(cluNames{i}).AllRasters(:,1),...
        s.(cluNames{i}).AllRasters(:,3),'k.','markersize',4)
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
    subplot(4,4,3)
    imagesc(s.(cluNames{i}).BinTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title({fileName;cluNames{i}},'fontweight','bold', 'Interpreter', 'none');

    %Plot binned response during general period
    subplot(4,4,7)
    imagesc(s.(cluNames{i}).BinGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Response (general)')
    %Plot peak response during tone period
    subplot(4,4,11)
    imagesc(s.(cluNames{i}).PeakMapTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Peak Response (tone)')
    %Plot peak response during general period
    subplot(4,4,15)
    imagesc(s.(cluNames{i}).PeakMapGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Peak Response (general)')

    %% Column 3
    %plot heatmap organized by frequency
    subplot(4,4,4)
    imagesc(s.(cluNames{i}).FrequencyHistograms(:,:))
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
    %plot out binned spikes (tone)
    subplot(4,4,8)
    hold on
    for cInd = 1:numDBs
        plot(s.(cluNames{i}).BinDiff(:,cInd,1),'Color',[cInd/numDBs 0 0])
        %find significant points, by p < 0.05
        findSigs = find(s.(cluNames{i}).BinSigVals(:,cInd,1)<0.05);
        plot(findSigs,s.(cluNames{i}).BinDiff(findSigs,cInd,1),'g*')
        %find significant points, by p < 0.01
        findSigs = find(s.(cluNames{i}).BinSigVals(:,cInd,1)<0.01);
        plot(findSigs,s.(cluNames{i}).BinDiff(findSigs,cInd,1),'go')
    end
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    ylabel('Binned Spikes/Fast Period')
    title('Tuning Curves Across Fast Period')

    %plot out binned spikes (tone)

    subplot(4,4,12)
    hold on
    for cInd = 1:numDBs
        plot(s.(cluNames{i}).BinDiff(:,cInd,2),'Color',[cInd/numDBs 0 0])
        %find significant points, by p < 0.05
        findSigs = find(s.(cluNames{i}).BinSigVals(:,cInd,2)<0.05);
        plot(findSigs,s.(cluNames{i}).BinDiff(findSigs,cInd,2),'g*')
        %find significant points, by p < 0.01
        findSigs = find(s.(cluNames{i}).BinSigVals(:,cInd,2)<0.01);
        plot(findSigs,s.(cluNames{i}).BinDiff(findSigs,cInd,2),'go')
    end
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    ylabel('Binned Spikes/Tone Period')
    title('Tuning Curves Across Tone Period')
    %plot out binned spikes (general)
    subplot(4,4,16)
    hold on
    for cInd = 1:numDBs
        plot(s.(cluNames{i}).BinDiff(:,cInd,3),'Color',[cInd/numDBs 0 0])
        %find significant points, by p < 0.05
        findSigs = find(s.(cluNames{i}).BinSigVals(:,cInd,3)<0.05);
        plot(findSigs,s.(cluNames{i}).BinDiff(findSigs,cInd,3),'g*')
        %find significant points, by p < 0.01
        findSigs = find(s.(cluNames{i}).BinSigVals(:,cInd,3)<0.01);
        plot(findSigs,s.(cluNames{i}).BinDiff(findSigs,cInd,3),'go')
    end
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    ylabel('Binned Spikes/Gen Period')
    title('Tuning Curves Across General Period')

    hold off
    spikeGraphName = strcat(fileName,cluNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
    close

end


