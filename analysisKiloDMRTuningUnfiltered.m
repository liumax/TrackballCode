%This code is meant to apply tuning analysis to selected kilosorted
%recordings. The goal here is to determine which units have cross
%correlogram relationships, and to examine their sorting. 



%for testing purposes
fileName = '190418ML190307G_RAudStrPen1Rec2_3546TuningDMR'
% fileName = '190123ML181105E_RAudStr3526pen1rec1tuningAndDMR';
% fileName = '190206ML181105F_RAudStr3633pen2rec1tuningDMR'
% fileName = '190205ML181105C_RAudStr3667pen2rec1tuningDMR';
% fileName = '190123_ML181105D_R_AudStr_3106_pen1_rec1_tuningAndDMR'
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
s.Parameters.PVLim = [0.00055 0.0006];
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

%     sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR400-6400-approx1kHz.mat';
%     % dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
%     ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimesCORRECTFreqSampling.mat';
    sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR400-6400-approx1kHz.mat';
    % dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
    ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimesCORRECTFreqSampling.mat';

%     ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimes.mat';
else 
%     sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR4kHz-16kHz_10mindsT30dsF5.mat';
%     % dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
%     ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-4000flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimes.mat';
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


prompt = {'What channels do you want to kill?'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'40'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
elimChans = str2num(answer{1});

% format short
disp('Parameters Set')

%% Set up file saving stuff. 

saveName = strcat(fileName,'Analysis','.mat');
fname = saveName;
pname = pwd;

%% Establishes folders
homeFolder = pwd;
% subFolders = genpath(homeFolder);
% addpath(subFolders)
% subFoldersCell = strsplit(subFolders,';')';

%% Extract files from kilosort. 

%load kilosort spikes
clustfolder ='kilosortoutput\';
spfolder = 'kilosortoutput\';
%get clusters
[cids, cgs] = readClusterGroupsCSV([clustfolder 'cluster_group.tsv']);
disp('Clusters Pulled')
phyGoodClusters = cids(cgs==2);
goodClu=phyGoodClusters;
numClu = numel(goodClu);

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

%190518 We need to do cleanup of data, not all ISIs are great. Lets remove
%anything with greater 


%This outputs a bunch of spikes and their designations. we need to
%appropriately sort these now. First lets generate the right names
cluNames = [];
numUnits = length(goodClu);
for i = 1:numUnits
    cluNames{i} = strcat('clu',num2str(goodClu(i)));
end
s.DesignationName = cluNames;
for i = 1:numUnits
    tarSpikes = find(spClu == goodClu(i));
    s.(cluNames{i}).SpikeTimes = spTimeSec(tarSpikes);
end

%% AT THIS TIME INTERRUPT TO GET DIO TIMES: THIS WILL ALLOW FOR SEPARATION OF DMR VS NON-DMR TIMES
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

%% Now return to waveform extraction. We can now pull first and last DIO times out too. 

dioLimiter = [dmrDIO(1) dmrDIO(end)];
spikeLimiter = [find(spTimeSec < dioLimiter(1),1,'last'),find(spTimeSec > dioLimiter(2),1,'first')];
%determine peak channel for given cluster. remember that cluster names are
%zero indexed!!

s.PeakChanVals = cluPkChan(goodClu+1);

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


%now lets extract the biggest individual template per cluster. This avoids
%issues that can arise from the averaging of multiple templates with
%different timing. 
for i = 1:length(goodClu)
    %determine templates associated with goodClu
    pullTemp = spTemp(spClu ==goodClu(i));
    %save number of merges here!
    numMerge = length(unique(pullTemp));
    mergeStore(i) = numMerge;
    
    %if we have more than one template, want to isolate to the correct one!
    if numMerge > 1
        [C ia ic] = unique(pullTemp);
        tempSpikeNums = hist(ic,[1:length(ia)]);
        [maxVal maxFind] = max(tempSpikeNums);
        tempStoreForWave(i) = C(maxFind);
    else
        tempStoreForWave(i) = unique(pullTemp);
    end
end

tempStoreForWave = double(tempStoreForWave);

%% Now pull from non-DMR periods!
%190514 Encountering error: if split clusters, then they are represented as
%multiple templates. Therefore, I think the thing we need to do is to
%basically go through and find the ones which are duplicates and remove
%them. Then go through and use templates for the ones that I can, and do a
%separate extraction for the others using their cluster identity. 
[n, bin] = histc(tempStoreForWave, unique(tempStoreForWave));
multiple = find(n > 1);
index    = find(ismember(bin, multiple));
revInd = [1:length(goodClu)];
revInd(index) = [];

tempStoreForWave(index) = [];
[B I] = sort(tempStoreForWave);
[B2 I2] = sort(I);
%extract only the clusters I want
truncSpTemp = spTemp;
truncSpTemp(spikeLimiter(1):spikeLimiter(2)) = [];
truncSpTime = spTime;
truncSpTime(spikeLimiter(1):spikeLimiter(2)) = [];
truncSpClu = spClu;
truncSpClu(spikeLimiter(1):spikeLimiter(2)) = [];
tester = ismember(truncSpTemp,tempStoreForWave);
newClu = truncSpTemp(tester);
newST = truncSpTime(tester);

tarPath = strcat(pwd,'\',fileName,'.phy.dat');
medWFsTemp = extractMedianWFs(newClu, newST, 30000,tarPath, 'int16', [64,nSamp], chanMap, 1);

medWFs(revInd,:,:) = medWFsTemp(I2,:,:);

%now we need to find the correct template numbers for the remaining
%values.While it would be ideal to find the next best cluster for merges,
%its too complicated. Instead, lets simply do those by their cluster
%values. 
remClu = goodClu(index);
tester2 = ismember(truncSpClu,remClu);
newClu = truncSpClu(tester2);
newST = truncSpTime(tester2);

medWFsSec = extractMedianWFs(newClu, newST, 30000,tarPath, 'int16', [64,nSamp], chanMap, 1);

medWFs(index,:,:) = medWFsSec;

medWFs(:,elimChans,:) = zeros;
% tarWave(elimChans,:) = zeros;

%% NOW DO DMR TIMES

%extract only the clusters I want
truncSpTemp = spTemp(spikeLimiter(1):spikeLimiter(2));
truncSpTime = spTime(spikeLimiter(1):spikeLimiter(2));
truncSpClu = spClu(spikeLimiter(1):spikeLimiter(2));
tester = ismember(truncSpTemp,tempStoreForWave);
newClu = truncSpTemp(tester);
newST = truncSpTime(tester);

tarPath = strcat(pwd,'\',fileName,'.phy.dat');
medWFsTemp = extractMedianWFs(newClu, newST, 30000,tarPath, 'int16', [64,nSamp], chanMap, 1);

medWFsDMR(revInd,:,:) = medWFsTemp(I2,:,:);

%now we need to find the correct template numbers for the remaining
%values.While it would be ideal to find the next best cluster for merges,
%its too complicated. Instead, lets simply do those by their cluster
%values. 
remClu = goodClu(index);
tester2 = ismember(truncSpClu,remClu);
newClu = truncSpClu(tester2);
newST = truncSpTime(tester2);

medWFsSec = extractMedianWFs(newClu, newST, 30000,tarPath, 'int16', [64,nSamp], chanMap, 1);

medWFsDMR(index,:,:) = medWFsSec;
medWFsDMR(:,elimChans,:) = zeros;

%% store the waveforms! also extract characteristics
for i = 1:numUnits
    tarWave = squeeze(medWFs(i,:,:));
    s.(cluNames{i}).medianWave = tarWave;
    tarWave2 = squeeze(medWFsDMR(i,:,:));
    s.(cluNames{i}).medianWaveDMR = tarWave2;
    %find minima for each wave. Note that waves are now downwards pointing
    for j = 1:64
        tarWave(j,:) = tarWave(j,:) - mean(tarWave(j,1:5));
    end
%     figure
%     subplot(2,1,1)
%     plot(tarWave')
%     subplot(2,1,2)
%     plot(tarWave)
%     tarWave(elimChans,:) = zeros; %190517 code to fix dead channel issues. 
%     testMin(i,:) = mean(tarWave(:,1:10));
    [pks inds] = min(tarWave(:,10:20)');
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
[sta_sig, ptd, siglevel] = ne_sig_sta_from_stim_obs_resp(sta, spikeArray, stimulus, 10, numLags, 95);

s.STAs = sta;
s.STASig = sta_sig;


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

s.RealCorrStore = realCorrStore;

s.DMRTimes = trueDMRtimes;
s.DMRStep = dmrStep;
s.DMRfaxis = faxis;

s.SpikeArray = spikeArray;
%now lets remove other stuff
rez = [];
% spikeArray = [];
% stimulus = [];
indivOut = [];



%% Plotting waveforms

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
    findTars = find(s.PeakChanVals <= 32*i & s.PeakChanVals >= 32*(i-1)+1);
    if findTars
        hFig = figure;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.07], [0.01 0.01]);
        set(hFig, 'Position', [10 10 1900 1000])
        for j = 1:length(findTars)
            subplot(1,length(findTars),j)
            hold on
            for k = 1:32
                if masterData(findTars(j),indPkTr) < s.Parameters.PVLim(1) %is FSI
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/min(min(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[1 0 0]);
                elseif masterData(findTars(j),indPkTr) > s.Parameters.PVLim(2) %is MSN
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/min(min(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[0 0 0]);
                else %unknown
                    plot(plotVect,s.(cluNames{findTars(j)}).medianWave(32*(i-1)+k,:)/min(min(s.(cluNames{findTars(j)}).medianWave)) - k + 1,'LineWidth',2,'Color',[.7 .7 .7]);
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

%% Calculating and plotting cross correlograms
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

% s.JitteredSpikes = shuffSpikes;

%now lets actually do the targeted cross correlograms! Try to only target FSIs and
%align MSN spikes to them. 
findFSIs =find(masterData(:,indPkTr) < s.Parameters.PVLim(1));
findMSNs = find(masterData(:,indPkTr) > s.Parameters.PVLim(2));

synchWind = 0.001; %time difference within which spikes are considered synchronous.
rasterWindow = [-0.01025 0.01025];
rasterVector = [-0.01:0.0005:0.01];

corrData = [];
%First, do cross correlograms, split by shanks
for i = 1:2
    findTars = find(s.PeakChanVals <= 32*i & s.PeakChanVals >= 32*(i-1)+1);
    if findTars
        disp(strcat('Performing Cross Corr for Shank-',num2str(i)))
        corrName = strcat('trueStoreShank',num2str(i));
        corrData.(corrName) = [];
        for k = 1:length(findTars)
            for j = 1:length(findTars)
                disp(strcat('Running xCorr Shank-',num2str(i),'Unit',num2str(cluNames{findTars(j)}),'-aligned to',num2str(cluNames{findTars(k)})))
                [rasters] = crosscorrelogram(s.(cluNames{findTars(k)}).SpikeTimes,s.(cluNames{findTars(j)}).SpikeTimes,rasterWindow);
                histStore = hist(rasters,rasterVector);
                rasters = [];
                corrData.(corrName)(k,j,:) = histStore;
            end
        end
    end
end

%make wide cross correlograms
rasterWindow = [-0.1025 0.1025];
rasterVector = [-0.1:0.005:0.1];
for i = 1:2
    findTars = find(s.PeakChanVals <= 32*i & s.PeakChanVals >= 32*(i-1)+1);
    if findTars
        disp(strcat('Performing Cross Corr for Shank-',num2str(i)))
        corrName = strcat('wideStoreShank',num2str(i));
        corrData.(corrName) = [];
        for k = 1:length(findTars)
            for j = 1:length(findTars)
                disp(strcat('Running xCorr Shank-',num2str(i),'Unit',num2str(cluNames{findTars(j)}),'-aligned to',num2str(cluNames{findTars(k)})))
                [rasters] = crosscorrelogram(s.(cluNames{findTars(k)}).SpikeTimes,s.(cluNames{findTars(j)}).SpikeTimes,rasterWindow);
                histStore = hist(rasters,rasterVector);
                rasters = [];
                corrData.(corrName)(k,j,:) = histStore;
            end
        end
    end
end

%now we want to go through and remove synchronous spikes and generate
%jittered xcorr
rasterWindow = [-0.01025 0.01025];
rasterVector = [-0.01:0.0005:0.01];
prctileBounds = [.5 99.5];
for i = 1:2
    findTars = find(s.PeakChanVals <= 32*i & s.PeakChanVals >= 32*(i-1)+1);
    if findTars
        disp(strcat('Performing Jittered xCorr for Shank-',num2str(i)))
        trueName = strcat('trueStoreShank',num2str(i));
        remName = strcat('RemStore',num2str(i));
        corrData.(remName) = [];
        shuffName = strcat('ShuffStore',num2str(i));
        corrData.(shuffName) = [];
        sigName = strcat('SigCross',num2str(i));
        corrData.(sigName) = [];
        for k = 1:length(findTars)
            for j = 1:length(findTars)
                if k ~= j
                    %now lets try and remove synchronous FSIs spikes. This way, we can
                    %eliminate any peak in the FSI aligned MSN spikes. 
                    diffArray = [];
                    subSpikes = s.(cluNames{findTars((k))}).SpikeTimes;
                    for n = length(s.(cluNames{findTars((j))}).SpikeTimes):-1:1
                        diffArray = subSpikes - s.(cluNames{findTars((j))}).SpikeTimes(n);
                        spikeFind = find(abs(diffArray) < synchWind);
                        subSpikes(spikeFind) = [];
                    end
                    [rasters] = crosscorrelogram(subSpikes,s.(cluNames{findTars((j))}).SpikeTimes,rasterWindow);
                    histStore = hist(rasters,rasterVector);
                    rasters = [];
                    corrData.(remName)((k),(j),:) = histStore;
                
                    %now lets run the jittered spikes
                    tempHist = [];
                    shuffTar = shuffSpikes{findTars((j))};
                    disp(strcat('Running Shuffle Cross Corr for Unit',num2str(cluNames{findTars(j)}),'-aligned to',num2str(cluNames{findTars((k))})))
                    for m = 1:shuffNum
    %                     disp(m)
                        [rasters] = crosscorrelogram(s.(cluNames{findTars((k))}).SpikeTimes,shuffTar(:,m),rasterWindow);
                        rasters(rasters < rasterWindow(1) | rasters > rasterWindow(2)) = [];
                        tempHist(:,m) = hist(rasters,rasterVector);
                    end
                    corrData.(shuffName){(k),(j)} = tempHist;
                    %now determine points of significance (exceeding percentile
                    %bounds)
                    Y = prctile(tempHist,prctileBounds,2);
                    trueVals = squeeze(corrData.(trueName)(k,j,:));
%                     figure
%                     hold on
%                     plot(Y(:,1))
%                     plot(Y(:,2))
%                     plot(trueVals,'r')
                    %now lets make sure anything lower gets picked up
                    sigVals = zeros(length(trueVals),1);
                    sigVals(trueVals - Y(:,1) < 0) = -1;
                    %now do higher
                    sigVals(trueVals - Y(:,2) > 0) = 1;
                    %now we want to clean this up. Remove all single change
                    %stretches, only keep prolonged changes. 
                    whileTrig = 0;
                    creepCount = 1;
                    minWidth = 2;
                    while whileTrig == 0
                        if creepCount > length(sigVals)
                            break
                        end
                        %check if sigVals == 0
                        if sigVals(creepCount) == 0
                            creepCount = creepCount + 1;
                        else
                            tmpVal = sigVals(creepCount);
                            storeCount = creepCount;
                            while whileTrig == 0
                                if creepCount > length(sigVals)
                                    break
                                end
                                if tmpVal == sigVals(creepCount) %this is still part of a continuous grouping
                                    creepCount = creepCount + 1;
                                elseif tmpVal ~= sigVals(creepCount) %no longer part of continuous block of numbers
                                    %determine length
                                    tmpLength = creepCount - storeCount;
                                    if tmpLength >= minWidth;
                                        disp('Continuous Sig Detected')
                                        break
                                    else
                                        disp('Insufficient Length, deleting...')
                                        sigVals(storeCount:creepCount - 1) = 0;
                                        break
                                    end
                                end
                            end
                        end
                    end
                    newStore{(k),(j)} = sigVals;
                    sigWarnPos((k),(j)) = length(find(sigVals == 1));
                    sigWarnNeg((k),(j)) = length(find(sigVals == -1));
                    
                end
            end
        end
        corrData.(sigName).SigVals = newStore;
        corrData.(sigName).PosWarn = sigWarnPos;
        corrData.(sigName).NegWarn = sigWarnNeg;
        newStore= [];
        sigWarnPos= [];
        sigWarnNeg= [];
    end
end


%plot out cross correlograms, zero eliminated cross corr, and jittered
%baseline with 99%ile bounds
for bigInd = 1:2
    findTars = find(s.PeakChanVals <= 32*bigInd & s.PeakChanVals >= 32*(bigInd-1)+1);
    if findTars
        hFig = figure;
        testData = s.PeakChanVals(findTars);
        [B,I] = sort(testData);
        rasterWindow = [-0.01025 0.01025];
        rasterVector = [-0.01:0.0005:0.01];
        subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.001], [0.02 0.02], [0.01 0.01]);
        set(hFig, 'Position', [10 10 1800 1000])
        corrName = strcat('trueStoreShank',num2str(bigInd));
        remName = strcat('RemStore',num2str(bigInd));
        shuffName = strcat('ShuffStore',num2str(bigInd));
        sigName = strcat('SigCross',num2str(bigInd));
        for i = 1:length(findTars)
            disp(strcat('unit',cluNames{findTars(I(i))}))
            for j = 1:length(findTars)
                histStore = squeeze(corrData.(corrName)(I(i),I(j),:));
                if j == i
                    histStore(21) = 0;
                end
                if ismember(findTars(I(i)),intersect(findFSIs,findTars)) && ismember(findTars(I(j)),intersect(findMSNs,findTars));
                    subplot(length(findTars),length(findTars),length(findTars)*(i-1)+j)
                    hold on
                    plot(rasterVector,mean(corrData.(shuffName){(I(i)),(I(j))}'),'Color',[0.7 0.7 0.7],'LineWidth',2)
                    Y = prctile(corrData.(shuffName){(I(i)),(I(j))},prctileBounds,2);
                    plot(rasterVector,Y(:,1),'Color',[0.7 0.7 0.7])
                    plot(rasterVector,Y(:,2),'Color',[0.7 0.7 0.7])
                    plot(rasterVector,smooth(corrData.(remName)((I(i)),(I(j)),:),3),'c','LineWidth',2)
                    plot(rasterVector,smooth(histStore,3),'r','LineWidth',2)
                    plot([0 0],[0 max(histStore)],'g')
                    if max(histStore) < 1
                        ylim([0 1])
                    else
                        ylim([0 max(histStore)])
                    end
                    if corrData.(sigName).PosWarn(I(i),I(j)) ~= 0 %if there are significantly positive xcorr points
                        %find the points
                        finder = find(corrData.(sigName).SigVals{I(i),I(j)} == 1);
                        smoothDat = smooth(histStore,3);
                        plot(rasterVector(finder),smoothDat(finder),'m*')
                    end
                    if corrData.(sigName).NegWarn(I(i),I(j)) ~= 0 %if there are significantly negative xcorr points
                        %find the points
                        finder = find(corrData.(sigName).SigVals{I(i),I(j)} == -1);
                        smoothDat = smooth(histStore,3);
                        plot(rasterVector(finder),smoothDat(finder),'g*')
                    end
                    
                    xlim([-0.01 0.01])
                    set(gca,'xtick',[])
                    set(gca,'ytick',[])
                    set(gca,'TickLength',[0 0])
                    rasters = [];
                elseif ismember(findTars(I(i)),intersect(findFSIs,findTars)) && ~ismember(findTars(I(j)),intersect(findMSNs,findTars));
                    subplot(length(findTars),length(findTars),length(findTars)*(i-1)+j)
                    hold on
                    plot(rasterVector,smooth(histStore,3),'r','LineWidth',2)
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
                    subplot(length(findTars),length(findTars),length(findTars)*(i-1)+j)
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
                if ~ismember(findTars(I(j)),[intersect(findFSIs,findTars);intersect(findMSNs,findTars)])
                    plot([0 0],[0 max(histStore)],'LineWidth',2,'Color',[0.7 0.7 0.7])
                    plot(rasterVector,smooth(histStore,3),'LineWidth',2,'Color',[0.7 0.7 0.7])
                end

                %now lets try and plot out titles only for top row.
                if i == 1
                    title(cluNames{findTars(I(j))})
                end

            end
        end
        
        spikeGraphName = strcat(fileName,'Shank-',num2str(bigInd),'CrossCorr');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end



%now lets try and separate out tone and non-tone period spikes
spikesNonTone = [];
spikesTone = [];
cluNonTone = [];
cluTone = [];
allSpikes = spTimeSec;
allClu = spClu;
counterTone = 1;
counterNonTone = 1; 
try
    toneTimes = s.TrialMatrix(:,1);
catch
    toneTimes = dioTimes;
end
for i = 1:length(toneTimes)
    if rem(i,500) == 0 
        disp('500 trials')
    end
    %subtract tone time
    toneSub = allSpikes - toneTimes(i);
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

%now lets do the cross correlograms
for i = 1:2
    findTars = find(s.PeakChanVals <= 32*i & s.PeakChanVals >= 32*(i-1)+1);
    if findTars
        disp(strcat('Performing Tone/NoneTone Cross Corr for Shank-',num2str(i)))
        nameNonTone = strcat('NonTone',num2str(i));
        nameTone = strcat('Tone',num2str(i));
        for k = 1:length(findTars)
            tarSpikesNonTone = spikesNonTone(cluNonTone == goodClu(findTars(k)));
            tarSpikesTone = spikesTone(cluTone == goodClu(findTars(k)));
            for j = 1:length(findTars)
                disp(strcat('Running xCorr Shank-',num2str(i),'Unit',num2str(cluNames{findTars(j)}),'-aligned to',num2str(cluNames{findTars(k)})))
                secSpikesNonTone = spikesNonTone(cluNonTone == goodClu(findTars(j)));
                secSpikesTone = spikesTone(cluTone == goodClu(findTars(j)));
                
                [rasters] = crosscorrelogram(tarSpikesNonTone,secSpikesNonTone,rasterWindow);
                histStore = hist(rasters,rasterVector);
                if j == k
                    histStore(21) = 0;
                end
                corrData.(nameNonTone)(k,j,:) = histStore;
                
                [rasters] = crosscorrelogram(tarSpikesTone,secSpikesTone,rasterWindow);
                histStore = hist(rasters,rasterVector);
                if j == k
                    histStore(21) = 0;
                end
                corrData.(nameTone)(k,j,:) = histStore;
            end
        end
    end
end



for bigInd = 1:2
    findTars = find(s.PeakChanVals <= 32*bigInd & s.PeakChanVals >= 32*(bigInd-1)+1);
    nameNonTone = strcat('NonTone',num2str(bigInd));
    nameTone = strcat('Tone',num2str(bigInd));
    if findTars
        hFig = figure;
        testData = s.PeakChanVals(findTars);
        [B,I] = sort(testData);
        rasterWindow = [-0.01025 0.01025];
        rasterVector = [-0.01:0.0005:0.01];
        subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.001], [0.02 0.02], [0.01 0.01]);
        set(hFig, 'Position', [10 10 1800 1000])
        nameTone = strcat('Tone',num2str(bigInd));
        nameNonTone = strcat('NonTone',num2str(bigInd));
        for i = 1:length(findTars)
            disp(strcat('unit',cluNames{findTars(I(i))}))
            for j = 1:length(findTars)
                subplot(length(findTars),length(findTars),length(findTars)*(i-1)+j)
                hold on
                plot(rasterVector,smooth(corrData.(nameNonTone)(I(i),I(j),:),3),'k','LineWidth',2)
                plot(rasterVector,smooth(corrData.(nameTone)(I(i),I(j),:),3),'m','LineWidth',2)
                
                %find max
                maxNonTone = max(smooth(corrData.(nameNonTone)(I(i),I(j),:),3));
                maxTone = max(smooth(corrData.(nameTone)(I(i),I(j),:),3));
                maxMax = max([maxNonTone,maxTone]);
                
                if max(maxMax) < 1
                    ylim([0 1])
                else
                    ylim([0 maxMax])
                end
                xlim([-0.01 0.01])
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                set(gca,'TickLength',[0 0])

                %grey out non-identified units
                if ~ismember(findTars(I(j)),[intersect(findFSIs,findTars);intersect(findMSNs,findTars)])
                    plot(rasterVector,smooth(corrData.(nameNonTone)(I(i),I(j),:),3),'LineWidth',2,'Color',[0.7 0.7 0.7])
                    plot(rasterVector,smooth(corrData.(nameTone)(I(i),I(j),:),3),'LineWidth',2,'Color',[0.7 0.7 0.7])
                end

                %now lets try and plot out titles only for top row.
                if i == 1
                    title(cluNames{I(j)})
                end

            end
        end
        
        spikeGraphName = strcat(fileName,'Shank-',num2str(bigInd),'ToneCrossCorr');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end

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
% oldInd = [1:length(dioTimes)];
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
%to prune out if recording was cut short. 
try
    targetDIOs = s.SoundData.ToneTimes;
catch
    targetDIOs = s.TrialMatrix(:,1);
end
subInd = subInd';
if length(subInd) > length(targetDIOs)
    subInd(length(targetDIOs)+1:end) = [];
end
tarWindow = [-0.2 0.4];
for bigInd = 1:2
    findTars = find(s.PeakChanVals <= 32*bigInd & s.PeakChanVals >= 32*(bigInd-1)+1);
    if findTars
        hFig = figure;
        testData = s.PeakChanVals(findTars);
        [B,I] = sort(testData);
        rasterWindow = [-0.01025 0.01025];
        rasterVector = [-0.01:0.0005:0.01];
        subplot = @(m,n,p) subtightplot (m, n, p, [0.001 0.001], [0.02 0.02], [0.01 0.01]);
        set(hFig, 'Position', [10 10 1800 1000])
        nameShank = strcat('trialShuff',num2str(bigInd));
        for i = 1:length(findTars)
            disp(strcat('unit',cluNames{findTars(I(i))}))
            for j = 1:length(findTars)
                histStore = zeros(1,length(rasterVector));
                for k = 1:length(subInd)
                    tarSpikes = s.(cluNames{I(i)}).SpikeTimes - targetDIOs(k);
                    tarSpikes(tarSpikes < tarWindow(1) | tarSpikes > tarWindow(2)) = [];
                    tarSpikesMatch = s.(cluNames{I(j)}).SpikeTimes - targetDIOs(subInd(k));
                    tarSpikesMatch(tarSpikesMatch < tarWindow(1) | tarSpikesMatch > tarWindow(2)) = [];
                    if length(tarSpikes)>0 & length(tarSpikesMatch) > 0
                        [rasters] = crosscorrelogram(tarSpikes,tarSpikesMatch,rasterWindow);
                        tempHist = hist(rasters,rasterVector);
                        histStore = histStore + tempHist;
                        corrData.(nameShank){i,j} = histStore;
                    end
                end
                
                subplot(length(findTars),length(findTars),length(findTars)*(i-1)+j)
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
                
                %grey out non-identified units
                if ~ismember(findTars(I(j)),[intersect(findFSIs,findTars);intersect(findMSNs,findTars)])
                    plot(rasterVector,histStore,'LineWidth',2,'Color',[0.7 0.7 0.7])
                end

                %now lets try and plot out titles only for top row.
                if i == 1
                    title(cluNames{I(j)})
                end

            end
        end
        
        spikeGraphName = strcat(fileName,'Shank-',num2str(bigInd),'TrialShuffCrossCorr');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end



%% Saving
save(fullfile(pname,fname),'s','masterData','masterHeader','corrData','-v7.3');
% save(fullfile(pname,strcat(fileName,'DMRData')),'sta','sta_sig','spikeArray','stimulus','faxis')
save(fullfile(pname,strcat(fileName,'DMRData')),'sta','spikeArray','stimulus','faxis','spikeHistVect')

 
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
    subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.02 0.02], [0.01 0.01]);
    % Column 1
    %plots average waveform
    subplot(1,8,1)
    hold on
    if s.PeakChanVals(i) <=32
        for j = 1:32
            plot(plotVect,s.(cluNames{i}).medianWaveDMR(j,:)/min(min(s.(cluNames{i}).medianWave))*2 - j + 1,'LineWidth',2,'Color',[0 1 0]);
            if masterData(i,indPkTr) < s.Parameters.PVLim(1) %is FSI
                plot(plotVect,s.(cluNames{i}).medianWave(j,:)/min(min(s.(cluNames{i}).medianWave))*2 - j + 1,'LineWidth',2,'Color',[1 0 0]);
            elseif masterData(i,indPkTr) > s.Parameters.PVLim(2) %is MSN
                plot(plotVect,s.(cluNames{i}).medianWave(j,:)/min(min(s.(cluNames{i}).medianWave))*2 - j + 1,'LineWidth',2,'Color',[0 0 0]);
            else %unknown
                plot(plotVect,s.(cluNames{i}).medianWave(j,:)/min(min(s.(cluNames{i}).medianWave))*2 - j + 1,'LineWidth',2,'Color',[.7 .7 .7]);
            end
        end
        ylim([-33 2])
        xlim([plotVect(1) plotVect(end)])
        title(strcat('Shank1-AvRate:',num2str(s.(cluNames{i}).AverageRate),'Tuning',num2str(numTuningSpikes(i)),'DMR',num2str(numDMRSpikes(i))))
    else
        for j = 1:32
            plot(plotVect,s.(cluNames{i}).medianWaveDMR(32+j,:)/min(min(s.(cluNames{i}).medianWave))*2 - j + 1,'LineWidth',2,'Color',[0 1 0]);
            if masterData(i,indPkTr) < s.Parameters.PVLim(1) %is FSI
                plot(plotVect,s.(cluNames{i}).medianWave(32+j,:)/min(min(s.(cluNames{i}).medianWave))*2 - j + 1,'LineWidth',2,'Color',[1 0 0]);
            elseif masterData(i,indPkTr) > s.Parameters.PVLim(2) %is MSN
                plot(plotVect,s.(cluNames{i}).medianWave(32+j,:)/min(min(s.(cluNames{i}).medianWave))*2 - j + 1,'LineWidth',2,'Color',[0 0 0]);
            else %unknown
                plot(plotVect,s.(cluNames{i}).medianWave(32+j,:)/min(min(s.(cluNames{i}).medianWave))*2 - j + 1,'LineWidth',2,'Color',[.7 .7 .7]);
            end
        end
        ylim([-33 2])
        xlim([plotVect(1) plotVect(end)])
        title(strcat('Shank2-AvRate:',num2str(s.(cluNames{i}).AverageRate),'Tuning',num2str(numTuningSpikes(i)),'DMR',num2str(numDMRSpikes(i))))
    end
    %next plot ISIs
    subplot(4,8,2)
    smallVect = [0:0.0005:0.025];
    tmpISI = hist(diff(s.(cluNames{i}).SpikeTimes),smallVect);
    bar(smallVect,tmpISI)
    xlim([smallVect(1) smallVect(end)])
    ylim([0 max(tmpISI(1:end-1))])
    rpvs = sum(tmpISI(1:5));
    title(['RPVs:',num2str(rpvs),'/',num2str(length(s.(cluNames{i}).SpikeTimes)),'aka',num2str(rpvs/length(s.(cluNames{i}).SpikeTimes)*100),'%'])
    
    subplot(4,8,10)
    bigVect = [0:0.001:0.1];
    tmpISI = hist(diff(s.(cluNames{i}).SpikeTimes),bigVect);
    bar(bigVect,tmpISI)
    xlim([bigVect(1) bigVect(end)])
    ylim([0 max(tmpISI(1:end-1))])
    
    
    %Column 2
    %plot FR and velocity
    subplot(4,4,2)
    hold on
    plot(newTimeVector,newVelVector/max(newVelVector),'b')
    plot([newTimeVector(1):s.Parameters.SpeedFiringBins:newTimeVector(end)],s.(cluNames{i}).SessionFiring/max(s.(cluNames{i}).SessionFiring),'r')
    xlim([newTimeVector(1),newTimeVector(end)])
    ylim([-0.1,1])
    title({fileName;cluNames{i}},'fontweight','bold', 'Interpreter', 'none');
    
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
    histLims = max(s.(cluNames{i}).AllHistograms + s.(cluNames{i}).HistogramStandardDeviation);
    ylim([0 histLims])
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

    % Column 3
    %plot heatmap tuning curves
    subplot(2,4,3)
    imagesc(s.(cluNames{i}).BinTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Resp (tone)')
    
    subplot(2,4,7)
    imagesc(s.(cluNames{i}).BinGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Resp (general)')
    
    
    % Column 4
    
    %plot DMR response
    subplot(2,4,4)
    imagesc(reshape(s.STAs(i,:),length(faxis),[]))
    colormap('parula')
    colorbar
    format short
    set(gca,'XTick',[0:20:numLags]);
    set(gca,'XTickLabel',[dmrStep*100:-20*dmrStep:dmrStep]);
    set(gca,'YTick',[1:10:length(faxis)]);
    set(gca,'YTickLabel',[faxis([1:10:end])]);
    title(strcat('DMR STA Based on',num2str(numDMRSpikes(i))))
    
    %plot out thresholded DMR
    subplot(2,4,8)
    imagesc(reshape(sta_sig(i,:),length(faxis),[]))
    colormap('parula')
    colorbar
    format short
    set(gca,'XTick',[0:20:numLags]);
    set(gca,'XTickLabel',[dmrStep*100:-20*dmrStep:dmrStep]);
    set(gca,'YTick',[1:10:length(faxis)]);
    set(gca,'YTickLabel',[faxis([1:10:end])]);
    title('Thresholded DMR a=0.95')
    
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



