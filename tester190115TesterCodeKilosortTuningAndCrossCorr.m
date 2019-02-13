%This code is meant to apply tuning analysis to selected kilosorted
%recordings. The goal here is to determine which units have cross
%correlogram relationships, and to examine their sorting. 



%for testing purposes
% fileName = '180718_ML180619B_L_AudStr_pen1_3000_fullTuning'
% fileName = '180718_ML180619C_R_AudStr_pen1_3000_fullTuning'
fileName = '180315_ML180306C_R17_3218mid1_fullTuning'

%% Constants

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

%% Set up file saving stuff. 

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
masterHeader = cell(10,1);
masterInd = 1;

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
dioTimes = dioTimes - timeFirst;

s.SoundData.ToneTimes = dioTimes;
s.SoundData.ToneTimeDiff = dioTimeDiff;



%pull D2 info (laser ID)
targetFolderSearch = dir;%pulls dir from DIO folder
targetFileNames = {targetFolderSearch.name}';%pulls names section
targetFileFinder = strfind(targetFileNames,'D2'); %examines names for D1
targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name

if length(targetFiles) == 0
    targetFolderSearch = dir;%pulls dir from DIO folder
    targetFileNames = {targetFolderSearch.name}';%pulls names section
    targetFileFinder = strfind(targetFileNames,'Din2'); %examines names for D1
    targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
    targetFiles = {targetFileNames{targetFileFinder}};%pulls out actual file name
end
D2FileName = targetFiles{1};

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

timeFinder = find(master(:,1) > spTimeSec(end));
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
%     
%     %now pull waveform data and isi data to categorize cells
%     waveForms = s.(cluNames{i}).AverageWaveForms;
%     
%     %find size of waveforms
%     waveSize = size(waveForms);
% 
%     %find maximum waveform by max peak size
%     maxWave = max(waveForms);
%     [maxVal maxInd] = max(maxWave);
% 
%     %chose the big wave, interpolate to fine degree
%     chosenWave = waveForms(:,maxInd);
%     interpVect = [1:0.1:40];
%     interpWave = interp1(1:40,chosenWave,interpVect,'spline');
%     
%     %now we need to find the peak. Find this starting at point 10. 
% 
%     [pkVal pkInd] = max(interpWave(100:end));
%     pkInd = pkInd + 100 - 1;
%     %now we need to find the minimum following the peak
% 
%     [troughVal troughInd] = min(interpWave(pkInd:end));
%     troughInd = troughInd + pkInd - 1;
% 
%     peakTrough = (troughInd - pkInd)/300000;
%     
%     %find ISIs
%     isiTimes = diff(spikeTimes);
%     isiCov = std(isiTimes)/mean(isiTimes);
%     
%     masterData(i,masterHolder) = peakTrough;
%     masterHeader{masterHolder} = 'PeakTrough(ms)';
%     masterHolder = masterHolder + 1;
%     
%     masterData(i,masterHolder) = isiCov;
%     masterHeader{masterHolder} = 'isiCov';
%     masterHolder = masterHolder + 1;
%     
%     
%     if peakTrough < s.Parameters.PVLim(1) & isiCov > s.Parameters.ChatLim
%         masterData(i,masterHolder) = 1; %pv cell
%     elseif peakTrough > s.Parameters.PVLim(2) & isiCov < s.Parameters.ChatLim & averageRate > 2
%         masterData(i,masterHolder) = 2; %ChAT Cell
%     elseif peakTrough > s.Parameters.PVLim(2) & isiCov > s.Parameters.ChatLim
%         masterData(i,masterHolder) = 0; %MSN
%     else
%         masterData(i,masterHolder) = NaN; %label as unknown
%     end
%     masterHeader{masterHolder} = 'CellType';
%     masterHolder = masterHolder + 1;
    
    
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

save(fullfile(pname,fname),'s','masterData','masterHeader');


%% Plotting

%first do cross correlogram
hFig = figure;
[B,I] = sort(s.PeakChanVals);
rasterWindow = [-0.01 0.01];
rasterVector = [-0.01:0.0005:0.01];
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.02 0.02], [0.01 0.01]);
set(hFig, 'Position', [10 80 1000 1000])
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
        subplot(numUnits,numUnits,numUnits*(i-1)+j)
        plot(rasterVector,smooth(histStore,3),'k','LineWidth',2)
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
        subplot(numUnits,numUnits,numUnits*(i-1)+j)
        plot(rasterVector,histStore,'k','LineWidth',2)
        max1 =max(histStore);
        hold on
        %plot out tone cross corr
        [rasters] = crosscorrelogram(tarSpikesTone,spikesTone(cluTone == clu(I(j))),rasterWindow);
        histStore = hist(rasters,rasterVector);
        if j == i
            histStore(21) = 0;
        end
        
        plot(rasterVector,histStore,'r','LineWidth',2)
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



% 
% 
% for i = 1:numUnits
%     hFig = figure;
%     set(hFig, 'Position', [10 80 1900 1000])
%     plot(s.(cluNames{i}).AllRasters(:,1),...
%         s.(cluNames{i}).AllRasters(:,3),'k.','markersize',4)
%     hold on
%     plot([0 0],[ylim],'b');
%     plot([toneDur toneDur],[ylim],'b');
%     rasterFreqLines = zeros(numFreqs,2);
%     if numDBs ==1
%         rasterFreqLines(:,1) = cumsum((matrixTrialNum'));
%     elseif numDBs > 1
%         rasterFreqLines(:,1) = cumsum(sum(matrixTrialNum'));
%     end
% 
%     rasterFreqLines(:,2) = uniqueFreqs;
%     %this generates green lines separating by Frequency
%     tempHold = 1;
%     for k = 1:size(uniqueFreqs,1)
%         plot(s.Parameters.RasterWindow,[tempHold+sum(matrixTrialNum(k,:)) tempHold+sum(matrixTrialNum(k,:))],'g','LineWidth',1)
%         tempHold = tempHold + sum(matrixTrialNum(k,:));
%     end
%     set(gca,'YTick',rasterFreqLines(:,1));
%     set(gca,'YTickLabel',rasterFreqLines(:,2));
%     set(gca,'Ydir','reverse')
%     ylim([0 totalTrialNum])
%     xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
% %     title('Descending = increase in amplitude and freq')
%     title({fileName;cluNames{i}},'fontweight','bold', 'Interpreter', 'none');
%     hold off
%     spikeGraphName = strcat(fileName,cluNames{i},'JustRaster');
%     savefig(hFig,spikeGraphName);
% 
%     %save as PDF with correct name
%     set(hFig,'Units','Inches');
%     pos = get(hFig,'Position');
%     set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(hFig,spikeGraphName,'-dpdf','-r0')
%     
%     close
% end
% 
% 
% for i = 1:numUnits
%     hFig = figure;
%     set(hFig, 'Position', [10 80 1900 1000])
%     imagesc(s.(cluNames{i}).BinTone')
%     colormap(parula)
%     colorbar
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     set(gca,'YTick',dbRange(:,2));
%     set(gca,'YTickLabel',dbRange(:,1));
% %     title('Descending = increase in amplitude and freq')
%     title({fileName;cluNames{i}},'fontweight','bold', 'Interpreter', 'none');
%     hold off
%     spikeGraphName = strcat(fileName,cluNames{i},'JustTuningPlot');
%     savefig(hFig,spikeGraphName);
% 
%     %save as PDF with correct name
%     set(hFig,'Units','Inches');
%     pos = get(hFig,'Position');
%     set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(hFig,spikeGraphName,'-dpdf','-r0')
%     
%     close
% end


for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])

    %% Column 1
    %plots average waveform
%     subplot(4,8,1)
%     hold on
%     plot(s.(cluNames{i}).AverageWaveForms,'LineWidth',2)
%     if ismember(i,findPVs)
%         title(strcat('PV AverageFiringRate:',num2str(s.(cluNames{i}).AverageRate)))
%     elseif ismember(i,findMSNs)
%         title(strcat('MSN AverageFiringRate:',num2str(s.(cluNames{i}).AverageRate)))
%     elseif ismember(i,findCHATs)
%         title(strcat('CHAT AverageFiringRate:',num2str(s.(cluNames{i}).AverageRate)))
%     else
%         title(strcat('UNK AverageFiringRate:',num2str(s.(cluNames{i}).AverageRate)))
%     end

    %plots ISI
%     subplot(4,8,2)
%     hist(s.(cluNames{i}).ISIGraph,1000)
%     histMax = max(hist(s.(cluNames{i}).ISIGraph,1000));
%     line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
%     xlim(s.Parameters.ClusterWindow)
%     title({strcat('ISI RPV %: ',num2str(s.(cluNames{i}).RPVPercent));...
%         strcat(num2str(s.(cluNames{i}).RPVNumber),'/',num2str(s.(cluNames{i}).TotalSpikeNumber))})

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

    %% Column 2
    %Plot binned response during tone period
    subplot(4,4,2)
    imagesc(s.(cluNames{i}).BinTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title({fileName;cluNames{i}},'fontweight','bold', 'Interpreter', 'none');

    %Plot binned response during general period
    subplot(4,4,6)
    imagesc(s.(cluNames{i}).BinGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean Binned Response (general)')
    %Plot peak response during tone period
    subplot(4,4,10)
    imagesc(s.(cluNames{i}).PeakMapTone')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Peak Response (tone)')
    %Plot peak response during general period
    subplot(4,4,14)
    imagesc(s.(cluNames{i}).PeakMapGen')
    colormap(parula)
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Peak Response (general)')

    %% Column 3

    %plot out binned spikes (tone)
    subplot(4,4,7)
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

    subplot(4,4,11)
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
    subplot(4,4,15)
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

    %% Column 4
    %plot velocity data
%     subplot(4,4,4)
%     hold on
%     plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'b')
%     plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],s.(cluNames{i}).SessionFiring/max(s.(cluNames{i}).SessionFiring),'r')
%     xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
%     ylim([-0.1,1])
%     if toggleROC == 1
%         title(strcat('Vel & Firing Rate. AUC:',num2str(s.(cluNames{i}).TrueAUC),'99%Range',num2str(prctile(s.(cluNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(cluNames{i}).ShuffleAUC,1))))
%     else
%         title('Vel & Firing Rate')
%     end


    %plot heatmap organized by frequency
    subplot(4,4,8)
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

%     if idToggle == 0
% 
%     elseif idToggle == 1
%         %plot laser raster
%         subplot(4,4,12)
%         plot(s.(cluNames{i}).LaserRasters(:,1),s.(cluNames{i}).LaserRasters(:,2),'k.')
%         xlim([s.Parameters.LaserWindow(1),s.Parameters.LaserWindow(2)])
%         title(strcat('Laser Raster, Response %:',num2str(sum(s.(cluNames{i}).LaserResps)/length(dio2Times))))
%         subplot(4,4,16)
%         plot(laserBinVect,s.(cluNames{i}).LaserHist)
%         xlim([s.Parameters.LaserWindow(1),s.Parameters.LaserWindow(2)])
%         title('Laser Histogram')
%     end



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

