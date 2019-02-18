



%% Constants and things you might want to tweak

%TOGGLES FOR ENABLING/DISABLING FEATURES
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.

%PARAMETERS FOR BASIC ARRANGEMENT OF DATA
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate

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
s.Parameters.LaserWindow = [-0.5 1]; %plotting window for rasters
s.Parameters.LaserBin = 0.01; %histogram bin size


%for dmr!
s.Parameters.STRFWin = [-0.2 0.03];
% s.Parameters.dmrTiming = 6;
% s.Parameters.expectedDur = 10; %expected time duration. 
%we have multiple datasets, so lets figure out which one is which. 
dateVal = str2num(fileName(1:6));
if dateVal < 190211
%     sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR10mindsT5dsF5.mat';
    sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR400-6400-approx1kHz.mat';
    % dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
    ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimesCORRECTFreqSampling.mat';

%     ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimes.mat';
else 
    sprFile = 'Z:\Max\dmrOutputFiles\extractedDMR4kHz-16kHz_10mindsT5dsF5.mat';
    % dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
    ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-4000flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimes.mat';
end
% s.Parameters.MinFreq = 4000;
% s.Parameters.MaxFreq = 64000;
% s.Parameters.DMRtimeDilation = 0.9995; %dilation of time. Multiply Trodes time by this to get DMR time.  
% s.Parameters.DMRtimeBin = 1/6000;
% newWin = round(s.Parameters.STRFWin/s.Parameters.DMRtimeBin);
% preDelay = 0.05;
% postDelay = 1/6;

%for looking at target window
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
masterData = zeros(numUnits,1);
masterHeader = cell(1,1);
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

%% QC of TTLs

%first step, is to isolate the big ITIs. This separates out the chunks 

findBig = find(dioTimeDiff > 10);
ttlChunks = zeros(length(findBig)+1,2);
for i = 1:length(ttlChunks)
    if i ==1
        ttlChunks(1,1) = 1;
        ttlChunks(1,2) = findBig(1);
    elseif i <= length(findBig)
        ttlChunks(i,1) = findBig(i-1) + 1;
        ttlChunks(i,2) = findBig(i);
    else
        ttlChunks(i,1) = findBig(i-1) + 1;
        ttlChunks(i,2) = length(dioTimes);
    end
end

%now lets determine if the length of any one segment is incorrect
blockLengths = ttlChunks(:,2) - ttlChunks(:,1)+1;
%load the data file.
load(ttlInfo)
lengthTrue = length(ttlOnsetTime);
for i = 1:length(blockLengths)
    if blockLengths(i) == lengthTrue
        disp(strcat('Block-',num2str(i),'-has correct number of TTLs'))
    else
        error(strcat('Block-',num2str(i),'-has incorrect number of TTLs'))
    end
end
%now lets go through more carefully in each block to look at the itis of
%these TTL pulses
devCounter = 1; %this counts number of deviating TTLs
devStore = [];
for i = 1:length(blockLengths)
    %display block number
    disp(strcat('Crawling through block-',num2str(i)))
    %now pull up appropriate ITIs
    blockITIs = diff(dioTimes(ttlChunks(i,1):ttlChunks(i,2)));
    %compare against true values with subtraction
    itiMod = ttlDiffTime - blockITIs;
    %find all remainders that are too big
    bigMod = find(abs(itiMod) > 0.001);
    if bigMod
        disp(strcat(num2str(length(bigMod)),'-deviations found, marking...'))
        for j = 1:length(bigMod)
            devStore(devCounter,1) = i; %store the chunk
            devStore(devCounter,2) = bigMod(i) - 1; %stores the TTLs before and after. 
            devStore(devCounter,3) = bigMod(i) +1;
            devStore(devCounter,4) = ttlChunks(i,1) - 1 + bigMod(i) - 1; %stores the TTLs before and after. 
            devStore(devCounter,5) = ttlChunks(i,1) - 1 + bigMod(i) + 1;
            
            devCounter = devCounter + 1;
        end
    else
        disp('No Deviations Found! :D')
    end
    
end

% the way in which I want to process things is by basically aligning spike
% times to teh start of each stimulus. I will then figure out when laser is
% on /off. 

%% Pull laser times
%find DIO folder and D1 file for analysis
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
if length(D2FileName) == 0
    [D2FileName] = functionFileFinder(subFoldersCell,'DIO','Din2');
end
D2FileName = D2FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%pulls out DIO up state onsets.
[dio2UpTimes,dio2TimeDiff] = functionBasicDIOCheck(DIO2Data,s.Parameters.trodesFS);
%pull down state shifts. 
findDown = find(DIO2Data.fields(2).data == 0);
findDown(1) =[];
dio2DownTimes = double(DIO2Data.fields(1).data(findDown))/s.Parameters.trodesFS;

if length(dio2UpTimes) ~= length(dio2DownTimes)
    error('LASER UP AND LASER DOWN TIMES MISMATCHED')
end

%% now we want to generate a series of times for laser and non-laser. 

load(sprFile)
%store stimulus
% s.SoundData.Stimulus = stimulus;

%downsample to 1 ms. 
if dateVal < 190211
else
    stimulus = stimulus(:,1:6:end);
end

tempVect = linspace(dmrTimes(1),dmrTimes(2),length(stimulus));
for i = 1:length(blockLengths)
    trueDMRtimes = interp1(ttlOnsetTime,dioTimes(ttlChunks(i,1):ttlChunks(i,2)),tempVect);
    dmrStore(:,i) = trueDMRtimes;
end
dmrStep = mean(mode(diff(dmrStore)));
newWin = round(s.Parameters.STRFWin/dmrStep);


%now lets create two sets of arrays, one for starts and stops of laser
%times, one for starts and stops of non-laser times. I will be using a 1
%second buffer period after each change. 

timesLaser = [];
timesNonLaser = [];

for i = 1:length(dio2UpTimes)-1
    timesLaser(i,:) = [dio2UpTimes(i)+1+s.Parameters.STRFWin(1),dio2DownTimes(i)];
    timesNonLaser(i,:) = [dio2DownTimes(i)+1+s.Parameters.STRFWin(1),dio2UpTimes(i+1)];
end

%now lets eliminate the points that arent directly in the stimulus. 
findNon = find(timesLaser(:,1) < dmrStore(1,1));
timesLaser(findNon,:) = [];
findNon = find(timesNonLaser(:,1) < dmrStore(1,1));
timesNonLaser(findNon,:) = [];

findNon = find(timesLaser(:,1) > dmrStore(end,1) & timesLaser(:,1) < dmrStore(1,2));
timesLaser(findNon,:) = [];
findNon = find(timesNonLaser(:,1) > dmrStore(end,1) & timesNonLaser(:,1) < dmrStore(1,2));
timesNonLaser(findNon,:) = [];

%NOW LETS CHECK!
if size(timesLaser) ~= size(timesNonLaser)
    error('ADJUSTED TIMES LASER AND TIMES NON LASER NOT MATCHED')
end

%% Now lets process the DMR. 

%going to go through each thing and perform STRFs
numStore= zeros(numUnits,4);
for i = 1:numUnits
    disp(strcat('Working on DMR For Unit ',num2str(i)))
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    %lets pull the stuff based on all spikes. we will pull in four phases:
    %laser and non-laser for first and second presentations. 
    
    %first laser
    indivStore = zeros(size(stimulus,1),-1*(newWin(1)-newWin(2)-1));
    %pull out spikes in the target window. 
    spikeFinder = [];
    spikeCounter = 1;
    findFirstEnd = find(timesLaser(:,1) < dmrStore(end,1),1,'last');
    for j = 1:findFirstEnd
        tempFind = find(spikeTimes > timesLaser(j,1) & spikeTimes < timesLaser(j,2));
        if tempFind
            spikeFinder(spikeCounter:spikeCounter+length(tempFind)-1) = tempFind;
            spikeCounter = spikeCounter + length(tempFind);
        end
    end
    
    for j = 1:length(spikeFinder)
        findTarget = find(dmrStore(:,1) - spikeTimes(spikeFinder(j)) > 0,1,'first');
        if findTarget
            takeChunk = stimulus(:,findTarget+newWin(1):findTarget+newWin(2));
            indivStore(:,:,j) = takeChunk;
            
        end
    end
    %also pull waveforms!
    averageWave = reshape(mean(s.(desigNames{i}).AllWaves(:,:,spikeFinder),3),[],1);
    s.(desigNames{i}).FirstLaserWave = averageWave;
    s.(desigNames{i}).FirstLaserAverage = mean(indivStore,3);
    s.(desigNames{i}).FirstLaserSpikes = length(spikeFinder);
    numStore(1,i) = length(spikeFinder);
    disp('Pulled First Laser Set')
    %first non laser
    indivStore = zeros(size(stimulus,1),-1*(newWin(1)-newWin(2)-1));
    %pull out spikes in the target window. 
    spikeFinder = [];
    spikeCounter = 1;
    findFirstEnd = find(timesNonLaser(:,1) < dmrStore(end,1),1,'last');
    for j = 1:findFirstEnd
        tempFind = find(spikeTimes > timesNonLaser(j,1) & spikeTimes < timesNonLaser(j,2));
        if tempFind
            spikeFinder(spikeCounter:spikeCounter+length(tempFind)-1) = tempFind;
            spikeCounter = spikeCounter + length(tempFind);
        end
    end
    for j = 1:length(spikeFinder)
        findTarget = find(dmrStore(:,1) - spikeTimes(spikeFinder(j)) > 0,1,'first');
        if findTarget
            takeChunk = stimulus(:,findTarget+newWin(1):findTarget+newWin(2));
            indivStore(:,:,j) = takeChunk;
        end
    end
    %also pull waveforms!
    averageWave = reshape(mean(s.(desigNames{i}).AllWaves(:,:,spikeFinder),3),[],1);
    s.(desigNames{i}).FirstNonLaserWave = averageWave;
    s.(desigNames{i}).FirstNonLaserAverage = mean(indivStore,3);
    s.(desigNames{i}).FirstNonLaserSpikes = length(spikeFinder);
    numStore(2,i) = length(spikeFinder);
    disp('Pulled First Non Laser Set')
    %now second laser
    indivStore = zeros(size(stimulus,1),-1*(newWin(1)-newWin(2)-1));
    %pull out spikes in the target window. 
    spikeFinder = [];
    spikeCounter = 1;
    findFirstEnd = find(timesLaser(:,1) > dmrStore(end,1),1,'first');
    for j = findFirstEnd:length(timesLaser)
        tempFind = find(spikeTimes > timesLaser(j,1) & spikeTimes < timesLaser(j,2));
        if tempFind
            spikeFinder(spikeCounter:spikeCounter+length(tempFind)-1) = tempFind;
            spikeCounter = spikeCounter + length(tempFind);
        end
    end
    for j = 1:length(spikeFinder)
        findTarget = find(dmrStore(:,2) - spikeTimes(spikeFinder(j)) > 0,1,'first');
        if findTarget
            takeChunk = stimulus(:,findTarget+newWin(1):findTarget+newWin(2));
            indivStore(:,:,j) = takeChunk;
        end
    end
    %also pull waveforms!
    averageWave = reshape(mean(s.(desigNames{i}).AllWaves(:,:,spikeFinder),3),[],1);
    s.(desigNames{i}).SecondLaserWave = averageWave;
    s.(desigNames{i}).SecondLaserAverage = mean(indivStore,3);
    s.(desigNames{i}).SecondLaserSpikes = length(spikeFinder);
    numStore(3,i) = length(spikeFinder);
    disp('Pulled Second Laser Set')
    %now second Non laser
    indivStore = zeros(size(stimulus,1),-1*(newWin(1)-newWin(2)-1));
    %pull out spikes in the target window. 
    spikeFinder = [];
    spikeCounter = 1;
    findFirstEnd = find(timesNonLaser(:,1) > dmrStore(end,1),1,'first');
    for j = findFirstEnd:length(timesNonLaser)
        tempFind = find(spikeTimes > timesNonLaser(j,1) & spikeTimes < timesNonLaser(j,2));
        if tempFind
            spikeFinder(spikeCounter:spikeCounter+length(tempFind)-1) = tempFind;
            spikeCounter = spikeCounter + length(tempFind);
        end
    end
    for j = 1:length(spikeFinder)
        findTarget = find(dmrStore(:,2) - spikeTimes(spikeFinder(j)) > 0,1,'first');
        if findTarget
            takeChunk = stimulus(:,findTarget+newWin(1):findTarget+newWin(2));
            indivStore(:,:,j) = takeChunk;
        end
    end
    %also pull waveforms!
    averageWave = reshape(mean(s.(desigNames{i}).AllWaves(:,:,spikeFinder),3),[],1);
    s.(desigNames{i}).SecondNonLaserWave = averageWave;
    s.(desigNames{i}).SecondNonLaserAverage = mean(indivStore,3);
    s.(desigNames{i}).SecondNonLaserSpikes = length(spikeFinder);
    numStore(4,i) = length(spikeFinder);
    disp('Pulled Second Non-Laser Set')
    
    %now we need to pull laser aligned histogram. for both laser onset and
    %offset.
    laserVect = [s.Parameters.LaserWindow(1)+s.Parameters.LaserBin/2:s.Parameters.LaserBin:s.Parameters.LaserWindow(end)-s.Parameters.LaserBin/2];
    laserRaster = functionBasicRaster(spikeTimes,dio2UpTimes,s.Parameters.LaserWindow);
    s.(desigNames{i}).LaserOnsetHist = hist(laserRaster(:,1),laserVect);
    laserRaster = functionBasicRaster(spikeTimes,dio2DownTimes,s.Parameters.LaserWindow);
    s.(desigNames{i}).LaserOffsetHist = hist(laserRaster(:,1),laserVect);
    disp('Pulled Laser Rasters')
end

%now determine cell type and session firing.
disp('Determining Cell Type')
% masterHolder = masterInd;
for i = 1:numUnits
    masterHolder = masterInd;
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
    elseif peakTrough > s.Parameters.PVLim(2) & isiCov > s.Parameters.ChatLim
        masterData(i,masterHolder) = 0; %MSN
    else
        masterData(i,masterHolder) = NaN; %label as unknown
    end
    masterHeader{masterHolder} = 'CellType';
    masterHolder = masterHolder + 1;
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
%     alignTimes = master(:,1);
    
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Distance(1,1)/s.Parameters.trodesFS:s.Parameters.SpeedFiringBins:s.RotaryData.Distance(end,1)/s.Parameters.trodesFS]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    s.(desigNames{i}).SessionFiring = sessionFiring;
end

masterInd = masterHolder;


%% Saving
save(fullfile(pname,fname),'s','masterData','masterHeader');


%% Plotting individual units
[indCellType] = functionCellStringFind(masterHeader,'CellType');
indCellType = indCellType(1);
%determine if there are units in each category
findPVs = find(masterData(:,indCellType) == 1);

findMSNs = find(masterData(:,indCellType) == 0);


for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])

    % Column 1: Average wave, ISI, response to laser on and off
    subplot(3,8,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    if ismember(i,findPVs)
        title('PUTATIVE FSI')
    elseif ismember(i,findMSNs)
        title('PUTATIVE MSN')
    else
        title('UNK TYPE')
    end
    %plots ISI
    subplot(3,8,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    
    %Plot out response to laser onset.
    subplot(3,4,5)
    plot(laserVect,s.(desigNames{i}).LaserOnsetHist)
    xlim(s.Parameters.LaserWindow)
    title('Hist Laser Onset')

    subplot(3,4,9)
    plot(laserVect,s.(desigNames{i}).LaserOffsetHist)
    xlim(s.Parameters.LaserWindow)
    title('Hist Laser Offset')
    
   
   
    %Column 2: waveforms and velocity vs FR
    %plots average waveform for first presentation
    subplot(3,4,2)
    hold on
    plot(s.(desigNames{i}).FirstLaserWave,'k','LineWidth',2)
    plot(s.(desigNames{i}).FirstNonLaserWave,'g','LineWidth',2)
    title({fileName;desigNames{i}},'fontweight','bold', 'Interpreter', 'none');
    
    
    %plots average waveform for first presentation
    subplot(3,4,6)
    hold on
    plot(s.(desigNames{i}).SecondLaserWave,'k','LineWidth',2)
    plot(s.(desigNames{i}).SecondNonLaserWave,'g','LineWidth',2)
    title('Waveforms For Second Presentation')

    %plot FR and velocity
    subplot(3,4,10)
    hold on
    plot(newTimeVector,newVelVector/max(newVelVector),'b')
    plot([newTimeVector(1):s.Parameters.SpeedFiringBins:newTimeVector(end)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
    xlim([newTimeVector(1),newTimeVector(end)])
    ylim([-0.1,1])
    title('FR (r) vs Velocity (b)')
    
    
    
    %Column 3: plot averages as well as subtraction
    subplot(3,4,3)
    newMean = (s.(desigNames{i}).FirstNonLaserAverage*numStore(2,i) + s.(desigNames{i}).SecondNonLaserAverage*numStore(4,i))/(numStore(2,i)+numStore(4,i));
    imagesc(newMean)
    colormap('parula')
    colorbar
    set(gca,'XTick',[0:60:length(s.(desigNames{i}).FirstNonLaserAverage)]);
    set(gca,'XTickLabel',[s.Parameters.STRFWin(1):0.01:s.Parameters.STRFWin(2)]);
    set(gca,'YTick',[1:10:40]);
    set(gca,'YTickLabel',[faxis([40:-10:1])]);
    title(strcat('AverageNonLaser from-',num2str(numStore(2,i)+numStore(4,i)),'-spikes'))
    
    subplot(3,4,7)
    newMean2 = (s.(desigNames{i}).FirstLaserAverage*numStore(1,i) + s.(desigNames{i}).SecondLaserAverage*numStore(3,i))/(numStore(1,i)+numStore(3,i));
    imagesc(newMean2)
    colormap('parula')
    colorbar
    set(gca,'XTick',[0:60:length(s.(desigNames{i}).SecondNonLaserAverage)]);
    set(gca,'XTickLabel',[s.Parameters.STRFWin(1):0.01:s.Parameters.STRFWin(2)]);
    set(gca,'YTick',[1:10:40]);
    set(gca,'YTickLabel',[faxis([40:-10:1])]);
    title(strcat('AverageLaser from-',num2str(numStore(1,i)+numStore(3,i)),'-spikes'))
    
    subplot(3,4,11)
    if sum(sum(abs(newMean))) > 0 & sum(sum(abs(newMean2))) > 0
        newMatrix = newMean2-newMean;
        imagesc(newMatrix)
        colormap('parula')
        colorbar
        set(gca,'XTick',[0:60:length(s.(desigNames{i}).SecondNonLaserAverage)]);
        set(gca,'XTickLabel',[s.Parameters.STRFWin(1):0.01:s.Parameters.STRFWin(2)]);
        set(gca,'YTick',[1:10:40]);
        set(gca,'YTickLabel',[faxis([40:-10:1])]);
        title('AVERAGE NON LASER RESPONSE')
    else
        title('INSUFFICIENT RESPONSES TO AVERAGE')
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


