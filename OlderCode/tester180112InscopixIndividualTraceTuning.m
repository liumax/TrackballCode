%This is code that will read and link MBED inputs with photometry, while
%recording from my big rig. 

% function [s] = analysisPhotometryTuning(fileName);
% fileName = '180105_ML171117A_Tuning';
% fileName = '180105_ML171117B_Tuning';
% fileName = '180105_ML171119A_Tuning';
% fileName = '180105_ML171119B_Tuning';
% fileName = '180201_ML171117A_InscopixTuning';
% fileName = '180201_ML171117B_InscopixTuning';
% fileName = '180205_ML171117A_InscopixTuning';
fileName = '180205_ML171117B_InscopixTuning'; %has missing tone TTL
% fileName = '180208_ML171117A_InscopixTuning';
% fileName = '180208_ML171117B_InscopixTuning';


% traceName = strcat(fileName(1:end-7),'_GMC_Crop_results.mat'); %for first set from 180105
traceName = strcat(fileName(1:16),'_GMC_Crop_results.mat'); % for second set.

%% Parameters

rasterWindow = [-3,5]; %raster window in seconds
thresh = 0.01;
locoTimeStep = 0.1;
photoToggle = 0;


%use diary function to save logfile of analysis. this is good for
%troubleshooting. 
diaryName = strcat(fileName,'ANALYSISLOGFILE');
diary(diaryName)

disp(fileName)

%store things in a big structure
s = struct;
s.Parameters.RasterWindow = rasterWindow;
s.Parameters.PeakThreshold = thresh;
s.Parameters.LocomotionTimeStep = locoTimeStep;

%% Extracts Sound Data from soundFile, including freq, db, reps.
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
expectedITI = soundData.Delays;

%store these in structured array
s.SoundData.UniqueFrequencies = uniqueFreqs;
s.SoundData.UniqueDBs = uniqueDBs;
s.SoundData.NumFreqs = numFreqs;
s.SoundData.NumDBs = numDBs;

toneDur = soundData.ToneDuration;
toneReps = soundData.ToneRepetitions;
totalTrialNum = length(soundData.Frequencies);


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
%% Next, lets pull the MBED stuff

%first, look for tmp file! first, we want to pull all files!
folderFiles = what;
folderFiles= folderFiles.mat;
%set TMP name
tmpName = strcat(fileName,'MBEDTMP.mat');

[findString] = functionCellStringFind(folderFiles,tmpName);
disp('LOOKING FOR MBED TMP FILE')
if findString %if there is a tmp file!
    disp('MBED TMP FILE FOUND! LOADING')
    load(folderFiles{findString})
else
    disp('NO MBED TMP FILE, EXTRACTING...')
    [trialStates, portStates, trialParams] = maxTrialVariablesLickingTask(fileName);
    disp('MBED Data Loaded!')
    save(tmpName,'trialStates','portStates','trialParams');
    disp('SAVED TMP FILE FOR MBED')
end
%now extract ports 1 (tone) and port 7 (inscopix sync)

mbedTone = [portStates.tStamps',portStates.inStates(:,1)];
mbedSync = [portStates.tStamps',portStates.inStates(:,7)];

[mbedTone] = functionSignalDuplicateElim(mbedTone,2);
[mbedSync] = functionSignalDuplicateElim(mbedSync,2);

s.MBED.ToneDelivery = mbedTone;
s.MBED.Sync = mbedSync;

frameTimes = mbedSync(mbedSync(:,2) == 1,1);
frameRate = 1/(mean(diff(frameTimes))/1000);
s.Framerate = frameRate;
%repair frameTimes!
meanITI = mean(diff(frameTimes));
whileTrig = 0;
counter = 1;
frameDiffs = diff(frameTimes);
while whileTrig == 0
    if counter == length(frameDiffs)
        break
    end
    if frameDiffs(counter) > meanITI*1.3
        disp(strcat('Missed Frame Found:',num2str(counter)))
        frameTimes(counter+2:end+1) = frameTimes(counter+1:end);
        frameTimes(counter+1) = frameTimes(counter) + meanITI;
        frameDiffs = diff(frameTimes);
        disp('Missing Frame Repaired')
    end
    counter = counter + 1;
end

toneTimes = mbedTone(mbedTone(:,2) == 1,1);
%CHECK FOR ISSUES BETWEEN ACTUAL SOUND DATA AND RECEIVED TONES
if length(toneTimes) == length(soundData.Delays)
    disp('Tone Times and Tuning Delay Number Matched')
else
    disp('Mismatch in Number of Tone Times and Tuning Delays')
    lengthDiff = length(toneTimes) - length(soundData.Delays);
    if lengthDiff>0
        %first, check to see if there are big ITIs
        diffTone = diff(toneTimes);
        maxDiff = max(expectedITI)*1000*1.1;
        diffFind = find(diffTone>maxDiff);
        toneTimes(diffFind) = [];
        disp('Removed Extra TTL')
        disp(strcat('MBED TTLs:',num2str(length(toneTimes))))
        disp(strcat('File TTLs:',num2str(length(soundData.Delays))))
        figure
        hold on
        plot(diff(toneTimes)/1000)
        plot(expectedITI,'r')
        title('Tone ITIs MBED Blue Expected R')
        xlabel('Tone Number')
        ylabel('Seconds')
    else
        error('Unplanned Error')
    end
end


%Now pull locomotor data, only if no TMP file
[locoData] = functionLocoTmp(fileName,portStates,locoTimeStep);
s.Locomotion = locoData;
disp('Locomotion Data Loaded')

%% Now pull up trace data!
% traceName = strcat(fileName(1:end-7),'_GMC_Crop_results.mat');
load(traceName);
traces = neuron_results.C;
numUnits = size(traces,1);

%% now lets make basic raster!

if length(traces) == length(frameTimes)
    disp('Trace and MBED Frames Equal')
elseif length(traces) ~= length(frameTimes)
    error('Mismatch in Trace and MBED Frames')
end

%now we need to convert things to inscopix frames
alignTimes = zeros(length(toneTimes),1);
for i = 1:length(alignTimes)
    finder = find(frameTimes - toneTimes(i) >=0,1,'first');
    alignTimes(i) = finder;
end

rasterPhotWindow = round(rasterWindow*frameRate);

maskDim = size(neuron_results.Cn);

for j = 1:numUnits
    %store mask
    disp('Storing mask')
    s.(strcat('unit',num2str(j))).Mask = reshape(neuron_results.A(:,j),maskDim(1),maskDim(2));
    
    %% generate point rasters too!
    disp('Generating Point Data...')
    traceData = smooth(traces(j,:),5);
    
    
    dSmoothDS = diff(traceData);
    ddSmoothDS= diff(dSmoothDS); %get double derivative

    shifter = zeros(100,3);
    shiftInd = 1;
    pSign = 1;


    for crawlInd = 1:length(dSmoothDS)
        %for the first point.
        if crawlInd == 1 & dSmoothDS(crawlInd) > 0
            pSign = 1;
        elseif crawlInd == 1 & dSmoothDS(crawlInd) < 0
            pSign = -1;
        elseif crawlInd == 1 & dSmoothDS(crawlInd) == 0
            whileTrig = 0;
            whilePlus = 1;
            while whileTrig == 0
                if dSmoothDS(crawlInd + whilePlus) >0
                    shifter(shiftInd,1) = crawlInd;
                    shifter(shiftInd,2) = 0;
                    shifter(shiftInd,3) = 1;
                    shiftInd = shiftInd + 1;
                    whileTrig = 1;
                elseif dSmoothDS(crawlInd + whilePlus) <0
                    shifter(shiftInd,1) = crawlInd;
                    shifter(shiftInd,2) = 0;
                    shifter(shiftInd,3) = -1;
                    shiftInd = shiftInd + 1;
                    whileTrig = 1;
                elseif dSmoothDS(crawlInd + whilePlus) == 0
                    whilePlus = whilePlus + 1;
                end
            end
            %now for later data points
            %in the case of positive value
        elseif crawlInd >1 & dSmoothDS(crawlInd) > 0
            if pSign == -1
    %             disp('Change Detected')
                shifter(shiftInd,1) = crawlInd;
                shifter(shiftInd,2) = pSign;
                shifter(shiftInd,3) = 1;
                shiftInd = shiftInd + 1;
                pSign = 1;
            elseif pSign == 0
    %             disp('Change Detected')
                shifter(shiftInd,1) = crawlInd;
                shifter(shiftInd,2) = pSign;
                shifter(shiftInd,3) = 1;
                shiftInd = shiftInd + 1;
                pSign = 1;
            end
            %in the case of negative value
        elseif crawlInd > 1 & dSmoothDS(crawlInd) < 0
            if pSign == 1
    %             disp('Change Detected')
                shifter(shiftInd,1) = crawlInd;
                shifter(shiftInd,2) = pSign;
                shifter(shiftInd,3) = -1;
                shiftInd = shiftInd + 1;
                pSign = -1;
            elseif pSign == 0
    %             disp('Change Detected')
                shifter(shiftInd,1) = crawlInd;
                shifter(shiftInd,2) = pSign;
                shifter(shiftInd,3) = -1;
                shiftInd = shiftInd + 1;
                pSign = -1;
            end
        elseif crawlInd > 1 & dSmoothDS(crawlInd) == 0
            pSign = 0;
        end
    end

    %now go through this stuff

    %classify peaks
    for crawlInd = 1:size(shifter,1)
        if shifter(crawlInd,2) == -1 & shifter(crawlInd,3) == 1
            shifter(crawlInd,4) = 1; %This is the dip
        elseif shifter(crawlInd,2) == 1 & shifter(crawlInd,3) == -1
            shifter(crawlInd,4) = 2; %This is the peak
        elseif shifter(crawlInd,2) == 0 & shifter(crawlInd,3) == 1
            shifter(crawlInd,4) = 1; % positive inflection
        elseif shifter(crawlInd,2) == 0 & shifter(crawlInd,3) == -1
            shifter(crawlInd,4) = 0; %negative inflection
        elseif shifter(crawlInd,2) == 1 & shifter(crawlInd,3) == 0
            shifter(crawlInd,4) = 0; %kind of a trough?
        elseif shifter(crawlInd,2) == -1 & shifter(crawlInd,3) == 0
            shifter(crawlInd,4) = 0; %kind of a peak?
        end
    end


    %now lets just remove all things that arent 1 or 2 (trough or peak)
    shifter(shifter(:,4) == 0,:) = [];
    
    whileTrig = 0;

    while whileTrig == 0;
        %check if first shifter value is 1. 
        if shifter(1,4) == 2;
            shifter(1,:) = [];
            disp('Removing Peak with No Trough')
        elseif shifter(1,4) ~= 2;
            whileTrig = 1;
        end
    end
    
    %now go through and eliminate any things where there is no state change (no
    %change in values

    whileTrig = 0;
    whileCounter = 2;
    prevValue = shifter(1,4);
    while whileTrig == 0;
        if whileCounter == size(shifter,1)
            whileTrig = 1;
        end
        currValue = shifter(whileCounter,4);
        if currValue == prevValue
            shifter(whileCounter,:) = [];
            disp('Cutting Duplicate')
        else
            prevValue = currValue;
            whileCounter = whileCounter + 1;
        end
    end

    %now remove the last negative peak, so there is only a positive peak at the
    %end
    if shifter(end,4) ~= 2
        shifter(end,:) = [];
        disp('Trimming Ends')
    end

    %now get peak and trough values.
    for crawlInd = 1:size(shifter,1)
        shifter(crawlInd,5) = traceData(shifter(crawlInd,1));
    end

    %stores peak information. first column is the height, second column is the
    %value of the peak itself, third column is the value of the trough, final
    %value is the time of the peak. 
    peakInds = find(shifter(:,4) == 2);
    troughInds = find(shifter(:,4) == 1);

    if length(peakInds) ~= length(troughInds)
        error('PEAKS AND TROUGHS NOT PRODUCING THE SAME NUMBERS')
    end
    peakVals = zeros(length(peakInds),5);
    peakVals(:,1) = shifter(peakInds,5) - shifter(peakInds-1,5);
    peakVals(:,2) = shifter(peakInds,5);
    peakVals(:,3) = shifter(peakInds-1,5);
    peakVals(:,4) = (shifter(peakInds,1));
    peakVals(:,5) = (shifter(troughInds,1));

    %now, sort peaks by size. 
    [Y I] = sort(peakVals(:,1));


    peakThresh = 1;
    targetPeaks = peakVals(peakVals(:,1) > peakThresh,:);
    targetPeaks(:,5) = targetPeaks(:,5)+3;
    s.(strcat('unit',num2str(j))).AllPeakTimes = targetPeaks(:,5);
    
    %% now lets also make a gaussian convolution
    gaussRange = [-5:1:5];
    norm = normpdf(gaussRange,0,1);
    convTrace = zeros(length(traces),1);
    for normInd = 1:length(targetPeaks(:,5))
        %check if within range
        if targetPeaks(normInd,5) + gaussRange(1) > 0 & targetPeaks(normInd,5) + gaussRange(end) < length(traces)
            convTrace(targetPeaks(normInd,5)-5:targetPeaks(normInd,5)+5) = convTrace(targetPeaks(normInd,5)-5:targetPeaks(normInd,5)+5) + norm';
        end
    end
    disp('Convolution Trace Made')
    s.(strcat('unit',num2str(j))).ConvTrace = convTrace;
    
    
    disp('Point Data Generated...')
    %% now make raster from this.
    [rasters] = functionBasicRaster(targetPeaks(:,5),alignTimes,rasterPhotWindow);
    s.(strcat('unit',num2str(j))).PointRaster = rasters;
    s.(strcat('unit',num2str(j))).PointRasterTime = rasters;
    s.(strcat('unit',num2str(j))).PointRasterTime(:,1) = rasters(:,1)/frameRate;
    
    
    for i = 1:length(alignTimes)
        photoRaster(:,i) = traces(j,alignTimes(i) + rasterPhotWindow(1):alignTimes(i) + rasterPhotWindow(2));
        photoRasterConv(:,i) =convTrace(alignTimes(i) + rasterPhotWindow(1):alignTimes(i) + rasterPhotWindow(2));
    end
    [photoRasterZ,baselineMean,baselineSTD] = functionZScore(photoRaster,rasterWindow(1)*-1,length(alignTimes));
    
    s.(strcat('unit',num2str(j))).PhotoAligned = photoRaster;
    s.(strcat('unit',num2str(j))).PhotoAlignedConv = photoRasterConv;
    s.(strcat('unit',num2str(j))).PhotoAlignedZ = photoRasterZ;
    s.(strcat('unit',num2str(j))).BaselineMean = baselineMean;
    s.(strcat('unit',num2str(j))).BaselineSTD = baselineSTD;
    desigNames{j} = strcat('unit',num2str(j));
    
    photoAverages = zeros(rasterPhotWindow(2)-rasterPhotWindow(1) + 1,numFreqs,numDBs);
    photoRasterStore = cell(numFreqs,numDBs);
    for ind = 1:numFreqs
        subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(ind),3));
        for k = 1:numDBs
            %find all trials of the particular setting
            targetFinder = find(trialMatrix(:,2) == uniqueFreqs(ind) & trialMatrix(:,3) == subUniqueDB(k));
            %pull and average these traces
            tempHolder = photoRasterZ(:,targetFinder);
            photoRasterStore{ind,k} = tempHolder;
            tempHolder = mean(tempHolder');
            photoAverages(:,ind,k) = tempHolder;
            %pull from convolved trace too
            tempHolder = photoRasterConv(:,targetFinder);
            photoRasterStoreConv{ind,k} = tempHolder;
            tempHolder = mean(tempHolder');
            photoAveragesConv(:,ind,k) = tempHolder;
        end
    end
    
    s.(strcat('unit',num2str(j))).PhotoAverages = photoAverages;
    s.(strcat('unit',num2str(j))).PhotoRasterStore = photoRasterStore;
    s.(strcat('unit',num2str(j))).PhotoAveragesConv = photoAveragesConv;
    s.(strcat('unit',num2str(j))).PhotoRasterStoreConv = photoRasterStoreConv;
    
end

s.DesignationName = desigNames;

rasterVelWindow = round(rasterWindow/locoTimeStep);
velRaster = zeros(rasterVelWindow(2) - rasterVelWindow(1) + 1,length(alignTimes));
for ind = 1:length(toneTimes)
    %find the time in the velocity trace
    alignTime = toneTimes(ind)/1000;
    velPoint = find(locoData.Velocity(:,1) - alignTime > 0,1,'first');
    if velPoint + rasterVelWindow(2) < length(locoData.Velocity(:,1)) & velPoint + rasterVelWindow(1) >= 1
        velRaster(:,ind) = locoData.Velocity(velPoint + rasterVelWindow(1):velPoint + rasterVelWindow(2),2); 
    else
        disp(ind)
        disp('Vel Rasters: Not matched up')
    end
end


s.Processed.VelRaster = velRaster;


%now I need to store different means

velAverages = zeros(rasterVelWindow(2) - rasterVelWindow(1) + 1,numFreqs,numDBs);
velRasterStore = cell(numFreqs,numDBs);
for ind = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(ind),3));
    for j = 1:numDBs
        %find all trials of the particular setting
        targetFinder = find(trialMatrix(:,2) == uniqueFreqs(ind) & trialMatrix(:,3) == subUniqueDB(j));
        %pull velocity traces and do the same
        tempVel = velRaster(:,targetFinder);
        velRasterStore{ind,j} = tempVel;
        tempVel = mean(tempVel');
        velAverages(:,ind,j) = tempVel;
        
    end
end
s.Processed.VelAverages = velAverages;
s.Processed.VelStore = velRasterStore;

%generate correct order for displaying things by freq/db
sortingCounter = 1;
for ind = 1:numFreqs
    subUniqueDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(ind),3));
    for j = 1:numDBs
        sortingFinder = find(trialMatrix(:,2) == uniqueFreqs(ind) & trialMatrix(:,3) == subUniqueDB(j));
        sortIndex(sortingFinder) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end

%now lets find the indices for various things to sort out by amplitude

%first, lets create an alternate trialMatrix with standardized DBs
standardDB = unique(trialMatrix(trialMatrix(:,2) == uniqueFreqs(1),3));
altMatrix = trialMatrix;
for ind = 2:numFreqs
    findDBs = unique(altMatrix(altMatrix(:,2) == uniqueFreqs(ind),3));
    for j = 1:length(findDBs)
        altMatrix(altMatrix(:,2) == uniqueFreqs(ind) & altMatrix(:,3) == findDBs(j),3) = standardDB(j);
    end
end

%now using the altMatrix, find DB related stuff!
dbSort = zeros(length(altMatrix),numDBs);
for ind = 1:numDBs
    sortingCounter = 1;
    sortIndex = zeros(length(altMatrix),1);
    for j = 1:numFreqs
        
        sortingFinder = find(altMatrix(:,2) == uniqueFreqs(j) & altMatrix(:,3) == uniqueDBs(ind));
        sortIndex(sortingFinder) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
    dbSort(:,ind) = sortIndex;
    sortIndex = [];
end

% for i = 1:numFreqs
%     for j = 1:numDBs
%         tempData = s.Processed.PhotoStore{i,j};
%         baseVal = min(tempData(130:145,:));
%         peakVal = max(tempData(145:192,:));
%         baseMean = mean(tempData(130:145,:));
%         peakMean = mean(tempData(145:192,:));
%         peakMag = peakVal - baseVal;
%         peakInt = peakMean - baseMean;
%         magStore(:,i,j) = peakMag;
%         intStore(:,i,j) = peakInt;
%     end
% end

% 
% 
% s.Processed.MagStore = magStore;
% s.Processed.IntStore = intStore;
s.SoundData.AltMatrix = altMatrix;
s.SoundData.DBSort = dbSort;


%% Plot everything
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

% create a set of ticks for labeling any display of octaves
if isfield(soundData,'WhiteNoise')
    if soundData.WhiteNoise == 1
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
        %add in the white noise
        octaveRange(2:end+1,:) = octaveRange(1:end,:);
        octaveRange(1,1) = 0;
        octaveRange(1,2) = 1;
    else soundData.WhiteNoise == 0;
        totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));
        
        %This then makes an array of the full octave steps I've made
        octaveRange = zeros(totalOctaves + 1,2);
        octaveRange(1,1) = uniqueFreqs(1);
        for ind = 1:totalOctaves
            octaveRange (ind+1,1) = octaveRange(ind,1)*2;
        end
        %next, I find the positions from uniqueFreqs that match octaveRange
        for ind = 1:size(octaveRange,1);
            octaveRange(ind,2) = find(round(uniqueFreqs) == octaveRange(ind,1));
        end
    end
else
    totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));
    
    %This then makes an array of the full octave steps I've made
    octaveRange = zeros(totalOctaves + 1,2);
    octaveRange(1,1) = uniqueFreqs(1);
    for ind = 1:totalOctaves
        octaveRange (ind+1,1) = octaveRange(ind,1)*2;
    end
    %next, I find the positions from uniqueFreqs that match octaveRange
    for ind = 1:size(octaveRange,1);
        octaveRange(ind,2) = find(round(uniqueFreqs) == octaveRange(ind,1));
    end
end

%create a set of points for the raster axis
rasterAxis = zeros(rasterWindow(2)-rasterWindow(1)+1,2);
rasterAxis(:,1) = [rasterWindow(1):rasterWindow(2)];
%find a second in photometry sample time
rasterAxis(:,2) = [0:size(photoRasterZ,1)/8:size(photoRasterZ,1)];
rasterAxis(1,2) = 1;


%create a set of points for the velRaster axis
velRasterAxis = zeros(rasterVelWindow(2)/10-rasterVelWindow(1)/10+1,2);
velRasterAxis(:,1) = [rasterVelWindow(1)/10:rasterVelWindow(2)/10];
%find a second in photometry sample time
velRasterAxis(:,2) = [0:length(velAverages)/8:length(velAverages)];
velRasterAxis(1,2) = 1;


%also create labels for the rasters. first column is frequencies, second is
%where lines are drawn, third is for placement of labels
rasterLabels = zeros(numFreqs,3);
rasterLabels(:,1) = uniqueFreqs;
rasterLabels(:,2) = [numDBs*soundData.ToneRepetitions:numDBs*soundData.ToneRepetitions:length(soundData.TrialMatrix)];
rasterLabels(:,3) = [numDBs*soundData.ToneRepetitions/2:numDBs*soundData.ToneRepetitions:length(soundData.TrialMatrix)-numDBs*soundData.ToneRepetitions/2];

dbRasterLabels = zeros(numFreqs,3);
dbRasterLabels(:,1) = uniqueFreqs;
dbRasterLabels(:,2) = [soundData.ToneRepetitions:soundData.ToneRepetitions:length(soundData.TrialMatrix)/numDBs];
dbRasterLabels(:,3) = [soundData.ToneRepetitions/2:soundData.ToneRepetitions:length(soundData.TrialMatrix)/numDBs-soundData.ToneRepetitions/2];


rasterTicks = cell(numFreqs,1);
for ind = 1:numFreqs
    rasterTicks{ind} = num2str(rasterLabels(ind,1));
end


%MAKE THE FIGURE
for plotInd = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    imagescLim = [min(min(min(s.(s.DesignationName{plotInd}).PhotoAverages))),max(max(max(s.(s.DesignationName{plotInd}).PhotoAverages)))];

    velMax = max(max(max(velAverages)));
    velMin = min(min(min(velAverages)));


    %plot heatmap of all responses
    subplot(2,4,1)
    imagesc(squeeze(s.(s.DesignationName{plotInd}).PhotoAverages(:,:)'))
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));


    %plot heatmap of all responses
    subplot(2,4,2)
    imagesc(squeeze(s.(s.DesignationName{plotInd}).PhotoAveragesConv(:,:)'))
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));


    %plot peak detection rasters
    subplot(2,4,3)
    hold on
    plot(s.(s.DesignationName{plotInd}).PointRaster(:,1)/frameRate,dbSort(s.(s.DesignationName{plotInd}).PointRaster(:,2)),'k.')

    for i = 1:numFreqs
        plot([rasterWindow(1) rasterWindow(2)],[dbRasterLabels(i,2) dbRasterLabels(i,2)],'g')
    end
    plot([0 0],[1 length(trialMatrix)],'b')
    plot([soundData.ToneDuration soundData.ToneDuration],[1 length(trialMatrix)],'b')
    title('RastersByFreq')
    set(gca,'YTick',dbRasterLabels(:,3))
    set(gca,'YTickLabel',rasterTicks)
    ylim([2 s.SoundData.ToneRepetitions*numFreqs])
    xlim(rasterWindow)
    set(gca,'Ydir','reverse')

    %plot out velocity changes
    subplot(2,4,4)
    imagesc(squeeze(velAverages(:,:))',[velMin velMax])
    colormap('parula')
    colorbar
    set(gca,'XTick',velRasterAxis(:,2));
    set(gca,'XTickLabel',velRasterAxis(:,1));
    title('Velocity Heatmap')

    %plot tuning curve
    subplot(2,1,2)
    tester = reshape(s.(s.DesignationName{plotInd}).PhotoAveragesConv,[],1);
    graphMax = max(tester);
    graphMin = min(tester);
    plot(smooth(tester,round(frameRate/2)))
    hold on
    spaceLength = length(s.(s.DesignationName{plotInd}).PhotoAveragesConv);
    for i = 1:numFreqs
        plot([spaceLength*i spaceLength*i],[graphMin graphMax],'k','LineWidth',2)
        plot([spaceLength*(i-1)+frameRate*rasterWindow(1)*-1 spaceLength*(i-1)+frameRate*rasterWindow(1)*-1],[graphMin graphMax],'r','LineWidth',1)
        plot([spaceLength*(i-1)+frameRate*(rasterWindow(1)*-1+1) spaceLength*(i-1)+frameRate*(rasterWindow(1)*-1+1)],[graphMin graphMax],'r','LineWidth',1)
    end
    xlim([1 length(tester)])
    ylim([graphMin graphMax])


    spikeGraphName = strcat(fileName,'unit',num2str(plotInd),'Figure');

    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
    close
end


saveName = strcat(fileName,'Analysis','.mat');
fname = saveName;
pname = pwd;

save(fullfile(pname,fname),'s');


diary off
% end


