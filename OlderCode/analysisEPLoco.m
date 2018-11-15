

function [s] = analysisEPLoco(fileName);

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
s.Parameters.RasterWindow = [-1 1];
s.Parameters.LocoSplit = 2;
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

histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins

%% sets up file saving stuff
saveName = strcat(fileName,'LightToneAnalysis','.mat');
fname = saveName;
pname = pwd;
%check to see if analysis file already exists.
% searchNames = what;
% searchNames = searchNames.mat;
% searchResult = strcmp(searchNames,saveName);
% if find(searchResult == 1) 
%     prevFile = 1;
%     disp('Previous Analysis Found!!')
% else
%     prevFile = 0;
% end

%% Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
if ispc
    subFoldersCell = strsplit(subFolders,';')';
elseif ismac
    subFoldersCell = strsplit(subFolders,':')';
end


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
%     tarFileName = uigetfile('*.mda',strcat('Select Shank-',num2str(i)));
        tarFileName = strcat('nt',num2str(i),'Sort.mda');
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
        
%         clipName = uigetfile('*.mda',strcat('Select Clips for-',num2str(i)));
        clipName = strcat('nt',num2str(i),'SortWave.mda');
        recordClips = readmda(clipName);
        %finally, pull spikes
        for j = 1:numUnits
            targetName = strcat('nt',num2str(i),'_',num2str(uniqueUnits(j)));
            desigNames{counter+j-1} = targetName;
            nameHold{counter - 1 + j} = targetName;
            targetSpikes = holder(2,holder(3,:) == uniqueUnits(j));
            
            
            %make crude autocorrelogram
            pruned = diff(targetSpikes);
            findDoubles = find(pruned < 0.0005);
            targetSpikes(findDoubles) = [];
            pruned = diff(targetSpikes);
            pruned(pruned>s.Parameters.ClusterWindow(2)) = [];
            s.(targetName).SpikeTimes = targetSpikes;
            s.(targetName).ISIGraph = hist(pruned,clusterVect);
            %store rpv violation number
            s.(targetName).RPVNum = length(find(pruned < s.Parameters.RPVTime));
%             %now store all spikes belonging to the target unit
%             s.(targetName).AllSpikes = -recordClips(:,:,holder(3,:) == uniqueUnits(j));
            %now pull template
            s.(targetName).TemplateSpike = recordClips(:,:,j);
        end
    end
    counter = counter + numUnits;
end

%store this information in combined manner
s.OrderData = unitData;


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
%find locomotor periods and split up into 2 second chunks?
s.Parameters.LocoSplit = 2;

runBoutLength = s.RotaryData.PosStartsEnds(:,2) - s.RotaryData.PosStartsEnds(:,1);
timeBin = mean(diff(s.RotaryData.Velocity(:,1)));
sampleRange = s.Parameters.LocoSplit/timeBin;
%find small bouts, remove
findBig = find(runBoutLength > sampleRange);
counter = 1;
locoPoints = [];
for i = 1:length(findBig)
    numFits = floor(runBoutLength(findBig(i))/sampleRange);
    locoPoints(counter:counter+numFits-1) = [s.RotaryData.PosStartsEnds(findBig(i),1):sampleRange:s.RotaryData.PosStartsEnds(findBig(i),1)+sampleRange*(numFits-1)+1];
    counter = counter + numFits;
end

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
%we need for whileCounter, which is the starting point of the search to
%equal the first zero value. 
whileCounter = find(magData < 200,1,'first');
% whileCounter = 1;
onsetStore = [];
while whileTrig == 0
%     whileInd
    %first, establish break for when exceeds length
    testFind = find(magData(whileCounter:end) > threshHiPass);
    if length(testFind) == 0
        break
    end
    %now lets do our iterative search. first, find the next point exceeding
    %the high threshold
    findNext = find(magData(whileCounter:end) > threshHiPass,1,'first');
%     disp('Finding Next')
    %next, find the earliest point before with a slope of less than 0
    %(bottom)
    findBottom = find(magDiff(1:whileCounter + findNext-2) < 0, 1, 'last');
%     disp('Finding Bottom')
    onsetStore(whileInd) = findBottom + 2;
    %store magnitude
    peakSizeStore(whileInd) = magData(findNext + whileCounter - 1) - magData(findBottom+2);
%     disp('Storing Mag')
    %now that we've stored the value, we need to find the next ok point. 
    findLim = find(magData(whileCounter + findNext:end) < threshReset,1,'first');
    if length(findLim) == 0
        break
    end
    whileCounter = whileCounter + findNext + findLim;
    wholeStore(whileInd,1) = findNext;
    wholeStore(whileInd,2) = findBottom;
    wholeStore(whileInd,3) = findLim;
    wholeStore(whileInd,4) = whileCounter;
%     disp('Storing whole')
    
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

s.EDR.OnsetTimes = edrOnsetTimes;
s.EDR.MagData = magData;
s.EDR.TimeVector = edrMagTimes;


%% Now process spiking data aligned to various points
%now that we have the desired alignment times, we can basically make simple
%rasters for every single thing. Tone related will be a bit more
%complicated...


for i = 1:numUnits
    %%Process general information
    masterHolder = masterInd;
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes';
    
    averageRate = length(spikeTimes)/(s.RotaryData.Velocity(end,1)-s.RotaryData.Velocity(1,1));
    
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    
    %make a fine grained one as well!
    fineSessFire = hist(spikeTimes,[s.RotaryData.Velocity(1,1):0.01:s.RotaryData.Velocity(end,1)]);
    fineSessFire(end) = 0; %this is to compensate for problems with spikes coming after the period
    fineSessFire(1) = 0; %this is to compensate for spikes coming before the tuning period. 
        
    %now pull waveform data and isi data to categorize cells
    waveForms = s.(desigNames{i}).TemplateSpike;
    
    %pull out biggest waveform
    maxWave = max(-waveForms');
    [maxVal findMax] = max(maxWave);
   
    %chose the big wave, interpolate to fine degree
    chosenWave = -waveForms(findMax,:);
    interpVect = [1:0.1:100];
    interpWave = interp1(1:100,chosenWave,interpVect,'spline');
    
    %now we need to find the peak. Find this starting at point 10. 

    [pkVal pkInd] = max(interpWave(300:end));
    pkInd = pkInd + 300 - 1;
    %now we need to find the minimum following the peak

    [troughVal troughInd] = min(interpWave(pkInd:end));
    troughInd = troughInd + pkInd - 1;

    peakTrough = (troughInd - pkInd)/300000;
    
    %now find spike width
    halfMax = pkVal/2;
    %Find front end of spike. 
    frontEnd = find(interpWave(1:pkInd)>halfMax,1,'first');
    %Find back end of spike. 
    backEnd = find(interpWave(pkInd:end)<halfMax,1,'first');
    backEnd = backEnd + pkInd - 1;
    
    spikeWidth = (backEnd - frontEnd)/300000;
    
    %find ISIs
    isiTimes = diff(spikeTimes);
    isiCov = std(isiTimes)/mean(isiTimes);
    
    masterData(i,masterHolder) = peakTrough;
    masterHeader{masterHolder} = 'PeakTrough(ms)';
    masterHolder = masterHolder + 1;
    
%     spikeWidth
    masterData(i,masterHolder) = spikeWidth;
    masterHeader{masterHolder} = 'SpikeWidth(ms)';
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
    
    
    s.(desigNames{i}).SessionFiring = sessionFiring;
    s.(desigNames{i}).FineSession = fineSessFire;
    s.(desigNames{i}).FineSessionVector = [s.RotaryData.Velocity(1,1):0.01:s.RotaryData.Velocity(end,1)];
    s.(desigNames{i}).HistBinVector = histBinVector;
    s.(desigNames{i}).AverageRate = averageRate;

    %% Now do rasters aligned to significant deviations based on EDR
    [rastersEDR] = functionBasicRaster(spikeTimes,edrOnsetTimes,s.Parameters.RasterWindow);
    [histCounts histCenters] = hist(rastersEDR(:,1),histBinVector);
    histEDR = histCounts'/length(edrOnsetTimes)/s.Parameters.histBin;
    s.(desigNames{i}).EDRRaster = rastersEDR;
    s.(desigNames{i}).EDRHist = histEDR;
    bigHistEDR(:,i) = histEDR; %for plotting purposes!
    
    %% Now do rasters aligned to locomotion bout onsets
    [rastersLocoStart] = functionBasicRaster(spikeTimes,s.RotaryData.Velocity(s.RotaryData.PosStartsEnds(:,1),1),s.Parameters.RasterWindow);
    [histCounts histCenters] = hist(rastersLocoStart(:,1),histBinVector);
    histLocoStart = histCounts'/length(edrOnsetTimes)/s.Parameters.histBin;
    s.(desigNames{i}).locoStartRaster = rastersLocoStart;
    s.(desigNames{i}).locoStartHist = histLocoStart;
    bigHistLocoStart(:,i) = histLocoStart; %for plotting purposes!
    
    %now stash values for the binned period before and after
    binningPeriod = 0.5;
    for j = 1:length(s.RotaryData.PosStartsEnds)
        targetSpikes = rastersLocoStart(rastersLocoStart(:,2) == j,1);
        s.(desigNames{i}).BinnedLocoStart(j,1) = length(find(targetSpikes > -binningPeriod & targetSpikes < 0));
        s.(desigNames{i}).BinnedLocoStart(j,2) = length(find(targetSpikes < binningPeriod & targetSpikes > 0));
    end
    
    %do stats, save in master. ########
    masterData(i,masterHolder) = signrank(s.(desigNames{i}).BinnedLocoStart(:,1),s.(desigNames{i}).BinnedLocoStart(:,2));
    masterHeader{masterHolder} = 'SignRankLocoStartSig';
    masterHolder = masterHolder + 1;
    %determine sign
    masterData(i,masterHolder) = (mean(s.(desigNames{i}).BinnedLocoStart(:,2)) - mean(s.(desigNames{i}).BinnedLocoStart(:,1)))/(mean(s.(desigNames{i}).BinnedLocoStart(:,2)) + mean(s.(desigNames{i}).BinnedLocoStart(:,1)));
    masterHeader{masterHolder} = 'SignRankLocoStartMod';
    masterHolder = masterHolder + 1;
    
    %and offsets
    [rastersLocoEnd] = functionBasicRaster(spikeTimes,s.RotaryData.Velocity(s.RotaryData.PosStartsEnds(:,2),1),s.Parameters.RasterWindow);
    [histCounts histCenters] = hist(rastersLocoEnd(:,1),histBinVector);
    histLocoEnd = histCounts'/length(edrOnsetTimes)/s.Parameters.histBin;
    s.(desigNames{i}).locoEndRaster = rastersLocoEnd;
    s.(desigNames{i}).locoEndHist = histLocoEnd;
    bigHistLocoEnd(:,i) = histLocoEnd; %for plotting purposes!
    %now stash values for the binned second before and after
    for j = 1:length(s.RotaryData.PosStartsEnds)
        targetSpikes = rastersLocoEnd(rastersLocoEnd(:,2) == j,1);
        s.(desigNames{i}).BinnedLocoEnd(j,1) = length(find(targetSpikes > -binningPeriod & targetSpikes < 0));
        s.(desigNames{i}).BinnedLocoEnd(j,2) = length(find(targetSpikes < binningPeriod & targetSpikes > 0));
    end
    
    %do stats, save in master. ########
    masterData(i,masterHolder) = signrank(s.(desigNames{i}).BinnedLocoEnd(:,1),s.(desigNames{i}).BinnedLocoEnd(:,2));
    masterHeader{masterHolder} = 'SignRankLocoEndSig';
    masterHolder = masterHolder + 1;
    
    %determine sign
    masterData(i,masterHolder) = (mean(s.(desigNames{i}).BinnedLocoEnd(:,2)) - mean(s.(desigNames{i}).BinnedLocoEnd(:,1)))/(mean(s.(desigNames{i}).BinnedLocoEnd(:,2)) + mean(s.(desigNames{i}).BinnedLocoEnd(:,1)));
    masterHeader{masterHolder} = 'SignRankLocoEndMod';
    masterHolder = masterHolder + 1;
    
    
    %% Now do locomotion overall correlation
    
    if toggleROC == 1
        [velOut] = functionLocomotionROC(spikeTimes,s.RotaryData.Velocity);
        s.(desigNames{i}).TrueAUC = velOut.TrueAUC;
        s.(desigNames{i}).ShuffleAUC = velOut.ShuffleAUC;
        masterData(i,masterHolder) = velOut.TrueAUC;
        masterHeader{masterHolder} = 'LocoAUC';
        masterHolder = masterHolder + 1;
        %now generate 99 to 1 percentile range. everything above considered
        %significant
        bigRange = [prctile(s.(desigNames{i}).ShuffleAUC,0.5) prctile(s.(desigNames{i}).ShuffleAUC,99.5)];
        if velOut.TrueAUC < bigRange(1)
            masterData(i,masterHolder) = -1;
            masterHeader{masterHolder} = 'AUCSignificance';
            masterHolder = masterHolder + 1;
        elseif velOut.TrueAUC > bigRange(2)
            masterData(i,masterHolder) = 1;
            masterHeader{masterHolder} = 'AUCSignificance';
            masterHolder = masterHolder + 1;
        else
            masterData(i,masterHolder) = 0;
            masterHeader{masterHolder} = 'AUCSignificance';
            masterHolder = masterHolder + 1;
        end
    else
        s.(desigNames{i}).TrueAUC = 0;
        s.(desigNames{i}).ShuffleAUC = zeros(1000,1);
        
        masterData(i,masterHolder) = 0;
        masterHeader{masterHolder} = 'AUCSignificance';
        masterHolder = masterHolder + 1;
    end
    
end

masterInd = masterHolder;

%now generate normalized overall histograms for light and movement
findZero = find(histBinVector < 0,1,'last');
for i = 1:numUnits
    bigHistEDRNorm(:,i) = bigHistEDR(:,i) / mean(bigHistEDR([1:findZero],i));
    bigHistLocoStartNorm(:,i) = bigHistLocoStart(:,i) / mean(bigHistLocoStart([1:findZero],i));
    bigHistLocoEndNorm(:,i) = bigHistLocoEnd(:,i) / mean(bigHistLocoEnd([1:findZero],i));
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

%now generate average of different aspects of velocity behavior. 
meanJump = nanmean(diff(edrMagTimes));
convertWindow = round(s.Parameters.RasterWindow/meanJump);
edrVect = [s.Parameters.RasterWindow(1):meanJump:s.Parameters.RasterWindow(2)];
% avEDRstore = zeros(length(edrOnsetTimes),length(edrVect));
avEDRstore = [];
for i = 1:length(edrOnsetTimes)
    %find target
    tarTime = find(edrMagTimes - edrOnsetTimes(i) > 0, 1, 'first');
    avEDRstore(i,:) = magData(tarTime + convertWindow(1):tarTime + convertWindow(2));
end

avEDR = mean(avEDRstore);

meanJump = s.Parameters.InterpolationStepRotary;
convertWindow = round(s.Parameters.RasterWindow/meanJump);
velVect = [s.Parameters.RasterWindow(1):meanJump:s.Parameters.RasterWindow(2)];
avVelStartStore = zeros(length(s.RotaryData.PosStartsEnds),length(velVect));
avVelEndStore = zeros(length(s.RotaryData.PosStartsEnds),length(velVect));
for i = 1:length(s.RotaryData.PosStartsEnds)
    avVelStartStore(i,:) = s.RotaryData.Velocity(s.RotaryData.PosStartsEnds(i,1) + convertWindow(1):s.RotaryData.PosStartsEnds(i,1) + convertWindow(2),2);
    avVelEndStore(i,:) = s.RotaryData.Velocity(s.RotaryData.PosStartsEnds(i,2) + convertWindow(1):s.RotaryData.PosStartsEnds(i,2) + convertWindow(2),2);
end

avVelStart = mean(avVelStartStore);
avVelEnd = mean(avVelEndStore);


%% Plotting!!

%% Want to plot by shanks, overall responses to light, sound, locomotion. Do heatmaps?

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
findFirst = find(masterData(:,1) == 1);
findSecond = find(masterData(:,1) == 2);

%plot out responses to EDR start
subplot(2,3,1)
imagesc(bigHistEDRNorm(:,findFirst)')
title('Shank1 Resp EDR')
colormap('parula')
colorbar

subplot(2,3,4)
imagesc(bigHistEDRNorm(:,findSecond)')
colormap('parula')
colorbar


%plot out responses to loco start
subplot(2,3,2)
imagesc(bigHistLocoStartNorm(:,findFirst)')
title('Shank1 Resp Loco')
colormap('parula')
colorbar

subplot(2,3,5)
imagesc(bigHistLocoStartNorm(:,findSecond)')
colormap('parula')
colorbar


%plot out responses to loco start
subplot(2,3,3)
imagesc(bigHistLocoEndNorm(:,findFirst)')
title('Shank1 Resp LocoEnd')
colormap('parula')
colorbar

subplot(2,3,6)
imagesc(bigHistLocoEndNorm(:,findSecond)')
colormap('parula')
colorbar

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
    normWave = -5*s.(desigNames{i}).TemplateSpike/(min(min(s.(desigNames{i}).TemplateSpike)));
    for j = 1:32
        plot(normWave(j,:)-j)
        tempMin = min(normWave(j,:)-j);
        tempMax = max(normWave(j,:)-j);
        if tempMin < fullMin
            fullMin = tempMin;
        end
        if tempMax > fullMax
            fullMax = tempMax;
        end
    end
    ylim([fullMin,fullMax])
    title({fileName;...
        strcat((desigNames{i}),'Shank',num2str(s.OrderData(i)),'Template')}, 'Interpreter', 'none');
    
    
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
%     X = s.(desigNames{i}).FineSession;
%     Fs = 100;
%     L = length(X);
%     Y = fft(X);
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = Fs*(0:(L/2))/L;
%     plot(f,smooth(P1,41)/max(smooth(P1,41)),'r') %180612, removing
%     general, going to put in just the temp. 
    posLocoFind = find(s.RotaryData.NewSimp == 0);
    X = s.(desigNames{i}).FineSession(posLocoFind);
    Fs = 100;
    L = length(X);
    Y = fft(X);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,smooth(P1,41)/max(smooth(P1,41)),'r') 
    plot(floco,smooth(P1floco,41)/max(smooth(P1floco,41)),'k')
    posLocoFind = find(s.RotaryData.NewSimp == 1);
    X = s.(desigNames{i}).FineSession(posLocoFind);
    Fs = 100;
    L = length(X);
    Y = fft(X);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,smooth(P1,41)/max(smooth(P1,41)),'g') 
    xlim([0.2 10])
    title({strcat('OverallRate:',num2str(s.(desigNames{i}).AverageRate));...
        'norm FFT of 10ms FR (r) FRloco (g) loco (k)'})
    
    
    %big velocity etc plot
    subplot(2,1,2)
    hold on
    %plot velocity
    plot(downsample(s.RotaryData.Velocity(:,1),2),downsample(s.RotaryData.Velocity(:,2)./max(s.RotaryData.Velocity(:,2)),2),'k')
    %plot EDR power
    plot(downsample(edrMagTimes,2),downsample(EDROut.PiezoPower(:,2)./max(EDROut.PiezoPower(:,2)),2),'b')
    %plot firing rate
    plot(downsample(s.(desigNames{i}).FineSessionVector,2),downsample(smooth(s.(desigNames{i}).FineSession,31)./max(smooth(s.(desigNames{i}).FineSession,31)),2),'r')
    %Consider looking at filtered version?
%     Wn = [0.5*2/100 5*2/100]; % %confirmed 180601 that this is correct
%     formula. do 2 x desired Hz /sampling rate
%     n = 1000; % 1000th order filter (slower? but 100-order was too low)
%     b = fir1(n, Wn); 
%     filtData = filtfilt(b,1,s.(desigNames{i}).FineSession);
    %restrict timeframe to EDR times.
    edrStart = edrMagTimes(find(~isnan(edrMagTimes),1,'first'));
    edrEnd = edrMagTimes(find(~isnan(edrMagTimes),1,'last'));
    xlim([edrStart edrEnd])
    if toggleROC == 1
        title(strcat('Normalized Vel (k), EDR power (b) and Firing Rate (r) 10ms bin 300ms smooth AUC:',num2str(s.(desigNames{i}).TrueAUC),'99%Range',num2str(prctile(s.(desigNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(desigNames{i}).ShuffleAUC,1))))
    else
        title('Normalized Vel (k), EDR power (b) and Firing Rate (r) 10ms bin 300ms smooth')
    end
    
    
    %Column 2: rasters for edr, histogram
    
    %EDR raster
    subplot(4,4,2)
    hold on
    rasterPlot(s.(desigNames{i}).EDRRaster(:,1),s.(desigNames{i}).EDRRaster(:,2))
    xlim(s.Parameters.RasterWindow)
    ylim([0 length(edrOnsetTimes)])
    title('EDR RASTER')
    
    subplot(4,4,6)
    plotyy(histBinVector,s.(desigNames{i}).EDRHist,edrVect,avEDR)
    xlim(s.Parameters.RasterWindow)
    title('Histogram of EDR Response')
    
    %column 3: histograms for loco onset/offset
    
%     %loco raster
%     subplot(4,4,3)
%     hold on
%     rasterPlot(s.(desigNames{i}).locoStartRaster(:,1),s.(desigNames{i}).locoStartRaster(:,2))
%     xlim(s.Parameters.RasterWindow)
%     ylim([0 length(s.RotaryData.PosStartsEnds)])
%     title('LocoStart RASTER')
    
    %loco Hist
    subplot(4,4,3)
    plotyy(histBinVector,s.(desigNames{i}).locoStartHist,velVect,avVelStart)
    xlim(s.Parameters.RasterWindow)
    title('Histogram of LocoStart Response')
    
%     %loco raster
%     subplot(4,4,4)
%     hold on
%     rasterPlot(s.(desigNames{i}).locoEndRaster(:,1),s.(desigNames{i}).locoEndRaster(:,2))
%     xlim(s.Parameters.RasterWindow)
%     ylim([0 length(s.RotaryData.PosStartsEnds)])
%     title('LocoEnd RASTER')
    
    %loco Hist
    subplot(4,4,7)
    plotyy(histBinVector,s.(desigNames{i}).locoEndHist,velVect,avVelEnd)
    xlim(s.Parameters.RasterWindow)
    title('Histogram of LocoEnd Response')
    
    %column 4: look at heatmap of responses during sustained locomotion. 
    subplot(2,4,4)
    %first, generate smoothed data of 10ms bin firing. 
    smoothFire = smooth(s.(desigNames{i}).FineSession,25);
    %next, fill in heatmap of this!
    heatStore = [];
    for j = 1:length(locoPoints)
        heatStore(:,j) = smoothFire(round(locoPoints(j)):round(locoPoints(j)) + 200);
    end
    imagesc(heatStore')
    colormap('parula')
    colorbar
    
    
    

    
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


