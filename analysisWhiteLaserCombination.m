

function [s] = analysisWhiteLaserCombination(fileName);
%% Constants and things you might want to tweak
%lets set some switches to toggle things on and off.
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.
toggleROCLoco = 0;


s.Parameters.RasterWindow = [-4 5]; %seconds for raster window. 
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.05; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate

%bins of interest
% s.Parameters.CalcBins = [-4,-1;-1,0;-1,2;0,0.2;2,5];
s.Parameters.PreBin = [-4,-1];
s.Parameters.PreToneLaserBin = [-1,0];
s.Parameters.PostBin = [2,5];
s.Parameters.LaserBin = [-1,2];
s.Parameters.ToneBin = [0,0.2];


%for duplicate elimination
s.Parameters.DownSampFactor = 10; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr

%for rotary encoder:
s.Parameters.InterpolationStepRotary = 0.01;
laserPeriod = [0 1]; %time window in seconds around laser onset that I want to analyze. This is if I want to select for trials with locomotion at a specifc time
locoThresh = 0.9; %threshold for time points in which the locomotion is active for trial to be considered locomotion trial
windowPref = [0 1.5]; %preference window in seconds. This is the period over which I determine whether locomotor starts are more or less common than expected by random chance
prefReps = 1000; %repetitions for bootstrapping to determine if locomotor starts are more common during stimulation
avLocoThresh = 0.1; %minimum speed for a trial to be considered a locomotion trial.

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%for selecting MSNs vs PV
s.Parameters.PVLim = 0.0005;

%% sets up file saving stuff
saveName = strcat(fileName,'FullTuningAnalysis','.mat');
fname = saveName;
pname = pwd;

%% Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%% Find and ID Matclust Files for Subsequent Analysis. Generates Structured Array for Data Storage
%pull matclust file names
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');
[paramFiles] = functionFileFinder(subFoldersCell,'matclust','param');
s.NumberTrodes = length(paramFiles)-length(matclustFiles);

%generate placeholder structure
% s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Parameters.RPVTime,...
    matclustFiles,s,s.Parameters.ClusterWindow,subFoldersCell);

%execute duplicate elimination
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

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%generate master array for 2-d storage of important values
master = zeros(numUnits,20);
masterHeader = cell(20,1);
masterInd = 1;

%determine PV vs MSN
for i = 1:numUnits
    tempInd = masterInd; %need to have this or masterInd will creep up with each for loop iteration.
    [funcOut,peakTrough] = functionPVMSNSeparator(s.(desigNames{i}).AverageWaveForms,s.Parameters.PVLim);
    master(i,tempInd) = funcOut; masterHeader{tempInd} = 'PV/MSN Desig'; tempInd = tempInd + 1;
    master(i,tempInd) = peakTrough; masterHeader{tempInd} = 'PeakTroughTimeinMS'; tempInd = tempInd + 1;
    master(i,tempInd) = s.(desigNames{i}).OverallFiringRate; masterHeader{tempInd} = 'OverallFiringRate'; tempInd = tempInd + 1;
    indChange = tempInd - masterInd;
end
masterInd = masterInd + indChange; %updates masterInd value.
disp('PV and MSN Designations Performed')
%calculate window sizes and binning
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes

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


%% Pull data from sound file
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;
toneDur = soundData.ToneDuration;
toneReps = soundData.ToneRepetitions;
laserLag = min(soundData.LaserLag);

%% Extract DIO information. Tuning curve should rely on just one DIO output, DIO1.

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);

%need to separate by trial type!
shortDIO = find(dioTimeDiff < laserLag*0.8);
shortDIODel = unique([shortDIO,shortDIO+1]);
dioTone = dioTimes;
dioTone(shortDIODel) = [];
%this isolates tones!
%now figure out which ones are laser
dioLaser = dioTimes(shortDIO+1);

%now find the tone laser DIO
medDIO = find(dioTimeDiff > laserLag * 0.8 & dioTimeDiff < laserLag * 1.2);
dioToneLaser1 = dioTimes(medDIO+1); %This generates times where the tone occurs during the tone laser
dioToneLaser2 = dioTimes(medDIO); %this pulls the laser time for tone laser trials


%now lets isolate tone only dio.
[C,ia,ib] = intersect(dioTone,dioToneLaser1);
dioToneOnly = dioTone;
dioToneOnly(ia) = [];


[C,ia,ib] = intersect(dioLaser,dioToneLaser2);
dioLaserOnly = dioLaser;
dioLaserOnly(ia) = [];


dioLaserOnlyFake = dioLaserOnly + laserLag; %This is a fake time point to generate rasters that are aligned to same times as tones.

%get indices for various things, as well as all the DIO values I care about
allDIO = unique([dioToneOnly;dioLaserOnly;dioToneLaser1]);
[C dioIndToneOnly ib] = intersect(allDIO,dioToneOnly);
[C dioIndLaserOnly ib] = intersect(allDIO,dioLaserOnly);
[C dioIndToneLaser ib] = intersect(allDIO,dioToneLaser1);

disp('Finished Extracting and Processing DIO data')


%% Extract data from rotary encoder.
[funcOut] = functionRotaryExtraction(s.Parameters.trodesFS,s.Parameters.InterpolationStepRotary,subFoldersCell);
s.RotaryData = funcOut;
disp('Velocity Data Extracted')
%rasterize this data
jumpsBack = round(s.Parameters.RasterWindow(1)/s.Parameters.InterpolationStepRotary);
jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);

velRaster = zeros(jumpsForward-jumpsBack+1,length(allDIO));
for i = 1:length(allDIO)
    %find the time from the velocity trace closest to the actual stim time
    targetInd = find(s.RotaryData.Velocity(:,1)-allDIO(i) > 0,1,'first');
    %pull appropriate velocity data
    velRaster(:,i) = s.RotaryData.Velocity([targetInd+jumpsBack:targetInd+jumpsForward],2);
end

velRasterTone = velRaster(:,dioIndToneOnly);
velRasterLaser = velRaster(:,dioIndLaserOnly);
velRasterToneLaser = velRaster(:,dioIndToneLaser);

%make average traces
averageVelTone = mean(velRasterTone,2);
averageVelLaser = mean(velRasterLaser,2);
averageVelToneLaser = mean(velRasterToneLaser,2);

%generate vectors for plotting things out later. 
velVector = [s.Parameters.RasterWindow(1):s.Parameters.InterpolationStepRotary:s.Parameters.RasterWindow(2)];
velDispVector = [s.Parameters.RasterWindow(1):1:s.Parameters.RasterWindow(2)];
velDispIndex = [1:round(1/s.Parameters.InterpolationStepRotary):(jumpsForward-jumpsBack+1)];
velZero = find(velVector >= 0,1,'first');

%now lets find locomotion based on average from the trial
avLoco = mean(velRaster);
avLocoPos = find(avLoco>avLocoThresh);
avLocoNeg = find(avLoco<=avLocoThresh);

avLocoPre = mean(velRaster(1:(s.Parameters.PreBin(2) - s.Parameters.RasterWindow(1))/s.Parameters.InterpolationStepRotary,:));
avLocoPrePos = find(avLocoPre>avLocoThresh);
avLocoPreNeg = find(avLocoPre<=avLocoThresh);

avLocoPost = mean(velRaster((s.Parameters.RasterWindow(2)-s.Parameters.PostBin(1))/s.Parameters.InterpolationStepRotary:end,:));
avLocoPostPos = find(avLocoPost>avLocoThresh);
avLocoPostNeg = find(avLocoPost<=avLocoThresh);
disp('Velocity Rasters and Averages Calculated')


%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db
disp('Starting Spike Processing')
for i = 1:numUnits  
    disp(strcat('Processing',desigNames{i}))
    tempInd = masterInd;%need to have this or masterInd will creep up with each for loop iteration.
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    s.(desigNames{i}).SessionFiring = sessionFiring;
    
    %I'm still keeping this separate individual code for pulling specific
    %trial types because it does make for a simple plotting of the rasters,
    %since I dont have to map the trial number in some funny way. 
    %now calculate rasters for each case
    [rastersTone] = functionBasicRaster(spikeTimes,dioToneOnly,s.Parameters.RasterWindow);
    [rastersLaser] = functionBasicRaster(spikeTimes,dioLaserOnlyFake,s.Parameters.RasterWindow);
    [rastersToneLaser] = functionBasicRaster(spikeTimes,dioToneLaser1,s.Parameters.RasterWindow);
    
    %store data into structured array
    s.(desigNames{i}).RasterToneOnly = rastersTone;
    s.(desigNames{i}).RasterLaserOnly = rastersLaser;
    s.(desigNames{i}).RasterToneLaser = rastersToneLaser;
    
    %generate histograms and store.
    [counts centers] = hist(rastersTone(:,1),histBinVector);
    s.(desigNames{i}).HistogramToneOnly = counts/length(dioToneOnly)/s.Parameters.histBin;
    
    [counts centers] = hist(rastersLaser(:,1),histBinVector);
    s.(desigNames{i}).HistogramLaserOnly = counts/length(dioLaserOnlyFake)/s.Parameters.histBin;
    
    [counts centers] = hist(rastersToneLaser(:,1),histBinVector);
    s.(desigNames{i}).HistogramToneLaser = counts/length(dioToneLaser1)/s.Parameters.histBin;
    
    
    %now pull rasters for all trials
    [rasters] = functionBasicRaster(spikeTimes,allDIO,s.Parameters.RasterWindow);
    %go through one by one and pull out values that I care about
    infoStore = zeros(length(allDIO),5);
    histStore = zeros(length(histBinVector),length(allDIO));
    for j = 1:length(allDIO)
        %find the spikes linked to the target trial
        selSpikes = rasters(rasters(:,2) == j,1);
        %first generate histogram and store!
        [counts centers] = hist(selSpikes,histBinVector);
        histStore(:,j) = counts;
        if length(selSpikes >0) %only take trials where there are spikes at all
            %fine pre spikes
            findSpikes = find(selSpikes > s.Parameters.PreBin(1) & selSpikes < s.Parameters.PreBin(2));
            infoStore(j,1) = length(findSpikes)/(s.Parameters.PreBin(2)-s.Parameters.PreBin(1));
            %find post spikes
            findSpikes = find(selSpikes > s.Parameters.PostBin(1) & selSpikes < s.Parameters.PostBin(2));
            infoStore(j,2) = length(findSpikes)/(s.Parameters.PostBin(2)-s.Parameters.PostBin(1));
            %find laser spikes
            findSpikes = find(selSpikes > s.Parameters.LaserBin(1) & selSpikes < s.Parameters.LaserBin(2));
            infoStore(j,3) = length(findSpikes)/(s.Parameters.LaserBin(2)-s.Parameters.LaserBin(1));
            %find tone spikes
            findSpikes = find(selSpikes > s.Parameters.ToneBin(1) & selSpikes < s.Parameters.ToneBin(2));
            infoStore(j,4) = length(findSpikes)/(s.Parameters.ToneBin(2)-s.Parameters.ToneBin(1));
            %find preToneLaser bins
            findSpikes = find(selSpikes > s.Parameters.PreToneLaserBin(1) & selSpikes < s.Parameters.PreToneLaserBin(2));
            infoStore(j,5) = length(findSpikes)/(s.Parameters.PreToneLaserBin(2)-s.Parameters.PreToneLaserBin(1));
        end
    end
    
    s.(desigNames{i}).TrialBinnedSpikes = infoStore;
    
    %now we want to go through and extract the meaningful time points
    %relative to specific trial types, and store in the master array.
    
    master(i,tempInd) = mean(infoStore(:,1)); masterHeader{tempInd} = 'PreAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(:,2)); masterHeader{tempInd} = 'PostAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(dioIndLaserOnly,3)); masterHeader{tempInd} = 'LaserOnlyAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(dioIndToneOnly,4)); masterHeader{tempInd} = 'ToneOnlyAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(dioIndToneLaser,4)); masterHeader{tempInd} = 'ToneLaserAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(avLocoPrePos,1)); masterHeader{tempInd} = 'PreAverageRunning'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(avLocoPreNeg,1)); masterHeader{tempInd} = 'PreAverageStationary'; tempInd = tempInd + 1;
    
    
    %now measure locomotion ROC if toggle is on
    if toggleROCLoco == 1
        [velOut] = functionLocomotionROC(spikeTimes,s.RotaryData.Velocity);
        s.(desigNames{i}).TrueAUC = velOut.TrueAUC;
        s.(desigNames{i}).ShuffleAUC = velOut.ShuffleAUC;
    else
        s.(desigNames{i}).TrueAUC = 0;
        s.(desigNames{i}).ShuffleAUC = zeros(1000,1);
    end
    
    %determine percentile of AUC value. Note that we are going to 99.9%
    %confidence, which is kind of cheating since I only do 1000 reps of
    %shuffling
    AUCLims(1) = prctile(s.(desigNames{i}).ShuffleAUC,0.05);
    AUCLims(2) = prctile(s.(desigNames{i}).ShuffleAUC,99.95);
    if s.(desigNames{i}).TrueAUC > AUCLims(2) | s.(desigNames{i}).TrueAUC < AUCLims(1)
        s.(desigNames{i}).AUCSig = 1;
    else
        s.(desigNames{i}).AUCSig = 0;
    end
    %store into master array. 
    master(i,tempInd) = s.(desigNames{i}).TrueAUC; masterHeader{tempInd} = 'LocoAUCScore'; tempInd = tempInd + 1;
    master(i,tempInd) = s.(desigNames{i}).AUCSig; masterHeader{tempInd} = 'LocoAUCSignificance'; tempInd = tempInd + 1;
    
    
    indChange = tempInd - masterInd; %serves as way to update master ind.
end
%display update about progress
disp('Finished Extracting Spike Data')

%update master index
masterInd = masterInd + indChange;

s.MasterSheet = master;
s.MasterDesigs = masterHeader;

%% Prepare summary figure information
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])

%plot AUC vs firing rate
subplot(3,3,1)
plot(master(:,11),master(:,4),'k.')
hold on
plot(master(master(:,12) == 1,11),master(master(:,12) == 1,4),'ro')
title('Scatter Plot of Baseline Rate vs LocoAUC')

%plot overall velocity trace. 
subplot(3,3,4)
hold on
plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2),'r')
maxVals = max(s.RotaryData.Velocity(:,2));
% for i = 1:length(dioIndToneOnly)
%     plot([allDIO(dioIndToneOnly(i)) allDIO(dioIndToneOnly(i))],[0 maxVals],'k')
% end
% for i = 1:length(dioIndLaserOnly)
%     plot([allDIO(dioIndLaserOnly(i)) allDIO(dioIndLaserOnly(i))],[0 maxVals],'b')
% end
% for i = 1:length(dioIndToneLaser)
%     plot([allDIO(dioIndToneLaser(i)) allDIO(dioIndToneLaser(i))],[0 maxVals],'g')
% end
xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
% title('Overall Vel With Tone(k) Laser(b) and ToneLaser(g) Times')
title('Overall Velocity Trace')

%plot rasterized velocities
subplot(3,3,7)
hold on
plot(velVector,averageVelTone,'k')
plot(velVector,averageVelLaser,'b')
plot(velVector,averageVelToneLaser,'g')
xlim([velVector(1) velVector(end)]);
title('Average Velocity for Tone(k) Laser(b) ToneLaser(g)')

%plot firing rate vs peak trough times
subplot(3,3,2)
plot(master(:,2),master(:,4),'k.')
hold on
plot(master(master(:,1) == 1,2),master(master(:,1) == 1,4),'ro')
title('Peak Trough vs Baseline Rate, PV in red')

%plot modulation index of pre vs laser for laser only trials
subplot(3,3,5)
%calculate modulation index, which is (laser - pre)/(pre + laser)
modInd1 = (master(:,6)-master(:,4))./(master(:,6)+master(:,4));
hist(modInd1,[-1:0.1:1])
xlim([-1 1]);
title('Modulation Index For Laser Only Trials')

%plot modulation index of tone responses
subplot(3,3,8)
modInd2 = (master(:,8)-master(:,7))./(master(:,8)+master(:,7));
hist(modInd2,[-1:0.1:1])
xlim([-1 1]);
title('Modulation Index For Tone Trials')

%plot out modulation indices by unit
subplot(3,3,3)
hold on
plot([0 0],[0 numUnits],'k')
for i = 1:numUnits
    plot([0 modInd1(i)],[i i],'b','LineWidth',2)
    plot([0 modInd2(i)],[i+0.2 i+0.2],'g','LineWidth',2)
end
ylim([0 numUnits+1])
title('Modulation Index Sorted By Unit')

%plot out by position, do one shank at a time. 
subplot(3,3,6)
hold on
plot([0 0],[0 -s.ShankLength],'k')
firstFind = find(posArray(:,2) == 1); %find units belonging to first shank
firstArray = posArray(firstFind,:);
for i = 1:length(firstFind)
    findOrder = find(s.SortedPeakWaveOrder == firstFind(i));
    plot([0 modInd1(findOrder)],[firstArray(i,1) firstArray(i,1)],'b')
    plot(modInd1(findOrder),firstArray(i,1),'b.')
    plot([0 modInd2(findOrder)],[firstArray(i,1) firstArray(i,1)],'g')
    plot(modInd2(findOrder),firstArray(i,1),'g.')
end
xlim([-1 1])
title('Shank 1 Sorted By Position')


subplot(3,3,9)
hold on
plot([0 0],[0 -s.ShankLength],'k')
secondFind = find(posArray(:,2) == 2); %find units belonging to first shank
secondArray = posArray(secondFind,:);
for i = 1:length(secondFind)
    findOrder = find(s.SortedPeakWaveOrder == secondFind(i));
    plot([0 modInd1(findOrder)],[secondArray(i,1) secondArray(i,1)],'b')
    plot(modInd1(findOrder),secondArray(i,1),'b.')
    plot([0 modInd2(findOrder)],[secondArray(i,1) secondArray(i,1)],'g')
    plot(modInd2(findOrder),secondArray(i,1),'g.')
end
xlim([-1 1])

title('Shank 2 Sorted By Position')


%% Now we need to plot out individual traces!
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    %plots average waveform
    subplot(3,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    title(strcat('OverallRate:',num2str(s.(desigNames{i}).OverallFiringRate)))
    %plots ISI
    subplot(3,6,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    
    %plot velocity and firing rate
    subplot(3,3,4)
    hold on
    plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'b')
    plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
    xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
    ylim([-0.1,1])
    title(strcat('Vel & Firing Rate. AUC:',num2str(s.(desigNames{i}).TrueAUC),'-99.9%Range',num2str(prctile(s.(desigNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(desigNames{i}).ShuffleAUC,1))))
    
    
    %plot histograms
    subplot(3,3,7)
    hold on
    plot(histBinVector,s.(desigNames{i}).HistogramToneOnly,'k','LineWidth',2)
    plot(histBinVector,s.(desigNames{i}).HistogramToneLaser,'g','LineWidth',2)
    plot(histBinVector,s.(desigNames{i}).HistogramLaserOnly,'b','LineWidth',2)

    plot([0 0],[ylim],'b');
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
%     title({fileName;desigNames{i}})
    title('Histograms Tone(k) Laser(b) ToneLaser(bg)')
    set(0, 'DefaulttextInterpreter', 'none')
    
    subplot(3,3,2)
    title({fileName;desigNames{i}})
    set(0, 'DefaulttextInterpreter', 'none')
    %plot out changes to response over time to laser only
    subplot(3,3,5)
    hold on
    plot(s.(desigNames{i}).TrialBinnedSpikes(dioIndLaserOnly,1),'ko')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(dioIndLaserOnly,1),11),'k.-')
    
    plot(s.(desigNames{i}).TrialBinnedSpikes(dioIndLaserOnly,3),'go')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(dioIndLaserOnly,3),11),'g.-')
    smoothTrace = smooth(s.(desigNames{i}).TrialBinnedSpikes(dioIndLaserOnly,3),11);
    %now we need to remap avLocoPos onto just the dioIndLaserOnly
    [C ia ib] = intersect(dioIndLaserOnly,avLocoPos); %want ia here.
    plot(ia,smoothTrace(ia),'bo')
    
    plot(s.(desigNames{i}).TrialBinnedSpikes(dioIndLaserOnly,2),'ro')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(dioIndLaserOnly,2),11),'r.-')
    title('TimeCourse Of Response to Laser Only Pre(k) Laser(g) Post(r)')
    
    %plot out changes over time for tone
    subplot(3,3,8)
    hold on
    plot(s.(desigNames{i}).TrialBinnedSpikes(dioIndToneOnly,4),'ko')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(dioIndToneOnly,4),11),'k.-')
    
    plot(s.(desigNames{i}).TrialBinnedSpikes(dioIndToneLaser,4),'go')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(dioIndToneLaser,4),11),'g.-')
    smoothTrace = smooth(s.(desigNames{i}).TrialBinnedSpikes(dioIndToneLaser,4),11);
    %now we need to remap avLocoPos onto just the dioIndLaserOnly
    [C ia ib] = intersect(dioIndLaserOnly,avLocoPos); %want ia here.
    plot(ia,smoothTrace(ia),'bo')
    title('TimeCourse of Response to Tone noLaser(k) and laser(g)')
    
    %plots rasters (chronological)
    subplot(3,3,3)
    plot(s.(desigNames{i}).RasterLaserOnly(:,1),...
        s.(desigNames{i}).RasterLaserOnly(:,2),'k.','markersize',5)
    hold on
    ylim([0 toneReps])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Laser Response')
    
    subplot(3,3,6)
    plot(s.(desigNames{i}).RasterToneOnly(:,1),...
        s.(desigNames{i}).RasterToneOnly(:,2),'k.','markersize',5)
    hold on
    ylim([0 toneReps])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Tone Response')
    
    subplot(3,3,9)
    plot(s.(desigNames{i}).RasterToneLaser(:,1),...
        s.(desigNames{i}).RasterToneLaser(:,2),'k.','markersize',5)
    hold on
    ylim([0 toneReps])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Tone Laser Response')
    
    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'threepeatAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end
%% Saving
save(fullfile(pname,fname),'s');

end