%This code is meant to process mountainsort CSV data. 

fileName = '171206_ML171117F_R10_2300_PVStimTwoLaser50Pulse';

%% Set parameters!
toggleROCLoco = 0;
%commenting out RPV code. meant to do this within mountainsort.
% s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
% s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
% s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.05; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info

%bins of interest, all in ratios
s.Parameters.RasterWindow = [-1 2]; %ratio relative to laser stim 
s.Parameters.PreBin = [-1,0];
s.Parameters.PostBin = [1,2];
s.Parameters.LaserBin = [0,1];
s.Parameters.ResPreBin = [-1,-0.5];
s.Parameters.ResPostBin = [1.5,2];
s.Parameters.ResLaserBin = [0.5,1];

%for rotary encoder:
s.Parameters.InterpolationStepRotary = 0.01;
avLocoThresh = 0.1; %minimum speed for a trial to be considered a locomotion trial.

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%for selecting MSNs vs PV
s.Parameters.PVLim = 0.0005;

%for laser related!
s.Parameters.TrialUnit = 50;
shankElectrodes = 16;
s.ShankLength = 200;




%% Pull out MDA stuffs!

%determine the number of shanks recorded
fileNames = dir('*.csv'); %find csv files
numNames = size(fileNames,1);
shankRecord = [0 0];
for i = 1:numNames
    nameStore{i} = fileNames(i).name;
    for j = 1:2
        shankName = strcat(fileName,'_nt',num2str(j),'firings');
        %check if shank exists
        nameTest = strfind(fileNames(i).name,shankName);
        if nameTest == 1
            shankRecord(j) = 1;
        end
    end
end
%find which shanks are represented
shankFind = find(shankRecord);

counter = 1;
unitData = zeros(2,3);

for i = shankFind
    %first, assign a fileName
    mountainName = strcat(fileName,'_nt',num2str(i),'firings');

    %lets pull the actual unit designations
    nameUnit = strcat(mountainName,'.curated_unitsInd.csv');
    tempHold = csvread(nameUnit);
    numUnits = length(tempHold);
    
    %lets pull the channels file
    nameCh = strcat(mountainName,'.curated_unitsCh.csv');
    channelDesig = csvread(nameCh);
    %this tells me which channel the waveform was biggest on
    
    %fill in unitData
    unitData(counter:counter + numUnits - 1,1) = i;
    unitData(counter:counter + numUnits - 1,2) = channelDesig; %store biggest channel
    unitData(counter:counter + numUnits - 1,3) = tempHold; %store channel name
    %finally, lets pull spikes
    nameSpikes = strcat(mountainName,'.curated_unitsSpikes.csv');
    rawSpikeTimes = csvread(nameSpikes);

    %pull out individual spikes
    for j = 1:numUnits
        targetName = strcat('nt',num2str(i),'_',num2str(tempHold(j)));
        nameHold{counter - 1 + j} = targetName;
        targetSpikes = rawSpikeTimes(rawSpikeTimes(:,2) == tempHold(j),1);
        s.(targetName).SpikeTimes = targetSpikes;
        s.(targetName).ISIGraph = 
        %since we are already here, why not also pull the position of
        %the electrode?
        fullMap(counter - 1 + j,:) = [i channelDesig(j)];
    end
    
    counter = counter + numUnits;
end
%now need to pull total number of units
numUnits = size(unitData,1);

%generate master array for 2-d storage of important values
master = zeros(numUnits,10);
masterHeader = cell(10,1);
masterInd = 1;

desigNames = nameHold;

%now we need to generate a real physical map. Need to add some variation
%for things found on the same channel, to increase spread.

%Go through each electrode and find matching sites. if exist, then store.
%if extras, then figure out spacing. 


cellCount = 1;
posArray = zeros(numUnits,2);
disp('Making Display Distances For Electrodes')
for j = shankFind
    disp(strcat('Starting Shank',num2str(j)))
    for i = 1:shankElectrodes
        %see if there are cells centered on given electrode
        cellFind = find(fullMap(:,2) == i & fullMap(:,1) == j);
        if length(cellFind) == 0
            disp(strcat('No Cells For Electrode',num2str(i)))
        elseif length(cellFind) == 1
            disp(strcat('Single Cell For Electrode',num2str(i)))
            posArray(cellCount,2) = j;
            posArray(cellCount,1) = fullMap(cellFind,2);
            cellCount = cellCount + 1;
        elseif length(cellFind) > 1
            disp(strcat('Multiple Cells For Electrode',num2str(i)))
            incAmount = 1/length(cellFind); %amount to increment different cells by
            for incInd = 1:length(cellFind)
                posArray(cellCount,2) = j;
                posArray(cellCount,1) = fullMap(cellFind(incInd),2) + (incInd-1)*incAmount;
                cellCount = cellCount + 1;
            end
        end
    end
end
%make things negative so they plot out correctly
posArray(:,1) = -1*posArray(:,1);
% posArray(posArray(:,2) == 2,1) = posArray(posArray(:,2) == 2,1)+s.ShankLength;

s.ElectrodePositionDisplayArray = posArray;

%store this into master
master(:,masterInd) = posArray(:,2); masterHeader{masterInd} = 'Shank Designation'; masterInd = masterInd + 1;
master(:,masterInd) = posArray(:,1); masterHeader{masterInd} = 'Electrode Plotting Distance'; masterInd = masterInd + 1;

s.DesignationNames = desigNames;
s.ProbeMap = fullMap;

%now we need to extract DIO binary files
%this is the filename for the .rec file. Please put it in the same dir as
%the config file if you decide to do post-hoc referencing etc.


%% Pull DIO datapoints
%this extracts timestamps
extractTimeBinaryFile(fileName)

%this will extract DIO data (inputs and outputs)
extractDioBinaryFiles(fileName);

homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
if length(D1FileName) == 0
    [D1FileName] = functionFileFinder(subFoldersCell,'DIO','Din2');
end
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
dioTimes = double(DIOData.fields(1).data)/s.Parameters.trodesFS;
%extracts states. States should be 0 or 1.
dioState = double(DIOData.fields(2).data);

%now do this for DIO going off
dioDiff = find(diff(dioState)==-1)+1; %finds all points where there is a positive diff. +1 fudge factor for diff()
dioHigh = find(dioState == 0); %finds all points at which state is 1

dioOff = intersect(dioDiff,dioHigh);
dioOffTimes = dioTimes(dioOff);

%now pull for DIO going on. 

dioDiff = find(diff(dioState)==1)+1; %finds all points where there is a positive diff. +1 fudge factor for diff()
dioHigh = find(dioState == 1); %finds all points at which state is 1

dioTrue = intersect(dioDiff,dioHigh); %finds intersect of dioDiff and dioHigh, which should only be onsets, not offsets or repeat DIO high states
%what i want is times, so extract dioTimes using dioTrue index.
dioTimes = dioTimes(dioTrue);

laserDur = mean(dioOffTimes - dioTimes);

%also calculate differences between DIO input times
dioTimeDiff = diff(dioTimes);

s.Timing.LaserTimes = dioTimes;
s.Timing.LaserDur = laserDur;

%check to see if this is expected number of laser powers?
numMultiple = length(dioTimes)/s.Parameters.TrialUnit;

if floor(numMultiple) == numMultiple
    disp('Checking Number of Trials: Whole Multiple')
    disp(num2str(numMultiple))
else
    error('Checking Number of Trials: Incomplete Multiple! Check File')
end

%now for an important step! We need to determine the duration of the laser
%pulse, since this will reset many things for later analysis.


rasterWindow = s.Parameters.RasterWindow * round(laserDur);
binPre = s.Parameters.PreBin * round(laserDur);
binPost = s.Parameters.PostBin * round(laserDur);
binLaser = s.Parameters.LaserBin * round(laserDur);
binPreRes = s.Parameters.ResPreBin * round(laserDur);
binPostRes = s.Parameters.ResPostBin  * round(laserDur);
binLaserRes = s.Parameters.ResLaserBin * round(laserDur);
disp('Windows Recalculated')


%% Extract data from rotary encoder.
[funcOut] = functionRotaryExtraction(s.Parameters.trodesFS,s.Parameters.InterpolationStepRotary,subFoldersCell);
s.RotaryData = funcOut;
disp('Velocity Data Extracted')
%rasterize this data
jumpsBack = round(rasterWindow(1)/s.Parameters.InterpolationStepRotary);
jumpsForward = round(rasterWindow(2)/s.Parameters.InterpolationStepRotary);

velRaster = zeros(jumpsForward-jumpsBack+1,length(dioTimes));
for i = 1:length(dioTimes)
    %find the time from the velocity trace closest to the actual stim time
    targetInd = find(s.RotaryData.Velocity(:,1)-dioTimes(i) > 0,1,'first');
    if targetInd + jumpsForward < length(s.RotaryData.Velocity) & targetInd - jumpsBack > 1
        %pull appropriate velocity data
        velRaster(:,i) = s.RotaryData.Velocity([targetInd+jumpsBack:targetInd+jumpsForward],2);
    end
end

averageVel = mean(velRaster,2);

%generate vectors for plotting things out later. 
velVector = [rasterWindow(1):s.Parameters.InterpolationStepRotary:rasterWindow(2)];
velDispVector = [rasterWindow(1):1:rasterWindow(2)];
velDispIndex = [1:round(1/s.Parameters.InterpolationStepRotary):(jumpsForward-jumpsBack+1)];
velZero = find(velVector >= 0,1,'first');

%now lets find locomotion based on average from the trial
avLoco = mean(velRaster);
avLocoPos = find(avLoco>avLocoThresh);
avLocoNeg = find(avLoco<=avLocoThresh);
disp('Velocity Rasters and Averages Calculated')

%calculate window sizes and binning
histBinNum = (rasterWindow(2)-rasterWindow(1))/s.Parameters.histBin;
histBinVector = [rasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:rasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins


%% Process spiking information: extract rasters and histograms
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
    [rastersLaser] = functionBasicRaster(spikeTimes,dioTimes,rasterWindow);

    %store data into structured array
    s.(desigNames{i}).RasterLaser = rastersLaser;

    %generate histograms and store.
    [counts centers] = hist(rastersLaser(:,1),histBinVector);
    %find zero point
    histZero = find(histBinVector<0,1,'last');
    s.(desigNames{i}).HistogramLaser = counts/length(dioTimes)/s.Parameters.histBin;
    
    %go through one by one and pull out values that I care about
    infoStore = zeros(length(dioTimes),6);
    histStore = zeros(length(histBinVector),length(dioTimes));
    for j = 1:length(dioTimes)
        %find the spikes linked to the target trial
        selSpikes = rastersLaser(rastersLaser(:,2) == j,1);
        %first generate histogram and store!
        [counts centers] = hist(selSpikes,histBinVector);
        histStore(:,j) = counts;
        if length(selSpikes >0) %only take trials where there are spikes at all
            %fine pre spikes
            findSpikes = find(selSpikes > binPre(1) & selSpikes < binPre(2));
            infoStore(j,1) = length(findSpikes)/(binPre(2)-binPre(1));
            %find post spikes
            findSpikes = find(selSpikes > binPost(1) & selSpikes < binPost(2));
            infoStore(j,2) = length(findSpikes)/(binPost(2)-binPost(1));
            %find laser spikes
            findSpikes = find(selSpikes > binLaser(1) & selSpikes < binLaser(2));
            infoStore(j,3) = length(findSpikes)/(binLaser(2)-binLaser(1));
            
            %now look at restricted bins (which cuts out the first second
            %from each bin)
            
            %find pre spikes
            findSpikes = find(selSpikes > binPreRes(1) & selSpikes < binPreRes(2));
            infoStore(j,4) = length(findSpikes)/(binPreRes(2)-binPreRes(1));
            %find post spikes
            findSpikes = find(selSpikes > binPostRes(1) & selSpikes < binPostRes(2));
            infoStore(j,5) = length(findSpikes)/(binPostRes(2)-binPostRes(1));
            %find laser spikes
            findSpikes = find(selSpikes > binLaserRes(1) & selSpikes < binLaserRes(2));
            infoStore(j,6) = length(findSpikes)/(binLaserRes(2)-binLaserRes(1));
        end
    end
    
    %use histStore to generate z-scored firing from unit
    baselineBins = reshape(histStore(1:histZero,:),1,[]);
    baselineMean = mean(baselineBins);
    baselineSTD = std(baselineBins);
    
    s.(desigNames{i}).zHistLaser = (s.(desigNames{i}).HistogramLaser - baselineMean)/baselineSTD;
    
    
    %now do some basic stats. lets pull a signRank (because data is
    %paired) on all the times!

    %first with the normal bins
    pPreLaser = signrank(infoStore(:,1),infoStore(:,3));
    pPostLaser = signrank(infoStore(:,2),infoStore(:,3));
    pPrePost = signrank(infoStore(:,1),infoStore(:,2));

    %now with the restricted bins
    pResPreLaser = signrank(infoStore(:,4),infoStore(:,6));
    pResPostLaser = signrank(infoStore(:,5),infoStore(:,6));
    pResPrePost = signrank(infoStore(:,4),infoStore(:,5));
    
    
    s.(desigNames{i}).TrialBinnedSpikes = infoStore;
    s.(desigNames{i}).TrialHists = histStore;
    
    %now we want to go through and extract the meaningful time points
    %relative to specific trial types, and store in the master array.
    
    master(i,tempInd) = mean(infoStore(:,1)); masterHeader{tempInd} = 'PreAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(:,2)); masterHeader{tempInd} = 'PostAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(:,3)); masterHeader{tempInd} = 'LaserAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(:,4)); masterHeader{tempInd} = 'PreResAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(:,5)); masterHeader{tempInd} = 'PostResAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = mean(infoStore(:,6)); masterHeader{tempInd} = 'LaserResAverage'; tempInd = tempInd + 1;
    
    %store pValues
    master(i,tempInd) = pPreLaser; masterHeader{tempInd} = 'pValPreLaser'; tempInd = tempInd + 1;
    master(i,tempInd) = pPostLaser; masterHeader{tempInd} = 'pValPostLaser'; tempInd = tempInd + 1;
    master(i,tempInd) = pPrePost; masterHeader{tempInd} = 'pValPrePost'; tempInd = tempInd + 1;
    master(i,tempInd) = pResPreLaser; masterHeader{tempInd} = 'pValResPreLaser'; tempInd = tempInd + 1;
    master(i,tempInd) = pResPostLaser; masterHeader{tempInd} = 'pValResPostLaser'; tempInd = tempInd + 1;
    master(i,tempInd) = pResPrePost; masterHeader{tempInd} = 'pValResPrePost'; tempInd = tempInd + 1;
    
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


%now lets extract data relative to the different kinds of trials
disp('Processing Information Regarding Different Laser Powers')
for i = 1:numUnits  
    disp(strcat('Processing',desigNames{i}))
    tempInd = masterInd;%need to have this or masterInd will creep up with each for loop iteration.
    %we basically now want to separate binned values and histograms into
    %bins corresponding to the number of repetitions. 
    
    %lets first make histograms of the different epocs
    laserHists = zeros(numMultiple,length(histBinVector));
    for j =1:numMultiple
        laserHists(j,:) = mean(s.(desigNames{i}).TrialHists(:,((j-1)*s.Parameters.TrialUnit+1:j*s.Parameters.TrialUnit))')/s.Parameters.histBin;
    end
    
    s.(desigNames{i}).LaserHistograms = laserHists;
    
    %now lets pull the time periods from infoStore
    
    laserSpikes = zeros(numMultiple,size(s.(desigNames{i}).TrialBinnedSpikes,2));
    for j = 1:numMultiple
        laserSpikes(j,:) = mean(s.(desigNames{i}).TrialBinnedSpikes(((j-1)*s.Parameters.TrialUnit+1:j*s.Parameters.TrialUnit),:));
    end
    
    s.(desigNames{i}).LaserBinnedSpikes = laserSpikes;
    
    master(i,tempInd) = laserSpikes(1,1); masterHeader{tempInd} = 'PreLowAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = laserSpikes(1,2); masterHeader{tempInd} = 'PostLowAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = laserSpikes(1,3); masterHeader{tempInd} = 'LaserLowAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = laserSpikes(2,1); masterHeader{tempInd} = 'PreHiAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = laserSpikes(2,2); masterHeader{tempInd} = 'PostHiAverage'; tempInd = tempInd + 1;
    master(i,tempInd) = laserSpikes(2,3); masterHeader{tempInd} = 'LaserHiAverage'; tempInd = tempInd + 1;
    
    indChange = tempInd - masterInd;
    
end

%update master index
masterInd = masterInd + indChange;

s.MasterSheet = master;
s.MasterDesigs = masterHeader;


%% Prepare summary figure information
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%find indices of master that I want for various things
% [indPVMSN] = functionCellStringFind(masterHeader,'PV/MSN Desig');
% [indPkTrough] = functionCellStringFind(masterHeader,'PeakTroughTimeinMS');
[indOverFire] = functionCellStringFind(masterHeader,'OverallFiringRate');
[indDistance] = functionCellStringFind(masterHeader,'Electrode Plotting Distance');
[indShankDes] = functionCellStringFind(masterHeader,'Shank Designation');
[indPreAverage] = functionCellStringFind(masterHeader,'PreAverage');
[indPostAverage] = functionCellStringFind(masterHeader,'PostAverage');
[indLaserAverage] = functionCellStringFind(masterHeader,'LaserAverage');
[indResPre] = functionCellStringFind(masterHeader,'PreResAverage');
[indResPost] = functionCellStringFind(masterHeader,'PostResAverage');
[indResLaser] = functionCellStringFind(masterHeader,'LaserResAverage');
[indLocoAUC] = functionCellStringFind(masterHeader,'LocoAUCScore');
[indLocoSig] = functionCellStringFind(masterHeader,'LocoAUCSignificance');
[indPValPreLaser] = functionCellStringFind(masterHeader,'pValPreLaser');
[indPValPostLaser] = functionCellStringFind(masterHeader,'pValPostLaser');
[indPValPrePost] = functionCellStringFind(masterHeader,'pValPrePost');
[indPValResPreLaser] = functionCellStringFind(masterHeader,'pValResPreLaser');
[indPValResPostLaser] = functionCellStringFind(masterHeader,'pValResPostLaser');
[indPValResPrePost] = functionCellStringFind(masterHeader,'pValResPrePost');
[indPreLowAverage] = functionCellStringFind(masterHeader,'PreLowAverage');
[indPostLowAverage] = functionCellStringFind(masterHeader,'PostLowAverage');
[indLaserLowAverage] = functionCellStringFind(masterHeader,'LaserLowAverage');
[indPreHiAverage] = functionCellStringFind(masterHeader,'PreHiAverage');
[indPostHiAverage] = functionCellStringFind(masterHeader,'PostHiAverage');
[indLaserHiAverage] = functionCellStringFind(masterHeader,'LaserHiAverage');

%% Plot Summary Figure

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])

%plot AUC vs firing rate
subplot(3,3,1)
plot(master(:,indLocoAUC),master(:,indPreAverage),'k.')
hold on
% plot(master(master(:,indLocoSig) == 1,indLocoAUC),master(master(:,indLocoSig) == 1,indPreAverage),'ro')
title('Scatter Plot of Baseline Rate vs LocoAUC')

%plot overall velocity trace. 
subplot(3,3,4)
hold on
plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2),'r')
maxVals = max(s.RotaryData.Velocity(:,2));
xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
% title('Overall Vel With Tone(k) Laser(b) and ToneLaser(g) Times')
title('Overall Velocity Trace')


%plot rasterized velocities
subplot(3,3,7)
hold on
plot(velVector,averageVel,'k')
xlim([velVector(1) velVector(end)]);
title('Average Velocity for Tone(k) Laser(b) ToneLaser(g)')


%plot firing rate vs peak trough times
% subplot(3,3,2)
% plot(master(:,indPkTrough),master(:,indPreAverage),'k.')
% hold on
% plot(master(master(:,indPVMSN) == 1,indPkTrough),master(master(:,indPVMSN) == 1,indPreAverage),'ro')
% title('Peak Trough vs Baseline Rate, PV in red')


%plot modulation index of pre vs laser for laser only trials
subplot(3,3,5)
%calculate modulation index, which is (laser - pre)/(pre + laser)
modInd1 = (master(:,indLaserLowAverage)-master(:,indPreLowAverage))./(master(:,indLaserLowAverage) + master(:,indPreLowAverage));
modInd2 = (master(:,indLaserHiAverage)-master(:,indPreHiAverage))./(master(:,indLaserHiAverage) + master(:,indPreHiAverage));
% modVals1 = hist(modInd1,[-1:0.1:1]);
% modVals2 = hist(modInd2,[-1:0.1:1]);
% bar([modVals1;modVals2])
hist([modInd1,modInd2],[-1:0.1:1]);
xlim([-1 1]);
title('Modulation Index For Laser Only Trials')
% 
% %plot p values!
% %First for normal bins
% subplot(3,6,15)
% hold on
% for i = 1:numUnits
%     plot(master(i,indPValPreLaser:indPValPrePost),'b.-')
% end
% set(gca,'XTick',[1:3]);
% set(gca,'XTickLabel',{'PreLaser','PostLaser','PrePost'});
% title('P Values for Comparisons of Timing')
% 
% 
% %now for restricted bins
% subplot(3,6,16)
% hold on
% for i = 1:numUnits
%     plot(master(i,indPValResPreLaser:indPValResPrePost),'b.-')
% end
% set(gca,'XTick',[1:3]);
% set(gca,'XTickLabel',{'PreLaser','PostLaser','PrePost'});
% title('P Values for Comparisons of Restricted Timing')

%plot out modulation indices by unit
%find PVs
% pvs = find(master(:,indPVMSN) == 1);
subplot(3,3,3)
hold on
plot([0 0],[0 numUnits],'k')
for i = 1:numUnits
    plot([0 modInd1(i)],[i i],'b','LineWidth',2)
%     if ismember(i,pvs)
%         plot(modInd1(i),i,'r.')
%     end
end
for i = 1:numUnits
    plot([0 modInd2(i)],[i i],'g--','LineWidth',2)
%     if ismember(i,pvs)
%         plot(modInd2(i),i,'r.')
%     end
end
ylim([0 numUnits+1])
xlim([-1 1])
title('Modulation Index Sorted By Unit')


%plot out by position, do one shank at a time. 
subplot(3,3,6)
hold on
plot([0 0],[min(posArray(:,1)) 0],'k')
firstFind = find(posArray(:,2) == 1); %find units belonging to first shank
% firstArray = posArray(firstFind,:);
%plot out mod-ind 1 (low laser)
for i = 1:length(firstFind)
    plot([0 modInd1(firstFind(i))],[posArray(firstFind(i),1) posArray(firstFind(i),1)],'b','LineWidth',2)
    plot(modInd1(firstFind(i)),posArray(firstFind(i),1),'b.')
%     if ismember(findOrder,pvs)
%         plot(modInd1(findOrder),firstArray(i,1),'r.')
%     end
end
for i = 1:length(firstFind)
    plot([0 modInd2(firstFind(i))],[posArray(firstFind(i),1) posArray(firstFind(i),1)],'g--','LineWidth',2)
    plot(modInd2(firstFind(i)),posArray(firstFind(i),1),'g.')
%     if ismember(findOrder,pvs)
%         plot(modInd2(findOrder),firstArray(i,1),'r.')
%     end
end
xlim([-1 1])
title('Shank 1 Sorted By Position')


subplot(3,3,9)
hold on
plot([0 0],[min(posArray(:,1)) 0],'k')
firstFind = find(posArray(:,2) == 2); %find units belonging to first shank
% firstArray = posArray(firstFind,:);
%plot out mod-ind 1 (low laser)
for i = 1:length(firstFind)
    plot([0 modInd1(firstFind(i))],[posArray(firstFind(i),1) posArray(firstFind(i),1)],'b','LineWidth',2)
    plot(modInd1(firstFind(i)),posArray(firstFind(i),1),'b.')
%     if ismember(findOrder,pvs)
%         plot(modInd1(findOrder),firstArray(i,1),'r.')
%     end
end
for i = 1:length(firstFind)
    plot([0 modInd2(firstFind(i))],[posArray(firstFind(i),1) posArray(firstFind(i),1)],'g--','LineWidth',2)
    plot(modInd2(firstFind(i)),posArray(firstFind(i),1),'g.')
%     if ismember(findOrder,pvs)
%         plot(modInd2(findOrder),firstArray(i,1),'r.')
%     end
end

xlim([-1 1])
title('Shank 2 Sorted By Position')


spikeGraphName = strcat(fileName,'LaserStimSummary');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')



%% Now we need to plot out individual traces!
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    %plots average waveform
%     subplot(3,6,1)
%     hold on
%     plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
%     title(strcat('OverallRate:',num2str(s.(desigNames{i}).OverallFiringRate)))
    %plots ISI
%     subplot(3,6,2)
%     hist(s.(desigNames{i}).ISIGraph,1000)
%     histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
%     line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
%     xlim(s.Parameters.ClusterWindow)
%     title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
%         strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
%     
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
    plot(histBinVector,s.(desigNames{i}).LaserHistograms(1,:),'k','LineWidth',2)
    plot(histBinVector,s.(desigNames{i}).LaserHistograms(2,:),'g','LineWidth',2)
    plot([0 0],[ylim],'b');
    xlim([rasterWindow(1) rasterWindow(2)])
%     title({fileName;desigNames{i}})
    title('Histograms Tone(k) Laser(b) ToneLaser(bg)')
    set(0, 'DefaulttextInterpreter', 'none')
    
    subplot(3,3,2)
    title({fileName;desigNames{i}})
    set(0, 'DefaulttextInterpreter', 'none')
    %plot out changes to response over time to laser only
    subplot(3,3,5)
    hold on
    plot(s.(desigNames{i}).TrialBinnedSpikes(:,1),'ko')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(:,1),11),'k.-')
    
    plot(s.(desigNames{i}).TrialBinnedSpikes(:,2),'bo')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(:,2),11),'b.-')
    
    plot(s.(desigNames{i}).TrialBinnedSpikes(:,3),'go')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(:,3),11),'g.-')
    smoothTrace = smooth(s.(desigNames{i}).TrialBinnedSpikes(:,3),11);
    plot(avLocoPos,smoothTrace(avLocoPos),'bo')
    
    xlim([0 length(dioTimes)])

    title('TimeCourse Of Response to Laser Only Pre(k) Laser(g) Post(r)')
    
    %plot out changes over time for restricted time bins
    subplot(3,3,8)
    hold on
    plot(s.(desigNames{i}).TrialBinnedSpikes(:,4),'ko')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(:,4),11),'k.-')
    
    plot(s.(desigNames{i}).TrialBinnedSpikes(:,5),'bo')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(:,5),11),'b.-')
    
    plot(s.(desigNames{i}).TrialBinnedSpikes(:,6),'go')
    plot(smooth(s.(desigNames{i}).TrialBinnedSpikes(:,6),11),'g.-')
    smoothTrace = smooth(s.(desigNames{i}).TrialBinnedSpikes(:,6),11);
    plot(avLocoPos,smoothTrace(avLocoPos),'bo')
    
    xlim([0 length(dioTimes)])
    
    title('TimeCourse of Restricted Response to Tone noLaser(k) and laser(g)')
    
    %plots rasters (chronological)
    subplot(3,3,3)
    plot(s.(desigNames{i}).RasterLaser(:,1),...
        s.(desigNames{i}).RasterLaser(:,2),'k.','markersize',5)
    hold on
    ylim([0 length(dioTimes)])
    xlim([rasterWindow(1) rasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Laser Response')
    
    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'LaserStimAnalysis');
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


