%This code is meant for the basic analysis of tuning curve data, examining
%and plotting basic properties of the tuning curve and displaying them as a
%figure. 

%%Inputs
%fileName: name of the target file, without file extension. 

%%Outputs
%s: structured array that has all the data from the analysis. s will
%contain the following: 
%
%n number of units, which are independent of clusters/trodes. This is to
%say that each unit will have its own named field. For example, 2 clusters
%on channel 10 will have two separate fields.

%DesignationArray/DesignationName: The array indicates which probe sites
%and clusters were used, designation names should be the names of all
%units.

%SoundData: all the sound data from the tuning curve, with addition of
%unique frequencies/dbs and number of frequencies/dbs.

function [s] = analysisDATStim(fileName);
%% Constants and things you might want to tweak
%lets set some switches to toggle things on and off.
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.
toggleROC = 1; %toggle for tuning on/off ROC analysis

s.Parameters.RasterWindow = [-4 6]; %seconds for raster window. NOT SCALED
s.Parameters.ToneWindow = [0 0.5];
s.Parameters.GenWindow = [0 1];
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.05; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate
s.Parameters.BaselineBin = [-4 0]; %time in seconds.
% s.Parameters.LFPWindow = [-1 1];
s.Parameters.PrePost = 3;%seconds that the code will look before and after zero point for binning spike responses
if s.Parameters.PrePost > s.Parameters.RasterWindow(2)
    disp('Resetting Pre/Post Window Based on Rasters')
    s.Parameters.PrePost = s.Parameters.RasterWindow(2);
elseif s.Parameters.PrePost > abs(s.Parameters.RasterWindow(1))
    disp('Resetting Pre/Post Window Based on Rasters')
    s.Parameters.PrePost = abs(s.Parameters.RasterWindow(1));
end

%for duplicate elimination
s.Parameters.DownSampFactor = 10; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr

%for rotary encoder:
s.Parameters.InterpolationStepRotary = 0.01;
laserPeriod = [0 1.5]; %time window in seconds around laser onset that I want to analyze. This is if I want to select for trials with locomotion at a specifc time
locoThresh = 0.9; %threshold for time points in which the locomotion is active for trial to be considered locomotion trial
windowPref = [0 1.5]; %preference window in seconds. This is the period over which I determine whether locomotor starts are more or less common than expected by random chance
prefReps = 1000; %repetitions for bootstrapping to determine if locomotor starts are more common during stimulation

%for edr
s.Parameters.EDRdownsamp = 20; %number of samples to downsample by. Smoothing is likely unnecessary
s.Parameters.EDRTimeCol = 1;
s.Parameters.EDRTTLCol = 3;
s.Parameters.EDRPiezoCol = 2;

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%for summary figure
s.Parameters.SumGraphSmooth = 0.200; %smoothing window in seconds
cLims = [-1,1]; %limits for color range for heatmap of z-scored firing rates
%other settings
format short

%% sets up file saving stuff
saveName = strcat(fileName,'DatStimAnalysis','.mat');
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
% s.NumberTrodes = length(paramFiles)-length(matclustFiles);
s.NumberTrodes = 8;
%generate placeholder structure
% s = struct;
%% Fill structure with correct substructures (units, not clusters/trodes) and then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Parameters.RPVTime,...
    matclustFiles,s,s.Parameters.ClusterWindow);

%% Execute duplicate elimination code. 
if toggleDuplicateElimination ==1
    if length(s.DesignationName) > 1
        disp('Now Selecting Based on xCORR')
        [s] = functionDuplicateElimination(s,s.Parameters.DownSampFactor,...
            s.Parameters.corrSlide,s.Parameters.ThresholdComparison,s.Parameters.trodesFS,s.Parameters.RPVTime,s.Parameters.ClusterWindow);
    end
else
    disp('NOT EXECUTING DUPLICATE ELIMINATION')
end


%% Pull variables from parameters for easier looking code.

numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;
rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes

%% Extract DIO information. Will be looking at DIO1, as this should be the output at the beginning of the stim. Can check DIO2 to validate.

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);
totalTrialNum = length(dioTimes);
s.AlignTimes = dioTimes;

%% Extract data from rotary encoder.
[s] = functionRotaryExtraction(s,s.Parameters.trodesFS,s.Parameters.InterpolationStepRotary,subFoldersCell);

%rasterize this data
jumpsBack = round(s.Parameters.RasterWindow(1)/s.Parameters.InterpolationStepRotary);
jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);
velRaster = zeros(jumpsForward-jumpsBack+1,totalTrialNum);
for i = 1:totalTrialNum
    %find the time from the velocity trace closest to the actual stim time
    targetInd = find(s.RotaryData.Velocity(:,1)-dioTimes(i) > 0,1,'first');
    %pull appropriate velocity data
    velRaster(:,i) = s.RotaryData.BinaryLocomotion([targetInd+jumpsBack:targetInd+jumpsForward]);
end

%make average trace:
averageVel = mean(velRaster,2);
velSTD = std(velRaster,0,2);
velVector = [s.Parameters.RasterWindow(1):s.Parameters.InterpolationStepRotary:s.Parameters.RasterWindow(2)];
velDispVector = [s.Parameters.RasterWindow(1):1:s.Parameters.RasterWindow(2)];
velDispIndex = [1:round(1/s.Parameters.InterpolationStepRotary):(jumpsForward-jumpsBack+1)];
velZero = find(velVector >= 0,1,'first');

%% Find locomotion trials
%find trials that have locomotion during the laser period.

%need to redefine laser period as time points on the velRaster system. 
firstPoint = find(velVector >= laserPeriod(1),1,'first');
secondPoint = find(velVector >= laserPeriod(2),1,'first');
laserPeriod = [firstPoint secondPoint];
locoThresh = (secondPoint-firstPoint)*locoThresh;

locoTrial = zeros(totalTrialNum,1);

for i = 1:totalTrialNum
    findTest = sum(velRaster(laserPeriod(1):laserPeriod(2),i));
    if findTest >locoThresh
        locoTrial(i) = 1;
    else
        locoTrial(i) = 0;
    end
end

findLoco = find(locoTrial == 1);

%% Find if there is a preference for locomotion start/stop during laser
%now see if there is a preference for locomotion during laser periods. Do
%this for both starts and stops.


%starts first
[X,Y] = meshgrid(s.RotaryData.LocoStarts,dioTimes);
testSub = X-Y;
testFind = find(testSub<windowPref(2) & testSub > windowPref(1));
timeRange = s.TimeFilterRange(2) - s.TimeFilterRange(1);
timeRatio = length(dioTimes)*(windowPref(2)-windowPref(1))/timeRange;
expNum = length(s.RotaryData.LocoStarts)*timeRatio;

s.RotaryData.PreferenceWindow = windowPref;
s.RotaryData.StartPreference = length(testFind);

%now ends
[X,Y] = meshgrid(s.RotaryData.LocoEnds,dioTimes);
testSub = X-Y;
testFind = find(testSub<windowPref(2) & testSub > windowPref(1));
timeRange = s.TimeFilterRange(2) - s.TimeFilterRange(1);
timeRatio = length(dioTimes)*(windowPref(2)-windowPref(1))/timeRange;
expNum = length(s.RotaryData.LocoStarts)*timeRatio;

s.RotaryData.EndPreference = length(testFind);
%generate distribution by simulating data by randomly pulling from velocity
%times

prefStore = zeros(prefReps,1);

for i = 1:prefReps
    timePoints = s.RotaryData.Velocity(randsample(length(s.RotaryData.Velocity),totalTrialNum),1);
    [X,Y] = meshgrid(s.RotaryData.LocoStarts,timePoints);
    testSub = X-Y;
    testFind = find(testSub<windowPref(2) & testSub > windowPref(1));
    prefStore(i) = length(testFind);
end
%store into structured array
s.RotaryData.PreferenceDistribution = prefStore;

%% Pull EDR Data
% %see if EDR data exists
% fileNames = dir(homeFolder);
% fileNames = {fileNames.name};
% targetFileFinder = strfind(fileNames,'.EDR'); %examines names for D1
% targetFileFinder = find(~cellfun(@isempty,targetFileFinder));%extracts index of correct file
% if ~isempty(targetFileFinder)
%     %flip toggle for graphing
%     edrToggle = 1;
%     %pull filename
%     targetFile = fileNames{targetFileFinder};%pulls out actual file name
%     %extract data
%     [s] = functionEDRPull(s,targetFile,s.Parameters.EDRTimeCol,s.Parameters.EDRTTLCol,s.Parameters.EDRPiezoCol);
%     %downsample for sake of space and time
%     s.EDR.NewData = downsample(s.EDR.FullData,s.Parameters.EDRdownsamp);
%     %confirm ttls line up
%     if length(s.EDR.TTLTimes) ~= totalTrialNum
%         error('Incorrect match of EDR TTLs and trial number')
%     end
%     %figure out raster window
%     edrTimeStep = mean(diff(s.EDR.FullData(:,s.Parameters.EDRTimeCol)));
%     jumpsBack = round(s.Parameters.RasterWindow(1)/(edrTimeStep*s.Parameters.EDRdownsamp));
%     jumpsForward = round(s.Parameters.RasterWindow(2)/(edrTimeStep*s.Parameters.EDRdownsamp));
%     edrRaster = zeros(jumpsForward-jumpsBack+1,totalTrialNum);
%     for i = 1:totalTrialNum
%         %find time closest to actual stim time
%         targetInd = find(s.EDR.NewData(:,s.Parameters.EDRTimeCol) - s.EDR.TTLTimes(i) > 0,1,'first');
%         edrRaster(:,i) = s.EDR.NewData([targetInd + jumpsBack:targetInd + jumpsForward],s.Parameters.EDRPiezoCol);
%     end
%     
%     %calculate overall mean
%     edrMean = mean(edrRaster');
%     edrAbsMean = mean(abs(edrRaster'));
%     edrVector = [s.Parameters.RasterWindow(1):edrTimeStep*s.Parameters.EDRdownsamp:s.Parameters.RasterWindow(2)];
%     edrZero = find(edrVector >= 0,1,'first');
%     
%     edrDispVector = [s.Parameters.RasterWindow(1):1:s.Parameters.RasterWindow(2)];
%     edrDispIndex = [1:round(1/(edrTimeStep*s.Parameters.EDRdownsamp)):(jumpsForward-jumpsBack+1)];
%     
% else
%     disp('NO EDR FILE FOUND')
%     edrToggle = 0;
% end


%% Process spiking information
for i = 1:numUnits
    %allocate a bunch of empty arrays.
    fullHistData = zeros(histBinNum,1);
    
    averageSpikeHolder = zeros(totalTrialNum,1);
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    alignTimes = dioTimes;
    
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    %calculates rasters based on spike information. 
    [rasters] = functionBasicRaster(spikeTimes,alignTimes,s.Parameters.RasterWindow);
    fullRasterData = rasters;
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    for k = 1:totalTrialNum
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > s.Parameters.BaselineBin(1) & rasters(:,1) < s.Parameters.BaselineBin(2)),1);
    end
    
    averageRate = mean(averageSpikeHolder/(s.Parameters.BaselineBin(2)-s.Parameters.BaselineBin(1)));
    averageSTD = std(averageSpikeHolder/(s.Parameters.BaselineBin(2)-s.Parameters.BaselineBin(1)));
    averageSTE = averageSTD/(sqrt(totalTrialNum-1));
    
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
    
    %bin spikes within pre/post period for each trial and store
    prepostStore = zeros(totalTrialNum,2);
    preInd = find(histBinVector > -s.Parameters.PrePost,1,'first');
    zeroInd = find(histBinVector > 0,1,'first');
    postInd = find(histBinVector > s.Parameters.PrePost,1,'first');
    for k = 1:totalTrialNum
        prePostStore(k,1) = sum(fullHistHolder(preInd:zeroInd,k));
        prePostStore(k,2) = sum(fullHistHolder(zeroInd:postInd,k));
    end
    
    [histCounts histCenters] = hist(rasters(:,1),histBinVector);
    fullHistData = histCounts'/totalTrialNum/s.Parameters.histBin;
    
        disp(strcat('Raster and Histogram Extraction Complete: Unit ',num2str(i),' Out of ',num2str(numUnits)))
    s.(desigNames{i}).AllRasters = fullRasterData;
    s.(desigNames{i}).AllHistograms = fullHistData;
    s.(desigNames{i}).IndividualHistograms = fullHistHolder; 
    s.(desigNames{i}).HistogramStandardDeviation = histSTE;
    s.(desigNames{i}).AverageRate = averageRate;
    s.(desigNames{i}).AverageSTD = averageSTD;
    s.(desigNames{i}).AverageSTE = averageSTE;
    s.(desigNames{i}).SessionFiring = sessionFiring;
    s.(desigNames{i}).HistBinVector = histBinVector;
    s.(desigNames{i}).PrePostBins = prePostStore;
    if toggleROC == 1
        targetName = desigNames{i};
        [s] = functionLocomotionROC(s,targetName);
    end
end

%% Put data into summary values for plotting as an initial summary figure
%first, lets get the average velocity trace
s.SumPlot.AveVel = averageVel;
s.SumPlot.AveVelSTD = velSTD;
s.SumPlot.AveVelVector = velVector;
%next, pull hist bin vector
s.SumPlot.HistBinVector = histBinVector;
%calculate the smoothing factor based on histogram bin size
smoothFact = round(s.Parameters.SumGraphSmooth/s.Parameters.histBin);

for i = 1:numUnits
    %get baseline firing rate
    s.SumPlot.BaselineRates(i) = s.(desigNames{i}).AverageRate;
    %check for AUC. If present, extract. 
    if isfield(s.(desigNames{i}),'TrueAUC');
        s.SumPlot.AUC(i) = s.(desigNames{i}).TrueAUC;
        s.SumPlot.AUCLims(i,1) = prctile(s.(desigNames{i}).ShuffleAUC,0.05);
        s.SumPlot.AUCLims(i,2) = prctile(s.(desigNames{i}).ShuffleAUC,99.95);
        if s.(desigNames{i}).TrueAUC > s.SumPlot.AUCLims(i,2) | s.(desigNames{i}).TrueAUC < s.SumPlot.AUCLims(i,1)
            s.SumPlot.AUCSig(i) = 1;
        else
            s.SumPlot.AUCSig(i) = 0;
        end
    end
    %extract overall histogram and peak normalize
    histMin = min(s.(desigNames{i}).AllHistograms);
    histMax = max(s.(desigNames{i}).AllHistograms);
    s.SumPlot.AllHists(i,:) = s.(desigNames{i}).AllHistograms;
    s.SumPlot.SmoothNormHist(i,:) = smooth((s.(desigNames{i}).AllHistograms - histMin)/(histMax-histMin),smoothFact);
    s.SumPlot.SmoothZHist(i,:) = smooth(zscore(s.(desigNames{i}).AllHistograms),smoothFact);
end
s.SumPlot.AverageFiringAllUnits = mean(s.SumPlot.AllHists);
s.SumPlot.AverageSTDAllUnits = std(s.SumPlot.AllHists);
[B,I] = sort(s.SumPlot.BaselineRates);


%calculate and plot LFP information
% [lfpStruct] = functionLFPaverage(master, s.Parameters.LFPWindow, s,homeFolder,fileName, 1, 1, 1, 1);
% s.LFP = lfpStruct;

%% Plotting
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
%first plot summary figure
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
%plot histogram of baseline rates
subplot(3,3,1)
hist(s.SumPlot.BaselineRates,50)
title('Histogram of Baseline Firing Rates')
%if there is auc, plot AUC vs firing rates
if toggleROC == 1
    subplot(3,3,4)
    plot(s.SumPlot.BaselineRates,s.SumPlot.AUC,'k.')
    %plot significant points
    hold on
    s.SumPlot.SignificantAUCs = find(s.SumPlot.AUCSig == 1);
    plot(s.SumPlot.BaselineRates(s.SumPlot.SignificantAUCs),s.SumPlot.AUC(s.SumPlot.SignificantAUCs),'ro')
    title('Scatter Plot of Baseline Rate vs AUC')
end
%plot normalized average velocity and average firing rate
subplot(3,3,2)
hold on
%plot normalized velocity with error lines
plot(s.SumPlot.AveVelVector,(s.SumPlot.AveVel-min(s.SumPlot.AveVel))/(max(s.SumPlot.AveVel)-min(s.SumPlot.AveVel)),'b','LineWidth',2)
plot(s.SumPlot.AveVelVector,(s.SumPlot.AveVel-s.SumPlot.AveVelSTD-min(s.SumPlot.AveVel))/(max(s.SumPlot.AveVel)-min(s.SumPlot.AveVel)),'b','LineWidth',1)
plot(s.SumPlot.AveVelVector,(s.SumPlot.AveVel+s.SumPlot.AveVelSTD-min(s.SumPlot.AveVel))/(max(s.SumPlot.AveVel)-min(s.SumPlot.AveVel)),'b','LineWidth',1)
%plot out average firing rate in the same manner
plot(s.SumPlot.HistBinVector,(s.SumPlot.AverageFiringAllUnits-min(s.SumPlot.AverageFiringAllUnits))/(max(s.SumPlot.AverageFiringAllUnits)-min(s.SumPlot.AverageFiringAllUnits)),'r','LineWidth',2)
plot(s.SumPlot.HistBinVector,(s.SumPlot.AverageFiringAllUnits-s.SumPlot.AverageSTDAllUnits-min(s.SumPlot.AverageFiringAllUnits))/(max(s.SumPlot.AverageFiringAllUnits)-min(s.SumPlot.AverageFiringAllUnits)),'r','LineWidth',1)
plot(s.SumPlot.HistBinVector,(s.SumPlot.AverageFiringAllUnits+s.SumPlot.AverageSTDAllUnits-min(s.SumPlot.AverageFiringAllUnits))/(max(s.SumPlot.AverageFiringAllUnits)-min(s.SumPlot.AverageFiringAllUnits)),'r','LineWidth',1)
title({strcat(fileName,desigNames{i});'Normalized Velocity (b) and Firing (r)'},'interpreter','none')
%plot out heatmaps of peak normalized firing rates
subplot(3,3,5)
imagesc(s.SumPlot.SmoothNormHist(I,:))
colormap(parula)
set(gca,'XTick',[0:round(length(histBinVector)/(s.Parameters.RasterWindow(2) - s.Parameters.RasterWindow(1))):length(histBinVector)])
set(gca,'XTickLabel',[s.Parameters.RasterWindow(1):1:s.Parameters.RasterWindow(2)])
title('Normalized Histograms Ordered by Firing Rate')
subplot(3,3,8)
imagesc(s.SumPlot.SmoothZHist(I,:),cLims)
colormap(parula)
colorbar
set(gca,'XTick',[0:round(length(histBinVector)/(s.Parameters.RasterWindow(2) - s.Parameters.RasterWindow(1))):length(histBinVector)])
set(gca,'XTickLabel',[s.Parameters.RasterWindow(1):1:s.Parameters.RasterWindow(2)])
title('Z Histograms Ordered by Firing Rate')
spikeGraphName = strcat(fileName,'DATStimSummaryFigure');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%next plot individual units
for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    %plots average waveform
    subplot(4,4,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    title({fileName;desigNames{i};strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate))})
    %plots ISI
    subplot(4,4,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    
    %plot velocity vs firing rate
    subplot(3,2,2)
    hold on
    plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'b')
    plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
    xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
    ylim([-0.1,1])
    if toggleROC == 1
        title(strcat('Vel & Firing Rate. AUC:',num2str(s.(desigNames{i}).TrueAUC),'99%Range',num2str(prctile(s.(desigNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(desigNames{i}).ShuffleAUC,1))))
    else
        title('Vel & Firing Rate')
    end
    % plot histogram.
    subplot(3,2,4)
    plot(histBinVector,s.(desigNames{i}).AllHistograms,'k','LineWidth',2)
    hold on
    plot(histBinVector,s.(desigNames{i}).AllHistograms - s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
    plot(histBinVector,s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
    %fill in laser onset
    plot([0 0],[ylim],'b');
    %if there are locomotion trials, plot these as well.
    if length(findLoco) > 5
        plot(histBinVector,mean(s.(desigNames{i}).IndividualHistograms(:,findLoco)'/s.Parameters.histBin),'r')
    end
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    title('Histogram')
    
    %plots rasters (chronological)
    subplot(3,2,6)
    plot(s.(desigNames{i}).AllRasters(:,1),...
        s.(desigNames{i}).AllRasters(:,2),'k.','markersize',5)
    hold on
    ylim([0 totalTrialNum])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    set(gca,'Ydir','reverse')
    plot([0 0],[ylim],'b');
    title({fileName;desigNames{i}},'fontweight','bold')
    set(0, 'DefaulttextInterpreter', 'none')
    
    %plot heatmap of firing
    subplot(4,2,3)
    imagesc(s.(desigNames{i}).IndividualHistograms')
    colormap(parula)
    set(gca,'XTick',[1:(1/s.Parameters.histBin):(size(s.(desigNames{i}).IndividualHistograms,2))]);
    set(gca,'XTickLabel',[s.Parameters.RasterWindow(1):1:s.Parameters.RasterWindow(2)]);
    colorbar
    title('Colormap of Firing')
    
    %plot heatmap of locomotion
    subplot(4,2,5)
    imagesc((velRaster'))
    colormap(parula)
    colorbar
    set(gca,'XTick',velDispIndex);
    set(gca,'XTickLabel',velDispVector);
    title('Colorized Velocity Trace Per Trial')
    
    %plot scatter of pre vs post bins
    subplot(4,2,7)
    plot(s.(desigNames{i}).PrePostBins(:,1),s.(desigNames{i}).PrePostBins(:,2),'k.')
    plotMax = max(max(s.(desigNames{i}).PrePostBins));
    hold on
    plot([0 plotMax],[0 plotMax],'b')
    title(strcat('Scatter Plot of Pre (x) and Post (y) Firing for',num2str(s.Parameters.PrePost),'Second Bins'))
    %plot edr data, if exists
%     if edrToggle == 1
%         subplot(4,2,7)
%         hold on
%         imagesc(flipud(edrRaster'),[-0.004 0.004])
%         xlim([0 size(edrRaster,1)])
%         ylim([0 size(edrRaster,2)])
%         colorbar
%         set(gca,'XTick',edrDispIndex);
%         set(gca,'XTickLabel',edrDispVector);
%         set(gca,'YTick',[1:10:totalTrialNum])
%         set(gca,'YTickLabel',[totalTrialNum:-10:1])
%         title('Colorized Piezo Data')
%     else
% %         subplot(4,2,7)
% %         hold on
% %         memFinder = ismember(s.(desigNames{i}).AllRasters(:,2),findLoco);
% %         memFinder = s.(desigNames{i}).AllRasters(memFinder,1);
% %         hist(memFinder,[s.Parameters.RasterWindow(1):0.1:s.Parameters.RasterWindow(2)])
% %         title('Histogram of Responses During Locomotion Trials')
%     end
    
    
    
    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'DATStimAnalysis');
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