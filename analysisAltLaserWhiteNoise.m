

function [s] = analysisAltLaserWhiteNoise(fileName);
%% Constants and things you might want to tweak
%lets set some switches to toggle things on and off.
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.

s.Parameters.RasterWindow = [-1 3]; %seconds for raster window. 
s.Parameters.ToneWindow = [0 0.5];
s.Parameters.GenWindow = [0 1];
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate
s.Parameters.zLimit = 3; %zlimit for calculating significant responses
s.Parameters.FirstSpikeWindow = [0 1];
s.Parameters.BaselineBin = [-1 0]; %ratio for bin from which baseline firing rate will be calculated
s.Parameters.LFPWindow = [-1 1];

%stuff for significance
s.Parameters.calcWindow = [0 2]; %defines period for looking for responses, based on toneDur
s.Parameters.zLimit = [0.05 0.01 0.001];
s.Parameters.numShuffle = 1000;
s.Parameters.minSpikes = 100; %minimum number of spikes to do spike shuffling
s.Parameters.minSigSpikes = 2; %minimum number of significant points to record a significant response.
s.Parameters.BaselineWindow = [-0.4 0]; %window for counting baseline spikes, in SECONDS. NOTE THIS IS DIFFERENT FROM RASTER WINDOW
s.Parameters.BaselineCalcBins = 1; %bin size in seconds if there are insufficient baseline spikes to calculate a baseline rate.
s.Parameters.PercentCutoff = 99.9;
s.Parameters.BaselineCutoff = 95;
s.Parameters.latBin = 0.001;
s.Parameters.ThresholdHz = 4; %minimum response in Hz to be counted as significant.

%for duplicate elimination
s.Parameters.DownSampFactor = 10; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr



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
    matclustFiles,s,s.Parameters.ClusterWindow);

if toggleDuplicateElimination ==1
    if length(s.DesignationName) > 1
        disp('Now Selecting Based on xCORR')
        [s] = functionDuplicateElimination(s,s.Parameters.DownSampFactor,...
            s.Parameters.corrSlide,s.Parameters.ThresholdComparison,s.Parameters.trodesFS,s.Parameters.RPVTime,s.Parameters.ClusterWindow);
    end
else
    disp('NOT EXECUTING DUPLICATE ELIMINATION')
end

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%calculate window sizes and binning
calcWindow = s.Parameters.calcWindow;
rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes

%% Pull data from sound file
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;
toneDur = soundData.ToneDuration;
toneReps = soundData.ToneRepetitions;
laserLag = min(soundData.LaserLag);

%recalculate some things!


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
shortDIODel = unique([shortDIO,shortDIO+1,shortDIO+2]);
toneDIO = dioTimes;
toneDIO(shortDIODel) = [];
%this isolates tones!
%now figure out which ones are laser
toneLaserDIO = dioTimes(shortDIO+2);
laserDIO = dioTimes(shortDIO);
allToneDIO = unique([toneDIO,toneLaserDIO]);

%check to see if numbers match up
if length(laserDIO) ~= length(toneLaserDIO)
    error('Laser and Tone Laser DIO mismatch')
elseif length(laserDIO) ~= length(toneDIO);
    error('Mismatch of Tones and Laser Presentations')
end

totalTrialNum = 2*toneReps;

%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db
for i = 1:numUnits
    %allocate a bunch of empty arrays.
    fullHistTone = zeros(histBinNum,1);
    fullHistLaserTone = zeros(histBinNum,1);
    
    averageSpikeHolder = zeros(totalTrialNum,1);
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    
    %calculates rasters based on spike information. 
    [rasters] = functionBasicRaster(spikeTimes,allToneDIO,s.Parameters.RasterWindow);
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    for k = 1:totalTrialNum
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > s.Parameters.BaselineBin(1) & rasters(:,1) < s.Parameters.BaselineBin(2)),1);
    end

    averageRate = mean(averageSpikeHolder/(s.Parameters.BaselineBin(2)-s.Parameters.BaselineBin(1)));
    averageSTD = std(averageSpikeHolder/(s.Parameters.BaselineBin(2)-s.Parameters.BaselineBin(1)));
    averageSTE = averageSTD/(sqrt(totalTrialNum-1));
        
    [rasters] = functionBasicRaster(spikeTimes,toneDIO,s.Parameters.RasterWindow);
    fullRasterDataTone = rasters;
    
    [rasters] = functionBasicRaster(spikeTimes,toneLaserDIO,s.Parameters.RasterWindow);
    fullRasterDataToneLaser = rasters;
    
    %calculates histograms of every single trial. 
    fullHistHolderTone = zeros(length(histBinVector),toneReps);
    for k =1:toneReps
        if size(fullRasterDataTone(fullRasterDataTone(:,2) == k,:),1) > 0
            [counts centers] = hist(fullRasterDataTone(fullRasterDataTone(:,2) == k,1),histBinVector);
            fullHistHolderTone(:,k) = counts;
        end
    end
    
    fullHistHolderToneLaser = zeros(length(histBinVector),toneReps);
    for k =1:toneReps
        if size(fullRasterDataToneLaser(fullRasterDataToneLaser(:,2) == k,:),1) > 0
            [counts centers] = hist(fullRasterDataToneLaser(fullRasterDataToneLaser(:,2) == k,1),histBinVector);
            fullHistHolderToneLaser(:,k) = counts;
        end
    end
    
    %calculate standard deviation for each bin
    histSTETone = (std(fullHistHolderTone'))';
    histSTEToneLaser = (std(fullHistHolderToneLaser'))';
    
    [histCounts histCenters] = hist(fullRasterDataTone(:,1),histBinVector);
    overallHistTone = histCounts'/toneReps/s.Parameters.histBin;
    
    [histCounts histCenters] = hist(fullRasterDataToneLaser(:,1),histBinVector);
    overallHistToneLaser = histCounts'/toneReps/s.Parameters.histBin;

    disp(strcat('Raster and Histogram Extraction Complete: Unit ',num2str(i),' Out of ',num2str(numUnits)))
    s.(desigNames{i}).ToneRasters = fullRasterDataTone;
    s.(desigNames{i}).ToneLaserRasters = fullRasterDataToneLaser;
    s.(desigNames{i}).OverallHistTone = overallHistTone;
    s.(desigNames{i}).OverallHistToneLaser = overallHistToneLaser;
    s.(desigNames{i}).IndividualHistogramsTone = fullHistHolderTone; 
    s.(desigNames{i}).IndividualHistogramsToneLaser = fullHistHolderToneLaser; 
    s.(desigNames{i}).HistogramStandardDeviationTone = histSTETone;
    s.(desigNames{i}).HistogramStandardDeviationToneLaser = histSTEToneLaser;
    s.(desigNames{i}).AverageRate = averageRate;
    s.(desigNames{i}).AverageSTD = averageSTD;
    s.(desigNames{i}).AverageSTE = averageSTE;
    s.(desigNames{i}).HistBinVector = histBinVector;
end

%calculate and plot LFP information
% [lfpStruct] = functionLFPaverage(master, s.Parameters.LFPWindow, s,homeFolder,fileName, 1, 1, 1, 1);
% s.LFP = lfpStruct;

%% Plotting
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    %plots average waveform
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    %plots ISI
    subplot(4,6,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    %plot histograms
    subplot(3,3,2)
    plot(histBinVector,s.(desigNames{i}).OverallHistTone,'k','LineWidth',1)
    hold on
    plot(histBinVector,s.(desigNames{i}).OverallHistTone - s.(desigNames{i}).HistogramStandardDeviationTone,'b','LineWidth',1)
    plot(histBinVector,s.(desigNames{i}).OverallHistTone + s.(desigNames{i}).HistogramStandardDeviationTone,'b','LineWidth',1)
    
    plot(histBinVector,s.(desigNames{i}).OverallHistToneLaser,'r','LineWidth',1)
    plot(histBinVector,s.(desigNames{i}).OverallHistToneLaser - s.(desigNames{i}).HistogramStandardDeviationToneLaser,'m','LineWidth',1)
    plot(histBinVector,s.(desigNames{i}).OverallHistToneLaser + s.(desigNames{i}).HistogramStandardDeviationToneLaser,'m','LineWidth',1)
    
    plot([0 0],[ylim],'b');
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    title('Histogram of Tone Response')
    
    %plots rasters (chronological)
    subplot(3,3,5)
    plot(s.(desigNames{i}).ToneRasters(:,1),...
        s.(desigNames{i}).ToneRasters(:,2),'k.','markersize',5)
    hold on
    ylim([0 toneReps])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Tone Response')
    
    subplot(3,3,6)
    plot(s.(desigNames{i}).ToneLaserRasters(:,1),...
        s.(desigNames{i}).ToneLaserRasters(:,2),'k.','markersize',5)
    hold on
    ylim([0 toneReps])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Tone Response With Laser')
    
    subplot(3,3,8)
    title({fileName;desigNames{i}},'fontweight','bold')
    set(0, 'DefaulttextInterpreter', 'none')
    
    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'LightAnalysis');
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