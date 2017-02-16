%This code is meant for the basic analysis of data from the stimulation
%protocol 

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

function [s] = analysisLaserLag(fileName);
%% Constants and things you might want to tweak

%lets set some switches to toggle things on and off.
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 1; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.


s.Parameters.RasterWindow = [-1 3]; %ratio for raster window. will be multiplied by toneDur
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate
s.Parameters.DefaultBins = 0.001;% bin size for calculating significant responses
s.Parameters.SmoothingBins = [0.01 0.001];%bins for smoothing
s.Parameters.CalcWindow = [0 2]; %window for calculating significant responses
s.Parameters.zLimit = 3; %zlimit for calculating significant responses
% s.Parameters.FirstSpikeWindow = [0 0.5 1 1.5]; %ratios! need to be multiplied by tone duration.
s.Parameters.FirstSpikeWindow = [0 1];
s.Parameters.ChosenSpikeBin = 1; %delineates which spike window I will graph.
s.Parameters.BaselineBin = [-1 0]; %ratio for bin from which baseline firing rate will be calculated
s.Parameters.LFPWindow = [-1 3];

%stuff for significance
s.Parameters.calcWindow = [0 2]; %defines period for looking for responses, based on toneDur
s.Parameters.zLimit = [0.05 0.01 0.001];
s.Parameters.numShuffle = 1000;
s.Parameters.firstSpikeWindow = [-1 0 0.5 1 1.5];%defines period for looking for first spike, based on toneDur
s.Parameters.chosenSpikeBin = 1; %spike bin selected in binSpike (in the event of multiple spike bins)
s.Parameters.minSpikes = 100; %minimum number of spikes to do spike shuffling
s.Parameters.minSigSpikes = 2; %minimum number of significant points to record a significant response.
s.Parameters.BaselineWindow = [-0.4 0]; %window for counting baseline spikes, in SECONDS. NOTE THIS IS DIFFERENT FROM RASTER WINDOW
s.Parameters.BaselineCalcBins = 1; %bin size in seconds if there are insufficient baseline spikes to calculate a baseline rate.



%for latency calculations
s.Parameters.ToneWindow = [0 1];
s.Parameters.GenWindow = [0 3];
s.Parameters.latBin = 0.001;
s.Parameters.PercentCutoff = 99.9;
s.Parameters.BaselineCutoff = 95;
s.Parameters.ThresholdHz = 4; %minimum response in Hz to be counted as significant.

%for duplicate elimination
s.Parameters.DownSampFactor = 3; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr




%% sets up file saving stuff
saveName = strcat(fileName,'LaserLagAnalysis','.mat');
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

disp('Spikes and Waveforms Allocated into Structured Array')

%eliminate duplicates.

if toggleDuplicateElimination ==1
    if length(s.DesignationName) > 1
        disp('Now Selecting Based on xCORR')
        [s] = functionDuplicateElimination(s,s.Parameters.DownSampFactor,...
            s.Parameters.corrSlide,s.Parameters.ThresholdComparison,s.Parameters.trodesFS);
    end
else
    disp('NOT EXECUTING DUPLICATE ELIMINATION')
end

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%% Extracts Sound Data from soundFile
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;

%pull important sound and trial information
targetFreq = unique(soundData.Frequencies);
targetDB = unique(soundData.dBs);

%store these in structured array
s.SoundData.TargetFreq = targetFreq;
s.SoundData.TargetDB = targetDB;

toneDur = soundData.ToneDuration;
toneReps = soundData.ToneRepetitions;
totalTrialNum = length(soundData.Frequencies);
numLags = length(soundData.LaserLags);
%Recalculate raster window and such
s.Parameters.RasterWindow = s.Parameters.RasterWindow * toneDur;
s.Parameters.FirstSpikeWindow = s.Parameters.FirstSpikeWindow * toneDur;
s.Parameters.BaselineBin = s.Parameters.BaselineBin * toneDur;
s.Parameters.ToneWindow = s.Parameters.ToneWindow * toneDur;
s.Parameters.GenWindow = s.Parameters.GenWindow * toneDur;
calcWindow = s.Parameters.calcWindow*toneDur;
rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes
s.Parameters.CalcWindow = s.Parameters.CalcWindow * toneDur;
s.Parameters.LFPWindow = s.Parameters.LFPWindow * toneDur;

%Sets up master array so that we have a size comparison for DIO
%information. 
master = zeros(totalTrialNum,6);

%master(:,1) is reserved for actual input times
%master(:,2) is frequency
master(:,2) = soundData.Frequencies;
%master(:,3) is dB
master(:,3) = soundData.dBs;
%master(:,4) is trial number (chronological)
master(:,4) = 1:1:totalTrialNum;
%master(:,6) is the designation of which lag to use.
master(:,6) = soundData.LaserDesig;
%master(:,5) is trial num, arranging trials in order from small dB to large
%dB, and low freq to high freq. frequency is larger category.
sortingCounter = 1;
for i = 1:length(soundData.LaserLags);
    sortFinder = find(master(:,6) == i);
    master(sortFinder,5) = sortingCounter:1:sortingCounter+size(sortFinder,1)-1;
    sortingCounter = sortingCounter + size(sortFinder,1);
end
%master(:,7) is reserved for the time of the laser pulse. IF NO laser
%pulse, then master(:,7) will equal 0. 

%% Extract DIO information. 

%extract DIO filenames
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};

%extract actual data
[DIO1Data] = readTrodesExtractedDataFile(D1FileName);
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%remove just time and state data. convert to double for easier processing
DIO1Data = double([DIO1Data.fields(1).data,DIO1Data.fields(2).data]);
DIO2Data = double([DIO2Data.fields(1).data,DIO2Data.fields(2).data]);

%pull all differences in state for DIO data (diff subtracts x(2)-x(1)). +1
%fudge factor adjusts for indexing. MUST BE DOUBLE!!
DIO1Diff = find(diff(DIO1Data(:,2))==1)+1;
DIO1High = find(DIO1Data(:,2) == 1);
%finds all points which are down to up states, extracts these times.
DIO1True = intersect(DIO1Diff,DIO1High);
DIO1True = DIO1Data(DIO1True,1)/s.Parameters.trodesFS;
%finds differences between time points
DIO1TrueDiff = diff(DIO1True);

%pull all differences in state for DIO data (diff subtracts x(2)-x(1)). +1
%fudge factor adjusts for indexing. MUST BE DOUBLE!!
DIO2Diff = find(diff(DIO2Data(:,2))==1)+1;
DIO2High = find(DIO2Data(:,2) == 1);
%finds all points which are down to up states, extracts these times.
DIO2True = intersect(DIO2Diff,DIO2High);
DIO2True = DIO2Data(DIO2True,1)/s.Parameters.trodesFS;
%finds differences between time points
DIO2TrueDiff = diff(DIO2True);

disp('DIO Signals Processed')

%now I need to assign the tone to the right TTLs

%exclude laser TTLs from tone TTLs. Laser TTLs are always doublets, so we
%can pick out the values from the doublets and remove those into a separate
%array
findLaser = find(DIO1TrueDiff< 0.8*min(soundData.LaserLags(soundData.LaserLags>0)));
laserArray = zeros(length(findLaser),2);
laserArray(:,1) = findLaser;
laserArray(:,2) = findLaser + 1;
laserArray = reshape(laserArray,[],1);

toneTTLs = DIO1True;
toneTTLs(laserArray) = [];
toneTTLs(diff(toneTTLs) < 0.7*min(soundData.Delays)) = [];

laserTTLs = DIO1True(findLaser);

%find laser end points. 
laserEndFind = find(DIO2TrueDiff > max(soundData.LaserLags));
laserEndTTLs = DIO2True(laserEndFind);
laserEndTTLs(end+1) = DIO2True(end);


%Check to make sure we have the appropriate number of toneTTLs and laser
%TTLs
corrLaserNum = length(soundData.LaserLags(soundData.LaserLags>0))*toneReps;
corrToneNum = length(soundData.LaserLags)*toneReps;

if length(toneTTLs) ~= corrToneNum
    disp(length(toneTTLs)-length(corrToneNum))
    error('Incorrect Number of Tones Detected.')
elseif length(laserTTLs) ~= corrLaserNum
    disp(length(laserTTLs)-length(corrLaserNum))
    error('Incorrect Number of Laser Deliveries Detected')
elseif length(laserEndTTLs) ~= length(laserTTLs)
    disp(length(laserEndTTLs)-length(laserTTLs))
    error('Incorrect Matchup of Laser Ends and Laser Starts')
end

s.TTLs.ToneTTLs = toneTTLs;
s.TTLs.LaserInitiationTTLs = laserTTLs;
s.TTLs.LaserEndTTLs = laserEndTTLs;

master(:,1) = toneTTLs;

%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db
for i = 1:numUnits
    %allocate a bunch of empty arrays.
    fullHistData = zeros(histBinNum,1);
    
    organizedRasters = cell(numLags,1); %allocates cell array for rasters
    organizedHist = zeros(numLags,histBinNum,1);
    
    freqSpecHist = zeros(histBinNum,1);
    
    firstSpikeTimeHolder = cell(numLags,1);
    firstSpikeStatsHolder = zeros(numLags,4,size(s.Parameters.FirstSpikeWindow,2)-1);
    
    binSpikeHolder = cell(numLags,1);
    binSpikeStatsHolder = zeros(numLags,2,size(s.Parameters.FirstSpikeWindow,2)-1);
    
    averageSpikeHolder = zeros(totalTrialNum,1);
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    alignTimes = master(:,1);
    
    %calculates rasters based on spike information. 
    [rasters] = functionBasicRaster(spikeTimes,alignTimes,s.Parameters.RasterWindow);
    rasters(:,3) = master(rasters(:,2),5); %adds information about frequency/amplitude
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
    
    [histCounts histCenters] = hist(rasters(:,1),histBinVector);
    fullHistData = histCounts'/totalTrialNum/s.Parameters.histBin;
    
    %calculate significant response for combined histogram of all responses
    %to all sounds.    
    [generalResponseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,length(master));

    disp(strcat('Baseline Spikes:',num2str(generalResponseHist.SpikeNumber),' Unit:',(desigNames{i})))
    s.BaselineSpikes = generalResponseHist.SpikeNumber;
    if generalResponseHist.Warning == 0 & generalResponseHist.SigSpike == 1
        s.SignificantSpikes(i) = 1;
    end
    
    %allocates empty array.
    organizedHist = zeros(length(soundData.LaserLags),length(histBinVector));
    organizedRasters = cell(length(soundData.LaserLags),1);
    responseHistHolder = cell(length(soundData.LaserLags),1);
    histErr = zeros(length(soundData.LaserLags),length(histBinVector));
    
    for k = 1:length(soundData.LaserLags)
        targetTrials = master(master(:,6) == k,4);
        findMatches = find(ismember(fullRasterData(:,2),targetTrials)); %uses these trial numbers to pull raster index
        targetRasters = fullRasterData(findMatches,:); %extracts target rasters from all rasters using the above index.
        organizedRasters{k} = targetRasters; %saves to organized rasters
        [histCounts,histCenters] = hist(targetRasters(:,1),histBinVector); %calculates histogram with defined bin size
        organizedHist(k,:) = histCounts/toneReps/s.Parameters.histBin; %saves histogram
        specHist = fullHistHolder(:,targetTrials);
        histErr(k,:) = std(specHist,0,2)/toneReps/s.Parameters.histBin;
        [firstSpikeTimes,firstSpikeStats,binSpikeTimes,binSpikeStats] = ...
            functionBasicFirstSpikeTiming(s.Parameters.FirstSpikeWindow,targetRasters,toneReps,2,targetTrials); %calculates information about first spike timing
        [latPeakBinOut] = functionLatPeakBinCalculation(s.Parameters.ToneWindow,s.Parameters.GenWindow,s.Parameters.RasterWindow,...
                targetRasters,length(targetTrials),2,targetTrials,s.Parameters.latBin,s.Parameters.histBin,s.Parameters.PercentCutoff,s.Parameters.BaselineCutoff);
        latStore{k} = latPeakBinOut;
        firstSpikeTimeHolder{k} = firstSpikeTimes; %saves first spike times
        firstSpikeStatsHolder(k,:,:) = firstSpikeStats; %saves statistics about first spikes
        binSpikeHolder{k} = binSpikeTimes; %binned spikes from the defined window.
        binSpikeStatsHolder(k,:,:) = binSpikeStats; %stats about those spikes
        [responseHist] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,toneReps);
        responseHistHolder{k} = responseHist;
    end

        disp(strcat('Raster and Histogram Extraction Complete: Unit ',num2str(i),' Out of ',num2str(numUnits)))
    s.(desigNames{i}).AllRasters = fullRasterData;
    s.(desigNames{i}).AllHistograms = fullHistData;
    s.(desigNames{i}).IndividualHistograms = fullHistHolder; 
    s.(desigNames{i}).HistogramStandardDeviation = histSTE;
    s.(desigNames{i}).LagRasters = organizedRasters;
    s.(desigNames{i}).LagHistograms = organizedHist;
    s.(desigNames{i}).LagHistogramErrors = histErr;
    s.(desigNames{i}).FirstSpikeTimes = firstSpikeTimeHolder;
    s.(desigNames{i}).FirstSpikeStats = firstSpikeStatsHolder;
    s.(desigNames{i}).BinSpikes = binSpikeHolder;
    s.(desigNames{i}).BinSpikeStats = binSpikeStatsHolder;
    s.(desigNames{i}).AverageRate = averageRate;
    s.(desigNames{i}).AverageSTD = averageSTD;
    s.(desigNames{i}).AverageSTE = averageSTE;
    s.(desigNames{i}).HistBinVector = histBinVector;
    s.(desigNames{i}).AllHistogramSig = generalResponseHist;
    s.(desigNames{i}).SpecHistogramSig = responseHistHolder;
    s.(desigNames{i}).LatencyData = latStore;
end

%calculate and plot LFP information
% [lfpStruct] = functionLFPaverage(master, s.Parameters.LFPWindow, s,homeFolder,fileName, uniqueFreqs, uniqueDBs, numFreqs, numDBs);
% s.LFP = lfpStruct;

%% Plotting
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.05 0.05], [0.05 0.01]);

if toggleTuneSelect == 1 %if you want tuning selection...
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    decisionTuning = zeros(numUnits,1);
    for i = 1:numUnits
        %plots average waveform
        subplot(4,6,1)
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
        %plot overall histogram
        subplot(4,3,4)
        plot(histBinVector,s.(desigNames{i}).AllHistograms,'k','LineWidth',2)
        hold on
        plot(histBinVector,s.(desigNames{i}).AllHistograms - s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
        plot(histBinVector,s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
        %plot significant values
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1,1),...
            'b*')
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1,1),...
            'c*')
        %plot negative values for first tuning
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
            'k*')
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Histogram')
        hold off
        %plot by lag
        subplot(4,3,7)
        hold on
        for k = 1:numLags
            plot(histBinVector,s.(desigNames{i}).LagHistograms(k,:),'LineWidth',2,'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
            plot(histBinVector,s.(desigNames{i}).LagHistograms(k,:)-s.(desigNames{i}).LagHistogramErrors(k,:),'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
            plot(histBinVector,s.(desigNames{i}).LagHistograms(k,:)+s.(desigNames{i}).LagHistogramErrors(k,:),'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
        end
        title('Histograms by Lag (later is Redder)')
        %plot zoomed in version
        subplot(4,3,10)
        hold on
        %figure out tone period
        lagHistVector = s.(desigNames{i}).LatencyData{1,1}.LatBinVector;
        timeZero = find(lagHistVector>0,1,'first');
        timeToneEnd = find(lagHistVector>s.SoundData.ToneDuration,1,'first');
        for k = 1:numLags
            plot(lagHistVector(timeZero:timeToneEnd),s.(desigNames{i}).LatencyData{k}.LatHist(timeZero:timeToneEnd),'LineWidth',2,'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
        end
        xlim([lagHistVector(timeZero) lagHistVector(timeToneEnd)])
        title('Tone Period Histograms by Lag (later is Redder)')
        %plots rasters (chronological)
        subplot(3,3,2)
        hold on
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,2),'k.','markersize',4)
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        title({fileName;desigNames{i}},'fontweight','bold')
        set(0, 'DefaulttextInterpreter', 'none')
        %plots rasters (organized by lag)
        subplot(3,3,5)
        hold on
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,3),'k.','markersize',4)
        %this generates green lines separating by Frequency
        for k = 1:numLags
            plot([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)],[toneReps*k toneReps*k],'g')
        end
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        rasterLagLines = zeros(numLags,2);
        rasterLagLines(:,1) = toneReps/2:toneReps:totalTrialNum;
        rasterLagLines(:,2) = soundData.LaserLags;

        set(gca,'YTick',rasterLagLines(:,1));
        set(gca,'YTickLabel',rasterLagLines(:,2));
        set(gca,'Ydir','reverse')
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Descending = increase in amplitude and freq')
        %plot heatmap organized by laser lag
        subplot(3,3,8)
        imagesc(s.(desigNames{i}).LagHistograms(:,:))
        colorbar
    %     set(gca,'YTick',octaveRange(:,2));
    %     set(gca,'YTickLabel',octaveRange(:,1));
    %     set(gca,'XTick',[1:10:size(histBinVector,2)]);
    %     set(gca,'XTickLabel',histBinVector(1:20:end));
        histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
        histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),toneDur);
        line([histBinZero histBinZero],[0 numLags],'LineWidth',3,'Color','green')
        line([histBinZero histBinZero],[0 numLags],'LineWidth',2,'Color','black')
        line([histBinTone histBinTone],[0 numLags],'LineWidth',3,'Color','green')
        line([histBinTone histBinTone],[0 numLags],'LineWidth',2,'Color','black')
    %         title('Heatmap by Frequency and Time Max')
        title('Frequency Arranged Heatmap')

        %plot timecourse of responses by binning
        subplot(3,3,3)
        hold on

        for k = 1:numLags
            plot(s.(desigNames{i}).LatencyData{k}.BinnedSpikesGen,'ko','MarkerEdgeColor',[(k-1)/numLags (numLags-k+1)/numLags 0])
            plot(smooth(s.(desigNames{i}).LatencyData{k}.BinnedSpikesGen,11),'LineWidth',2,'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
        end
        xlim([0 length(s.(desigNames{i}).LatencyData{k}.BinnedSpikesGen)])
        title('Responses to Different Lags Over Time')

        hold off
        spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
        
        %ask for input! 
        promptCounter = 1; %This is used to run the while loop.
        whileCounter = 0; %this is the counter that gets updated to exit the loop

        while whileCounter ~= promptCounter
            try
                prompt = 'Is this unit tuned? (y/n)';
                str = input(prompt,'s');
                if str~='n' & str~='y'
                    error
                else
                    whileCounter = 1;
                end
            catch
            end
        end
        if strfind(str,'y')
            decisionTuning(i) = 1;
        elseif strfind(str,'n')
            decisionTuning(i) = 0;
        end

        %clear figure.
        clf
    end
    s.TuningDecision = decisionTuning;
else
    for i = 1:numUnits
        hFig = figure;
        set(hFig, 'Position', [10 80 1240 850])
        %plots average waveform
        subplot(4,6,1)
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
        %plot overall histogram
        subplot(4,3,4)
        plot(histBinVector,s.(desigNames{i}).AllHistograms,'k','LineWidth',2)
        hold on
        plot(histBinVector,s.(desigNames{i}).AllHistograms - s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
        plot(histBinVector,s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
        %plot significant values
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1,1),...
            'b*')
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1,1),...
            'c*')
        %plot negative values for first tuning
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
            'k*')
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Histogram')
        hold off
        %plot by lag
        subplot(4,3,7)
        hold on
        for k = 1:numLags
            plot(histBinVector,s.(desigNames{i}).LagHistograms(k,:),'LineWidth',2,'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
            plot(histBinVector,s.(desigNames{i}).LagHistograms(k,:)-s.(desigNames{i}).LagHistogramErrors(k,:),'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
            plot(histBinVector,s.(desigNames{i}).LagHistograms(k,:)+s.(desigNames{i}).LagHistogramErrors(k,:),'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
        end
        title('Histograms by Lag (later is Redder)')
        %plot zoomed in version
        subplot(4,3,10)
        hold on
        %figure out tone period
        lagHistVector = s.(desigNames{i}).LatencyData{1,1}.LatBinVector;
        timeZero = find(lagHistVector>0,1,'first');
        timeToneEnd = find(lagHistVector>s.SoundData.ToneDuration,1,'first');
        for k = 1:numLags
            plot(lagHistVector(timeZero:timeToneEnd),s.(desigNames{i}).LatencyData{k}.LatHist(timeZero:timeToneEnd),'LineWidth',2,'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
        end
        xlim([lagHistVector(timeZero) lagHistVector(timeToneEnd)])
        title('Tone Period Histograms by Lag (later is Redder)')
        %plots rasters (chronological)
        subplot(3,3,2)
        hold on
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,2),'k.','markersize',4)
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        title({fileName;desigNames{i}},'fontweight','bold')
        set(0, 'DefaulttextInterpreter', 'none')
        %plots rasters (organized by lag)
        subplot(3,3,5)
        hold on
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,3),'k.','markersize',4)
        %this generates green lines separating by Frequency
        for k = 1:numLags
            plot([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)],[toneReps*k toneReps*k],'g')
        end
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        rasterLagLines = zeros(numLags,2);
        rasterLagLines(:,1) = toneReps/2:toneReps:totalTrialNum;
        rasterLagLines(:,2) = soundData.LaserLags;

        set(gca,'YTick',rasterLagLines(:,1));
        set(gca,'YTickLabel',rasterLagLines(:,2));
        set(gca,'Ydir','reverse')
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Descending = increase in amplitude and freq')
        %plot heatmap organized by laser lag
        subplot(3,3,8)
        imagesc(s.(desigNames{i}).LagHistograms(:,:))
        colorbar
    %     set(gca,'YTick',octaveRange(:,2));
    %     set(gca,'YTickLabel',octaveRange(:,1));
    %     set(gca,'XTick',[1:10:size(histBinVector,2)]);
    %     set(gca,'XTickLabel',histBinVector(1:20:end));
        histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
        histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),toneDur);
        line([histBinZero histBinZero],[0 numLags],'LineWidth',3,'Color','green')
        line([histBinZero histBinZero],[0 numLags],'LineWidth',2,'Color','black')
        line([histBinTone histBinTone],[0 numLags],'LineWidth',3,'Color','green')
        line([histBinTone histBinTone],[0 numLags],'LineWidth',2,'Color','black')
    %         title('Heatmap by Frequency and Time Max')
        title('Frequency Arranged Heatmap')

        %plot timecourse of responses by binning
        subplot(3,3,3)
        hold on

        for k = 1:numLags
            plot(s.(desigNames{i}).LatencyData{k}.BinnedSpikesGen,'ko','MarkerEdgeColor',[(k-1)/numLags (numLags-k+1)/numLags 0])
            plot(smooth(s.(desigNames{i}).LatencyData{k}.BinnedSpikesGen,11),'LineWidth',2,'Color',[(k-1)/numLags (numLags-k+1)/numLags 0])
        end
        xlim([0 length(s.(desigNames{i}).LatencyData{k}.BinnedSpikesGen)])
        title('Responses to Different Lags Over Time')

        hold off
        spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end

    
%% Saving
save(fullfile(pname,fname),'s');

end