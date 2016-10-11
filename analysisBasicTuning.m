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

function [s] = analysisBasicTuning(fileName);
%% Constants and things you might want to tweak
rasterWindow = [-1 3]; %ratio for raster window. will be multiplied by toneDur
rpvTime = 0.001; %time limit in seconds for consideration as an RPV
clusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
histBin = 0.005; %histogram bin size in seconds
sampleRate = 30000;%trodes sampling rate
defaultBins = 0.001;% bin size for calculating significant responses
smoothingBins = [0.01 0.001];%bins for smoothing
calcWindow = [0 2]; %window for calculating significant responses
zLimit = 3; %zlimit for calculating significant responses
% firstSpikeWindow = [0 0.5 1 1.5]; %ratios! need to be multiplied by tone duration.
firstSpikeWindow = [0 1];
chosenSpikeBin = 1; %delineates which spike window I will graph.
baselineBin = [-1 0]; %ratio for bin from which baseline firing rate will be calculated
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
%generate placeholder structure
s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(rpvTime,...
    matclustFiles,s,clusterWindow);


%% Extracts Sound Data from soundFile, including freq, db, reps.
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;

%pull important sound and trial information
uniqueFreqs = unique(soundData.Frequencies);
uniqueDBs = unique(soundData.dBs);
numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);
%store these in structured array
s.SoundData.UniqueFrequencies = uniqueFreqs;
s.SoundData.UniqueDBs = uniqueDBs;
s.SoundData.NumFreqs = numFreqs;
s.SoundData.NumDBs = numDBs;

toneDur = soundData.ToneDuration;
toneReps = soundData.ToneRepetitions;
totalTrialNum = length(soundData.Frequencies);

%Recalculate raster window and such
rasterWindow = rasterWindow * toneDur;
firstSpikeWindow = firstSpikeWindow * toneDur;
baselineBin = baselineBin * toneDur;
lfpWindow = rasterWindow;
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes
calcWindow = calcWindow * toneDur;

% Generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for i =1:numFreqs
    freqNameHolder{i} = num2str(uniqueFreqs(i));
end

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

% Does the same for dBs. 
dbSteps = uniqueDBs(2) - uniqueDBs(1);
totalDBs = (uniqueDBs(end) - uniqueDBs(1))/dbSteps;
dbRange = zeros(totalDBs + 1,2);
dbRange(:,1) = uniqueDBs(1):dbSteps:uniqueDBs(end);
for i = 1:size(dbRange,1)
    dbRange(i,2) = find(uniqueDBs == dbRange(i,1));
end

%Sets up master array so that we have a size comparison for DIO
%information. 
master = zeros(size(soundFile.soundData.Frequencies,1),5);

%master(:,1) is reserved for actual input times
%master(:,2) is frequency
master(:,2) = soundData.Frequencies;
%master(:,3) is dB
master(:,3) = soundData.dBs;
%master(:,4) is trial number (chronological)
master(:,4) = 1:1:totalTrialNum;
%master(:,5) is trial num, arranging trials in order from small dB to large
%dB, and low freq to high freq. frequency is larger category.
sortingCounter = 1;
for i = 1:numFreqs
    for j = 1:numDBs
        sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == uniqueDBs(j));
        master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end

%% Extract DIO information. Tuning curve should rely on just one DIO output, DIO1.

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,sampleRate);

%insert to master. check for errors
if length(dioTimes) ~= length(master)
    length(dioTimes)
    length(master)
    error('dioTimes and master mismatched')
elseif length(dioTimes) == length(master)
    master(:,1) = dioTimes;
end

numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db
for i = 1:numUnits
    %allocate a bunch of empty arrays.
    fullHistData = zeros(histBinNum,1);
    
    organizedRasters = cell(numFreqs,numDBs,1); %allocates cell array for rasters
    organizedHist = zeros(numFreqs,numDBs,histBinNum,1);
    
    freqSpecHist = zeros(numFreqs,histBinNum,1);
    
    firstSpikeTimeHolder = cell(numFreqs,numDBs,1);
    firstSpikeStatsHolder = zeros(numFreqs,numDBs,4,size(firstSpikeWindow,2)-1);
    
    binSpikeHolder = cell(numFreqs,numDBs,1);
    binSpikeStatsHolder = zeros(numFreqs,numDBs,2,size(firstSpikeWindow,2)-1);
    
    averageSpikeHolder = zeros(totalTrialNum,1);
    
    sigResp = cell(1,1);
    sigRespGraph = zeros(numFreqs,numDBs,3);
    fullRespGraph = zeros(1,3);
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    alignTimes = master(:,1);
    
    %calculates rasters based on spike information. 
    [rasters] = functionBasicRaster(spikeTimes,alignTimes,rasterWindow);
    rasters(:,3) = master(rasters(:,2),5); %adds information about frequency/amplitude
    fullRasterData = rasters;
    %pulls all spikes from a specific trial in the baseline period, holds
    %the number of these spikes for calculation of average rate and std.
    for k = 1:totalTrialNum
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > baselineBin(1) & rasters(:,1) < baselineBin(2)),1);
    end
    
    averageRate = mean(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
    averageSTD = std(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
    
    [histCounts histCenters] = hist(rasters(:,1),histBinVector);
    fullHistData = histCounts'/totalTrialNum/histBin;
    
    %calculate significant response for combined histogram of all responses
    %to all sounds.
    inputRaster = rasters(rasters(:,1) > calcWindow(1) & rasters(:,1) < calcWindow(2),1);
    [respStore] = ...
    functionBasicResponseSignificance(smoothingBins,defaultBins,calcWindow,...
    averageRate,averageSTD,inputRaster,zLimit,totalTrialNum);
    fullResp = respStore;
    %if there is a significant response, stores important values: start,
    %duration, and peak size.
    if ~isempty(respStore{1})
        fullRespGraph(1) = respStore{1}(1,1);
        fullRespGraph(2) = respStore{1}(1,2)-respStore{1}(1,1);
        fullRespGraph(3) = respStore{1}(1,5);
    end
    
    %allocates empty array.
    sigResp = cell(numFreqs,numDBs);
    
    for k = 1:numFreqs
        for l = 1:numDBs
            targetTrials = master(master(:,2) == uniqueFreqs(k) & master(:,3) == uniqueDBs(l),4); %finds the trial number for all trials of given frequency and amplitude
            findMatches = find(ismember(fullRasterData(:,2),targetTrials)); %uses these trial numbers to pull raster index
            targetRasters = fullRasterData(findMatches,:); %extracts target rasters from all rasters using the above index.
            organizedRasters{k,l} = targetRasters; %saves to organized rasters
            [histCounts,histCenters] = hist(targetRasters(:,1),histBinVector); %calculates histogram with defined bin size
            organizedHist(k,l,:) = histCounts/toneReps/histBin; %saves histogram
            [firstSpikeTimes,firstSpikeStats,binSpikeTimes,binSpikeStats] = ...
                functionBasicFirstSpikeTiming(firstSpikeWindow,targetRasters,toneReps,2,targetTrials); %calculates information about first spike timing
            firstSpikeTimeHolder{k,l} = firstSpikeTimes; %saves first spike times
            firstSpikeStatsHolder(k,l,:,:) = firstSpikeStats; %saves statistics about first spikes
            binSpikeHolder{k,l} = binSpikeTimes; %binned spikes from the defined window.
            binSpikeStatsHolder(k,l,:,:) = binSpikeStats; %stats about those spikes
            inputRaster = targetRasters(targetRasters(:,1) > calcWindow(1) & targetRasters(:,1) < calcWindow(2),1);
            [respStore] = ...
            functionBasicResponseSignificance(smoothingBins,defaultBins,calcWindow,...
            averageRate,averageSTD,inputRaster,zLimit,toneReps); %calculates significance of responses to specific frequencies and dBs
            sigResp{k,l} = respStore;
            if ~isempty(respStore{1}) %if there is significant response, stores this for later display.
                sigRespGraph(k,l,1) = respStore{1}(1,1);
                sigRespGraph(k,l,2) = respStore{1}(1,2)-respStore{1}(1,1);
                sigRespGraph(k,l,3) = respStore{1}(1,5);
            end
        end
        freqSpecHist(k,:) = mean(squeeze(organizedHist(k,:,:)));
    end
        disp(strcat('Raster and Histogram Extraction Complete: Unit ',num2str(i),' Out of ',num2str(numUnits)))
    s.(desigNames{i}).AllRasters = fullRasterData;
    s.(desigNames{i}).AllHistograms = fullHistData;
    s.(desigNames{i}).FreqDBRasters = organizedRasters;
    s.(desigNames{i}).FreqDBHistograms = organizedHist;
    s.(desigNames{i}).FirstSpikeTimes = firstSpikeTimeHolder;
    s.(desigNames{i}).FirstSpikeStats = firstSpikeStatsHolder;
    s.(desigNames{i}).BinSpikes = binSpikeHolder;
    s.(desigNames{i}).BinSpikeStats = binSpikeStatsHolder;
    s.(desigNames{i}).FrequencyHistograms = freqSpecHist;
    s.(desigNames{i}).AverageRate = averageRate;
    s.(desigNames{i}).AverageSTD = averageSTD;
    s.(desigNames{i}).ResponseStats = sigResp;
    s.(desigNames{i}).ResponseStatsGraph = sigRespGraph;
    s.(desigNames{i}).FullResponseStats = fullResp;
    s.(desigNames{i}).FullResponseGraphs = fullRespGraph;
end


% Plotting

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 10 1280 1000])
    %plots average waveform
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms(:,2),'LineWidth',2)
    plot(s.(desigNames{i}).AverageWaveForms(:,1),'r','LineWidth',1)
    plot(s.(desigNames{i}).AverageWaveForms(:,3),'r','LineWidth',1)
    title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    %plots ISI
    subplot(4,6,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(clusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    %plots first spike latency
    subplot(4,3,4)
    imagesc(s.(desigNames{i}).FirstSpikeStats(:,:,1,chosenSpikeBin)')
    colormap hot
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Mean First Spike Latency')
    %plots heatmap of binned spikes to the chosen spike timing window.
    subplot(4,3,7)
    imagesc(squeeze(s.(desigNames{i}).BinSpikeStats(:,:,1,chosenSpikeBin))')
    colormap hot
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Binned Response')
    %plots heatmaps of response reliability in chosen bin 
    subplot(4,3,10)
    imagesc(squeeze(s.(desigNames{i}).FirstSpikeStats(:,:,3,chosenSpikeBin))')
    colormap hot
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Probability of Response')
    %plots rasters (chronological)
    subplot(3,3,2)
    plot(s.(desigNames{i}).AllRasters(:,1),...
        s.(desigNames{i}).AllRasters(:,2),'k.','markersize',4)
    hold on
    ylim([0 totalTrialNum])
    xlim([rasterWindow(1) rasterWindow(2)])
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    title({fileName;desigNames{i}})
    set(0, 'DefaulttextInterpreter', 'none')
    %plots rasters (frequency and amplitude organized)
    subplot(3,3,5)
    plot(s.(desigNames{i}).AllRasters(:,1),...
        s.(desigNames{i}).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    rasterFreqLines = zeros(numFreqs,2);
    rasterFreqLines(:,1) = toneReps*size(uniqueDBs,1)/2:toneReps*size(uniqueDBs,1):totalTrialNum;
    rasterFreqLines(:,2) = uniqueFreqs;
    %this generates green lines separating by Frequency
    for k = 1:size(uniqueFreqs,1)
        plot(rasterWindow,[toneReps*numDBs*k toneReps*numDBs*k],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 totalTrialNum])
    xlim([rasterWindow(1) rasterWindow(2)])
    title('Descending = increase in amplitude and freq')
    %plot heatmap organized by frequency
    subplot(3,3,8)
    imagesc(s.(desigNames{i}).FrequencyHistograms(:,:))
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
    % plot histogram.
    subplot(4,3,3)
    plot(histBinVector,s.(desigNames{i}).AllHistograms,'k','LineWidth',2)
    hold on
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    if s.(desigNames{i}).FullResponseGraphs(1) ~= 0
        plot([s.(desigNames{i}).FullResponseGraphs(1)/1000 ...
            s.(desigNames{i}).FullResponseGraphs(1)/1000],[ylim],'r');
        plot([(s.(desigNames{i}).FullResponseGraphs(1) + ...
            s.(desigNames{i}).FullResponseGraphs(2))/1000 
            (s.(desigNames{i}).FullResponseGraphs(1) + ...
            s.(desigNames{i}).FullResponseGraphs(2))/1000],[ylim],'r');
    end
    xlim([rasterWindow(1) rasterWindow(2)])
    title('Histogram')
    %plot information on onset time (significance calculation)
    subplot(4,3,6)
    imagesc(squeeze(s.(desigNames{i}).ResponseStatsGraph(:,:,1))')
    colormap hot
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Response Latency (Zlimit)')
    %plot information about response duration
    subplot(4,3,9)
    imagesc(squeeze(s.(desigNames{i}).ResponseStatsGraph(:,:,2))')
    colormap hot
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Response Duration (Zlimit)')
    hold off
    %plot information about response peak magnitude
    subplot(4,3,12)
    imagesc(squeeze(s.(desigNames{i}).ResponseStatsGraph(:,:,3))')
    colormap hot
    colorbar
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title('Response Peak (Zlimit)')
    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end

save(fullfile(pname,fname),'s');

end