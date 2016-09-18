

function [] = analysisBasicTuning(fileName);
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
baselineBin = [-2 0]; %ratio for bin from which baseline firing rate will be calculated
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
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');

%extracts matclust file names and removes periods which allow structured array formation.
truncatedNames = matclustFiles;
numTrodes = length(truncatedNames);
for i = 1:length(truncatedNames)
    truncatedNames{i} = truncatedNames{i}(16:find(truncatedNames{i} == '.')-1);
end

%generates structured array for storage of data
matclustStruct = struct;
for i = 1:length(truncatedNames);
    matclustStruct.(truncatedNames{i}) = [];
end
matclustStruct.NumberTrodes = numTrodes;

%% Extracts Sound Data from soundFile, including freq, db, reps.
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
matclustStruct.SoundData = soundData;

%pull important sound and trial information
uniqueFreqs = unique(soundData.Frequencies);
uniqueDBs = unique(soundData.dBs);
numFreqs = length(uniqueFreqs);
numDBs = length(uniqueDBs);

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

%% Extract data from matclust files
for i = 1:numTrodes
    %this should extract spikes, put them into matcluststruct. Will also
    %calculate RPVs and overall firing rate, as well as average waveform.
    [matclustStruct, clusterSizer] = functionSpikeWaveExtraction(rpvTime,...
        i,matclustFiles,matclustStruct,truncatedNames,clusterWindow);
    disp(strcat('Spike Extraction Complete: nTrode ',num2str(i)))
end

%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db
for i = 1:numTrodes
    clusterSizer = matclustStruct.(truncatedNames{i}).Clusters; %pulls number of clusters
    
    fullRasterData = cell(clusterSizer,1); %allocates matrix space for storage of raster data
    fullHistData = zeros(histBinNum,clusterSizer);
    
    organizedRasters = cell(numFreqs,numDBs,clusterSizer); %allocates cell array for rasters
    organizedHist = zeros(numFreqs,numDBs,histBinNum,clusterSizer);
    
    freqSpecHist = zeros(numFreqs,histBinNum,clusterSizer);
    
    firstSpikeTimeHolder = cell(numFreqs,numDBs,clusterSizer);
    firstSpikeStatsHolder = zeros(numFreqs,numDBs,clusterSizer,4,size(firstSpikeWindow,2)-1);
    
    binSpikeHolder = cell(numFreqs,numDBs,clusterSizer);
    binSpikeStatsHolder = zeros(numFreqs,numDBs,clusterSizer,2,size(firstSpikeWindow,2)-1);
    
    averageRate = zeros(clusterSizer,1);
    averageSTD = zeros(clusterSizer,1);
    averageSpikeHolder = zeros(totalTrialNum,1);
    
    sigResp = cell(clusterSizer,1);
    sigRespGraph = zeros(clusterSizer,numFreqs,numDBs,3);
    fullResp = cell(clusterSizer,1);
    fullRespGraph = zeros(clusterSizer,3);
    
    for j = 1:clusterSizer
        spikeTimes = matclustStruct.(truncatedNames{i}).SpikeTimes{j};
        alignTimes = master(:,1);
        
        [rasters] = functionBasicRaster(spikeTimes,alignTimes,rasterWindow);
        rasters(:,3) = master(rasters(:,2),5); %adds information about frequency/amplitude
        fullRasterData{j} = rasters;
        
        for k = 1:totalTrialNum
            averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > baselineBin(1) & rasters(:,1) < baselineBin(2)),1);
        end
        averageRate(j) = mean(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
        averageSTD(j) = std(averageSpikeHolder/(baselineBin(2)-baselineBin(1)));
        
        [histCounts histCenters] = hist(rasters(:,1),histBinVector);
        fullHistData(:,j) = histCounts'/totalTrialNum/histBin;
        
        %calculate significant response for whole thing.
        inputRaster = rasters(rasters(:,1) > calcWindow(1) & rasters(:,1) < calcWindow(2),1);
        [respStore] = ...
        functionBasicResponseSignificance(smoothingBins,defaultBins,calcWindow,...
        averageRate(j),averageSTD(j),inputRaster,zLimit,totalTrialNum);
        fullResp{j} = respStore;
        if ~isempty(respStore{1})
            fullRespGraph(j,1) = respStore{1}(1,1);
            fullRespGraph(j,2) = respStore{1}(1,2)-respStore{1}(1,1);
            fullRespGraph(j,3) = respStore{1}(1,5);
        end
        
        sigResp{j} = cell(numFreqs,numDBs);
        
        for k = 1:numFreqs
            for l = 1:numDBs
                targetTrials = master(master(:,2) == uniqueFreqs(k) & master(:,3) == uniqueDBs(l),4); %identifies files with the correct trial numbers
                findMatches = find(ismember(fullRasterData{j}(:,2),targetTrials));
                targetRasters = fullRasterData{j}(findMatches,:);
                organizedRasters{k,l,j} = targetRasters; %saves to organized rasters
                [histCounts,histCenters] = hist(targetRasters(:,1),histBinVector);
                organizedHist(k,l,:,j) = histCounts/toneReps/histBin;
                [firstSpikeTimes,firstSpikeStats,binSpikeTimes,binSpikeStats] = ...
                    functionBasicFirstSpikeTiming(firstSpikeWindow,targetRasters,toneReps,2,targetTrials);
                firstSpikeTimeHolder{k,l,j} = firstSpikeTimes;
                firstSpikeStatsHolder(k,l,j,:,:) = firstSpikeStats;
                binSpikeHolder{k,l,j} = binSpikeTimes;
                binSpikeStatsHolder(k,l,j,:,:) = binSpikeStats;
                inputRaster = targetRasters(targetRasters(:,1) > calcWindow(1) & targetRasters(:,1) < calcWindow(2),1);
                [respStore] = ...
                functionBasicResponseSignificance(smoothingBins,defaultBins,calcWindow,...
                averageRate(j),averageSTD(j),inputRaster,zLimit,toneReps);
                sigResp{j}{k,l} = respStore;
                if ~isempty(respStore{1})
                    sigRespGraph(j,k,l,1) = respStore{1}(1,1);
                    sigRespGraph(j,k,l,2) = respStore{1}(1,2)-respStore{1}(1,1);
                    sigRespGraph(j,k,l,3) = respStore{1}(1,5);
                end
            end
            freqSpecHist(k,:,j) = mean(squeeze(organizedHist(k,:,:,j)));
        end
        disp(strcat('Raster and Histogram Extraction Complete: nTrode ',num2str(i),' Cluster ',num2str(j)))
    end
    matclustStruct.(truncatedNames{i}).AllRasters = fullRasterData;
    matclustStruct.(truncatedNames{i}).AllHistograms = fullHistData;
    matclustStruct.(truncatedNames{i}).FreqDBRasters = organizedRasters;
    matclustStruct.(truncatedNames{i}).FreqDBHistograms = organizedHist;
    matclustStruct.(truncatedNames{i}).FirstSpikeTimes = firstSpikeTimeHolder;
    matclustStruct.(truncatedNames{i}).FirstSpikeStats = firstSpikeStatsHolder;
    matclustStruct.(truncatedNames{i}).BinSpikes = binSpikeHolder;
    matclustStruct.(truncatedNames{i}).BinSpikeStats = binSpikeStatsHolder;
    matclustStruct.(truncatedNames{i}).FrequencyHistograms = freqSpecHist;
    matclustStruct.(truncatedNames{i}).AverageRate = averageRate;
    matclustStruct.(truncatedNames{i}).AverageSTD = averageSTD;
    matclustStruct.(truncatedNames{i}).ResponseStats = sigResp;
    matclustStruct.(truncatedNames{i}).ResponseStatsGraph = sigRespGraph;
    matclustStruct.(truncatedNames{i}).FullResponseStats = fullResp;
    matclustStruct.(truncatedNames{i}).FullResponseGraphs = fullRespGraph;
end

%% Plotting

for i = 1:numTrodes
    for j = 1:matclustStruct.(truncatedNames{i}).Clusters
        hFig = figure;
        set(hFig, 'Position', [10 10 1280 1000])
        %plots average waveform
        subplot(4,6,1)
        hold on
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,2),'LineWidth',2)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,1),'r','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,3),'r','LineWidth',1)
        title(strcat('AverageFiringRate:',num2str(matclustStruct.(truncatedNames{i}).AverageRate(j))))
        %plots ISI
        subplot(4,6,2)
        hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000)
        histMax = max(hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000));
        line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
        xlim(clusterWindow)
        title({strcat('ISI RPV %: ',num2str(matclustStruct.(truncatedNames{i}).RPVs(j)));...
            strcat(num2str(matclustStruct.(truncatedNames{i}).RPVNumber(j)),'/',num2str(matclustStruct.(truncatedNames{i}).TotalSpikeNumber(j)))})
        %plots first spike latency
        subplot(4,3,4)
        imagesc(matclustStruct.(truncatedNames{i}).FirstSpikeStats(:,:,j,1,chosenSpikeBin)')
        colormap hot
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Mean First Spike Latency')
        %plots heatmap of binned spikes to the chosen spike timing window.
        subplot(4,3,7)
        imagesc(squeeze(matclustStruct.(truncatedNames{i}).BinSpikeStats(:,:,j,1,chosenSpikeBin))')
        colormap hot
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Frequency Response')
        %plots heatmaps of response reliability in chosen bin 
        subplot(4,3,10)
        imagesc(squeeze(matclustStruct.(truncatedNames{i}).FirstSpikeStats(:,:,j,3,chosenSpikeBin))')
        colormap hot
        colorbar
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title('Probability of Response')
        %plots rasters (chronological)
        subplot(3,3,2)
        plot(matclustStruct.(truncatedNames{i}).AllRasters{j}(:,1),...
            matclustStruct.(truncatedNames{i}).AllRasters{j}(:,2),'k.','markersize',4)
        hold on
        ylim([0 totalTrialNum])
        xlim([rasterWindow(1) rasterWindow(2)])
        plot([0 0],[ylim],'r');
        plot([toneDur toneDur],[ylim],'r');
        title({fileName;strcat(truncatedNames{i},' Cluster ',num2str(j))})
        set(0, 'DefaulttextInterpreter', 'none')
        %plots rasters (frequency and amplitude organized)
        subplot(3,3,5)
        plot(matclustStruct.(truncatedNames{i}).AllRasters{j}(:,1),...
            matclustStruct.(truncatedNames{i}).AllRasters{j}(:,3),'k.','markersize',4)
        hold on
        plot([0 0],[ylim],'r');
        plot([toneDur toneDur],[ylim],'r');
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
        imagesc(matclustStruct.(truncatedNames{i}).FrequencyHistograms(:,:,j))
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
        subplot(2,3,3)
        plot(histBinVector,matclustStruct.(truncatedNames{i}).AllHistograms(:,j),'k','LineWidth',2)
        hold on
        plot([0 0],[ylim],'r');
        plot([toneDur toneDur],[ylim],'r');
        xlim([rasterWindow(1) rasterWindow(2)])
        title('Histogram')
        %plot information about responses?
        hold off
        
        spikeGraphName = strcat(fileName,truncatedNames{i},' Cluster ',num2str(j),'SpikeAnalysis');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end

save(fullfile(pname,fname),'matclustStruct');

end