function [] = functionTuningAnalysisMod(fileName);

%this should be basic code to use cluster data to pick out spike time
%points and then align them to auditory stimuli. 

%This needs the following in the same folder: matclust file of picked
%spikes, matlab file with audio order
%and log file from MBED. 

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
    truncatedNames{i} = truncatedNames{i}(1:find(truncatedNames{i} == '.')-1);
end

trodesDesignation = cell(size(truncatedNames));
for i = 1:length(truncatedNames)
    trodesDesignation{i} = truncatedNames{i}(17:end);
end

%generates structured array for storage of data
matclustStruct = struct;
for i = 1:length(truncatedNames);
    matclustStruct.(truncatedNames{i}) = [];
end

%% Extracts Sound Data from soundFile, including freq, db, reps.
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);

%extracts frequency information.
master = zeros(size(soundFile.soundData.Frequencies,1),5);
master(:,2) = soundFile.soundData.Frequencies;
uniqueFreqs = unique(master(:,2));
master(:,3) = soundFile.soundData.dBs;
uniqueDBs = unique(master(:,3));

numFreqs = size(uniqueFreqs,1);
numDBs = size(uniqueDBs,1);

%% Store Info into structured array.
matclustStruct.UniqueFreqs = uniqueFreqs;
matclustStruct.UniqueDBs = uniqueDBs;
matclustStruct.SoundTimes = master(:,1);
matclustStruct.Frequencies = master(:,2);
matclustStruct.dBs = master(:,3);
%also stores parameters for rep number and tone duration.
matclustStruct.ToneReps = soundFile.soundData.ToneRepetitions;
matclustStruct.ToneDur = soundFile.soundData.ToneDuration;

%% Set and/or Generate Raster and Histogram Parameters, store in Structured Array
rasterWindow = [-matclustStruct.ToneDur,matclustStruct.ToneDur*3];
lfpWindow = [-matclustStruct.ToneDur,matclustStruct.ToneDur*3];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
%These are for refractory period violations and looking at spike ITIs
rpvTime = 0.0013; %limit to be considered an RPV.
clusterWindow = [-0.01,0.03]; %this is hardcoded for consistency
%clims1 is a value that sets the limits for heatmaps for displaying firing.
clims1 = [-1 1]; %limits for the range of heatmaps for firing. Adjust if reach saturation. Currently based on log10
%Use raster parameters to set histogram settings.
histBin = 0.005; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes

%store values in structured array.
matclustStruct.RasterAxis = rasterAxis;
matclustStruct.HistogramAxis = histBinVector;
matclustStruct.RasterLimits = rasterWindow;
%%

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%extracts port states and times of changes
dioState = double(DIOData.fields(2).data);
dioTime = double(DIOData.fields(1).data);

inTimes = dioTime(dioState == 1)/30000;
inTimesSpacing = diff(inTimes);
dioState = [];
dioTime = [];
if size(inTimes,1) ~= size(master,1)
    if size(inTimes,1) == size(master,1) + 1
        inTimes(end) = [];
        master(:,1) = inTimes;
    elseif inTimesSpacing(size(master,1)) > 3*mean(inTimesSpacing)
        master(:,1) = inTimes(1:size(master,1));
        disp('HOLY SHIT YOUR TTL PULSES TO DIO 1 ARE FUCKED, BUT ADJUSTED')
    else
        disp('I DONT KNOW WHATS GOING ON')
        pause
    end
elseif size(inTimes,1) == size(master,1)
    master(:,1) = inTimes;
end

inTimes = [];

%% Finds the number of octaves, makes array of octave steps. This will be used for imagesc graphing applications
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

%% Does the same for dBs. 
dbSteps = uniqueDBs(2) - uniqueDBs(1);
totalDBs = (uniqueDBs(end) - uniqueDBs(1))/dbSteps;
dbRange = zeros(totalDBs + 1,2);
dbRange(:,1) = uniqueDBs(1):dbSteps:uniqueDBs(end);
for i = 1:size(dbRange,1)
    dbRange(i,2) = find(uniqueDBs == dbRange(i,1));
end


%% Generates index that goes from low freq to high freq, and within each freq, goes from low to high amplitude
master(:,4) = 1:1:size(master,1);
master(:,5) = zeros;

sortingCounter = 1;

for i = 1:numFreqs
    for j = 1:numDBs
        sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == uniqueDBs(j));
        master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end
%now, master(:,5) is the index if I want to sort rasters by frequency and
%amplitude

%% Generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for i =1:numFreqs
    freqNameHolder{i} = num2str(uniqueFreqs(i));
end


%%
%here I need to calculate LFPs!
% [s] = functionLFPaverage(master, lfpWindow, matclustStruct,homeFolder,fileName, uniqueFreqs, uniqueDBs, numFreqs, numDBs);
% matclustStruct.LFPData = s;
%%
%extracts times of clustered spikes.
for i = 1:numTrodes
    %extract all spike times total. this is to calculate ISIs for individual
%clusters
    matclustFile = open(matclustFiles{i});
    %this extracts indexes for cluster components
    clusterSizer= size(matclustFile.clustattrib.clustersOn,1); %pulls number of clusters
    matclustStruct.(truncatedNames{i}).ClusterNumber =  clusterSizer; %stores this number into structured array
    
    %preps series of arrays
    clusterSpikes = cell(clusterSizer,1);
    clusterWaves = cell(clusterSizer,1);
    clearedSpikes = cell(clusterSizer,1);
    rpvViolationPercent = zeros(clusterSizer,1);
    %calculates inter-spike interval, saves this information
    for j = 1:clusterSizer
        %this line pulls the actual indices of spikes for the cluster
        clusterSpikes{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
        %this subtracts all adjacent spikes
        diffSpikes = diff(matclustFile.clustdata.params(clusterSpikes{j},1));
        %this calculates all spikes within the refractory period violation
        %period
        rpvViolationPercent(j) = size(diffSpikes(diffSpikes<rpvTime),1)/size(clusterSpikes{j},1)*100;
        %this then windows out spikes, to only plot the small ISIs
        clearedSpikes{j} = diffSpikes(diffSpikes<clusterWindow(2)); %removes long pauses
        diffSpikes = [];
    end
    matclustStruct.(truncatedNames{i}).ISIData = clearedSpikes;
    matclustStruct.(truncatedNames{i}).RPVs = rpvViolationPercent;
    clearedSpikes = [];
    
    %prepares cluster indices. Finds actual points of index per cluster.
    %Then replaces with real times. 
    clusterIndex = cell(clusterSizer,1);
    for j = 1:clusterSizer
        clusterIndex{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
        clusterIndex{j} = matclustFile.clustdata.params(clusterIndex{j},1);
    end
    %puts clusterIndex into structured array
    matclustStruct.(truncatedNames{i}).SpikeTimes = clusterIndex;
    
    %calculates average firing rate
    totalTime = matclustFile.clustdata.datarange(2,1)-matclustFile.clustdata.datarange(1,1);
    averageFiring = zeros(clusterSizer,1);
    for j = 1:clusterSizer
        averageFiring(j) = size(clusterIndex{j},1)/totalTime;
    end
    matclustStruct.(truncatedNames{i}).AverageFiringRate = averageFiring;
    
    %pull out average waveform and standard error
    targetWaveName = strcat('waves',truncatedNames{i}(15:end),'.mat');
    waveLoader = open(targetWaveName);
    waveLoader =squeeze(waveLoader.waves);
    waveHolder = cell(clusterSizer,1);
    averageWaveHolder = zeros(size(waveLoader,1),clusterSizer,3);
    for j = 1:clusterSizer
        waveHolder{j} = waveLoader(:,clusterSpikes{j});
        averageWaveHolder(:,j,2) = mean(waveHolder{j},2);
        averageWaveHolder(:,j,1) = mean(waveHolder{j},2)-std(waveHolder{j},0,2)/sqrt(size(waveHolder{j},2));
        averageWaveHolder(:,j,3) = mean(waveHolder{j},2)+std(waveHolder{j},0,2)/sqrt(size(waveHolder{j},2));
    end
    for j = 1:clusterSizer
        waveHolder{j} = waveLoader(:,clusterSpikes{j});
        averageWaveHolder(:,j,2) = mean(waveHolder{j},2);
        averageWaveHolder(:,j,1) = mean(waveHolder{j},2)-std(waveHolder{j},0,2);
        averageWaveHolder(:,j,3) = mean(waveHolder{j},2)+std(waveHolder{j},0,2);
    end
    matclustStruct.(truncatedNames{i}).AverageWaveForms = averageWaveHolder;
    matclustStruct.(truncatedNames{i}).AllWaveForms = waveHolder;
    
    masterToneRaster = [];
    masterToneHist = [];
    %generate rasters and histograms, ignores frequency information
    for j = 1:clusterSizer
        rasterHolder = 1;
        for k = 1:size(master,1)
            toneHolder = clusterIndex{j}(clusterIndex{j}>master(k,1)+rasterWindow(1) & clusterIndex{j}<master(k,1)+rasterWindow(2));
            toneHolder = toneHolder - master(k,1);
            %This will generate a histogram based on each individual trace,
            %which can be used to generate standard deviations.
            [counts,centers] = hist(toneHolder,histBinVector);
            %below code necessary to prevent bugs with row vectors
            countSize = size(counts);
            centerSize = size(centers);
            if countSize(1)>countSize(2)
                counts = counts';
            end
            if centerSize(1)>centerSize(2)
                centers = centers';
            end
            indivToneHist{j,k} = [counts'*(1/histBin),centers'];
            counts = [];
            centers = [];
            %fills in large raster plot. holder updated position.
            toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = k;
            toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
            toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,3) = master(k,2); %stores frequency
            toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,4) = master(k,3); %stores amplitude
            toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,5) = master(k,5); %stores freq/amp index
            rasterHolder = rasterHolder + length(toneHolder);
            toneRaster(toneRaster(:,1) == 0,:) = [];
            toneHolder = [];
        end
        masterToneRaster{j} = toneRaster;
        toneRaster = zeros(100000,2);
        [counts,centers] = hist(masterToneRaster{j}(:,2),histBinVector);
        countSize = size(counts);
        centerSize = size(centers);
        if countSize(1)>countSize(2)
            counts = counts';
        end
        if centerSize(1)>centerSize(2)
            centers = centers';
        end
        masterToneHist{j} = [counts'*(1/histBin)/length(master(:,1)),centers'];
        counts = [];
        centers = [];
    end
    matclustStruct.(truncatedNames{i}).Rasters = masterToneRaster;
    matclustStruct.(truncatedNames{i}).Histogram = masterToneHist;
    
    %code to make histograms for each set of frequencies and amplitudes.
    %Extracts data from rasters, generates new histograms for each
    %frequency/amplitude pair.
    respReliab = cell(clusterSizer,1);
    
    for k = 1:clusterSizer
        percentResponse = zeros(size(matclustStruct.UniqueFreqs,1),size(matclustStruct.UniqueDBs,1));
        freqHistHolder = matclustStruct.(truncatedNames{i}).Rasters{k}(:,2:4);
        for j=1:size(matclustStruct.UniqueFreqs,1);
            for l = 1:size(matclustStruct.UniqueDBs,1);
                aveFreqHolder = zeros(histBinNum,1);
                %finds all trials of specific frequency and amplitude
                trialNumHolder = intersect(find(matclustStruct.Frequencies == matclustStruct.UniqueFreqs(j)),...
            find(matclustStruct.dBs == matclustStruct.UniqueDBs(l)));
        %generates placeholder for responses per trial
                responseHolder = zeros(size(trialNumHolder,1),2);
                responseHolder(:,1) = trialNumHolder;
                %this goes through the individual trials and finds if there
                %are any responses during the tone period
                for m = 1:matclustStruct.ToneReps;
                    indivResp = size(find(matclustStruct.(truncatedNames{i}).Rasters{k}(:,2)>0 & ...
            matclustStruct.(truncatedNames{i}).Rasters{k}(:,2)<matclustStruct.ToneDur &...
            matclustStruct.(truncatedNames{i}).Rasters{k}(:,1)==responseHolder(m)),1);
                    responseHolder(m,2) = indivResp;
                end
                %uses these values to calculate response reliability,
                %stores this as a percentage
                percentResponse(j,l) = (matclustStruct.ToneReps - size(find(responseHolder(:,2) == 0),1))...
                    /matclustStruct.ToneReps;
                %this pulls all times from all trials for the given dB and
                %frequency, and puts them all together 
                spikeTimeHolder = matclustStruct.(truncatedNames{i}).Rasters{k}...
                    (matclustStruct.(truncatedNames{i}).Rasters{k}(:,3) == ...
                matclustStruct.UniqueFreqs(j) & ...
                matclustStruct.(truncatedNames{i}).Rasters{k}(:,4) ...
                == matclustStruct.UniqueDBs(l),2);
            %converts this data into the form of a histogram for storage
                [counts centers] = hist(spikeTimeHolder,histBinVector);
                countSize = size(counts);
                centerSize = size(centers);
                if countSize(1)>countSize(2)
                    counts = counts';
                end
                if centerSize(1)>centerSize(2)
                    centers = centers';
                end
                freqDBHist{j,l} = [counts'*(1/histBin)*(1/matclustStruct.ToneReps),centers'];%This should make it such taht units are in Hz/trial overall. That is, the average response per trial.
                aveFreqHolder(:,1) = aveFreqHolder(:,1) + freqDBHist{j,l}(:,1)/size(matclustStruct.UniqueDBs,1); %160529 adjusted so that this doesnt produce overinflated values due to repetitions of multiple dbs.
                counts = [];
                centers = [];
            end
            averageFreqResp(:,j) = aveFreqHolder;
            aveFreqHolder = [];
        end
        masterFreqDBHist{k} = freqDBHist;
        respReliab{k} = percentResponse;
        aveFreqRespMaster{k} = averageFreqResp;
    end
    
    matclustStruct.(truncatedNames{i}).FreqDBSpecificHist = masterFreqDBHist;
    matclustStruct.(truncatedNames{i}).ResponseReliability = respReliab;
    matclustStruct.(truncatedNames{i}).AverageFrequencyHistogram = aveFreqRespMaster;
    
    stdHolder = zeros(length(histBinVector),length(master(:,1)));
    steHolder = zeros(length(histBinVector),clusterSizer);

    for j = 1:clusterSizer
        for k = 1:length(master(:,1))
            stdHolder(:,k) = indivToneHist{j,k}(:,1);
        end
        steHolder(:,j) = std(stdHolder,0,2)/sqrt(length(master(:,1)));
    end

    %generates lines representing standard error.
    stePlotter = zeros(length(histBinVector),clusterSizer,2);
    for j = 1:clusterSizer
        stePlotter(:,j,1) = masterToneHist{j}(:,1)-steHolder(:,j);
        stePlotter(:,j,2) = masterToneHist{j}(:,1)+steHolder(:,j);
    end
    matclustStruct.(truncatedNames{i}).StandardError = steHolder;
    matclustStruct.(truncatedNames{i}).StandardErrorPlotting = stePlotter;
    
    %calculates responses binned per tone frequency
    processRaster = cell(clusterSizer,1);
    
    %this pulls out correct frequency and correct timing, then finds which
    %indices are correct for both.
    for j = 1:clusterSizer
        processRaster{j} = zeros(length(uniqueDBs),length(uniqueFreqs));
        for k = 1:length(uniqueFreqs)
            for m = 1:length(uniqueDBs)
            corrFreq = find(masterToneRaster{j}(:,3) == uniqueFreqs(k));
            corrDBs = find(masterToneRaster{j}(:,4) == uniqueDBs(m));
            corrTime = find(masterToneRaster{j}(:,2) > 0 & masterToneRaster{j}(:,2) < matclustStruct.ToneDur);
            processRaster{j}(m,k) = length(intersect(intersect(corrFreq,corrTime),corrDBs));  
            end
        end
    end
    
    matclustStruct.(truncatedNames{i}).FrequencyResponse = processRaster;
    matclustStruct.(truncatedNames{i}).Frequencies = uniqueFreqs;
    
    masterToneRaster = [];
    masterToneHist = [];
    steHolder = [];
    stePlotter = [];
    indivToneHist = [];
    processRaster = [];
end


%% Graphing!
for i = 1:numTrodes
    for j = 1:matclustStruct.(truncatedNames{i}).ClusterNumber
        hFig = figure;
        set(hFig, 'Position', [10 10 1280 1000])
        %plots average waveform
        subplot(4,3,1)
        hold on
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,2),'LineWidth',2)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,1),'r','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,3),'r','LineWidth',1)
        title(strcat('AverageFiringRate:',num2str(matclustStruct.(truncatedNames{i}).AverageFiringRate(j))))
        %plots ISI
        subplot(4,3,4)
        hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000)
        histMax = max(hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000));
        line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
        xlim(clusterWindow)
        title(strcat('ISI RPV %: ',num2str(matclustStruct.(truncatedNames{i}).RPVs(j))))
        
         %plots histogram
        subplot(2,3,4)
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Histogram{j}(:,1),'k','LineWidth',2)
        hold on
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).StandardErrorPlotting(:,j,1),'b')
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).StandardErrorPlotting(:,j,2),'b')
        plot([0 0],[ylim],'r');
        plot([matclustStruct.ToneDur matclustStruct.ToneDur],[ylim],'r');
        xlim([matclustStruct.RasterLimits(1) matclustStruct.RasterLimits(2)])
        title('Histogram')
        
        %plots simple rasters
        subplot(2,3,2)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,1),'k.','markersize',4)
        hold on
        ylim([0 size(matclustStruct.SoundTimes,1)])
        xlim([matclustStruct.RasterLimits(1) matclustStruct.RasterLimits(2)])
        plot([0 0],[ylim],'r');
        plot([matclustStruct.ToneDur matclustStruct.ToneDur],[ylim],'r');
        
        title(strcat(truncatedNames{i},' Cluster ',num2str(j)))
        
        %plots rasters organized by frequency and amp
        subplot(2,3,5)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,5),'k.','markersize',4)
        hold on
        plot([0 0],[ylim],'r');
        plot([matclustStruct.ToneDur matclustStruct.ToneDur],[ylim],'r');
        %removed blue lines to be able to see rasters better
%         for k = 1:size(matclustStruct.UniqueFreqs,1)*size(matclustStruct.UniqueDBs,1)
%             plot(matclustStruct.RasterLimits,...
%                 [soundFile.soundData.ToneRepetitions*k soundFile.soundData.ToneRepetitions*k])
%         end
        rasterFreqLines = zeros(numFreqs,2);
        rasterFreqLines(:,1) = soundFile.soundData.ToneRepetitions*size(uniqueDBs,1)/2:soundFile.soundData.ToneRepetitions*size(uniqueDBs,1):size(matclustStruct.SoundTimes,1);
        rasterFreqLines(:,2) = uniqueFreqs;
        %this generates green lines separating by Frequency
        for k = 1:size(uniqueFreqs,1)
            plot(matclustStruct.RasterLimits,...
                [soundFile.soundData.ToneRepetitions*size(uniqueDBs,1)*k soundFile.soundData.ToneRepetitions*size(uniqueDBs,1)*k],...
                'g','LineWidth',2)
        end
        set(gca,'YTick',rasterFreqLines(:,1));
        set(gca,'YTickLabel',rasterFreqLines(:,2));
        ylim([0 size(matclustStruct.SoundTimes,1)])
        xlim([matclustStruct.RasterLimits(1) matclustStruct.RasterLimits(2)])
        title('Sorted Ascending')
        
        %plots heatmap. This uses a log10 scaling for change in firing
        %rates. This way, no change is essentially zero on the imagesc
        %scale, rather than having a linear scale where green actually
        %represents a fairly large increase in response, and there is
        %little room for inhibition
        subplot(4,3,3)
        imagesc(log10(matclustStruct.(truncatedNames{i}).FrequencyResponse{j}/(matclustStruct.ToneReps*matclustStruct.ToneDur*matclustStruct.(truncatedNames{i}).AverageFiringRate(j))),clims1)
        colormap hot
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title(strcat('Normalized Frequency Response.Max',...
            num2str(max(max(matclustStruct.(truncatedNames{i}).FrequencyResponse{j}))/matclustStruct.(truncatedNames{i}).AverageFiringRate(j)),...
            'Min',num2str(min(min(matclustStruct.(truncatedNames{i}).FrequencyResponse{j}))/matclustStruct.(truncatedNames{i}).AverageFiringRate(j))))
        
        %plots reliability of response in heat map
        subplot(4,3,6)
        clims = [0,1];
        imagesc(matclustStruct.UniqueFreqs,...
            matclustStruct.UniqueDBs,...
            matclustStruct.(truncatedNames{i}).ResponseReliability{j}',clims)
        colormap hot
%         colormap hot
        set(gca,'XTick',octaveRange(:,2));
        set(gca,'XTickLabel',octaveRange(:,1));
        set(gca,'YTick',dbRange(:,2));
        set(gca,'YTickLabel',dbRange(:,1));
        title(strcat('ResponseReliability.Max',...
            num2str(max(max(matclustStruct.(truncatedNames{i}).ResponseReliability{j}))),...
            'Min',num2str(min(min(matclustStruct.(truncatedNames{i}).ResponseReliability{j})))))
        
        %plots heatmap by frequencies
        subplot(2,3,6)
        x = matclustStruct.(truncatedNames{i}).AverageFrequencyHistogram{j}/matclustStruct.(truncatedNames{i}).AverageFiringRate(j);
        imagesc(log10(x'),clims1)
        set(gca,'YTick',octaveRange(:,2));
        set(gca,'YTickLabel',octaveRange(:,1));
        set(gca,'XTick',[1:10:size(histBinVector,2)]);
        set(gca,'XTickLabel',histBinVector(1:10:end));
        histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
        histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),matclustStruct.ToneDur);
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
%         title('Heatmap by Frequency and Time Max')
        title(strcat('Normalized Heatmap by F and T.Max',...
            num2str(max(max(matclustStruct.(truncatedNames{i}).AverageFrequencyHistogram{j}))/matclustStruct.(truncatedNames{i}).AverageFiringRate(j)),...
            'Min',num2str(min(min(matclustStruct.(truncatedNames{i}).AverageFrequencyHistogram{j}))/matclustStruct.(truncatedNames{i}).AverageFiringRate(j))))
        hold off
        %save as matlab figure with correct name (fileName+LFP)
        spikeGraphName = strcat(trodesDesignation{i},' Cluster ',num2str(j),'SpikeAnalysis');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end
%%
rmpath(subFolders) %removes folders from the path. This is to reduce chances for confusion with multiple runs of code

%saves matclustStruct
save(fullfile(pname,fname),'matclustStruct');
