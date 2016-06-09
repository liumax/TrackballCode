%this should be basic code to use cluster data to pick out spike time
%points and then align them to auditory stimuli. 

%This needs the following in the same folder: matclust file of picked
%spikes, matlab file with audio order
%and log file from MBED. 


fileName = '160405_ML160218B_R17_2559_fullTune2';
%sets up file saving stuff
saveName = strcat(fileName,'FullTuningAnalysis','.mat');
[fname pname] = uiputfile(saveName);

%parameters I can play with
rasterWindow = [-0.5,0.5];
clusterWindow = [0,0.05];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];

histBin = 0.010; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes.

%
%Establishes folders and extracts files!
currFolder = pwd;
subFolders = genpath(currFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%find DIO folder and D1 file for analysis
dioFinder = strfind(subFoldersCell,'DIO');%finds DIO folder
dioFinder = find(~cellfun(@isempty,dioFinder)); %determines empty cells, finds index for target cell
dioFolderName = subFoldersCell{dioFinder}; %pulls out the name for the DIO folder
dioFolderSearch = dir(dioFolderName);%pulls dir from DIO folder
dioFileNames = {dioFolderSearch.name}';%pulls names section
D1FileFinder = strfind(dioFileNames,'D1'); %examines names for D1
D1FileFinder = find(~cellfun(@isempty,D1FileFinder));%extracts index of correct file
D1FileName = dioFileNames{D1FileFinder};%pulls out actual file name

%%
%extracts matclust file names
files = dir(fullfile(pwd,'*.mat'));
files = {files.name};
matclustFiles = cell(0);

fileHolder = 1;

for i = 1:length(files)
    if strfind(files{i},'matclust')==1
        matclustFiles{fileHolder} = files{i};
        fileHolder = fileHolder + 1;
    end
end

%removes periods which allow structured array formation.
truncatedNames = matclustFiles;
numTrodes = length(truncatedNames);
for i = 1:length(truncatedNames)
    truncatedNames{i} = truncatedNames{i}(1:find(truncatedNames{i} == '.')-1);
end

%generates structured array for storage of data
matclustStruct = struct;
for i = 1:length(truncatedNames);
    matclustStruct.(truncatedNames{i}) = [];
end

%%
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);
%extracts port states and times of changes
dioState = double(DIOData.fields(2).data);
dioTime = double(DIOData.fields(1).data);

inTimes = dioTime(dioState == 1)/30000;
inTimes = inTimes(1:2:end);
dioState = [];
dioTime = [];
master = zeros(size(inTimes,1),5);
master(:,1) = inTimes;
inTimes = [];
%% pulls out sound data array
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);

%extracts frequency information.
master(:,2) = soundFile.soundData.Frequencies;
uniqueFreqs = unique(master(:,2));
master(:,3) = soundFile.soundData.dBs;
uniqueDBs = unique(master(:,3));

%stores info on total frequencies and magnitudes
matclustStruct.UniqueFreqs = uniqueFreqs;
matclustStruct.UniqueDBs = uniqueDBs;
matclustStruct.RasterAxis = rasterAxis;
matclustStruct.HistogramAxis = histBinVector;
matclustStruct.RasterLimits = rasterWindow;
matclustStruct.SoundTimes = master(:,1);
matclustStruct.Frequencies = master(:,2);
matclustStruct.dBs = master(:,3);
%also stores parameters for rep number and tone duration.
matclustStruct.ToneReps = soundFile.soundData.ToneRepetitions;
matclustStruct.ToneDur = soundFile.soundData.ToneDuration;

%here, need to find new index that goes from low freq to high freq, and
%within each frequency, goes from low amp to high amp

master(:,4) = 1:1:size(master,1);
master(:,5) = zeros;

sortingCounter = 1;

for i = 1:size(uniqueFreqs,1)
    for j = 1:size(uniqueDBs,1)
        sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == uniqueDBs(j));
        master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end
%now, master(:,5) is the index if I want to sort rasters by frequency and
%amplitude
%%

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
    %calculates inter-spike interval, saves this information
    for j = 1:clusterSizer
        %this line pulls the actual indices of spikes for the cluster
        clusterSpikes{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
        %this subtracts all adjacent spikes
        diffSpikes = diff(matclustFile.clustdata.params(clusterSpikes{j},1));
        %this then windows that out, to only plot the small ISIs
        clearedSpikes{j} = diffSpikes(diffSpikes<clusterWindow(2)); %removes long pauses
        diffSpikes = [];
    end
    matclustStruct.(truncatedNames{i}).ISIData = clearedSpikes;
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
    totalTime = matclustFile.clustdata.datarange(2,1);
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
                freqDBHist{j,l} = [counts'*(1/histBin)*(1/matclustStruct.ToneReps),centers'];
                aveFreqHolder(:,1) = aveFreqHolder(:,1) + freqDBHist{j,l}(:,1);
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
%%
%%Graphing!

for i = 1:numTrodes
    for j = 1:matclustStruct.(truncatedNames{i}).ClusterNumber
        hFig = figure
        set(hFig, 'Position', [100 100 1000 800])
        %plots average waveform
        subplot(4,4,1)
        hold on
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,2),'LineWidth',2)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,1),'r','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,3),'r','LineWidth',1)
        title(strcat('AverageFiringRate:',num2str(matclustStruct.(truncatedNames{i}).AverageFiringRate(j))))
        %plots ISI
        subplot(4,4,5)
        hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000)
        xlim([0 0.05])
        title('ISI')
        
        %plots heatmap
        subplot(4,4,9)
        imagesc(matclustStruct.UniqueFreqs,...
            matclustStruct.UniqueDBs,...
            matclustStruct.(truncatedNames{i}).FrequencyResponse{j})
        colormap hot
        title(strcat('Frequency Response.Max',...
            num2str(max(max(matclustStruct.(truncatedNames{i}).FrequencyResponse{j}))),...
            'Min',num2str(min(min(matclustStruct.(truncatedNames{i}).FrequencyResponse{j})))))
        
        %plots reliability of response in heat map
        subplot(4,4,13)
        imagesc(matclustStruct.UniqueFreqs,...
            matclustStruct.UniqueDBs,...
            matclustStruct.(truncatedNames{i}).ResponseReliability{j}')
        colormap hot
        title(strcat('ResponseReliability.Max',...
            num2str(max(max(matclustStruct.(truncatedNames{i}).ResponseReliability{j}))),...
            'Min',num2str(min(min(matclustStruct.(truncatedNames{i}).ResponseReliability{j})))))

        %plots simple rasters
        subplot(2,4,2)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,1),'k.')
        ylim([0 size(matclustStruct.SoundTimes,1)])
        title(strcat(truncatedNames{i},' Cluster ',num2str(j)))
        %plots rasters organized by frequency and amp
        subplot(2,4,6)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,5),'k.')
        hold on
        for k = 1:size(matclustStruct.UniqueFreqs,1)*size(matclustStruct.UniqueDBs,1)
            plot(matclustStruct.RasterLimits,...
                [soundFile.soundData.ToneRepetitions*k soundFile.soundData.ToneRepetitions*k])
        end
        for k = 1:size(uniqueFreqs,1)
            plot(matclustStruct.RasterLimits,...
                [soundFile.soundData.ToneRepetitions*size(uniqueDBs,1)*k soundFile.soundData.ToneRepetitions*size(uniqueDBs,1)*k],...
                'k','LineWidth',2)
        end
        ylim([0 size(matclustStruct.SoundTimes,1)])
        title('Sorted Ascending')
        %plots histogram
        subplot(2,4,3)
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Histogram{j}(:,1),'k','LineWidth',2)
        hold on
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).StandardErrorPlotting(:,j,1),'b')
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).StandardErrorPlotting(:,j,2),'b')
        title('Histogram')
        %plots histograms by frequencies
        subplot(2,4,7)
        x = matclustStruct.(truncatedNames{i}).AverageFrequencyHistogram{j};
        y = x/(max(max(x)));
        for k = 1:size(matclustStruct.UniqueFreqs,1)
            y(:,k) = y(:,k) + (k-1)*0.3;
        end
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            y,'LineWidth',1)
        ylim([0 size(matclustStruct.UniqueFreqs,1)*0.3+0.5])
        title('Histogram By Frequency Ascending')
    end
end
%%
rmpath(subFolders) %removes folders from the path. This is to reduce chances for confusion with multiple runs of code
clearvars -except matclustStruct fname pname
%saves matclustStruct
save(fullfile(pname,fname),'matclustStruct');
% clear






















