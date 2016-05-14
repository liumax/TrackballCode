%this should be basic code to use cluster data to pick out spike time
%points and then align them to auditory stimuli. 

%This needs the following in the same folder: matclust file of picked
%spikes, matlab file with audio order
%and log file from MBED. 


fileName = '160513_ML160410E_R17_3000_fullTuning';
%sets up file saving stuff
saveName = strcat(fileName,'FullTuningAnalysis','.mat');
[fname pname] = uiputfile(saveName);

%parameters I can play with
rasterWindow = [-0.3,0.4];
clusterWindow = [0,0.05];
lfpWindow = [-0.1,0.2];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];

histBin = 0.005; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes.

%Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%find matclust folder and matclust files for analysis
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');


%%
%extracts matclust file names and removes periods which allow structured array formation.
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

%% pulls out sound data array
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

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);
%extracts port states and times of changes
dioState = double(DIOData.fields(2).data);
dioTime = double(DIOData.fields(1).data);

inTimes = dioTime(dioState == 1)/30000;
dioState = [];
dioTime = [];
if size(inTimes,1) ~= size(master,1)
    if size(inTimes,1) == size(master,1) + 1
        inTimes(end) = [];
    else
        disp('HOLY SHIT YOUR TTL PULSES TO DIO 1 ARE FUCKED')
    end
end
master(:,1) = inTimes;
inTimes = [];


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

for i = 1:numFreqs
    for j = 1:numDBs
        sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == uniqueDBs(j));
        master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end
%now, master(:,5) is the index if I want to sort rasters by frequency and
%amplitude
%%
%now I will extract LFP information!
disp('Starting LFP Analysis')
%this code picks out the LFP folder and moves to that folder
lfpFinder =dir;
lfpFinder = {lfpFinder.name};
lfpIndex = strfind(lfpFinder,'LFP');
lfpIndex = find(not(cellfun('isempty', lfpIndex)));
lfpFolder = lfpFinder{lfpIndex};
newDir = strcat(pwd,'\',lfpFolder);
cd(newDir)
%this code pulls the file names of the LFP files, and reorders them in
%natural order (1,2,10, not 1,10,2)
lfpFinder = dir;
lfpFinder = {lfpFinder.name};
lfpIndex = strfind(lfpFinder,'LFP');
lfpIndex = find(not(cellfun('isempty', lfpIndex)));
lfpFiles = cell(size(lfpIndex,2),1);
lfpFiles= lfpFinder(lfpIndex);
[cs index] = sort_nat(lfpFiles);
lfpFiles = cs;
numLFPs = size(lfpFiles,2);
%this generates the set of colors I want to use
colorArray = zeros(60,3);
colorArray(1,:) = [0,0,1];
for i = 1:29
    colorArray(i+1,:) = [0,(i)/29,(29-i)/29];
end
for i = 1:30
    colorArray(i+30,:) = [(i)/30,(30-i)/30,0];
end 

spacer = size(colorArray,1)/size(uniqueFreqs,1);
spacerArray = 1:1:size(uniqueFreqs,1);
spacerArray = round(spacerArray*spacer);
plotColors = colorArray(spacerArray,:);

%this cell array will hold lfp information organized by nTrode
lfpMaster = cell(numLFPs,1);
%this next code extracts all LFP traces for all trials. 
% These are stored into the cell array
for i = 1:numLFPs
    lfp = readTrodesExtractedDataFile(lfpFiles{i});
    %counts number of LFP samples
    lfpSamples = size(lfp.fields.data,1);
    %makes LFP time points based on # of samples and decimation
    %adjust time to actual time (to match master)
    lfpTimes = (((0:1:lfpSamples-1)*lfp.decimation)+lfp.first_timestamp)'/30000;
    lfpSignals = lfp.fields.data;
    %calculates number of samples in the viewing window
    viewSamples = round((lfpWindow(2)-lfpWindow(1))*lfp.clockrate/lfp.decimation);
    %holds all LFP traces
    lfpHolder = zeros(size(master,1),viewSamples);
    for j = 1:size(master,1)
        finder = find(lfpTimes>master(j,1)+lfpWindow(1),1);
        lfpHolder(j,:) = lfpSignals(finder:finder+viewSamples-1);
    end
    lfpMaster{i} = lfpHolder;
    lfp = [];
    lfpSamples = [];
    lfpTimes = [];
    lfpSignals = [];
    finder = [];
    lfpHolder = [];
    disp(i)
end
disp('Done with Intial LFP Storage')
%This generates a 4D array for storage of plotting data
lfpPlotHolder = zeros(size(lfpFiles,2),size(uniqueDBs,1),size(uniqueFreqs,1),viewSamples);

for i = 1:numLFPs
    tempHolder = lfpMaster{i};
    for j = 1:numDBs
        for k = 1:numFreqs
            traceHolder = tempHolder(intersect(find(master(:,3) == uniqueDBs(j)),find(master(:,2) == uniqueFreqs(k))),:);
            meanHolder = mean(traceHolder);
            lfpPlotHolder(i,j,k,:) = meanHolder;
        end
    end
end
%removes lfpMaster to clear memory
lfpMaster = [];
disp('Done with Storing All Plot Data')
%makes array of LFP means by frequency
lfpMeans = zeros(viewSamples,numLFPs,numFreqs);

for i = 1:numLFPs
    for j = 1:numFreqs
        lfpMeans(:,i,j) = mean(squeeze(lfpPlotHolder(i,:,j,:)));
    end
end
disp('Done Calulating Means')
%calculates absolute min and max for averages, this is to set graph bounds
minLFP = min(min(min(lfpMeans)));
maxLFP = max(max(max(lfpMeans)));

%calculates the zero point based on number of samples
totalLFPWindow = lfpWindow(2)-lfpWindow(1);
lfpZero = abs(lfpWindow(1))*viewSamples/totalLFPWindow;
toneEnd = abs(matclustStruct.ToneDur)*viewSamples/totalLFPWindow;

%saves to structured array!
matclustStruct.LFP.Means = lfpMeans;
matclustStruct.LFP.Window = lfpWindow;
matclustStruct.LFP.ZeroTime = lfpZero;
matclustStruct.LFP.ToneEndTime = toneEnd;

%generates figure with one plot per nTrode, and frequency coded by color. 
h = figure;
set(h, 'Position', [100 100 1000 1000])
%generates text box with information!
mTextBox = uicontrol('style','text');
descr = {'LFP By Frequency';
    'File:';
    fileName;
    'Window (s):';
    lfpWindow;
    'Tone Duration(s):';
    matclustStruct.ToneDur};
set(mTextBox,'String',descr);
set(mTextBox,'Units','normalized');
set(mTextBox,'Position',[0.5,0.74,0.3,0.2])

%generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for i =1:numFreqs
    freqNameHolder{i} = num2str(uniqueFreqs(i));
end

for i = 1:numLFPs
    subplot(numLFPs,2,1+(2*(i-1)))
    plot([lfpZero lfpZero],[minLFP maxLFP],'k');
    hold on
    plot([lfpZero+toneEnd lfpZero+toneEnd],[minLFP maxLFP],'k')
    hold on
    for j = 1:numFreqs
        plot(lfpMeans(:,i,j),'color',plotColors(j,:));
    end
    
    xlim([0 viewSamples])
    ylim([minLFP maxLFP])
    set(gca, 'XTick', [], 'YTick', [])
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1 1.4]) % stretch its width and height
end

%generates legend and places correctly in figure.
hL = legend(freqNameHolder);
set(hL,'Position', [0.5 0.4 0.3 0.2],'Units','normalized');

%returns to original directory
cd(homeFolder)

%save as matlab figure with correct name (fileName+LFP)
lfpName = strcat(fileName,'LFPGraph');
savefig(h,lfpName);

%save as PDF with correct name
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,lfpName,'-dpdf','-r0')

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

%%Graphing!

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
        xlim([0 0.05])
        title('ISI')
        
        %plots heatmap
        subplot(4,3,7)
        imagesc(matclustStruct.UniqueFreqs,...
            matclustStruct.UniqueDBs,...
            matclustStruct.(truncatedNames{i}).FrequencyResponse{j})
        colormap hot
        title(strcat('Frequency Response.Max',...
            num2str(max(max(matclustStruct.(truncatedNames{i}).FrequencyResponse{j}))),...
            'Min',num2str(min(min(matclustStruct.(truncatedNames{i}).FrequencyResponse{j})))))
        
        %plots reliability of response in heat map
        subplot(4,3,10)
        imagesc(matclustStruct.UniqueFreqs,...
            matclustStruct.UniqueDBs,...
            matclustStruct.(truncatedNames{i}).ResponseReliability{j}')
        colormap hot
        title(strcat('ResponseReliability.Max',...
            num2str(max(max(matclustStruct.(truncatedNames{i}).ResponseReliability{j}))),...
            'Min',num2str(min(min(matclustStruct.(truncatedNames{i}).ResponseReliability{j})))))

        %plots simple rasters
        subplot(2,3,2)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,1),'k.')
        hold on
        plot([0 0],[ylim],'r');
        plot([matclustStruct.ToneDur matclustStruct.ToneDur],[ylim],'r');
        ylim([0 size(matclustStruct.SoundTimes,1)])
        title(strcat(truncatedNames{i},' Cluster ',num2str(j)))
        %plots rasters organized by frequency and amp
        subplot(2,3,5)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,5),'k.')
        hold on
        plot([0 0],[ylim],'r');
        plot([matclustStruct.ToneDur matclustStruct.ToneDur],[ylim],'r');
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
        subplot(2,3,3)
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Histogram{j}(:,1),'k','LineWidth',2)
        hold on
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).StandardErrorPlotting(:,j,1),'b')
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).StandardErrorPlotting(:,j,2),'b')
        plot([0 0],[ylim],'r');
        plot([matclustStruct.ToneDur matclustStruct.ToneDur],[ylim],'r');
        title('Histogram')
        %plots histograms by frequencies
        subplot(2,3,6)
        x = matclustStruct.(truncatedNames{i}).AverageFrequencyHistogram{j};
        y = x/(max(max(x)));
        for k = 1:size(matclustStruct.UniqueFreqs,1)
            y(:,k) = y(:,k) + (k-1)*0.3;
        end
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            y,'LineWidth',1)
        hold on
        plot([0 0],[0 k*0.3],'k');
        plot([matclustStruct.ToneDur matclustStruct.ToneDur],[0 k*0.3],'k');
        ylim([0 size(matclustStruct.UniqueFreqs,1)*0.3+0.5])
        title('Histogram By Frequency Ascending')
        hL = legend(freqNameHolder);
        set(hL,'Position', [0.9 0.2 0.1 0.2],'Units','normalized');
        %save as matlab figure with correct name (fileName+LFP)
        spikeGraphName = strcat(truncatedNames{i},' Cluster ',num2str(j),'SpikeAnalysis');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end

rmpath(subFolders) %removes folders from the path. This is to reduce chances for confusion with multiple runs of code
clearvars -except matclustStruct fname pname
%saves matclustStruct
save(fullfile(pname,fname),'matclustStruct');
% clear




















