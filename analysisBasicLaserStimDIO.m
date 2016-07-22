
function [] = analysisBasicLaserStimDIO(fileName);
%this should be basic code to use cluster data to pick out spike time
%points and then align them to auditory stimuli. 

%This needs the following in the same folder: matclust file of picked
%spikes, matlab file with audio order
%and log file from MBED. 

% clear

% dirName = 'C:\TrodesRecordings\160203_ML150108A_R12_2600\160203_ML150108A_R12_2600_toneFinder.matclust';
% fileName = '160511_ML160410D_R17_2402_fullTuningFineGrain';

inputPort = 2;
rasterWindow = [-0.5,0.5];
clusterWindow = [0,0.05];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
tuningWindow = [0,0.1]; %window over which responses are integrated for calculation of tuning!

histBin = 0.025; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes.
%%
%Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';
%%
%finds matclust folder, extracts matclust file names
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');

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
%finds DIO folder, extracts D2 specifically
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};
%extracts DIO stuffs!

[DIOData] = readTrodesExtractedDataFile(D2FileName);
%extracts port states and times of changes
dioState = double(DIOData.fields(2).data);
dioTime = double(DIOData.fields(1).data);

inTimes = dioTime(dioState == 1)/30000;
master(:,1) = inTimes;

histBin = 0.001; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes.
%%


matclustStruct.SoundTimes = master(:,1);

%%
%extracts times of clustered spikes.
for i = 1:numTrodes
    matclustFile = open(matclustFiles{i});
    %this extracts indexes for cluster components
    clusterSizer= size(matclustFile.clustattrib.clustersOn,1);
    matclustStruct.(truncatedNames{i}).ClusterNumber = clusterSizer;
    %this will find which clusters actually have info

    clusterIndex = cell(clusterSizer,1);

    for j = 1:clusterSizer
        clusterIndex{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
    end
    %replaces indices with real times
    for j = 1:clusterSizer
        clusterIndex{j} = matclustFile.clustdata.params(clusterIndex{j},1);
    end
    
    %puts clusterIndex into structured array
    matclustStruct.(truncatedNames{i}).SpikeIndices = clusterIndex;
    
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
        for k = 1:length(inTimes)
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
    
    masterToneRaster = [];
    masterToneHist = [];
    steHolder = [];
    stePlotter = [];
    indivToneHist = [];
end

%%Graphing!

for i = 1:numTrodes
    for j = 1:matclustStruct.(truncatedNames{i}).ClusterNumber
        figure
        %plots average waveform
        subplot(3,2,1)
        hold on
        j
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,2),'LineWidth',2)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,1),'r','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,3),'r','LineWidth',1)
        title(strcat('AverageFiringRate:',num2str(matclustStruct.(truncatedNames{i}).AverageFiringRate(j))))
        %plots ISI
        subplot(3,2,2)
        hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000)
        xlim([0 0.05])
        title('ISI')
        %plots rasters
        subplot(3,1,2)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,1),'k.')
        title(strcat(truncatedNames{i},' Cluster ',num2str(j)))
        %plots histogram
        subplot(3,1,3)
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Histogram{j}(:,1))
        title('Histogram')
    end
end

% clearvars -except matclustStruct fname pname
%saves matclustStruct
pname = pwd;
fname = strcat(fileName,'LaserIDAnalysis');
save(fullfile(pname,fname),'matclustStruct');

end















