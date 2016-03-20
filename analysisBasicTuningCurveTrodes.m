%this should be basic code to use cluster data to pick out spike time
%points and then align them to auditory stimuli. 

%This needs the following in the same folder: matclust file of picked
%spikes, matlab file with audio order
%and log file from MBED. 

% dirName = 'C:\TrodesRecordings\160203_ML150108A_R12_2600\160203_ML150108A_R12_2600_toneFinder.matclust';
soundName = '160225_ML160218A_L12_2500_toneFinder';
mbedName = '160225_ML160218A_L12_2500_toneFinder';

inputPort = 2;
rasterWindow = [-0.5,0.5];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
tuningWindow = [0,0.1]; %window over which responses are integrated for calculation of tuning!

histBin = 0.025; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes.
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
% matclustName = 'matclust_param_nt1';


% matclustName = strcat(matclustName,'.mat');
soundName = strcat(soundName,'.mat');
mbedName = strcat(mbedName,'.txt');

% matclustFile = open(matclustName);
soundFile = open(soundName);
%%
%extracts port states!
[portStates] = maxTrialVariableNoTask(mbedName);

%Extracts times of audio inputs.
inTimes = find(diff([0;portStates.inStates(:,2)])==1);
inTimes = portStates.tStamps(inTimes)';
master(:,1) = inTimes/1000;

%extracts frequency information.
master(:,2) = soundFile.soundData.Frequencies;
uniqueFreqs = unique(master(:,2));

matclustStruct.SoundTimes = master(:,1);
matclustStruct.Frequencies = master(:,2);

%%
%extracts times of clustered spikes.
for i = 1:numTrodes
    matclustFile = open(matclustFiles{i});
    %this extracts indexes for cluster components
    clusterSizer= length(matclustFile.clustattrib.clusters);
    matclustStruct.(truncatedNames{i}).ClusterNumber = length(matclustFile.clustattrib.clusters);

    clusterIndex = cell(clusterSizer,1);

    for j = 1:clusterSizer
        clusterIndex{j} = matclustFile.clustattrib.clusters{1,j}.index;
    end
    %replaces indices with real times
    for j = 1:clusterSizer
        clusterIndex{j} = matclustFile.clustdata.params(clusterIndex{j},1);
    end
    
    %puts clusterIndex into structured array
    matclustStruct.(truncatedNames{i}).SpikeIndices = clusterIndex;
    
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
            toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,3) = master(k,2);
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
    
    %calculates responses binned per tone frequency
    processRaster = cell(clusterSizer,1);
    
    %this pulls out correct frequency and correct timing, then finds which
    %indices are correct for both.
    for j = 1:clusterSizer
        processRaster{j} = zeros(length(uniqueFreqs),1);
        for k = 1:length(uniqueFreqs)
            corrFreq = find(masterToneRaster{j}(:,3) == uniqueFreqs(k));
            corrTime = find(masterToneRaster{j}(:,2) > tuningWindow(1) & masterToneRaster{j}(:,2) < tuningWindow(2));
            processRaster{j}(k) = length(intersect(corrFreq,corrTime));
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
        figure
        
        %plots rasters
        subplot(3,1,1)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,1),'k.')
        title(strcat(truncatedNames{i},' Cluster ',num2str(j)))
        %plots histogram
        subplot(3,1,2)
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Histogram{j}(:,1))
        title('Histogram')
        subplot(3,1,3)
        plot(matclustStruct.(truncatedNames{i}).FrequencyResponse{j})
        title('Frequency Response')
        set(gca,'XTickLabel',uniqueFreqs)
    end
end





















