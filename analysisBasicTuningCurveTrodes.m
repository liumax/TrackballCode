%this should be basic code to use cluster data to pick out spike time
%points and then align them to auditory stimuli. 

%This needs the following in the same folder: matclust file of picked
%spikes, matlab file with audio order
%and log file from MBED. 

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
matclustName = 'matclust_param_nt16';
soundName = '160129_ML150108A_L17_2300_toneFinder';
mbedName = '160129_ML150108A_L17_2300_toneFinder';

matclustName = strcat(matclustName,'.mat');
soundName = strcat(soundName,'.mat');
mbedName = strcat(mbedName,'.txt');

matclustFile = open(matclustName);
soundFile = open(soundName);
%%
%extracts port states!
[portStates] = maxTrialVariableNoTask(mbedName);

%Extracts times of audio inputs.
inTimes = find(portStates.inStates(:,inputPort)==1);
inTimes = portStates.tStamps(inTimes)';
master(:,1) = inTimes/1000;

%extracts frequency information.
master(:,2) = soundFile.soundData.Frequencies;

%%
%extracts times of clustered spikes.

%this extracts indexes for cluster components
clusterSizer = length(matclustFile.clustattrib.clusters);

clusterIndex = cell(clusterSizer,1);

for i = 1:clusterSizer
    clusterIndex{i} = matclustFile.clustattrib.clusters{1,i}.index;
end
%replaces indices with real times
for i = 1:clusterSizer
    clusterIndex{i} = matclustFile.clustdata.params(clusterIndex{i},1);
end

%%
rasterHolder = 1;
%generate rasters and histograms, ignores frequency information
for j = 1:clusterSizer
    rasterHolder = 1;
    for i = 1:length(inTimes)
        toneHolder = clusterIndex{j}(clusterIndex{j}>master(i,1)+rasterWindow(1) & clusterIndex{j}<master(i,1)+rasterWindow(2));
        toneHolder = toneHolder - master(i,1);
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
        indivToneHist{j,i} = [counts'*(1/histBin),centers'];
        counts = [];
        centers = [];
        %fills in large raster plot. holder updated position.
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = i;
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

stdHolder = zeros(length(histBinVector),length(master(:,1)));
steHolder = zeros(length(histBinVector),clusterSizer);

for i = 1:clusterSizer
    for j = 1:length(master(:,1))
        stdHolder(:,j) = indivToneHist{i,j}(:,1);
    end
    steHolder(:,i) = std(stdHolder,0,2)/sqrt(length(master(:,1)));
end

%generates lines representing standard error.
stePlotter = zeros(length(histBinVector),clusterSizer,2);
for i = 1:clusterSizer
    stePlotter(:,i,1) = masterToneHist{i}(:,1)-steHolder(:,i);
    stePlotter(:,i,2) = masterToneHist{i}(:,1)+steHolder(:,i);
end
%%
%This section is for generation of a tuning curve plot. Sadly, no amplitude
%data appears to be available for this :D

rasterHolder = 1;
toneRaster = zeros(100000,4);
freqToneRaster = cell(clusterSizer,1);

for j = 1:clusterSizer
    for i = 1:length(inTimes)
        toneHolder = clusterIndex{j}(clusterIndex{j}>master(i,1)+tuningWindow(1) & clusterIndex{j}<master(i,1)+tuningWindow(2));
        toneHolder = toneHolder - master(i,1);
        %fills in large raster plot. holder updated position.
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = i;
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
        toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,3) = master(i,2);
%         toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,4) = dBs(i,1);
        rasterHolder = rasterHolder + length(toneHolder);
        toneRaster(toneRaster(:,1) == 0,:) = [];
        toneHolder = [];
    end
    freqToneRaster{j} = toneRaster;
    toneRaster = zeros(100000,2);
end

processRaster = cell(clusterSizer,1);
uniqueFreqs = unique(master(:,2));

for i = 1:clusterSizer
    processRaster{i} = zeros(length(uniqueFreqs),1);
    for j = 1:length(uniqueFreqs)
            processRaster{i}(j) = length(find(freqToneRaster{i}(:,3) == uniqueFreqs(j)));
    end
end






















