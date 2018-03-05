load('170810_ML170621B_R10_2250_3mW2sOn4sOffLaserStimAnalysis.mat')

%first find simultaneous spikes
lagWindow = 0.001;
crossWindow = [-lagWindow,lagWindow];
lags = [-lagWindow:0.001:lagWindow];


spikeTimes1 = s.nt5cluster1.SpikeTimes;
spikeTimes2 = s.nt3cluster2.SpikeTimes;
infoStore = zeros(100000,1);
histCounter = 1;
for spikeInd = 1:length(spikeTimes1)
    %subtract spike time from first train out of all of second train
    subSpikes = spikeTimes2 - spikeTimes1(spikeInd);
    %find indices of the targeted shit
    spikeIndices = find(subSpikes <= lagWindow & subSpikes >= -lagWindow);
    infoStore(histCounter:histCounter + length(spikeIndices) - 1) = spikeIndices;
    histCounter = histCounter + length(spikeIndices);
end
infoStore(histCounter:end) = [];

masterSpikes = infoStore;
%now apply these to dataset.


lagWindow = 0.05;
crossWindow = [-lagWindow,lagWindow];
lags = [-lagWindow:0.001:lagWindow];
spikeTimes1 = masterSpikes;
for i = 1:length(s.DesignationName)
    spikeTimes2 = s.(s.DesignationName{i}).SpikeTimes;
    histStore = ones(100000,1);
    histCounter = 1;
    for spikeInd = 1:length(spikeTimes1)
        %subtract spike time from first train out of all of second train
        subSpikes = spikeTimes2 - spikeTimes1(spikeInd);
        %remove things outside the window of interest
        subSpikes(subSpikes>crossWindow(2) | subSpikes<crossWindow(1)) = [];
        histStore(histCounter:(histCounter + length(subSpikes)-1)) = subSpikes;
        histCounter = histCounter + length(subSpikes);
    end
    histStore(histCounter:end) = [];
    figure
    hist(histStore,lags)
    xlim(crossWindow)
    title((s.DesignationName{i}))
end

