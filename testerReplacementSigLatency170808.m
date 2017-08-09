
%Things I need:

%variables: rasterWindow, toneWindow, generalWindow, 
%numTrials,trialColumn(for rasters), timeColumn(for rasters),
%baselineCutoff, percentCutoff (for latency, should be 99 and 95%
%respectively), percentileRange for significance. Probably can use baseline
%and percentCutoff interchangeably with this, since I'm going to be using
%99 and 95%. 


%inputs: baseline histogram values across all trials. want 1 ms bins
%smoothed over 10 ms. 

%rasters for specific target stimulus






%% Process basic binning of data

%generate binned values for tone, general, and pre windows. This goes
%through each trial to pull individual trial responses, will allow
%calculation of stats at later point if necessary. Uses equivalent windows
%for pre tone and pre gen so that stats will match up more properly.
binSpikeTone = zeros(numTrials,1);
binSpikeGen = zeros(numTrials,1);
binSpikePreTone = zeros(numTrials,1);
binSpikePreGen = zeros(numTrials,1);

for latCount = 1:numTrials
    targetSpikes = targetRasters(targetRasters(:,trialColumn) == trialNumbers(latCount),1);
    if ~isempty(targetSpikes)
        genFinder = targetSpikes(targetSpikes>generalWindow(1) & targetSpikes<generalWindow(2));
        if ~isempty(genFinder)
            binSpikeGen(latCount) = length(genFinder);
            toneFinder = find(genFinder > toneWindow(1) & genFinder < toneWindow(2));
            if ~isempty(toneFinder)
                binSpikeTone(latCount) = length(toneFinder);
            end
        end
        preFinder = targetSpikes(targetSpikes>-1*toneWindow(2) & targetSpikes<toneWindow(1));
        if ~isempty(preFinder)
            binSpikePreTone(latCount) = length(preFinder);
        end
        preFinder = targetSpikes(targetSpikes>-1*generalWindow(2) & targetSpikes<generalWindow(1));
        if ~isempty(preFinder)
            binSpikePreGen(latCount) = length(preFinder);
        end
    end
end

%save information.
output.Binned.ToneIndiv = binSpikeTone;
output.Binned.GenIndiv = binSpikeGen;
output.Binned.PreToneIndiv = binSpikePreTone;
output.Binned.PreGenIndiv = binSpikePreGen;
output.Binned.ToneAv = mean(binSpikeTone);
output.Binned.GenAv = mean(binSpikeGen);
output.Binned.PreToneAv = mean(binSpikePreTone);
output.Binned.PreGenAv = mean(binSpikePreGen);


%% Process spiking data from rasters to histograms
%make a histogram with 1ms bins, use this for everything to follow
hist1msVector = [rasterWindow(1) + 0.0005:0.001:rasterWindow(2)-0.0005];
zeroPoint = find(hist1msVector>0,1,'first');
toneEnd = find(hist1msVector>toneWindow(2),1,'first');
genEnd = find(hist1msVector>generalWindow(2),1,'first');

hist1ms = hist(targetRasters(:,timeColumn),hist1msVector)/numTrials/0.001;

%for latency calculations, I'm going to want 3 ms bins.

hist3ms = smooth(hist1ms,3);

%for calculating significance, I think that 5-10 ms bins should be good. I
%believe that I settled on 5 ms bins for good reason, so we'll stick to
%those
hist5ms = smooth(hist1ms,5);

output.Histogram1ms = hist1ms;
output.Histogram3ms = hist3ms;
output.Histogram5ms = hist5ms;
output.HistVector = hist1msVector;

%% Detect Peak Values
%Here we will use 5ms smoothing for detection of peaks.
%Do this for the tone period, post tone period, and the general period. I
%do the post tone period here because I wish to determine whether there are
%large peaks in the post, which could be masked by larger tones in the
%general period. 

[peakToneVal,peakToneTime] = max(hist5ms(zeroPoint:toneEnd));
[peakGenVal,peakGenTime] = max(hist5ms(zeroPoint:genEnd));
[peakPostVal,peakPostTime] = max(hist5ms(toneEnd:genEnd));

output.Peaks.ToneVal = peakToneVal;
output.Peaks.ToneTime = peakToneTime; %since things are stored at 1ms, this should produce the timing in ms. 
output.Peaks.GenVal = peakGenVal;
output.Peaks.GenTime = peakGenTime;
output.Peaks.PostVal = peakPostVal;
output.Peaks.PostTime = peakPostTime;

%% Look at Excitatory Responses

%going to use 5ms histogram for tracking excitatory responses

%first, set benchmarks for significance
percentileRange = 100.-(zLimit*100);
for respInd = 1:length(percentileRange)
    valueRange(respInd) = prctile(hist5ms(1:zeroPoint),percentileRange(respInd));
end

%find values above these specific cutoffs
posSigHold = zeros(length(hist5ms),length(percentileRange));
for respInd = 1:length(percentileRange)
    sigFinder = find(hist5ms(zeroPoint:genEnd)> valueRange(respInd));
    posSigHold(sigFinder,respInd) = 1;
end





















