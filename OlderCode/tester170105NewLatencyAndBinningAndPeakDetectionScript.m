
%This function is meant to calculate the latency of a significant response,
%the binned number of spikes in a given window, and the peak response
%following a tone. 

function [latPeakBinOut] = functionLatPeakBinCalculation(toneWindow,generalWindow,fullWindow,alignedSpikes,numTrials,trialColumn,trialNumbers,latBins,peakBins,percentCutoff,baselineCutoff);

%% Calculate binned spikes occurring during tone or general period
binSpikeTone = zeros(length(numTrials),1);
binSpikeGen = zeros(length(numTrials),1);

for latCount = 1:numTrials
    targetSpikes = alignedSpikes(alignedSpikes(:,trialColumn) == trialNumbers(latCount),1);
    if ~isempty(targetSpikes)
        genFinder = targetSpikes(targetSpikes>generalWindow(1) & targetSpikes<generalWindow(2));
        if ~isempty(genFinder)
            binSpikeGen(latCount) = length(genFinder);
            toneFinder = find(genFinder > toneWindow(1) & genFinder < toneWindow(2));
            if ~isempty(toneFinder)
                binSpikeTone(latCount) = length(toneFinder);
            end
        end
    end
end

%% Calculate probability of a spike falling in the tone or general period
zeroFindTone = length(find(binSpikeTone == 0));
zeroFindGen = length(find(binSpikeGen == 0));
probSpikeTone = (numTrials - zeroFindTone)/numTrials;
probSpikeGen = (numTrials - zeroFindGen)/numTrials;

%% Calculate peak response, using peakBins, which will be larger than latBins
%First, calculate histogram
peakBinVector = [fullWindow(1)+ peakBins/2:peakBins:fullWindow(2) - peakBins/2];
peakHist = hist(alignedSpikes(:,1),peakBinVector);
%then, find indices for windows I care about
[~,genStart] = min(abs((peakBinVector - generalWindow(1))));
[~,genEnd] = min(abs((peakBinVector - generalWindow(2))));
[~,toneStart] = min(abs((peakBinVector - toneWindow(1))));
[~,toneEnd] = min(abs((peakBinVector - toneWindow(2))));
%finally, find max values during this window
[peakGenVal,peakGenTime]= max(peakHist(genStart:genEnd));
[peakToneVal,peakToneTime]= max(peakHist(toneStart:toneEnd));

%supplement indices for time with the correct offsets
peakGenTime = peakBinVector(peakGenTime + genStart - 1);
peakToneTime = peakBinVector(peakToneTime + toneStart - 1);

%% Now calculate latency of response! XD
%first, calculate histogram
latBinVector = [fullWindow(1)+ latBins/2:latBins:fullWindow(2) - latBins/2];
latHist = smooth(hist(alignedSpikes(:,1),latBinVector),3); %smooth over three bins to reduce craziness

%recalculate indices for windows I care about
[~,genStart] = min(abs((latBinVector - generalWindow(1))));
[~,genEnd] = min(abs((latBinVector - generalWindow(2))));


%calculate percentile cutoffs
percentLimit = prctile(latHist(1:genStart),percentCutoff);
baseLimit = prctile(latHist(1:genStart),baselineCutoff);

%calculate bins above cutoff. 
sigBins = find(latHist(genStart:genEnd)>percentLimit);

%find number of bins that are significant
respToneLength = length(sigBins);
%find if these bins are consecutive
respToneDiff = diff(sigBins);
%find first of consecutive string
respToneFirst = find(respToneDiff == 1,1,'first');

%now we check if the response satisfies two criteria: more than one
%significant bin, and consecutive responses
if respToneLength > 1 && ~isempty(respToneFirst)
    latSwitch = 0; %switch to get out of while loop
    latInd = 1; %index that increments to move backwards in time
    while latSwitch == 0
        respVal = latHist(genStart + sigBins(respToneFirst) - latInd); %calculates histogram value 
        if respVal > baseLimit
            latInd = latInd + 1;
        elseif respVal <= baseLimit
            latInd = latInd - 1;
            latSwitch = 1;
        end
        if latInd >= sigBins(respToneFirst)
            latSwitch = 1;
        end
    end
    if latInd <= sigBins(respToneFirst)
        toneRespLat = latBinVector(genStart + sigBins(respToneFirst)-latInd);
    else
        toneRespLat = 0;
    end
else
    toneRespLat = 0;
end

latPeakBinOut = struct;
latPeakBinOut.ResponseLatency = toneRespLat;
latPeakBinOut.ResponseFirstSig = latBinVector(genStart + sigBins(respToneFirst)-1);

latPeakBinOut.BinnedSpikesTone = binSpikeTone;
latPeakBinOut.BinnedSpikesGen = binSpikeGen;

latPeakBinOut.ProbSpikeTone = probSpikeTone;
latPeakBinOut.ProbSpikeGen = probSpikeGen;

latPeakBinOut.PeakRespTone = peakToneVal;
latPeakBinOut.PeakRespToneTime = peakToneTime;
latPeakBinOut.PeakRespGen = peakGenVal;
latPeakBinOut.PeakRespGenTime = peakGenTime;

end