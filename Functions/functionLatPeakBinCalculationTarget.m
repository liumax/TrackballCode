
%181216 This function is meant to calculate the latency of a significant response,
%the binned number of spikes in a given window, and the peak response
%following a tone. This is only to look at a targeted window that was calculated beforehand.  

%Inputs:
%targetWindow: two element vector delineating beginning and end of tone response in seconds

%fullWindow: two element vector representing the bounds of the raster, in
%seconds. This should correspond to raster window from main analysis
%scripts.

%alignedSpikes: n x m matrix, with n spikes, and spikes in column 1. Column
%2 should contain information about trial number, to allow for separation
%by trial number for calculation of response binning

%numTrials: scalar input for the number of trials to examine

%trialColumn: column in alignedSpikes with the information about which
%trial number

%trialNumbers: numbers of the trials of interest

%latBins: bin size for calculation of latency. Default should be 0.001

%peakBins: bin size for calculation of peak response. Default should be
%0.005

%percentCutoff: percentile cutoff for calculation as a significant
%response. This is important for calculation of latency. Default is 99.9

%baselineCutoff: percentile cutoff for tracing back response initiation.
%This is lower than percentCutoff, and is used as an arbitrary bound to
%determine when responses are no longer significant, thereby deriving the
%first response bin. Default is 95.

%Outputs:

%latPeakBinOut: structured array with multiple components. 

% Concerning latency, we have: 
% ResponseLatency (latency in ms), and 
% ResponseFirstSig, which is the first bin reaching significance (in ms)

% Concerning binned responses, we have:
% BinnedSpikesTone: which bins all spikes in the tone period, and
% BinnedSpikesGen: which bins all spikes in the general period.

% Concerning Probability of response, we have:
% ProbSpikeTone: overall probability of having any spikes in the tone period
% ProbSpikeGen: overall probability of having any spikes in the general period

% Concerning Peak Responses, we have:
% PeakRespTone: peak response during tone period
% PeakRespToneTime: time of the peak response, in ms
% PeakRespGen: peak response during general period
% PeakRespGenTime: time of general peak response, in ms.


function [latPeakBinOut] = functionLatPeakBinCalculationTarget(targetWindow,fullWindow,alignedSpikes,numTrials,trialColumn,trialNumbers,latBins,percentCutoff,baselineCutoff);
smoothWindow = 9;
%% Calculate binned spikes occurring during target period
binSpikeTar = zeros(numTrials,1);

for latCount = 1:numTrials
    targetSpikes = alignedSpikes(alignedSpikes(:,trialColumn) == trialNumbers(latCount),1);
    if ~isempty(targetSpikes)
        tarFinder = targetSpikes(targetSpikes>targetWindow(1) & targetSpikes<targetWindow(2));
        if ~isempty(tarFinder)
            binSpikeTar(latCount) = length(tarFinder);
        end
    end
end

binSpikeTarBase = zeros(numTrials,1);

for latCount = 1:numTrials
    targetSpikes = alignedSpikes(alignedSpikes(:,trialColumn) == trialNumbers(latCount),1);
    if ~isempty(targetSpikes)
        tarFinder = targetSpikes(targetSpikes<-1*targetWindow(1) & targetSpikes>-1*targetWindow(2));
        if ~isempty(tarFinder)
            binSpikeTarBase(latCount) = length(tarFinder);
        end
    end
end

%% Calculate probability of a spike falling in the tone or general period

zeroFindTar = length(find(binSpikeTar == 0));

probSpikeTar = (numTrials - zeroFindTar)/numTrials;

%% Calculate peak response, using latBins
%First, calculate histogram
peakBinVector = [fullWindow(1)+ latBins/2:latBins:fullWindow(2) - latBins/2];
peakHist = smooth(hist(alignedSpikes(:,1),peakBinVector)/latBins/numTrials,smoothWindow);
% peakHist = smooth(hist(alignedSpikes(:,1),peakBinVector)/latBins/numTrials,smoothWindow);
%then, find indices for target window
[~,tarStart] = min(abs((peakBinVector - generalWindow(1))));
[~,tarEnd] = min(abs((peakBinVector - generalWindow(2))));
%finally, find max values during this window
[peakTarVal,peakTarTime]= max(peakHist(tarStart:tarEnd));
[peakOverallVal,peakOverallTime]= max(peakHist);
[peakFastVal,peakFastTime]= max(peakHist(fastStart:fastEnd));
%compensate for bin size and number of trials
peakTarVal = peakTarVal/latBins/numTrials;
peakOverallVal = peakOverallVal/latBins/numTrials;
%supplement indices for time with the correct offsets
peakTarTime = peakBinVector(peakTarTime + tarStart - 1);
peakOverallTime = peakBinVector(peakOverallTime);


%% Now calculate latency of response! XD
%first, calculate histogram
latBinVector = [fullWindow(1)+ latBins/2:latBins:fullWindow(2) - latBins/2];
latHist = smooth(hist(alignedSpikes(:,1),latBinVector),3); %smooth over three bins to reduce craziness

%recalculate indices for windows I care about
tarStart = find(latBinVector > generalWindow(1),1,'first');
tarEnd = find(latBinVector < generalWindow(2),1,'last');

%calculate percentile cutoffs
percentLimit = prctile(latHist(1:tarStart),percentCutoff);
baseLimit = prctile(latHist(1:tarStart),baselineCutoff);

%calculate bins above cutoff. 
sigBins = find(latHist(tarStart:tarEnd)>percentLimit);

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
        respVal = latHist(tarStart + sigBins(respToneFirst) - latInd); %calculates histogram value 
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
        toneRespLat = latBinVector(tarStart + sigBins(respToneFirst)-latInd);
    else
        toneRespLat = 0;
    end
else
    toneRespLat = 0;
end

latPeakBinOut = struct;

latPeakBinOut.LatHist = latHist;
latPeakBinOut.LatBinVector = latBinVector;

latPeakBinOut.ResponseLatency = toneRespLat;
latPeakBinOut.ResponseFirstSig = latBinVector(tarStart + sigBins(respToneFirst)-1);

latPeakBinOut.BinnedSpikesTar = binSpikeTar;
latPeakBinOut.BinnedSpikesTarBase = binSpikeTarBase;

latPeakBinOut.BinDiff = mean(binSpikeTar) - mean(binSpikeTarBase);

latPeakBinOut.BinSigVals = signrank(binSpikeTarBase,binSpikeTar);

latPeakBinOut.ProbSpikeTar = probSpikeTar;

latPeakBinOut.PeakTarVal = peakTarVal;
latPeakBinOut.PeakTarTime = peakTarTime;
latPeakBinOut.PeakOverallVal = peakOverallVal;
latPeakBinOut.PeakOverallTime = peakOverallTime;

end