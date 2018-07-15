
%This is a function to calculate statistically significant responses
%following a tone. This utilizes spike shuffling of the baseline period,
%which generates a series of histogram values. These histogram values
%provide a distribution against which any future bins can be compared. 

%Inputs:
%s: main structured array %baselineWindow: window that is being used, in time (seconds).
%calcWindow: window in real time (seconds), that you want to read out
%inputRaster: raster input, should be a vector of spike times
%trialNum: number of trials being processed. This helps to normalize things
%properly.
 

%Outputs
%sigResp: output with response information. 

%sigResp.Histogram: n x m, where n is number of bins, m is 1 +
%length(zLimit). First column is actual histogram data. Columns 2 to m are
%zeros for non significant values, 1s for significant values. Different
%columns represent different levels of significance

%sigResp.Centers: n x 1 vector, should indicate the time periods for the
%bins above

%sigResp.Warning: if there are insufficient spikes in the baseline period,
%this warning will be triggered (producing a 1). 

% sigResp.SpikeNumber: this will be the number of spikes detected in the
% baseline period. 


function [sigResp] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,trialNum,minSpikes,sigCalcBin,baselineWindow,zLimit,minSigSpikes,smoothWindow);

baselineSpikes = functionBasicRaster(spikeTimes,alignTimes,baselineWindow);
inputRaster = functionBasicRaster(spikeTimes,alignTimes,[calcWindow(1)-(round(smoothWindow/2)*sigCalcBin) calcWindow(2)+(round(smoothWindow/2)*sigCalcBin)]);
inputRaster = inputRaster(:,1);
baselineVector = [baselineWindow(1) + sigCalcBin/2:sigCalcBin:baselineWindow(2)];

%% First step is to compute an overall histogram of baseline bins.
baselineHist = smooth(hist(baselineSpikes(:,1),baselineVector)/sigCalcBin/trialNum,smoothWindow);
%These bins will serve as the distribution against which we will
%compare. Currently, with 0.4 seconds/0.005 sec bins, this should
%produce 80 bins.
%170425 adjustment: instead, moving to 1ms bins with 10ms smoothing. 

%% Based on  distribution, generate positive percentile values
%generate a percentile range from the data. 
percentileRange = 100.-(zLimit*100);

%find the values from the distribution that are associated with the given
%percentiles
valueRange = zeros(length(percentileRange),1);
%     simValueRange = zeros(length(percentileRange),1);

for respInd = 1:length(percentileRange)
    valueRange(respInd) = prctile(baselineHist,percentileRange(respInd));
%         simValueRange(respInd) = prctile(baselineSim,percentileRange(respInd));
end

%% Now we use the  distribution and percentiles to calculate positive significant responses
%generate vector for bin centers for histogram
targetHistVector = [calcWindow(1)-(round(smoothWindow/2)*sigCalcBin) + sigCalcBin/2: sigCalcBin:calcWindow(2)+(round(smoothWindow/2)*sigCalcBin)]; %180302 added in code to have histogram overlap over both edges so that smoothing is more accurate. 
%calculate histogram
targetHist = hist(inputRaster,targetHistVector);
targetHist = smooth(reshape(targetHist,[],1)/sigCalcBin/trialNum,smoothWindow);
% length(targetHist)
%remove edges after smoothing
targetHist(1:(round(smoothWindow/2))) = [];
targetHist(end-round(smoothWindow/2)+1:end) = [];
targetHistVector(1:(round(smoothWindow/2))) = [];
targetHistVector(end-round(smoothWindow/2)+1:end) = [];
% length(targetHist)

%allocate space for information about significance
targetHist(:,2:length(percentileRange)+1) = zeros;
%find all values above specific thresholds
for respInd = 1:length(percentileRange)
    sigFinder = find(targetHist(:,1)> valueRange(respInd));
    targetHist(sigFinder,respInd+1) = 1;
end

%% Now look for negative responses
negRange = zLimit*100;
negValues = zeros(length(negRange),1);
%     simNegValues = zeros(length(negRange),1);
for respInd = 1:length(percentileRange)
    negValues(respInd) = prctile(baselineHist,negRange(respInd));
%         simNegValues(respInd) = prctile(baselineSim,negRange(respInd));
end
%allocate space
targetHist(:,end+1:end+length(percentileRange)) = zeros;
%find all values below specific thresholds
for respInd = 1:length(percentileRange)
    sigFinder = find(targetHist(:,1)< negValues(respInd));
    targetHist(sigFinder,length(percentileRange)+respInd+1) = 1;
end

if length(find(targetHist(:,4)>0))>minSigSpikes; %180714 adjusted so that this can work on just the 0.001 percentile values. 
    sigSpikePos = 1;
end
if length(find(targetHist(:,6)>0))>minSigSpikes; %180714 negatives still will be at 0.01 level. 
    sigSpikeNeg = 1;
end
if  length(find(targetHist(:,4)>0))<=minSigSpikes;
    sigSpikePos = 0;
end
if  length(find(targetHist(:,6)>0))<=minSigSpikes;
    sigSpikeNeg = 0;
end


%% save data!
sigResp.Histogram = targetHist;
sigResp.BaselineHist = baselineHist;
sigResp.MeanBaseline = mean(baselineHist);
sigResp.STDBaseline = std(mean(baselineHist));
sigResp.Centers = targetHistVector;
if length(baselineSpikes) > minSpikes
    sigResp.Warning = 0; %warning indicates baseline was too low for anything.
else
    sigResp.Warning = 1; %warning indicates baseline was too low for anything.
end
sigResp.SpikeNumber = length(baselineSpikes);
sigResp.SigSpikePos = sigSpikePos;
sigResp.SigSpikeNeg = sigSpikeNeg;
sigResp.MaxResp = max(targetHist(:,1));

end