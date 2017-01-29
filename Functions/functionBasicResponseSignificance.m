
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


function [sigResp] = functionBasicResponseSignificance(s,calcWindow,spikeTimes,alignTimes,trialNum);
thresholdHz = s.Parameters.ThresholdHz;
baselineSpikes = functionBasicRaster(spikeTimes,alignTimes,s.Parameters.BaselineWindow);
inputRaster = functionBasicRaster(spikeTimes,alignTimes,calcWindow);
inputRaster = inputRaster(:,1);
baselineVector = [s.Parameters.BaselineWindow(1) + s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.BaselineWindow(2)];

histBin = s.Parameters.histBin;

if length(baselineSpikes) > s.Parameters.minSpikes
    %This is the case where there are sufficient baseline spikes to go
    %ahead directly
    %% First step is to compute an overall histogram of baseline bins.
    
    baselineHist = hist(baselineSpikes(:,1),baselineVector)/histBin/trialNum;
    %These bins will serve as the distribution against which we will
    %compare. Currently, with 0.4 seconds/0.005 sec bins, this should
    %produce 80 bins.
    
    %% Based on  distribution, generate positive percentile values
    %generate a percentile range from the data. 
    percentileRange = 100.-(s.Parameters.zLimit*100);

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
    targetHistVector = [calcWindow(1) + histBin/2: histBin:calcWindow(2)];
    %calculate histogram
    targetHist = hist(inputRaster,targetHistVector);
    targetHist = reshape(targetHist,[],1)/histBin/trialNum;
    %allocate space for information about significance
    targetHist(:,2:length(percentileRange)+1) = zeros;
    %find all values above specific thresholds
    for respInd = 1:length(percentileRange)
        sigFinder = find(targetHist(:,1)> valueRange(respInd));
        targetHist(sigFinder,respInd+1) = 1;
    end
    
    %% Now look for negative responses
    negRange = s.Parameters.zLimit*100;
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
    
    if length(find(targetHist(:,2:end)>0))>s.Parameters.minSigSpikes;
        sigSpike = 1;
    else
        sigSpike = 0;
    end
    
    %% save data!
    sigResp.Histogram = targetHist;
    sigResp.BaselineHist = baselineHist;
    sigResp.MeanBaseline = mean(baselineHist);
    sigResp.Centers = targetHistVector;
    sigResp.Warning = 0; %warning indicates baseline was too low for anything.
    sigResp.SpikeNumber = length(baselineSpikes);
    sigResp.SigSpike = sigSpike;
    sigResp.MaxResp = max(targetHist(:,1));
elseif length(baselineSpikes) <= s.Parameters.minSpikes
    %This is the case where there are insufficient baseline spikes. I think
    %my approach here will be to still use the baseline histogram for
    %percentiles, but then also use raster data to get a firm idea of total
    %number of spikes. 

    baselineHist = hist(baselineSpikes(:,1),baselineVector)/histBin/trialNum;
    %These bins will serve as the distribution against which we will
    %compare. Currently, with 0.4 seconds/0.005 sec bins, this should
    %produce 80 bins.
    
    %% Based on  distribution, generate positive percentile values
    %generate a percentile range from the data. 
    percentileRange = 100.-(s.Parameters.zLimit*100);

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
    targetHistVector = [calcWindow(1) + histBin/2: histBin:calcWindow(2)];
    %calculate histogram
    targetHist = hist(inputRaster,targetHistVector);
    targetHist = reshape(targetHist,[],1)/histBin/trialNum;
    %allocate space for information about significance
    targetHist(:,2:length(percentileRange)+1) = zeros;
    %find all values above specific thresholds
    for respInd = 1:length(percentileRange)
        sigFinder = find(targetHist(:,1)> valueRange(respInd));
        numBinSpikes = targetHist(sigFinder,1);
        sigFinder(numBinSpikes<thresholdHz) = [];
        targetHist(sigFinder,respInd+1) = 1;
    end
    
    %% Now look for negative responses
    negRange = s.Parameters.zLimit*100;
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
    
    if length(find(targetHist(:,2:end)>0))>s.Parameters.minSigSpikes;
        sigSpike = 1;
    else
        sigSpike = 0;
    end
    
    %% save data!
    sigResp.Histogram = targetHist;
    sigResp.BaselineHist = baselineHist;
    sigResp.MeanBaseline = mean(baselineHist);
    sigResp.Centers = targetHistVector;
    sigResp.Warning = 1; %warning indicates baseline was too low for anything.
    sigResp.SpikeNumber = length(baselineSpikes);
    sigResp.SigSpike = sigSpike;
    sigResp.MaxResp = max(targetHist(:,1));

end
end