
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


function [sigResp] = functionBasicResponseSignificance(s,calcWindow,baselineSpikes,inputRaster,trialNum,rasterWindow);

if length(baselineSpikes) > s.Parameters.minSpikes
    %% First step is to calculate a bootstrapped baseline distribution
    %Calculate spike ISIs
    negDiff = diff(baselineSpikes);
    %pull number of times to shuffle spikes
    numShuffle = s.Parameters.numShuffle;
    %pull number of spikes
    numSpikes = length(baselineSpikes);
%     baselineWindow = s.Parameters.BaselineWindow;

    %allocate array for storage of data
    simDataHolder = zeros(numShuffle,numSpikes);
    simDataHolder(:,1) = baselineSpikes(1); %seed with first spike time

    for respInd = 1:numShuffle
        %generate random permutation of ISIs
        randInd = randperm(numSpikes-1);
        %add these to the data!
        simDataHolder(respInd,2:end) = negDiff(randInd);
    end
    %sum to generate appropriate values (adding ISIs to times)
    simDataHolder = cumsum(simDataHolder,2);

    %process so that produce histogram with appropriate bin size
    histBin = s.Parameters.histBin;
    baseHistVector = [rasterWindow(1)+histBin/2:histBin:0];
    simHist = hist(simDataHolder',baseHistVector);

    %reshape and normalize data so that reporting firing rate in Hz. 
    baselineSim = reshape(simHist,[],1)/histBin/trialNum;


    %% Based on bootstrapped distribution, generate positive percentile values
    %generate a percentile range from the bootstrapped data. 
    percentileRange = 100.-(s.Parameters.zLimit*100);

    %find the values from the distribution that are associated with the given
    %percentiles
    valueRange = zeros(length(percentileRange),1);

    for respInd = 1:length(percentileRange)
        valueRange(respInd) = prctile(baselineSim,percentileRange(respInd));
    end

    %% Now we use the bootstrapped distribution and percentiles to calculate positive significant responses
    %eliminate data points outside of desired calculation window
    targetData = inputRaster(inputRaster>calcWindow(1) & inputRaster<calcWindow(2));
    %generate vector for bin centers for histogram
    targetHistVector = [calcWindow(1) + histBin/2: histBin:calcWindow(2)];
    %calculate histogram
    targetHist = hist(targetData,targetHistVector);
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
    for respInd = 1:length(percentileRange)
        negValues(respInd) = prctile(baselineSim,negRange(respInd));
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
    sigResp.Centers = targetHistVector;
    sigResp.Warning = 0; %warning indicates baseline was too low for shuffling.
    sigResp.SpikeNumber = length(baselineSpikes);
    sigResp.SigSpike = sigSpike;
else
    %% If too few spikes, generate histogram with no significant values.
    histBin = s.Parameters.histBin;
    %eliminate data points outside of desired calculation window
    targetData = inputRaster(inputRaster>calcWindow(1) & inputRaster<calcWindow(2));
    %generate vector for bin centers for histogram
    targetHistVector = [calcWindow(1) + histBin/2:histBin:calcWindow(2)];
    %calculate histogram
    targetHist = hist(targetData,targetHistVector);
    %allocate space for information about significance
    targetHist(:,2:length(s.Parameters.zLimit)+1) = zeros;
    
    %save data!
    sigResp.Histogram = targetHist;
    sigResp.Centers = targetHistVector;
    sigResp.Warning = 1; %warning indicates baseline was too low for shuffling.
    sigResp.SpikeNumber = length(baselineSpikes);
%     disp('Insufficient Spikes for Shuffling. Not performing analysis')
end
end