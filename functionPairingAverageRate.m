%This function is meant to calculate average firing rate off of a denoted
%period of time, which will be given by a 2 element array. spikeTimes
%should be the structured array containing spike times from the baseline
%period. Since time is in ms, there is a 1000 fudge factor. 

function [masterStruct] = functionPairingAverageRate(masterStruct,truncatedNames,...
    spikeStruct,spikeTimes);
averageFiringRate = cell(size(truncatedNames,2),1);
for i = 1:size(truncatedNames,2)
    clusters = spikeStruct.(truncatedNames{i}).Clusters;
    averageFireHolder = zeros(clusters,1);
    for j = 1:clusters
        numSpikes = size(spikeTimes{j},1);
        initialTime = (masterStruct.TimePeriods.Baseline(2) - masterStruct.TimePeriods.Baseline(1))*1000;
        averageFireHolder(j) = numSpikes/initialTime;
    end
    averageFiringRate{i} = averageFireHolder;
end

masterStruct.AverageFiringRates = averageFiringRate;

end