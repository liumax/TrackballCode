%This function is meant to calculate average firing rate off of a denoted
%period of time, which will be given by a 2 element array. spikeTimes
%should be the structured array containing spike times from the baseline
%period. 

function [masterStruct] = functionPairingAverageRate(masterStruct,truncatedNames);

for i = 1:size(truncatedNames,2)
    spikeTimes = masterStruct.(truncatedNames{i}).BaselineSpikes;
    clusters = masterStruct.(truncatedNames{i}).Clusters;
    averageFireHolder = zeros(clusters,1);
    for j = 1:clusters
        numSpikes = size(spikeTimes{j},1);
        initialTime = (masterStruct.TimePeriods.Baseline(2) - masterStruct.TimePeriods.Baseline(1));
        averageFireHolder(j) = numSpikes/initialTime;
    end
    masterStruct.(truncatedNames{i}).AverageFiringRates = averageFireHolder;
end

end