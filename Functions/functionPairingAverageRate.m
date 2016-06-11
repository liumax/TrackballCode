%This function is meant to calculate average firing rate off of a denoted
%period of time, which will be given by a 2 element array. spikeTimes
%should be the structured array containing spike times from the baseline
%period. 

function [masterStruct] = functionPairingAverageRate(masterStruct,truncatedNames,...
    spikeNames,names);
%pulls average firing rate from baseline period.
for i = 1:size(truncatedNames,2)
    spikeTimes = masterStruct.(truncatedNames{i}).BaselineSpikes;
    clusters = masterStruct.(truncatedNames{i}).Clusters;
    averageFireHolder = zeros(clusters,1);
    for j = 1:clusters
        numSpikes = size(spikeTimes{j},1);
        timePeriod = (masterStruct.TimePeriods.Baseline(2) - masterStruct.TimePeriods.Baseline(1));
        averageFireHolder(j) = numSpikes/timePeriod;
    end
    masterStruct.(truncatedNames{i}).AverageFiringRates = averageFireHolder;
end
%pulls average firing rate from all other periods.
for i = 1:size(truncatedNames,2)
    clusters = masterStruct.(truncatedNames{i}).Clusters;
    averageFireHolder = zeros(clusters,1);
    for j = 1:size(spikeNames,1)
        spikeTimes = masterStruct.(truncatedNames{i}).(spikeNames{j});
        for k = 1:clusters
            numSpikes = size(spikeTimes{k},1);
            timePeriod = (masterStruct.TimePeriods.(names{j})(2) - masterStruct.TimePeriods.(names{j})(1));
            averageFireHolder(k) = numSpikes/timePeriod;
        end
        masterStruct.(truncatedNames{i}).Averages.(names{j}) = averageFireHolder;
    end
end

end