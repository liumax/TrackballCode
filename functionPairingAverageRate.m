%This function is meant to calculate average firing rate off of a denoted
%period of time, which will be given by a 2 element array. Assumed sampling
%rate is 30,000 Hz for recording data.

function [] = functionPairingAverageRate(masterStruct,truncatedNames,...
    spikeStruct,timesBaseline);
averageFiringRate = cell(size(truncatedNames,2),1);
for i = 1:size(truncatedNames,2)
    spikes = spikeStruct.(truncatedNames{i}).SpikeTimes;
    clusters = spikeStruct.(truncatedNames{i}).Clusters;
    averageFireHolder = zeros(clusters,1);
    for j = 1:clusters
        modSpikes = spikes{j};
        modSpikes = modSpikes(modSpikes<timesBaseline(2));
        initialTime = timesBaseline(2) - timesBaseline(1);
        averageFireHolder(j) = size(modSpikes,1)/(initialTime/30000);
    end
    averageFiringRate{i} = averageFireHolder;
end

masterStruct.AverageFiringRates = averageFiringRate;

end