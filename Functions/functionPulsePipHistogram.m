
%this code extracts the timing of the first spike (mean time of first
%spike), probability of first spike within a certain period, and timing of
%spikes in a set number of bins. want to include this to existing code.

function [matclustStruct] = functionPulsePipHistogram(rasterData,pipReps,blockNum,histBinVector,histBinNum,histBin,matclustStruct,j,truncatedNames);

pipHistograms = zeros(round(histBinNum),pipReps);

for counter1 = 1:pipReps
    histHolder = rasterData(rasterData(:,5) == (counter1),:);
    [counts,centers] = hist(histHolder(:,2),histBinVector);
    %below code necessary to prevent bugs with row vectors
    countSize = size(counts);
    centerSize = size(centers);
    if countSize(1)>countSize(2)
        counts = counts';
    end
    if centerSize(1)>centerSize(2)
        centers = centers';
    end
    pipHistograms(:,counter1) = counts/histBin/blockNum;
end

matclustStruct.(truncatedNames{i}).PipHistograms{j} = pipHistograms


end