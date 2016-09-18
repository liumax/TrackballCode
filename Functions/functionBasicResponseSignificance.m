
%This is a function to calculate statistically significant responses
%following a tone. 
%Inputs:
%smoothingBins: an nx1 vector, with the size of intended bins for analysis
%defaultBins: single scalar, should be 0.001 seconds
%calcWindow: window over which you look for responses to the tone. Should
%be two element vector
%averageSpikes: mean firing rate based on baseline period
%averageSTD: standard deviation of firing rate
%inputRaster: selected raster. Should be nx1, with only timestamps
%zLimit: limit above which things are considered significant. should be
%scalar = 3
%trialNum = number of trials in the selected inputRaster

%Outputs
%respStore: output with response information. 



function [respStore] = ...
    functionBasicResponseSignificance(smoothingBins,defaultBins,calcWindow,...
    averageRate,averageSTD,inputRaster,zLimit,trialNum);

%calculates histogram centers based on calc window (2 element vector) and
%defaultBins (which should be 0.001)
defaultBinCenters = [calcWindow(1)+defaultBins/2:defaultBins:calcWindow(2)-defaultBins/2];
%calculates histogram based on inputRasters and defaultBinCenters
rawHist = hist(inputRaster,defaultBinCenters)
%normalizes based on bin size and number of repetitions
zHist = (rawHist/defaultBins/trialNum - averageRate) / averageSTD;
singleStandard = ((1/defaultBins/trialNum) - averageRate) / averageSTD; %calculates the zscore of a single spike
if singleStandard > zLimit
    zLimit = singleStandard;
end
smoothedHistHolder = zeros(length(defaultBinCenters),length(smoothingBins));


respStore = cell(length(smoothingBins),1);

for respCounter = 1:size(smoothingBins,2)
    %first smooth. if statement is for unsmoothed curve
    smoothSpan = smoothingBins(respCounter) / defaultBins;
    if smoothSpan == 1;
        smoothedHist = zHist';
    else
        smoothedHist = smooth(zHist,smoothSpan,'lowess');
    end
    smoothedHistHolder(:,respCounter) = smoothedHist;
    %next, find if any threshold crossings
    if ~isempty(find(smoothedHist > zLimit,1))
        %pull all values above zLimit
        zSignalIndex = find(smoothedHist > zLimit);
        zSigDiff = diff(zSignalIndex);
        bigDiffFinder = [0;find(zSigDiff > 1);length(zSignalIndex)];
        respNum = length(bigDiffFinder)-1;

        respStore{respCounter} = zeros(respNum,7);
        %stores information about responses. 1 and 2 are onset and offset
        %of responses, 3 is peak of response in z score, 4 is location of
        %said peak, 5 and 6 are the same, but for the raw histogram,
        for storeCounter = 1:respNum
            respStore{respCounter}(storeCounter,1) = zSignalIndex(bigDiffFinder(storeCounter)+1);
            respStore{respCounter}(storeCounter,2) = zSignalIndex(bigDiffFinder(storeCounter+1));
            [M I] = max(smoothedHist(zSignalIndex(bigDiffFinder(storeCounter)+1):zSignalIndex(bigDiffFinder(storeCounter+1))));
            respStore{respCounter}(storeCounter,3) = M;
            respStore{respCounter}(storeCounter,4) = I + zSignalIndex(bigDiffFinder(storeCounter)+1);
            [M I] = max(rawHist(zSignalIndex(bigDiffFinder(storeCounter)+1):zSignalIndex(bigDiffFinder(storeCounter+1))));
            respStore{respCounter}(storeCounter,5) = M;
            respStore{respCounter}(storeCounter,6) = I + zSignalIndex(bigDiffFinder(storeCounter)+1);
            respStore{respCounter}(storeCounter,7) = zLimit;
        end
    end
end


end