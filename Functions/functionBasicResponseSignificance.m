
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
rawHist = hist(inputRaster,defaultBinCenters);
%normalizes based on bin size and number of repetitions
zHist = (rawHist/defaultBins/trialNum - averageRate) / averageSTD;
singleStandard = ((1/defaultBins/trialNum) - averageRate) / averageSTD; %calculates the zscore of a single spike
if singleStandard > zLimit
    zLimit = zeros(length(smoothingBins),1);
    exampleSmooth = zeros(50,1);
    exampleSmooth(25) = 1;
    for zCount = 1:length(smoothingBins)
        zLimit(zCount) = max(smooth(exampleSmooth,smoothingBins(zCount)/defaultBins,'lowess'))*singleStandard*1.2;
    end
else
    oneHolder = ones(length(smoothingBins),1);
    zLimit = zLimit*oneHolder;
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
    if ~isempty(find(smoothedHist > zLimit(respCounter),1))
        %pull all values above zLimit
        zSignalIndex = find(smoothedHist > zLimit(respCounter));
        zSigDiff = diff(zSignalIndex);
        pointEnd = [zSignalIndex(zSigDiff>1);zSignalIndex(end)];
        pointStart = [zSignalIndex(1);zSignalIndex(find(zSigDiff>1)+1)];
        pointsCheck = length(pointEnd) - length(pointStart);
        if pointsCheck ~= 0
            pointsCheck
            error('Incorrect number of starts and stops for response significance')
        end
        pointDiff = pointEnd - pointStart;
        pointEnd(pointDiff <= 1) = [];
        pointStart(pointDiff <= 1) = [];
        respNum = length(pointEnd);
        
        respStore{respCounter} = zeros(respNum,7);
        for storeCounter = 1:respNum
            respStore{respCounter}(storeCounter,1) = pointStart(storeCounter);
            respStore{respCounter}(storeCounter,2) = pointEnd(storeCounter);
            [M I] = max(smoothedHist(pointStart(storeCounter):pointEnd(storeCounter)));
            respStore{respCounter}(storeCounter,3) = M;
            respStore{respCounter}(storeCounter,4) = I + pointStart(storeCounter);
            [M I] = max(rawHist(pointStart(storeCounter):pointEnd(storeCounter)));
            respStore{respCounter}(storeCounter,5) = M;
            respStore{respCounter}(storeCounter,6) = I + pointStart(storeCounter);
            respStore{respCounter}(storeCounter,7) = zLimit(respCounter);
        end
%         bigDiffFinder = [0;find(zSigDiff > 1);length(zSignalIndex)];
% %         bigDiffCheck = find(diff(bigDiffFinder)==1)+1;
% %         bigDiffFinder((bigDiffCheck == 1)+1) = [];
%         respNum = length(bigDiffFinder)-1;
% %PROBELM: DOES NOT DETECT OFFSETS FROM THRESHOLD CROSSINGS< SO THIS WILL
% %MISREPRESENT DURATION OF RESPONSES. ALSO DOES NOT ELIMINATE SINGLE POINT
% %RESPONSES< WHICH SHOULD BE GOTTEN RID OF
%         respStore{respCounter} = zeros(respNum,7);
%         %stores information about responses. 1 and 2 are onset and offset
%         %of responses, 3 is peak of response in z score, 4 is location of
%         %said peak, 5 and 6 are the same, but for the raw histogram,
%         for storeCounter = 1:respNum
%             respStore{respCounter}(storeCounter,1) = zSignalIndex(bigDiffFinder(storeCounter)+1);
%             respStore{respCounter}(storeCounter,2) = zSignalIndex(bigDiffFinder(storeCounter+1));
%             [M I] = max(smoothedHist(zSignalIndex(bigDiffFinder(storeCounter)+1):zSignalIndex(bigDiffFinder(storeCounter+1))));
%             respStore{respCounter}(storeCounter,3) = M;
%             respStore{respCounter}(storeCounter,4) = I + zSignalIndex(bigDiffFinder(storeCounter)+1);
%             [M I] = max(rawHist(zSignalIndex(bigDiffFinder(storeCounter)+1):zSignalIndex(bigDiffFinder(storeCounter+1))));
%             respStore{respCounter}(storeCounter,5) = M;
%             respStore{respCounter}(storeCounter,6) = I + zSignalIndex(bigDiffFinder(storeCounter)+1);
%             respStore{respCounter}(storeCounter,7) = zLimit(respCounter);
%         end
    end
end


end