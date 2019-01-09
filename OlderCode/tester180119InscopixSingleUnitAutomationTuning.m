%This code is meant to automate response detection for inscopix single unit
%data. 


%set parameters
rasterWindow = [-3 3];
binSize = 0.5;
alphaVal = 0.001;


%pull sound data
numTrials = length(s.SoundData.Delays);
numFreqs = length(s.SoundData.UniqueFrequencies);
uniqueFreqs = s.SoundData.UniqueFrequencies;

%determine number of units
numUnits = length(s.DesignationName);

for unitInd = 1:numUnits
    tester = s.(s.DesignationName{unitInd}).PointRasterTime;
    %make raster into hist with correct bin size
    histVector = [rasterWindow(1)+binSize/2:binSize:rasterWindow(2)-binSize/2];
    for subInd = 1:numTrials
        %check if any spikes
        testHolder = tester(tester(:,2) == subInd,1);
        testHolder(testHolder>rasterWindow(2)) = [];
        testHist = hist(testHolder,histVector);
        histHolder(:,subInd) = testHist;
    end

    baselineBins = [1 find(histVector<0,1,'last')];
    %grab baseline bins
    baseStore = histHolder(baselineBins(1):baselineBins(2),:);
    baseStore = reshape(baseStore,[],1);
    %this produces a single column vector of baseline bins with which we can
    %compare tone related bins. 

    %now store bins based on frequencies

    onsetBin = find(histVector > 0,1,'first');
    offsetBin = find(histVector > 1,1,'first');

    for subInd = 1:numFreqs
        %find target trials
        findFreqs = find(s.SoundData.Frequencies == uniqueFreqs(subInd));
        targetSet = histHolder(:,findFreqs);
        freqStore{subInd} = targetSet;
        freqBins(subInd,1) =mean(targetSet(onsetBin,:));
        freqBins(subInd,2) =mean(targetSet(offsetBin,:));
    end


    %now lets zscore freqBins
    zfreqBins = (freqBins-mean(baseStore))/std(baseStore);
    bigPreHolder(unitInd,:) = zfreqBins(:,1);
    bigPostHolder(unitInd,:) = zfreqBins(:,2);
    n.(s.DesignationName{unitInd}).zScoreOnsetOffset = zfreqBins;
end

























