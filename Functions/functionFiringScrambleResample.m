

%alignWindows should be an N x 2 array, with the limits of the window in
%each row, and different rows representing different windows

function [meanStore,stdStore] = functionFiringScrambleResample(firingTimes,alignTimes,alignWindows,scrambleReps);


%first things first, get the spike ISIs
spikeISI = diff(firingTimes);

numWindows = size(alignWindows,1);

meanStore = zeros(numWindows,scrambleReps);
% medianStore = zeros(numWindows,scrambleReps);
stdStore = zeros(numWindows,scrambleReps);

for shuffInd = 1:scrambleReps
    %make random index of length equal to spikeISI
    shuffleTimes = randperm(length(spikeISI));
    shuffleTimes = spikeISI(shuffleTimes);
    %generate the new scrambled spike train
    shuffSpikeTimes = zeros(length(firingTimes),1);
    shuffSpikeTimes(1) = firingTimes(1)+rand;
    shuffSpikeTimes(2:end) = cumsum(shuffleTimes)+shuffSpikeTimes(1);
    
    %now we want to pull out the windows that we care about
    %first, set up blank spaces for the things!
    
    valStore = zeros(numWindows,length(alignTimes));
    for alignInd = 1:length(alignTimes)
        subTimes = shuffSpikeTimes - alignTimes(alignInd);
        for winInd = 1:numWindows
            valStore(winInd,alignInd) = length(find(subTimes >= alignWindows(winInd,1) & subTimes <= alignWindows(winInd,2)));
        end
    end
    
    meanStore(:,shuffInd) = mean(valStore,2);
%     medianStore(:,shuffInd) = median(valStore,2);
    stdStore(:,shuffInd) = std(valStore,0,2);
    
end











































end
