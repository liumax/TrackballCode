
%this code extracts the timing of the first spike (mean time of first
%spike), probability of first spike within a certain period, and timing of
%spikes in a set number of bins. want to include this to existing code.

function [s] = functionPulseFirstSpikeTiming(totalTrialNum,rasterData,pipPos,pipReps,toneDur,blockNum);

spikeProbWindow = [0 0.5*toneDur toneDur 1.5*toneDur]; %window over which spike responses will count towards probability of response

%finds all trials and records first spike time (if it exists)
firstSpikeTimes = zeros(totalTrialNum,4);

for counter1 = 1:totalTrialNum
    if ~isempty(find(rasterData(:,1) == counter1 & rasterData(:,2) >0,1,'first'))
        findFirstSpike = find(rasterData(:,1) == counter1 & rasterData(:,2) >0,1,'first');
        firstSpikeTimes(counter1,1) = rasterData(findFirstSpike,1); %pulls trial number
        firstSpikeTimes(counter1,2) = rasterData(findFirstSpike,2); %pulls timing
    end
end

firstSpikeTimes(firstSpikeTimes(:,1) == 0,:) = [];
firstSpikeTimes(:,3)= fix((firstSpikeTimes(:,1)+pipReps-1)/pipReps); %determines which rep this is part of. adding 14 allows for this to work properly
firstSpikeTimes(:,4)= rem(firstSpikeTimes(:,1),pipReps); %determines which pip in sequence
firstSpikeTimes(firstSpikeTimes(:,4) == 0,4) = pipReps; %fixed glitch with rem

firstSpikeStats = zeros(length(pipPos),5);
%here, output dimension 1 is average first spike, 2 is std of that first
%spike, and 3:5 are probability of first spike in three time bins.

%finds timing of the first spike in three time bins, averaged across all
%repetitions of tone presentation
for counter1 = 1:size(pipPos,2)
    if ~isempty(find(firstSpikeTimes(:,4) == pipPos(counter1),1))
        firstSpikeStats(counter1,1) = mean(firstSpikeTimes(find(firstSpikeTimes(:,4) == pipPos(counter1)),2));
        firstSpikeStats(counter1,2) = std(firstSpikeTimes(find(firstSpikeTimes(:,4) == pipPos(counter1)),2));
    end
    if ~isempty(find(firstSpikeTimes(:,4) == pipPos(counter1) & firstSpikeTimes(:,2) < spikeProbWindow(2) & firstSpikeTimes(:,2) > spikeProbWindow(1),1))
        firstSpikeStats(counter1,3) = size(find(firstSpikeTimes(:,4) == pipPos(counter1) & firstSpikeTimes(:,2) < spikeProbWindow(2) & firstSpikeTimes(:,2) > spikeProbWindow(1)),1)/blockNum;
    end
    if ~isempty(find(firstSpikeTimes(:,4) == pipPos(counter1) & firstSpikeTimes(:,2) < spikeProbWindow(3) & firstSpikeTimes(:,2) > spikeProbWindow(2),1))
        firstSpikeStats(counter1,4) = size(find(firstSpikeTimes(:,4) == pipPos(counter1) & firstSpikeTimes(:,2) < spikeProbWindow(3) & firstSpikeTimes(:,2) > spikeProbWindow(2)),1)/blockNum;
    end
    if ~isempty(find(firstSpikeTimes(:,4) == pipPos(counter1) & firstSpikeTimes(:,2) < spikeProbWindow(4) & firstSpikeTimes(:,2) > spikeProbWindow(3),1))
        firstSpikeStats(counter1,5) = size(find(firstSpikeTimes(:,4) == pipPos(counter1) & firstSpikeTimes(:,2) < spikeProbWindow(4) & firstSpikeTimes(:,2) > spikeProbWindow(3)),1)/blockNum;
    end
end

firstSpikeMatrix = zeros(pipReps,blockNum);
for counter1 = 1:size(firstSpikeTimes,1)
    firstSpikeMatrix(firstSpikeTimes(counter1,4),firstSpikeTimes(counter1,3)) = firstSpikeTimes(counter1,2);
end

%% looking at spikes in general. this way i can look at multiple windows?

toneSpikeTimes = zeros(totalTrialNum,5);

%pulls out all the spikes per trial within defined bins.
for counter1 = 1:totalTrialNum
    if ~isempty(find(rasterData(:,1) == counter1 & rasterData(:,2) > 0& rasterData(:,2) < spikeProbWindow(4),1))
        findToneSpikes = find(rasterData(:,1) == counter1 & rasterData(:,2) >0 & rasterData(:,2) < spikeProbWindow(4));
        toneSpikes = rasterData(findToneSpikes,2);
        toneSpikeTimes(counter1,1) = size(toneSpikes(toneSpikes<spikeProbWindow(2)),1);
        toneSpikeTimes(counter1,2) = size(toneSpikes(toneSpikes<spikeProbWindow(3)&toneSpikes>spikeProbWindow(2)),1);
        toneSpikeTimes(counter1,3) = size(toneSpikes(toneSpikes<spikeProbWindow(4)&toneSpikes>spikeProbWindow(3)),1);
        toneSpikeTimes(counter1,4) = counter1;
        toneSpikeTimes(counter1,5) = fix((toneSpikeTimes(counter1,4)+pipReps - 1)/pipReps); %determines which rep this is part of
        toneSpikeTimes(counter1,6) = rem(toneSpikeTimes(counter1,4),pipReps); %determines which pip in sequence
    end
end

%clears empty stuff
toneSpikeTimes(toneSpikeTimes(:,4) == 0,:) = [];
toneSpikeTimes(toneSpikeTimes(:,6) == 0,6) = 15;

toneSpikeStats = zeros(length(pipPos),3);
for counter1 =1:size(pipPos,2)
    if ~isempty(find(toneSpikeTimes(:,6) == pipPos(counter1),1))
        toneSpikeStats(counter1,1) = size(find(toneSpikeTimes(:,6) == pipPos(counter1) & toneSpikeTimes(:,1) > 0),1)/blockNum;
        toneSpikeStats(counter1,2) = size(find(toneSpikeTimes(:,6) == pipPos(counter1) & toneSpikeTimes(:,2) > 0),1)/blockNum;
        toneSpikeStats(counter1,3) = size(find(toneSpikeTimes(:,6) == pipPos(counter1) & toneSpikeTimes(:,3) > 0),1)/blockNum;
    end
end




s = struct;
s.firstSpikeTimes = firstSpikeTimes;
s.firstSpikeStats = firstSpikeStats;
s.toneSpikesTimes = toneSpikeTimes;
s.toneSpikeStats = toneSpikeStats;
s.firstSpikeMatrix = firstSpikeMatrix;


end