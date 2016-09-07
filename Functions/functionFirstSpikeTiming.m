
%this code extracts the timing of the first spike (mean time of first
%spike), probability of first spike within a certain period, and timing of
%spikes in a set number of bins. want to include this to existing code.

function [s] = functionFirstSpikeTiming(totalTrialNum,rasterData,uniqueFreqs,uniqueDBs,toneDur,toneReps);

spikeProbWindow = [0 0.5*toneDur toneDur 1.5*toneDur]; %window over which spike responses will count towards probability of response

%finds all trials and records first spike time (if it exists)
firstSpikeTimes = zeros(totalTrialNum,3);

for counter1 = 1:totalTrialNum
    if ~isempty(find(rasterData(:,1) == counter1 & rasterData(:,2) >0,1,'first'))
        findFirstSpike = find(rasterData(:,1) == counter1 & rasterData(:,2) >0,1,'first');
        firstSpikeTimes(counter1,1) = rasterData(findFirstSpike,2);
        firstSpikeTimes(counter1,2) = rasterData(findFirstSpike,3);
        firstSpikeTimes(counter1,3) = rasterData(findFirstSpike,4);
    end
end

firstSpikeTimes(firstSpikeTimes(:,3) == 0,:) = [];

firstSpikeStats = zeros(length(uniqueFreqs),length(uniqueDBs),5);
%finds timing of the first spike in three time bins
for counter1 = 1:size(uniqueFreqs,1)
    for counter2 = 1:size(uniqueDBs,1)
        if ~isempty(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) & firstSpikeTimes(:,3) == uniqueDBs(counter2),1))
            firstSpikeStats(counter1,counter2,1) = mean(firstSpikeTimes(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) & firstSpikeTimes(:,3) == uniqueDBs(counter2)),1));
            firstSpikeStats(counter1,counter2,2) = std(firstSpikeTimes(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) & firstSpikeTimes(:,3) == uniqueDBs(counter2)),1));
        end
        if ~isempty(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) & firstSpikeTimes(:,3) == uniqueDBs(counter2) & firstSpikeTimes(:,1) < spikeProbWindow(2) & firstSpikeTimes(:,1) > spikeProbWindow(1),1))
            firstSpikeStats(counter1,counter2,3) = size(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) &...
            firstSpikeTimes(:,3) == uniqueDBs(counter2) & firstSpikeTimes(:,1) < spikeProbWindow(2) & firstSpikeTimes(:,1) > spikeProbWindow(1)),1)/toneReps;
        end
        if ~isempty(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) & firstSpikeTimes(:,3) == uniqueDBs(counter2) & firstSpikeTimes(:,1) < spikeProbWindow(3) & firstSpikeTimes(:,1) > spikeProbWindow(2),1))
            firstSpikeStats(counter1,counter2,4) = size(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) &...
            firstSpikeTimes(:,3) == uniqueDBs(counter2) & firstSpikeTimes(:,1) < spikeProbWindow(3) & firstSpikeTimes(:,1) > spikeProbWindow(2)),1)/toneReps;
        end
        if ~isempty(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) & firstSpikeTimes(:,3) == uniqueDBs(counter2) & firstSpikeTimes(:,1) < spikeProbWindow(4) & firstSpikeTimes(:,1) > spikeProbWindow(3),1))
            firstSpikeStats(counter1,counter2,5) = size(find(firstSpikeTimes(:,2) == uniqueFreqs(counter1) &...
            firstSpikeTimes(:,3) == uniqueDBs(counter2) & firstSpikeTimes(:,1) < spikeProbWindow(4) & firstSpikeTimes(:,1) > spikeProbWindow(3)),1)/toneReps;
        end
        
    end
end

%% looking at spikes in general. this way i can look at multiple windows?

toneSpikeTimes = zeros(totalTrialNum,5);

%pulls out all the spikes per trial within defined bins.
for counter1 = 1:totalTrialNum
    if ~isempty(find(rasterData(:,1) == counter1 & rasterData(:,2) >0& rasterData(:,2) < spikeProbWindow(4),1))
        findToneSpikes = find(rasterData(:,1) == counter1 & rasterData(:,2) >0 & rasterData(:,2) < spikeProbWindow(4));
        toneSpikes = rasterData(findToneSpikes,2);
        toneSpikeTimes(counter1,1) = size(toneSpikes(toneSpikes<spikeProbWindow(2)),1);
        toneSpikeTimes(counter1,2) = size(toneSpikes(toneSpikes<spikeProbWindow(3)&toneSpikes>spikeProbWindow(2)),1);
        toneSpikeTimes(counter1,3) = size(toneSpikes(toneSpikes<spikeProbWindow(4)&toneSpikes>spikeProbWindow(3)),1);
        toneSpikeTimes(counter1,4) = rasterData(findToneSpikes(1),3);
        toneSpikeTimes(counter1,5) = rasterData(findToneSpikes(1),4);
    end
end

%clears empty stuff
toneSpikeTimes(toneSpikeTimes(:,5) == 0,:) = [];

toneSpikeStats = zeros(length(uniqueFreqs),length(uniqueDBs),3);
for counter1 =1:size(uniqueFreqs,1)
    for counter2 = 1:size(uniqueDBs,1)
        if ~isempty(find(toneSpikeTimes(:,4) == uniqueFreqs(counter1)&toneSpikeTimes(:,5) == uniqueDBs(counter2),1))
            toneSpikeStats(counter1,counter2,1) = size(find(toneSpikeTimes(:,4) == uniqueFreqs(counter1) & toneSpikeTimes(:,5) == uniqueDBs(counter2) & toneSpikeTimes(:,1) > 0),1)/toneReps;
            toneSpikeStats(counter1,counter2,2) = size(find(toneSpikeTimes(:,4) == uniqueFreqs(counter1) & toneSpikeTimes(:,5) == uniqueDBs(counter2) & toneSpikeTimes(:,2) > 0),1)/toneReps;
            toneSpikeStats(counter1,counter2,3) = size(find(toneSpikeTimes(:,4) == uniqueFreqs(counter1) & toneSpikeTimes(:,5) == uniqueDBs(counter2) & toneSpikeTimes(:,3) > 0),1)/toneReps;
        end
    end
end


s = struct;
s.firstSpikeTimes = firstSpikeTimes;
s.firstSpikeStats = firstSpikeStats;
s.toneSpikesTimes = toneSpikeTimes;
s.toneSpikeStats = toneSpikeStats;



end