
%This function is meant to extract the first spike timing of a set of
%aligned spikes. This code will also extract the number of spikes per set
%bin, defined by the timingWindows. 
%The inputs are: 
%timingWindows: a vector of minimum of 2 values, no maximum. should be
%monotonically increasing values. 
%alignedSpikes: an m x n matrix, with m spikes, n parameters. first column
%should always be time.
%numTrials: number of trials performed, so as to calculate probability of
%spiking. Should be a scalar
%trialColumn: column in alignedSpikes that contains trial numbers. This is
%to allow the code to find the correct trials. 
%trialNumbers: actual numbers of the trials in question

%Outputs: 
%firstSpikeTimes: an m x n matrix, with m trials, n tuning windows. Zeros
%indicate no spike. 
%firstSpikeStats: a 4 x n maxtrix, with n tuning windows. First row is
%mean, second row is standard deviation, third row is probability of
%response, fourth row is number of trials with spiking responses.
%binSpike is spikes binned over the entire period, rather than the first
%spike itself.

function [firstSpikeTimes,firstSpikeStats,binSpikeTimes,binSpikeStats] = functionBasicFirstSpikeTiming(timingWindows,alignedSpikes,numTrials,trialColumn,trialNumbers);
%computes the actual number of windows
windowLength = (length(timingWindows)-1);
%allocates a matrix for spikes
firstSpikeTimes = zeros(numTrials,windowLength);
binSpikeTimes = zeros(numTrials,windowLength);
%pulls spikes in double for loop
for counter1 = 1:windowLength
    for counter2 = 1:numTrials
        %if statement is to ensure that the find function will not return
        %an empty array. 
        if ~isempty(find(alignedSpikes(:,trialColumn) == trialNumbers(counter2) & ...
                alignedSpikes(:,1) < timingWindows(counter1+1) & ...
                alignedSpikes(:,1) > timingWindows(counter1)))
            %save the first spike found.
            firstSpikeTimes(counter2,counter1) = alignedSpikes(find(alignedSpikes(:,trialColumn) == trialNumbers(counter2) & ...
                alignedSpikes(:,1) < timingWindows(counter1+1) & ...
                alignedSpikes(:,1) > timingWindows(counter1),1,'first'),1);
            %saves the number of spikes in the specific timing window. 
            binSpikeTimes(counter2,counter1) = size(find(alignedSpikes(:,trialColumn) == trialNumbers(counter2) & ...
                alignedSpikes(:,1) < timingWindows(counter1+1) & ...
                alignedSpikes(:,1) > timingWindows(counter1)),1);
        end
    end
end
%now calculate summary statistics: mean, std, probability of response, and number of responses.
firstSpikeStats = zeros(4,windowLength);
binSpikeStats = zeros(2,windowLength);
%allocates a placeholder for realSpikes, which will hold all spikes, and
%sort out the zeros.
realSpikes = zeros(1,numTrials);
for counter1 = 1:windowLength
    %imports first spike data from firstSpikeTimes
    realSpikes = firstSpikeTimes(:,counter1);
    %weeds out zeros, which indicate a spike was not fired.
    realSpikes(realSpikes == 0) = [];
    if ~isempty(realSpikes)
        %calculates and stores mean first spike timing
        firstSpikeStats(1,counter1) = mean(realSpikes);
        binSpikeStats(1,counter1) = mean(binSpikeTimes(:,counter1));
        %calculates and stores standard deviation
        firstSpikeStats(2,counter1) = std(realSpikes);
        binSpikeStats(2,counter1) = std(binSpikeTimes(:,counter1));
        %calculate and stores the probability of spiking in the window at all
        firstSpikeStats(3,counter1) = length(realSpikes)/numTrials;
        %saves number of spike responses
        firstSpikeStats(4,counter1) = length(realSpikes);
    end
end

end