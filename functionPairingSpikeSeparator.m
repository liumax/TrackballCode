%This code is meant to use time periods delineated by signal pulses to
%separate out chunks of spike times. This will save these spike files in
%the masterStruct file.

function [masterStruct] = functionPairingSpikeSeparator(masterStruct,...
   truncatedNames);
%First, I want to pull all the time periods.
timesBaseline = masterStruct.TimePeriods.Baseline;
timesTuningFirst = masterStruct.TimePeriods.TuningFirst;
timesPresentationFirst = masterStruct.TimePeriods.PresentationFirst;
timesPairing = masterStruct.TimePeriods.Pairing;
timesPresentationSecond = masterStruct.TimePeriods.PresentationSecond;
timesTuningSecond = masterStruct.TimePeriods.TuningSecond;

%stores spikes as separated spike times. 
for i = 1:size(truncatedNames,2);
    for j = 1:masterStruct.(truncatedNames{i}).Clusters
        masterStruct.(truncatedNames{i}).BaselineSpikes{j} = ...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}(...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}>timesBaseline(1) &...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}<timesBaseline(2));
        masterStruct.(truncatedNames{i}).TuningFirstSpikes{j} = ...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}(...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}>timesTuningFirst(1) &...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}<timesTuningFirst(2));
        masterStruct.(truncatedNames{i}).TuningSecondSpikes{j} = ...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}(...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}>timesTuningSecond(1) &...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}<timesTuningSecond(2));
        masterStruct.(truncatedNames{i}).PresentFirstSpikes{j} = ...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}(...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}>timesPresentationFirst(1) &...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}<timesPresentationFirst(2));
        masterStruct.(truncatedNames{i}).PresentSecondSpikes{j} = ...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}(...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}>timesPresentationSecond(1) &...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}<timesPresentationSecond(2));
        masterStruct.(truncatedNames{i}).PairingSpikes{j} = ...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}(...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}>timesPairing(1) &...
            masterStruct.(truncatedNames{i}).SpikeTimes{j}<timesPairing(2));
    end
end


end