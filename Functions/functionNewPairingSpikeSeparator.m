%This code is meant to use time periods delineated by signal pulses to
%separate out chunks of spike times. This will save these spike files in
%the master structure

%Inputs:
%s: master structure
%desigNames: names of units (1 x n)
%spikeNames: names for spike periods
%soundNames: names for sound periods

%Outputs: 
%s: master structure. Will include all spike data into the master structure
%under each unit, under the title 'SpikeTimes'

function [s] = functionNewPairingSpikeSeparator(s,...
   desigNames,spikeNames,soundNames);
%First, I want to pull all the time periods.
timesBaseline = s.TimePeriods.Baseline;
timesTuningFirst = s.TimePeriods.(names{1});
timesTuningSecond = s.TimePeriods.(names{2});
timesPairing = s.TimePeriods.(names{3});

%stores spikes as separated spike times. 
for i = 1:size(desigNames,2);
    s.(desigNames{i}).BaselineSpikes = ...
        s.(desigNames{i}).SpikeTimes(...
        s.(desigNames{i}).SpikeTimes>timesBaseline(1) &...
        s.(desigNames{i}).SpikeTimes<timesBaseline(2));
    s.(desigNames{i}).(spikeNames{1}) = ...
        s.(desigNames{i}).SpikeTimes(...
        s.(desigNames{i}).SpikeTimes>timesTuningFirst(1) &...
        s.(desigNames{i}).SpikeTimes<timesTuningFirst(2));
    s.(desigNames{i}).(spikeNames{2}) = ...
        s.(desigNames{i}).SpikeTimes(...
        s.(desigNames{i}).SpikeTimes>timesTuningSecond(1) &...
        s.(desigNames{i}).SpikeTimes<timesTuningSecond(2));
    s.(desigNames{i}).(spikeNames{3}) = ...
        s.(desigNames{i}).SpikeTimes(...
        s.(desigNames{i}).SpikeTimes>timesPairing(1) &...
        s.(desigNames{i}).SpikeTimes<timesPairing(2));
end


end
