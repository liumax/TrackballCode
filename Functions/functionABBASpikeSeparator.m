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

function [s] = functionABBASpikeSeparator(s,...
   desigNames,spikeNames,soundNames);
%First, I want to pull all the time periods.
numDivs = length(spikeNames);

timesTuningFirst = s.TimePeriods.(soundNames{1});
timesTuningSecond = s.TimePeriods.(soundNames{2});
timesPairing = s.TimePeriods.(soundNames{3});

%stores spikes as separated spike times. 
for i = 1:size(desigNames,2);
    for j = 1:numDivs
        s.(desigNames{i}).(spikeNames{j}) = ...
            s.(desigNames{i}).SpikeTimes(...
            s.(desigNames{i}).SpikeTimes>s.TimePeriods.(soundNames{j})(1) &...
            s.(desigNames{i}).SpikeTimes<s.TimePeriods.(soundNames{j})(2));
    end
end


end
