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

function [s] = functionPairingSpikeSeparator(s,...
   desigNames,spikeNames,soundNames);

%stores spikes as separated spike times. 
for i = 1:size(desigNames,2);
    for j = 1:length(spikeNames)
        s.(desigNames{i}).(spikeNames{j}) = ...
            s.(desigNames{i}).SpikeTimes(...
            s.(desigNames{i}).SpikeTimes>s.TimePeriods.(soundNames{j})(1) &...
            s.(desigNames{i}).SpikeTimes<s.TimePeriods.(soundNames{j})(2));
    end
end


end
