%This function is meant to calculate overall firing rates of multiple
%periods of the task. 

%Inputs: 
%s: master structure
%desigNames: names of units
% spikeNames: names of field holding spikes in desigNames
% soundNames: names of sound periods.

%Outputs:
%s: master structure. Inserts information into s.unitName.BaselineSpikes
%and s.UnitName.OverallFiringRates.

function [s] = functionABBAOverallRate(s,desigNames,...
    spikeNames,soundNames);
%pulls average firing rate from baseline period.
for i = 1:size(desigNames,2)
    for j = 1:length(spikeNames)
        spikeTimes = s.(desigNames{i}).(spikeNames{j});
        numSpikes = size(spikeTimes,1);
        timePeriod = (s.TimePeriods.(soundNames{j})(2) - s.TimePeriods.(soundNames{j})(1));
        overallFireHolder(j) = numSpikes/timePeriod;
    end
    s.(desigNames{i}).OverallFiringRates = overallFireHolder;
end

end