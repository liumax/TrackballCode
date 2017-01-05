%% This code is meant to produce averaged LFP traces relative to a targeted event. 

%% Inputs:
% lfpTimes: in real time (seconds)
%lfpTrace: LFP signal, should be vector of same length as lfpTimes
%lfpWindow: window of analysis in seconds
%lfpSampleRate: lfp samples per second
%alignTimes: times I want to align to, in seconds. 

%% Outputs
%lfpAverage: average lfp trace
%lfpSTD: standard deviation.


function [lfpAverage,lfpSTD] = functionLFPpsth(lfpTimes,lfpTrace,lfpWindow,lfpSampleRate,alignTimes);

%calculate the number of LFP samples in the targeted window
sampWindow = round((lfpWindow(2)-lfpWindow(1))*lfpSampleRate);
%allocate array
lfpHolder = zeros(sampWindow,length(alignTimes));
%fill array with values
for lfpInd = 1:length(alignTimes)
    finder = find(lfpTimes>(alignTimes(lfpInd)+lfpWindow(1)),1,'first');
    lfpHolder(:,lfpInd) = lfpTrace(finder:finder+sampWindow-1);
end
%calculate average
lfpAverage = mean(lfpHolder,2);
lfpSTD = std(lfpHolder,0,2);

end





















