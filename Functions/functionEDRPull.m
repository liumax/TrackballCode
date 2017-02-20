
%This code is meant to pull EDR files, and extract them based on a set of
%TTL pulses. This system is reliant on only using TTL output 1, which is
%being fed into the NIDAQ system directly, for timing and syncing the two
%signals

%Inputs:

%s: structured array for data storage
%fileName: full name of file, with file extension
%colTime: column for time variable
%colTTL: column for TTL inputs
%colPiezo: column for piezo inputs

function [s] = functionEDRPull(s,fileName,colTime,colTTL,colPiezo);

%parameters
% colTime = 1;
% colTTL = 3;
% colPiezo = 2;
threshTTL = 0.4; %threshold used to pull TTL data.

[data, h] = import_edr(fileName);
%this outputs two data: data, which is a nxm array of n inputs and m
%timepoints, with the first column as time, following columns as additional
%information (in this case, 2 is piezo, 3 is TTL)

%h is a structured array with information about settings for the recording
%session. Discard h for now

%pull TTL data
edrTTLs = data(:,colTTL);
%find all points above the threshold
ttlFinder = find(edrTTLs>threshTTL);
%find onsets of pulses
ttlFinderDiff = diff(ttlFinder);
onsetFinder = [1;find(ttlFinderDiff > 1)+1];
%finds associated NIDAQ time. NIDAQ time should always start at zero
ttlTimes = data(ttlFinder(onsetFinder),colTime); %ttl Times in ms
%pull indices for TTL onsets
ttlInds = ttlFinder(onsetFinder);

s.EDR.FullData = data;
s.EDR.TTLTimes = ttlTimes;


end







