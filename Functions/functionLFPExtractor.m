% Lets write a function that goes, extracts the LFP file, and outputs the
% times and data in real time (seconds).

%% Inputs
%targetName: name for the target LFP file


%% Outputs
%lfpTrace: actual LFP trace
%lfpTimes: corresponding time points

function [lfpTrace lfpTimes lfpSampleRate] = functionLFPExtractor(targetName);

%pull lfp structure
lfp = readTrodesExtractedDataFile(targetName);
%extract trace
lfpTrace = lfp.fields.data;
%extract parameters from trace
lengthTrace = length(lfpTrace);
sampleRate = lfp.clockrate;
lfpDecimation = lfp.decimation;
%calculate time step for each LFP sample
timeStep = 1/sampleRate*lfpDecimation;
%calculate first and last timepoints
firstTime = lfp.first_timestamp/sampleRate;
endTime = lfp.first_timestamp/sampleRate+((lengthTrace-1)*timeStep);
%determine lfpTimes that correspond to lfpTrace
lfpTimes = firstTime:timeStep:endTime;

%calculate lfp sample rate
lfpSampleRate = sampleRate/lfpDecimation;
end
