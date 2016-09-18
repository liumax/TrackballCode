%This function is meant to be a basic checker for DIO information. This
%will specifically detect onsets for DIO inputs.
%Inputs:
%DIOData: should be a structured array, output of
%readTrodesExtractedDataFile. 
%sampleRate: should be scalar, 30,000

%Outputs:
%dioTimes: true times of DIO onsets
%dioTimeDiff: differences between dioTimes
function [dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,sampleRate);

%extracts dioTimes. Converts to double to avoid issues downstream. Divides
%by sample rate to give actual time, not sample number
dioTimes = double(DIOData.fields(1).data)/sampleRate;
%extracts states. States should be 0 or 1.
dioState = double(DIOData.fields(2).data);

dioDiff = find(diff(dioState)==1)+1; %finds all points where there is a positive diff. +1 fudge factor for diff()
dioHigh = find(dioState == 1); %finds all points at which state is 1

dioTrue = intersect(dioDiff,dioHigh); %finds intersect of dioDiff and dioHigh, which should only be onsets, not offsets or repeat DIO high states
%what i want is times, so extract dioTimes using dioTrue index.
dioTimes = dioTimes(dioTrue);
%also calculate differences between DIO input times
dioTimeDiff = diff(dioTimes);

end