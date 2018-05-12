
%This code is meant to pull EDR files, and extract them based on a set of
%TTL pulses. This system is reliant on only using TTL output 1, which is
%being fed into the NIDAQ system directly, for timing and syncing the two
%signals

%Inputs:

%s: structured array for data storage
%fileName: full name of file, with file extension
%colTime: column for time variable
%colTTL: column for TTL inputs (can be multiple)
%colPiezo: column for piezo inputs

function [EDROut] = functionEDRPull(fileName,colTime,colPiezo,colTTL);

%parameters
% colTime = 1;
% colTTL = 3;
% colPiezo = 2;
threshTTL = 100; %threshold used to pull TTL data.

[data, h] = import_edr(fileName);
%this outputs two data: data, which is a nxm array of n inputs and m
%timepoints, with the first column as time, following columns as additional
%information (in this case, 2 is piezo, 3 is TTL)

%h is a structured array with information about settings for the recording
%session. Discard h for now

%pull out discrete TTLs and save them in 
for edrInd = 1:length(colTTL)
    %pull TTL data
    edrTTLs = [];
    edrTTLs = diff(data(:,colTTL(edrInd)));
    %find all points above the threshold, fix the fudge factor
    ttlFinder = find(edrTTLs>threshTTL);
    ttlFinder = ttlFinder+1;
    %save onset times
    onsetSaver{edrInd} = data(ttlFinder,1);
end

%now lets look at the piezo signal and do some filtering. Looks like most
%of the locomotor signal is at 5-50 Hz. 
Wn = [0.005 0.04]; % this establishes a bandpass filter. based on data from spectrogram 180511
n = 1000; % 1000th order filter (slower? but 100-order was too low)
b = fir1(n, Wn); 
filteredData = filtfilt(b,1,data(:,colPiezo));

%try pulling spectrogram as well. 

[sig,w,t] = spectrogram(data(:,colPiezo),128,120,128,1e3); %this makes spectrograms using 128 sample time bins. Produces 65 frequency bins

desiredW = [2:3];

for edrInd = 1:length(desiredW)
    exampleData = sig(desiredW(edrInd),:);
    realData = real(exampleData);
    imData = imag(exampleData);
    mag(edrInd,:) = realData.^2 + imData.^2;
end

meanMag = mean(mag);

EDROut = struct;
EDROut.TTLs = onsetSaver;
EDROut.Piezo = filteredData;
EDROut.PiezoPower = [t',meanMag'];
EDROut.Time = data(:,1);
EDROut.RawPiezo = data(:,colPiezo);
end







