
%Outputs: filtSig1(470 channel), filtSig2(405 channel), normSig (dF/F
%channel), spectTimes (time points)

%Inputs: data (the raw data structure from TDT). 

function [filtSig1,filtSig2,normSig,spectTimes] = functionPhotometryRawExtraction(data);

rx = data.streams.Fi1r.data(3,:);
Fs = data.streams.Fi1r.fs; % Sampling rate
disp('Loaded Raw Data')


freqRange = 100:5:400; % Frequencies to calculate spectrogram in Hz
winSize = 0.04; % Window size for spectrogram (sec)
spectSample = 0.005; % Step size for spectrogram (sec)
inclFreqWin = 3; % Number of frequency bins to average (on either side of peak freq)
filtCut = 300; % Cut off frequency for low pass filter of data

% Convert spectrogram window size and overlap from time to samples
spectWindow = 2.^nextpow2(Fs .* winSize);
spectOverlap = ceil(spectWindow - (spectWindow .* (spectSample ./ winSize)));
disp(['Calculating spectrum using window size ', num2str(spectWindow ./ Fs)])

% Create low pass filter fot final data
lpFilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',300,...
    'PassbandRipple',0.01, 'SampleRate',Fs);

% Calculate spectrogram
tic
[spectVals,spectFreqs,spectTimes]=spectrogram(rx,spectWindow,spectOverlap,freqRange,Fs);
spectAmpVals = double(abs(spectVals));
toc

% Find the two carrier frequencies
tic
avgFreqAmps = mean(spectAmpVals,2);
locs =findpeaks(double(avgFreqAmps),max(avgFreqAmps./10));
locs = locs.loc;

if length(locs) > 1;
    % Calculate signal at each frequency band
    sig1 = mean(abs(spectVals((locs(1)-inclFreqWin):(locs(1)+inclFreqWin),:)),1);
    sig2 = mean(abs(spectVals((locs(2)-inclFreqWin):(locs(2)+inclFreqWin),:)),1);
else
    %this is in case you shut off a laser for some reason
    disp('PhotometryRawExtraction: Failed to Extract Two Peaks')
    freqs = [211 330]; %these are hard coded in. 211 is for 470, 330 for 405
    [mins locs(1)] = min(abs(spectFreqs - freqs(1)));
    [mins locs(2)] = min(abs(spectFreqs - freqs(2)));
    
    sig1 = mean(abs(spectVals((locs(1)-inclFreqWin):(locs(1)+inclFreqWin),:)),1);
    sig2 = mean(abs(spectVals((locs(2)-inclFreqWin):(locs(2)+inclFreqWin),:)),1);
end

% Low pass filter the signals
filtSig1 = filtfilt(lpFilt,double(sig1));
filtSig2 = filtfilt(lpFilt,double(sig2));

% Convert signal to dF/F
normSig = (filtSig1 ./ filtSig2) ./ (mean(filtSig1 ./ filtSig2));



end
