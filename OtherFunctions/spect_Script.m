clear variables

load('CD170307A_170410_01.mat')

rx = data.streams.Fi1r.data(3,:);%1:100000); % take first 100k points
Fs = data.streams.Fi1r.fs; % Sampling rate

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
[pks,locs]=findpeaks(double(avgFreqAmps),'minpeakheight',max(avgFreqAmps./10));

% Calculate signal at each frequency band
sig1 = mean(abs(spectVals((locs(1)-inclFreqWin):(locs(1)+inclFreqWin),:)),1);
sig2 = mean(abs(spectVals((locs(2)-inclFreqWin):(locs(2)+inclFreqWin),:)),1);

% Low pass filter the signals
filtSig1 = filtfilt(lpFilt,double(sig1));
filtSig2 = filtfilt(lpFilt,double(sig2));

% Convert signal to dF/F
normSig = (filtSig1 ./ filtSig2) ./ (mean(filtSig1 ./ filtSig2));

% Read in TDT data
tdtTs = (1:length(data.streams.x70G.data))./data.streams.x70G.fs;
normSigTdt = (data.streams.x70G.data./data.streams.x05G.data) ./ nanmean(data.streams.x70G.data./data.streams.x05G.data);

% Create figure to plot
plotFig = figure('color','w');
imAx = subplot(3,3,1:3,'parent',plotFig);
sigAx = subplot(3,3,4:6,'parent',plotFig); hold(sigAx,'on');
normAx = subplot(3,3,7:9,'parent',plotFig); hold(normAx,'on');

% Plot spectrogram image
imagesc('XData',spectTimes,'YData',spectFreqs,'CData',spectAmpVals,'parent',imAx);

% Plot 405 and 470 signals (unfiltered)
plot(spectTimes,sig1,'color',[0.5 0.5 1],'linewidth',0.5,'parent',sigAx);
plot(spectTimes,sig2,'color',[0.5 1 0.5],'linewidth',0.5,'parent',sigAx);
% Plot filtered signals
plot(spectTimes,filtSig1,'color',[0 0 0.7],'linewidth',2,'parent',sigAx);
plot(spectTimes,filtSig2,'color',[0 0.7 0],'linewidth',2,'parent',sigAx);

% Plot TDT signals
plot(spectTimes,normSig,'color',[0 0 0.7],'linewidth',2,'parent',normAx); 
plot(tdtTs,normSigTdt,'color',[0 0.7 0],'linewidth',2,'parent',normAx);
legend(normAx,{'SFO Calc','TDT Calc'},'location','northwest');
% plot(spectTimes,lpFiltData,'color','r','linewidth',2,'parent',sigAx);
% plot(spectFreqs,avgFreqAmps,'parent',freqAx);
toc