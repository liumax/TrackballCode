%plotting function that plots out information for simple tuning curves.

function [s] = functionPairingMasterPlot(i,s,...
    desigNames,params,fileName,soundNames,figName);

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

histBinVector = [(params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur) + ...
params.histBin/2:params.histBin:(params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur)-params.histBin/2];

spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
%% plots average waveform
subplot(4,6,1)
hold on
plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
title({fileName;desigNames{i};strcat('AverageFiringRate:',num2str(mean(s.(desigNames{i}).OverallFiringRates)))});
set(0, 'DefaulttextInterpreter', 'none')
%% plots ISI
subplot(4,6,2)
hist(s.(desigNames{i}).ISIGraph,1000)
histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
line([params.rpvTime params.rpvTime],[0 histMax],'LineWidth',1,'Color','red')
xlim(params.clusterWindow)
title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
    strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))});

%% plots histogram
subplot(4,3,4)
%plots histogram from first tuning.
plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistograms,'k','LineWidth',1)
hold on
plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistograms,'r','LineWidth',1)
%plot significant positive values for first tuning
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Centers(...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
    'k*')
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Centers(...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
    'ko')
%plot negative values for first tuning
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Centers(...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
    'b*')
%plot significant positive values for second tuning
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Centers(...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
    'r*')
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Centers(...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
    'ro')
%plot negative values for second tuning
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Centers(...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
    'm*')

%draws in tone!
plot([0 0],[ylim],'r');
plot([s.SoundData.(soundNames{1}).ToneDur ...
    s.SoundData.(soundNames{1}).ToneDur],[ylim],'r');
%sets xlimits to avoid awkward graphs
xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur...
    params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
title('Histogram K before B after')
%% plot rasters organized by freq for first tuning
subplot(4,3,7)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllRasters(:,1),...
    s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllRasters(:,3),'k.','markersize',4)
hold on
plot([0 0],[ylim],'b');
plot([s.SoundData.(soundNames{1}).ToneDur s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
rasterFreqLines1 = zeros(size(s.SoundData.(soundNames{1}).UniqueFreqs,1),2);
rasterFreqLines1(:,1) = s.SoundData.(soundNames{1}).ToneReps*...
    size(s.SoundData.(soundNames{1}).UniqueDBs,1):s.SoundData.(soundNames{1}).ToneReps*...
    size(s.SoundData.(soundNames{1}).UniqueDBs,1):length(s.SoundData.(soundNames{1}).Frequencies);
rasterFreqLines1(:,2) = s.SoundData.(soundNames{1}).UniqueFreqs;
for k = 1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)
    plot(params.rasterWindow*s.SoundData.(soundNames{1}).ToneDur,[rasterFreqLines1(k,1) rasterFreqLines1(k,1)],'g','LineWidth',1)
end
set(gca,'YTick',rasterFreqLines1(:,1)-s.SoundData.(soundNames{1}).ToneReps/2*size(s.SoundData.(soundNames{1}).UniqueDBs,1));
set(gca,'YTickLabel',rasterFreqLines1(:,2));
set(gca,'Ydir','reverse')
ylim([0 length(s.SoundData.(soundNames{1}).Frequencies)])
xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
title({'Pre-Pairing Rasters';'dB Increase in Descending Order'})
%% plot rasters organized by freq for second tuning
subplot(4,3,10)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllRasters(:,1),...
    s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllRasters(:,3),'k.','markersize',4)
hold on
plot([0 0],[ylim],'b');
plot([s.SoundData.(soundNames{2}).ToneDur s.SoundData.(soundNames{2}).ToneDur],[ylim],'b');
rasterFreqLines2 = zeros(size(s.SoundData.(soundNames{2}).UniqueFreqs,1),2);
rasterFreqLines2(:,1) = s.SoundData.(soundNames{2}).ToneReps*...
    size(s.SoundData.(soundNames{2}).UniqueDBs,1):s.SoundData.(soundNames{2}).ToneReps*...
    size(s.SoundData.(soundNames{2}).UniqueDBs,1):length(s.SoundData.(soundNames{2}).Frequencies);
rasterFreqLines2(:,2) = s.SoundData.(soundNames{2}).UniqueFreqs;

for k = 1:size(s.SoundData.(soundNames{2}).UniqueFreqs,1)
    plot(params.rasterWindow*s.SoundData.(soundNames{2}).ToneDur,[rasterFreqLines2(k,1) rasterFreqLines2(k,1)],'g','LineWidth',1)
end
set(gca,'YTick',rasterFreqLines2(:,1)-s.SoundData.(soundNames{2}).ToneReps/2*size(s.SoundData.(soundNames{2}).UniqueDBs,1));
set(gca,'YTickLabel',rasterFreqLines2(:,2));
set(gca,'Ydir','reverse')
ylim([0 length(s.SoundData.(soundNames{2}).Frequencies)])
xlim([params.rasterWindow(1)*s.SoundData.(soundNames{2}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{2}).ToneDur])
title({'Post-Pairing Rasters';'dB Increase in Descending Order'})

%% Now plot tuning data based on binned/peak data
%first calculate standard errors
freqLength = length(s.SoundData.(soundNames{1}).UniqueFreqs);
dbLength = length(s.SoundData.(soundNames{1}).UniqueDBs);
toneReps = s.SoundData.(soundNames{1}).ToneReps;

preDataTone = zeros(freqLength,dbLength,toneReps);
postDataTone = zeros(freqLength,dbLength,toneReps);

preDataGen = zeros(freqLength,dbLength,toneReps);
postDataGen = zeros(freqLength,dbLength,toneReps);

for freqInd = 1:freqLength
    for dbInd = 1:dbLength
        preDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
        postDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
        preDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
        postDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
    end
end

preToneSTD = std(preDataTone,1,3)/sqrt(toneReps); 
postToneSTD = std(postDataTone,1,3)/sqrt(toneReps); 
preGenSTD = std(preDataGen,1,3)/sqrt(toneReps); 
postGenSTD = std(postDataGen,1,3)/sqrt(toneReps);

%do the same for compensated binned values
preDataToneComp = zeros(freqLength,dbLength,toneReps);
postDataToneComp = zeros(freqLength,dbLength,toneReps);

preDataGenComp = zeros(freqLength,dbLength,toneReps);
postDataGenComp = zeros(freqLength,dbLength,toneReps);

for freqInd = 1:freqLength
    for dbInd = 1:dbLength
        preDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
        postDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
        preDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
        postDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
    end
end
preToneCompSTD = std(preDataToneComp,1,3)/sqrt(toneReps); 
postToneCompSTD = std(postDataToneComp,1,3)/sqrt(toneReps); 
preGenCompSTD = std(preDataGenComp,1,3)/sqrt(toneReps); 
postGenCompSTD = std(postDataGenComp,1,3)/sqrt(toneReps);

%first plot binned responses, first to tone
%first pairing
subplot(4,3,2)
hold on
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreTone,'k','LineWidth',2)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreTone+preToneSTD,'k','LineWidth',1)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreTone-preToneSTD,'k','LineWidth',1)

plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreTone,'r','LineWidth',2)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreTone+postToneSTD,'r','LineWidth',1)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreTone-postToneSTD,'r','LineWidth',1)

%plot compensated values
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreToneComp,'k:','LineWidth',4)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreToneComp+preToneCompSTD,'k:','LineWidth',1)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreToneComp-preToneCompSTD,'k:','LineWidth',1)

plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreToneComp,'r:','LineWidth',4)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreToneComp+postToneCompSTD,'r:','LineWidth',1)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreToneComp-postToneCompSTD,'r:','LineWidth',1)

%check if timing is string or num
titleCheck = isstr(s.SoundData.(soundNames{5}).OptoStimDelay);
if titleCheck == 1
    title({'Binned K R (Tone)';strcat(num2str(s.SoundData.(soundNames{5}).TargetFreq/1000),'kHz',s.SoundData.(soundNames{5}).OptoStimDelay,'Timing')})
elseif titleCheck == 0
    title({'Binned K R (Tone)';strcat(num2str(s.SoundData.(soundNames{5}).TargetFreq/1000),'kHz',num2str(s.SoundData.(soundNames{5}).OptoStimDelay),'Timing')})
end
xlim([1 freqLength])
set(gca,'XTick',[1:2:freqLength]);
set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);

%now for general range
subplot(4,3,3)
hold on
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGen,'k','LineWidth',2)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGen+preGenSTD,'k','LineWidth',1)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGen-preGenSTD,'k','LineWidth',1)

plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGen,'r','LineWidth',2)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGen+postGenSTD,'r','LineWidth',1)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGen-postGenSTD,'r','LineWidth',1)

plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGenComp,'k:','LineWidth',2)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGenComp+preGenCompSTD,'k:','LineWidth',1)
plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGenComp-preGenCompSTD,'k:','LineWidth',1)

plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGenComp,'r:','LineWidth',2)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGenComp+postGenCompSTD,'r:','LineWidth',1)
plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGenComp-postGenCompSTD,'r:','LineWidth',1)


title('Binned K R (Gen)')
xlim([1 freqLength])
set(gca,'XTick',[1:2:freqLength]);
set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);

%% Plot information from pairing style presentations

histBinVector = [(params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration) + ...
params.histBin/2:params.histBin:(params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration)-params.histBin/2];

%First, plot general histograms of each period, for control, then
%target
subplot(4,3,5)
hold on
plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{3},'ControlAnalysis')).Histograms,'k','LineWidth',1)
plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{4},'ControlAnalysis')).Histograms,'r','LineWidth',1)
plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{5},'ControlAnalysis')).Histograms,'g','LineWidth',1)
%draws in tone!
plot([0 0],[ylim],'r');
plot([s.SoundData.(soundNames{3}).ToneDuration ...
    s.SoundData.(soundNames{3}).ToneDuration],[ylim],'r');
%sets xlimits to avoid awkward graphs
xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration...
    params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
title(strcat('Control',num2str(s.SoundData.(soundNames{5}).ControlFreq),'kHz',num2str(s.SoundData.(soundNames{5}).ControlDB),' dB K G R'))
%Now for target
subplot(4,3,6)
hold on
plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{3},'TargetAnalysis')).Histograms,'k','LineWidth',1)
plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{4},'TargetAnalysis')).Histograms,'r','LineWidth',1)
plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{5},'TargetAnalysis')).Histograms,'g','LineWidth',1)
%draws in tone!
plot([0 0],[ylim],'r');
plot([s.SoundData.(soundNames{3}).ToneDuration ...
    s.SoundData.(soundNames{3}).ToneDuration],[ylim],'r');
%sets xlimits to avoid awkward graphs
xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration...
    params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
title(strcat('Control',num2str(s.SoundData.(soundNames{5}).TargetFreq),'kHz',num2str(s.SoundData.(soundNames{5}).TargetDB),' dB K G R'))

%Now plot rasters, in chronological order
%first plot control
subplot(4,3,8)
hold on
%need to link up rasters!
rastersPre = s.(desigNames{i}).(strcat(soundNames{3},'ControlAnalysis')).Rasters;
rastersPost = s.(desigNames{i}).(strcat(soundNames{4},'ControlAnalysis')).Rasters;
rastersPairing = s.(desigNames{i}).(strcat(soundNames{5},'ControlAnalysis')).Rasters;
%find number of presentations
rasterNumPre = s.SoundData.(soundNames{3}).ToneRepetitions;
rasterNumPair = s.SoundData.(soundNames{5}).ToneRepetitions;
rasterNumPost = s.SoundData.(soundNames{4}).ToneRepetitions;
%perform compensation for trial number
rastersPairing(:,2) = rastersPairing(:,2) + rasterNumPre;
rastersPost(:,2) = rastersPost(:,2) + (rasterNumPre + rasterNumPair);
%plot
plot(rastersPre(:,1),rastersPre(:,2),'bo','markersize',4)
plot(rastersPre(:,1),rastersPre(:,2),'k.','markersize',4)
plot(rastersPost(:,1),rastersPost(:,2),'ro','markersize',4)
plot(rastersPost(:,1),rastersPost(:,2),'k.','markersize',4)
plot(rastersPairing(:,1),rastersPairing(:,2),'go','markersize',4)
plot(rastersPairing(:,1),rastersPairing(:,2),'k.','markersize',4)

plot([0 0],[ylim],'b');
plot([s.SoundData.(soundNames{3}).ToneDuration s.SoundData.(soundNames{3}).ToneDuration],[ylim],'b');

ylim([0 rasterNumPre+rasterNumPair+rasterNumPost])
xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
title({'Control Rasters';'K G R'})

%now plot Target
subplot(4,3,9)
hold on
%need to link up rasters!
rastersPre = s.(desigNames{i}).(strcat(soundNames{3},'TargetAnalysis')).Rasters;
rastersPost = s.(desigNames{i}).(strcat(soundNames{4},'TargetAnalysis')).Rasters;
rastersPairing = s.(desigNames{i}).(strcat(soundNames{5},'TargetAnalysis')).Rasters;
%find number of presentations
rasterNumPre = s.SoundData.(soundNames{3}).ToneRepetitions;
rasterNumPair = s.SoundData.(soundNames{5}).ToneRepetitions;
rasterNumPost = s.SoundData.(soundNames{4}).ToneRepetitions;
%perform compensation for trial number
rastersPairing(:,2) = rastersPairing(:,2) + rasterNumPre;
rastersPost(:,2) = rastersPost(:,2) + (rasterNumPre + rasterNumPair);
%plot
plot(rastersPre(:,1),rastersPre(:,2),'bo','markersize',4)
plot(rastersPre(:,1),rastersPre(:,2),'k.','markersize',4)
plot(rastersPost(:,1),rastersPost(:,2),'ro','markersize',4)
plot(rastersPost(:,1),rastersPost(:,2),'k.','markersize',4)
plot(rastersPairing(:,1),rastersPairing(:,2),'go','markersize',4)
plot(rastersPairing(:,1),rastersPairing(:,2),'k.','markersize',4)

plot([0 0],[ylim],'b');
plot([s.SoundData.(soundNames{3}).ToneDuration s.SoundData.(soundNames{3}).ToneDuration],[ylim],'b');

ylim([0 rasterNumPre+rasterNumPair+rasterNumPost])
xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
title({'Target Rasters';'K G R'})

%% NEED TO PLOT TIME COURSE
subplot(4,3,11)
hold on
%use the trial information from the rasters for the right designations
%for now, going to use general period. however, can consider switching
%this in the future.
binnedSpikes = zeros(rasterNumPre+rasterNumPair+rasterNumPost,2);
binnedSpikes(:,1) = [1:1:rasterNumPre+rasterNumPair+rasterNumPost];
binnedSpikesComp = zeros(rasterNumPre+rasterNumPair+rasterNumPost,2);
binnedSpikesComp(:,1) = [1:1:rasterNumPre+rasterNumPair+rasterNumPost];
binnedSpikeCounter = 1;
for k = [3,5,4]
    findBin = s.(desigNames{i}).(strcat(soundNames{k},'ControlAnalysis')).LatBinPeakCalcs.BinnedSpikesGen;
    findBinComp = s.(desigNames{i}).(strcat(soundNames{k},'ControlAnalysis')).BinGenComp;
    binnedSpikes(binnedSpikeCounter:binnedSpikeCounter + length(findBin)-1,2) = findBin;
    binnedSpikesComp(binnedSpikeCounter:binnedSpikeCounter + length(findBin)-1,2) = findBinComp;
    binnedSpikeCounter = binnedSpikeCounter + length(findBin);
end
plot(binnedSpikes(:,1),binnedSpikes(:,2),':ko')
plot(binnedSpikes(:,1),smooth(binnedSpikes(:,2),11),'-k.')
%plot in important time points. 
plot([rasterNumPre rasterNumPre],[ylim],'r')
plot([rasterNumPre+rasterNumPair rasterNumPre+rasterNumPair],[ylim],'r')
xlim([1 rasterNumPre+rasterNumPair+rasterNumPost])

title('Control Time Course')

subplot(4,3,12)
hold on
%use the trial information from the rasters for the right designations
%for now, going to use general period. however, can consider switching
%this in the future.
binnedSpikes = zeros(rasterNumPre+rasterNumPair+rasterNumPost,2);
binnedSpikes(:,1) = [1:1:rasterNumPre+rasterNumPair+rasterNumPost];
binnedSpikesComp = zeros(rasterNumPre+rasterNumPair+rasterNumPost,2);
binnedSpikesComp(:,1) = [1:1:rasterNumPre+rasterNumPair+rasterNumPost];
binnedSpikeCounter = 1;
for k = [3,5,4]
    findBin = s.(desigNames{i}).(strcat(soundNames{k},'TargetAnalysis')).LatBinPeakCalcs.BinnedSpikesGen;
    findBinComp = s.(desigNames{i}).(strcat(soundNames{k},'TargetAnalysis')).BinGenComp;
    binnedSpikes(binnedSpikeCounter:binnedSpikeCounter + length(findBin)-1,2) = findBin;
    binnedSpikesComp(binnedSpikeCounter:binnedSpikeCounter + length(findBin)-1,2) = findBinComp;
    binnedSpikeCounter = binnedSpikeCounter + length(findBin);
end
plot(binnedSpikes(:,1),binnedSpikes(:,2),':ko')
plot(binnedSpikes(:,1),smooth(binnedSpikes(:,2),11),'-k.')
%plot in important time points. 
plot([rasterNumPre rasterNumPre],[ylim],'r')
plot([rasterNumPre+rasterNumPair rasterNumPre+rasterNumPair],[ylim],'r')
xlim([1 rasterNumPre+rasterNumPair+rasterNumPost])

title('Target Time Course')

%%
hold off
%save as matlab figure with correct name (fileName+LFP)
spikeGraphName = strcat(fileName,desigNames{i},'PairingAnalysis');

%save as PDF with correct name
set(figName,'Units','Inches');
pos = get(figName,'Position');
set(figName,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(figName,spikeGraphName,'-dpdf','-r0')

end