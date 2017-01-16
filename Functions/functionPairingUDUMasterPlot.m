%plotting function that plots out information for simple tuning curves.

function [s] = functionPairingUDUMasterPlot(numUnits,s,...
    desigNames,params,fileName,soundNames);

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

histBinVector = [(params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur) + ...
    params.histBin/2:params.histBin:(params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur)-params.histBin/2];

for i=1:numUnits
    %% FIRST FIGURE
    hFig = figure;
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    set(hFig, 'Position', [10 80 1280 1000])
    %plots average waveform
    subplot(4,8,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    title({fileName;desigNames{i};strcat('AverageFiringRate:',num2str(s.(desigNames{i}).BaselineFiringRate))});
    set(0, 'DefaulttextInterpreter', 'none')
    %plots ISI
    subplot(4,8,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([params.rpvTime params.rpvTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(params.clusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))});
    
    %% set 1 hist plots histogram from set 1
    subplot(4,4,2)
    %plots histogram from first tuning.
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistograms,'k','LineWidth',1)
    hold on
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistograms,'r','LineWidth',1)
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistograms,'g','LineWidth',1)
    %draws in tone!
    plot([0 0],[ylim],'r');
    plot([s.SoundData.(soundNames{1}).ToneDur ...
        s.SoundData.(soundNames{1}).ToneDur],[ylim],'r');
    %sets xlimits to avoid awkward graphs
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur...
        params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    %plot significant positive values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'ko')
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'k*')
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
        'ro')
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'r*')
    %plot negative values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'm*')
    %plot significant positive values for f/u first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'go')
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'g*')
    %plot negative values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'c*')
    title('Set 1 K before R After G NextPre')
    
    %% set 2 hist %plots histogram from set 2
    subplot(4,4,3)
    %plots histogram from first tuning.
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistograms,'k','LineWidth',1)
    hold on
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistograms,'r','LineWidth',1)
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistograms,'g','LineWidth',1)
    %draws in tone!
    plot([0 0],[ylim],'r');
    plot([s.SoundData.(soundNames{1}).ToneDur ...
        s.SoundData.(soundNames{1}).ToneDur],[ylim],'r');
    %sets xlimits to avoid awkward graphs
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur...
        params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    %plot significant positive values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'ko')
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'k*')
    %plot negative values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'b*')
    %plot significant positive values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'ro')
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'r*')
    %plot negative values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'm*')
    %plot significant positive values for f/u first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'go')
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'g*')
    %plot negative values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'c*')
    title('Set 2 K before B After G NextPre')
    
    %% set 3 hist %plots histogram from set 3
    subplot(4,4,4)
    %plots histogram from first tuning.
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistograms,'k','LineWidth',1)
    hold on
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistograms,'r','LineWidth',1)
    %draws in tone!
    plot([0 0],[ylim],'r');
    plot([s.SoundData.(soundNames{1}).ToneDur ...
        s.SoundData.(soundNames{1}).ToneDur],[ylim],'r');
    %sets xlimits to avoid awkward graphs
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur...
        params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    %plot significant positive values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'ko')
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'k*')
    %plot negative values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'k*')
    %plot significant positive values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'ro')
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'r*')
    %plot negative values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'm*')
    title('Set 3 K before B After')
    
    %% Set 1 Rasters plot rasters organized by frequency and dB
    subplot(4,3,4)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{1}).ToneDur s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
    rasterFreqLines = zeros(size(s.SoundData.(soundNames{1}).UniqueFreqs,1),2);
    rasterFreqLines(:,1) = s.SoundData.(soundNames{1}).ToneReps*...
        size(s.SoundData.(soundNames{1}).UniqueDBs,1):s.SoundData.(soundNames{1}).ToneReps*...
        size(s.SoundData.(soundNames{1}).UniqueDBs,1):length(s.SoundData.(soundNames{1}).Frequencies);
    rasterFreqLines(:,2) = s.SoundData.(soundNames{1}).UniqueFreqs/1000;
    for k = 1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{1}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1)-(s.SoundData.(soundNames{1}).ToneReps/2));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 length(s.SoundData.(soundNames{1}).Frequencies)])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    title({'Set 1 Pre-Pairing Rasters';'dB Increase in Descending Order'})
    
    subplot(4,3,7)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{2}).ToneDur s.SoundData.(soundNames{2}).ToneDur],[ylim],'b');
    for k = 1:size(s.SoundData.(soundNames{2}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{2}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1)-(s.SoundData.(soundNames{1}).ToneReps/2));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 length(s.SoundData.(soundNames{2}).Frequencies)])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{2}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{2}).ToneDur])
    title({'Set 1 Post-Pairing Rasters';'dB Increase in Descending Order'})
    
    %plot rasters from pairing
    subplot(4,3,10)
    plot(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).AllRasters(:,2),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{3}).ToneDuration s.SoundData.(soundNames{3}).ToneDuration],[ylim],'b');
    ylim([0 s.SoundData.(soundNames{3}).ToneRepetitions])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
    title('Set 1 Pairing Rasters')
    
    %% Set 2 Rasters plot rasters organized by frequency and dB
    subplot(4,3,5)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{1}).ToneDur s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
    for k = 1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{1}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1)-(s.SoundData.(soundNames{1}).ToneReps/2));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 length(s.SoundData.(soundNames{1}).Frequencies)])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    title({'Set 2 Pre-Pairing Rasters';'dB Increase in Descending Order'})
    
    subplot(4,3,8)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{2}).ToneDur s.SoundData.(soundNames{2}).ToneDur],[ylim],'b');
    for k = 1:size(s.SoundData.(soundNames{2}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{2}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1)-(s.SoundData.(soundNames{1}).ToneReps/2));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 length(s.SoundData.(soundNames{2}).Frequencies)])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{2}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{2}).ToneDur])
    title({'Set 2 Post-Pairing Rasters';'dB Increase in Descending Order'})
    
    %plot rasters from pairing
    subplot(4,3,11)
    plot(s.(desigNames{i}).(strcat(soundNames{6},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{6},'Analysis')).AllRasters(:,2),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{3}).ToneDuration s.SoundData.(soundNames{3}).ToneDuration],[ylim],'b');
    ylim([0 s.SoundData.(soundNames{3}).ToneRepetitions])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
    title('Set 2 Pairing Rasters')
    
    %% Set 3 Rasters plot rasters organized by frequency and dB
    subplot(4,3,6)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{1}).ToneDur s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
    for k = 1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{1}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1)-(s.SoundData.(soundNames{1}).ToneReps/2));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 length(s.SoundData.(soundNames{1}).Frequencies)])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    title({'Set 3 Pre-Pairing Rasters';'dB Increase in Descending Order'})
    
    subplot(4,3,9)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{2}).ToneDur s.SoundData.(soundNames{2}).ToneDur],[ylim],'b');
    for k = 1:size(s.SoundData.(soundNames{2}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{2}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1)-(s.SoundData.(soundNames{1}).ToneReps/2));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 length(s.SoundData.(soundNames{2}).Frequencies)])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{2}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{2}).ToneDur])
    title({'Set 3 Post-Pairing Rasters';'dB Increase in Descending Order'})
    
    %plot rasters from pairing
    subplot(4,3,12)
    plot(s.(desigNames{i}).(strcat(soundNames{9},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{9},'Analysis')).AllRasters(:,2),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{3}).ToneDuration s.SoundData.(soundNames{3}).ToneDuration],[ylim],'b');
    ylim([0 s.SoundData.(soundNames{3}).ToneRepetitions])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
    title('Set 3 Pairing Rasters')
    
    %%
    hold off
    %save as matlab figure with correct name (fileName+LFP)
    spikeGraphName = strcat(fileName,desigNames{i},'PairingAnalysis');

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
    %% Second figure
    hFig = figure;
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis2');
    set(hFig, 'Position', [10 80 1280 1000])
    
    %% First Pairing
    %170116 Rank sum code has some weird glitches. Instead plot standard
    %error for visualization purposes
    %first pull all the data
    freqLength = length(s.SoundData.(soundNames{1}).UniqueFreqs);
    dbLength = length(s.SoundData.(soundNames{1}).UniqueDBs);
    toneReps = s.SoundData.(soundNames{1}).ToneReps;
    
    preDataTone = zeros(freqLength,dbLength,toneReps);
    postDataTone = zeros(freqLength,dbLength,toneReps);
    postPostDataTone = zeros(freqLength,dbLength,toneReps);
    
    
    preDataGen = zeros(freqLength,dbLength,toneReps);
    postDataGen = zeros(freqLength,dbLength,toneReps);
    postPostDataGen = zeros(freqLength,dbLength,toneReps);

    
    for freqInd = 1:freqLength
        for dbInd = 1:dbLength
            preDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
            postDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
            postPostDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
            
            preDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
            postDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
            postPostDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
            
        end
    end
    
    preToneSTD = std(preDataTone,1,3)/sqrt(toneReps); 
    postToneSTD = std(postDataTone,1,3)/sqrt(toneReps); 
    postPostToneSTD = std(postPostDataTone,1,3)/sqrt(toneReps); 
    
    preGenSTD = std(preDataGen,1,3)/sqrt(toneReps); 
    postGenSTD = std(postDataGen,1,3)/sqrt(toneReps);
    postPostGenSTD = std(postPostDataGen,1,3)/sqrt(toneReps);
    
    %do the same for compensated binned values
    preDataToneComp = zeros(freqLength,dbLength,toneReps);
    postDataToneComp = zeros(freqLength,dbLength,toneReps);
    postPostDataToneComp = zeros(freqLength,dbLength,toneReps);
    
   
    
    preDataGenComp = zeros(freqLength,dbLength,toneReps);
    postDataGenComp = zeros(freqLength,dbLength,toneReps);
    postPostDataGenComp = zeros(freqLength,dbLength,toneReps);
    
    
    
    for freqInd = 1:freqLength
        for dbInd = 1:dbLength
            preDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
            postDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
            postPostDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
            
            preDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
            postDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
            postPostDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
            
        end
    end
    
    preToneCompSTD = std(preDataToneComp,1,3)/sqrt(toneReps); 
    postToneCompSTD = std(postDataToneComp,1,3)/sqrt(toneReps); 
    postPostToneCompSTD = std(postPostDataToneComp,1,3)/sqrt(toneReps); 
    
    preGenCompSTD = std(preDataGenComp,1,3)/sqrt(toneReps); 
    postGenCompSTD = std(postDataGenComp,1,3)/sqrt(toneReps);
    postPostGenCompSTD = std(postPostDataGenComp,1,3)/sqrt(toneReps);
    
    %first plot binned responses, first to tone
    %first pairing
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreTone,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreTone+preToneSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreTone-preToneSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreTone,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreTone+postToneSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreTone-postToneSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreTone,'g','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreTone+postPostToneSTD,'g','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreTone-postPostToneSTD,'g','LineWidth',1)

    %check if timing is string or num
    titleCheck = isstr(s.SoundData.(soundNames{3}).OptoStimDelay);
    if titleCheck == 1
        title({'Binned K R G (Tone)';strcat(num2str(s.SoundData.(soundNames{3}).Frequency/1000),'kHz',s.SoundData.(soundNames{3}).OptoStimDelay,'Timing')})
    elseif titleCheck == 0
        title({'Binned K R G (Tone)';strcat(num2str(s.SoundData.(soundNames{3}).Frequency/1000),'kHz',num2str(s.SoundData.(soundNames{3}).OptoStimDelay),'Timing')})
    end
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %first pairing, with compensation for background
    subplot(4,6,2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreToneComp,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreToneComp+preToneCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreToneComp-preToneCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreToneComp,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreToneComp+postToneCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreToneComp-postToneCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreToneComp,'g','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreToneComp+postPostToneCompSTD,'g','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreToneComp-postPostToneCompSTD,'g','LineWidth',1)

    
    xlim([1 freqLength])
    title('Binned K R G (Tone) BACK')
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %now for general range
    subplot(4,6,7)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGen,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGen+preGenSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGen-preGenSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGen,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGen+postGenSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGen-postGenSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGen,'g','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGen+postPostGenSTD,'g','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGen-postPostGenSTD,'g','LineWidth',1)
    
    title('Binned K R G (Gen)')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    %first pairing, with compensation for background
    subplot(4,6,8)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGenComp,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGenComp+preGenCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinStoreGenComp-preGenCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGenComp,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGenComp+postGenCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinStoreGenComp-postGenCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGenComp,'g','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGenComp+postPostGenCompSTD,'g','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGenComp-postPostGenCompSTD,'g','LineWidth',1)
    
    
    title('Binned K R G (Gen) BACK')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %now plot peak values, all for general
    subplot(4,6,13)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).PeakStoreGen,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).PeakStoreGen,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).PeakStoreGen,'g','LineWidth',2)
    title('Peak K R G (Tone)')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    %plot compensated version
    subplot(4,6,14)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).PeakStoreGenComp,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).PeakStoreGenComp,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).PeakStoreGenComp,'g','LineWidth',2)
    title('Peak K R G (Tone) BACK')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %finally, plot latency
    subplot(4,3,10)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).LatencyStore,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).LatencyStore,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).LatencyStore,'g','LineWidth',2)
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    title('Latency K R G (Gen)')
    xlim([1 freqLength])
    
    
    %% Second Pairing
    %lets calculate significance using a rank sum!
    %first pull all the data
    freqLength = length(s.SoundData.(soundNames{1}).UniqueFreqs);
    dbLength = length(s.SoundData.(soundNames{1}).UniqueDBs);
    toneReps = s.SoundData.(soundNames{1}).ToneReps;
    
    preDataTone = zeros(freqLength,dbLength,toneReps);
    postDataTone = zeros(freqLength,dbLength,toneReps);
    postPostDataTone = zeros(freqLength,dbLength,toneReps);

    
    preDataGen = zeros(freqLength,dbLength,toneReps);
    postDataGen = zeros(freqLength,dbLength,toneReps);
    postPostDataGen = zeros(freqLength,dbLength,toneReps);
    
    for freqInd = 1:freqLength
        for dbInd = 1:dbLength
            preDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
            postDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
            postPostDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
            
            preDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
            postDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
            postPostDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
            
        end
    end
    
    preToneSTD = std(preDataTone,1,3)/sqrt(toneReps); 
    postToneSTD = std(postDataTone,1,3)/sqrt(toneReps); 
    postPostToneSTD = std(postPostDataTone,1,3)/sqrt(toneReps); 
    preGenSTD = std(preDataGen,1,3)/sqrt(toneReps); 
    postPostGenSTD = std(postPostDataGen,1,3)/sqrt(toneReps);
    
    %do the same for compensated binned values
    preDataToneComp = zeros(freqLength,dbLength,toneReps);
    postDataToneComp = zeros(freqLength,dbLength,toneReps);
    postPostDataToneComp = zeros(freqLength,dbLength,toneReps);
    

    
    preDataGenComp = zeros(freqLength,dbLength,toneReps);
    postDataGenComp = zeros(freqLength,dbLength,toneReps);
    postPostDataGenComp = zeros(freqLength,dbLength,toneReps);
    

    
    for freqInd = 1:freqLength
        for dbInd = 1:dbLength
            preDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
            postDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
            postPostDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
            
            preDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
            postDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
            postPostDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
           
        end
    end
    preToneCompSTD = std(preDataToneComp,1,3)/sqrt(toneReps); 
    postToneCompSTD = std(postDataToneComp,1,3)/sqrt(toneReps); 
    postPostToneCompSTD = std(postPostDataToneComp,1,3)/sqrt(toneReps); 
    preGenCompSTD = std(preDataGenComp,1,3)/sqrt(toneReps); 
    postGenCompSTD = std(postDataGenComp,1,3)/sqrt(toneReps);
    postPostGenCompSTD = std(postPostDataGenComp,1,3)/sqrt(toneReps);
    
    %first plot binned responses, first to tone
    %first pairing
    subplot(4,6,3)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreTone,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreTone + preToneSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreTone - preToneSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreTone,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreTone + postToneSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreTone - postToneSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreTone,'g','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreTone + postPostToneSTD,'g','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreTone - postPostToneSTD,'g','LineWidth',1)
    %plot significant values as asterisks?
    sigVal1 = find(squeeze(rankSumTone(:,:,1))<0.05);
    sigVal2 = find(squeeze(rankSumTone(:,:,2))<0.05);
    plot(sigVal1,s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreTone(sigVal1),'k*')
    plot(sigVal2,s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreTone(sigVal2),'r*')
    %check if timing is string or num
    titleCheck = isstr(s.SoundData.(soundNames{6}).OptoStimDelay);
    if titleCheck == 1
        title({'Binned K R G (Tone)';strcat(num2str(s.SoundData.(soundNames{6}).Frequency/1000),'kHz',s.SoundData.(soundNames{6}).OptoStimDelay,'Timing')})
    elseif titleCheck == 0
        title({'Binned K R G (Tone)';strcat(num2str(s.SoundData.(soundNames{6}).Frequency/1000),'kHz',num2str(s.SoundData.(soundNames{6}).OptoStimDelay),'Timing')})
    end
    
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %first pairing, with compensation for background
    subplot(4,6,4)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreToneComp,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreToneComp + preToneCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreToneComp - preToneCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreToneComp,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreToneComp + postToneCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreToneComp - postToneCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreToneComp,'g','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreToneComp + postPostToneCompSTD,'g','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreToneComp - postPostToneCompSTD,'g','LineWidth',1)
    
    
    xlim([1 freqLength])
    title('Binned K R G (Tone) BACK')
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %now for general range
    subplot(4,6,9)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGen,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGen + preGenSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGen - preGenSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreGen,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreGen + postGenSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreGen - postGenSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGen,'g','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGen + postPostGenSTD,'g','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGen - postPostGenSTD,'g','LineWidth',1)
    

    title('Binned K R G (Gen)')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    %first pairing, with compensation for background, general window
    subplot(4,6,10)
    hold on
    
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGenComp,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGenComp + preGenCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).BinStoreGenComp - preGenCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreGenComp,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreGenComp + postGenCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).BinStoreGenComp - postGenCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGenComp,'g','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGenComp + postPostGenCompSTD,'g','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGenComp - postPostGenCompSTD,'g','LineWidth',1)
    
    title('Binned K R G (Gen) BACK')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %now plot peak values, all for general
    subplot(4,6,15)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).PeakStoreGen,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).PeakStoreGen,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).PeakStoreGen,'g','LineWidth',2)
    title('Peak K R G (Tone)')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    %plot compensated version
    subplot(4,6,16)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).PeakStoreGenComp,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).PeakStoreGenComp,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).PeakStoreGenComp,'g','LineWidth',2)
    title('Peak K R G (Tone) BACK')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %finally, plot latency
    subplot(4,3,11)
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).LatencyStore,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).LatencyStore,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).LatencyStore,'g','LineWidth',2)
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    title('Latency K R G (Gen)')
    xlim([1 freqLength])
    
    %% Third Pairing
    %lets calculate significance using a rank sum!
    %first pull all the data
    freqLength = length(s.SoundData.(soundNames{1}).UniqueFreqs);
    dbLength = length(s.SoundData.(soundNames{1}).UniqueDBs);
    toneReps = s.SoundData.(soundNames{1}).ToneReps;
    
    preDataTone = zeros(freqLength,dbLength,toneReps);
    postDataTone = zeros(freqLength,dbLength,toneReps);

    
    preDataGen = zeros(freqLength,dbLength,toneReps);
    postDataGen = zeros(freqLength,dbLength,toneReps);

    
    for freqInd = 1:freqLength
        for dbInd = 1:dbLength
            preDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;
            postDataTone(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesTone;

            
            preDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;
            postDataGen(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGen;

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
            preDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;
            postDataToneComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesToneComp;

            
            preDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;
            postDataGenComp(freqInd,dbInd,:) = s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).IndividualLatbinPeakCalcs{freqInd,dbInd}.BinnedSpikesGenComp;

        end
    end
    
    preToneCompSTD = std(preDataToneComp,1,3)/sqrt(toneReps); 
    postToneCompSTD = std(postDataToneComp,1,3)/sqrt(toneReps); 
    preGenCompSTD = std(preDataGenComp,1,3)/sqrt(toneReps); 
    postGenCompSTD = std(postDataGenComp,1,3)/sqrt(toneReps);
    
    %first plot binned responses, first to tone
    %first pairing
    subplot(4,6,5)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreTone,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreTone+preToneSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreTone-preToneSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreTone,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreTone+postToneSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreTone-postToneSTD,'r','LineWidth',1)
    
    
    
    %check if timing is string or num
    titleCheck = isstr(s.SoundData.(soundNames{9}).OptoStimDelay);
    if titleCheck == 1
        title({'Binned K R G (Tone)';strcat(num2str(s.SoundData.(soundNames{9}).Frequency/1000),'kHz',s.SoundData.(soundNames{9}).OptoStimDelay,'Timing')})
    elseif titleCheck == 0
        title({'Binned K R G (Tone)';strcat(num2str(s.SoundData.(soundNames{9}).Frequency/1000),'kHz',num2str(s.SoundData.(soundNames{9}).OptoStimDelay),'Timing')})
    end
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %first pairing, with compensation for background
    subplot(4,6,6)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreToneComp,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreToneComp+preToneCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreToneComp-preToneCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreToneComp,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreToneComp+postToneCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreToneComp-postToneCompSTD,'r','LineWidth',1)

    xlim([1 freqLength])
    title('Binned K R (Tone) BACK')
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %now for general range
    subplot(4,6,11)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGen,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGen+preGenSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGen-preGenSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreGen,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreGen+postGenSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreGen-postGenSTD,'r','LineWidth',1)

    title('Binned K R (Gen)')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %first pairing, with compensation for background
    subplot(4,6,12)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGenComp,'k','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGenComp+preGenCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).BinStoreGenComp-preGenCompSTD,'k','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreGenComp,'r','LineWidth',2)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreGenComp+postGenCompSTD,'r','LineWidth',1)
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).BinStoreGenComp-postGenCompSTD,'r','LineWidth',1)
    title('Binned K R (Gen) BACK')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %now plot peak values, all for general
    subplot(4,6,17)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).PeakStoreGen,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).PeakStoreGen,'r','LineWidth',2)

    title('Peak K R (Tone)')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    %plot compensated version
    subplot(4,6,18)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).PeakStoreGenComp,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).PeakStoreGenComp,'r','LineWidth',2)

    title('Peak K R (Tone) BACK')
    xlim([1 freqLength])
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    
    %finally, plot latency
    subplot(4,3,12)
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).LatencyStore,'k','LineWidth',2)
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).LatencyStore,'r','LineWidth',2)
    set(gca,'XTick',[1:2:freqLength]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs(1:2:end)/1000);
    title('Latency K R G (Gen)')
    xlim([1 freqLength])
    
    hold off
    %save as matlab figure with correct name (fileName+LFP)

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
    
    
end
end