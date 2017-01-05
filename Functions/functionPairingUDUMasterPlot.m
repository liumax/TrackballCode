%plotting function that plots out information for simple tuning curves.

function [s] = functionPairingUDUMasterPlot(numUnits,s,...
    desigNames,params,fileName,soundNames);

histBinVector = [(params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur) + ...
    params.histBin/2:params.histBin:(params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur)-params.histBin/2];

for i=1:numUnits
    hFig = figure;
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    set(hFig, 'Position', [5 5 1280 1000])
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
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistograms,'b','LineWidth',1)
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
        'b*')
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'c*')
    %plot negative values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'k*')
    %plot significant positive values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'm*')
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'g*')
    %plot negative values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'r*')
    title('Set 1 K before B After')
    
    %% set 2 hist %plots histogram from set 2
    subplot(4,4,3)
    %plots histogram from first tuning.
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistograms,'k','LineWidth',1)
    hold on
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistograms,'b','LineWidth',1)
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
        'b*')
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'c*')
    %plot negative values for first tuning
    plot(s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{4},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'k*')
    %plot significant positive values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,3) == 1,1),...
        'm*')
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'g*')
    %plot negative values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{5},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'r*')
    title('Set 2 K before B After')
    
    %% set 3 hist %plots histogram from set 3
    subplot(4,4,4)
    %plots histogram from first tuning.
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistograms,'k','LineWidth',1)
    hold on
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistograms,'b','LineWidth',1)
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
        'b*')
    plot(s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{7},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'c*')
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
        'm*')
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,4) == 1,1),...
        'g*')
    %plot negative values for second tuning
    plot(s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Centers(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(...
        s.(desigNames{i}).(strcat(soundNames{8},'Analysis')).AllHistogramSig.Histogram(:,6) == 1,1),...
        'r*')
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
    rasterFreqLines(:,2) = s.SoundData.(soundNames{1}).UniqueFreqs;
    for k = 1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{1}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1));
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
    set(gca,'YTick',rasterFreqLines(:,1));
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
    rasterFreqLines = zeros(size(s.SoundData.(soundNames{1}).UniqueFreqs,1),2);
    rasterFreqLines(:,1) = s.SoundData.(soundNames{1}).ToneReps*...
        size(s.SoundData.(soundNames{1}).UniqueDBs,1):s.SoundData.(soundNames{1}).ToneReps*...
        size(s.SoundData.(soundNames{1}).UniqueDBs,1):length(s.SoundData.(soundNames{1}).Frequencies);
    rasterFreqLines(:,2) = s.SoundData.(soundNames{1}).UniqueFreqs;
    for k = 1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{1}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1));
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
    set(gca,'YTick',rasterFreqLines(:,1));
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
    rasterFreqLines = zeros(size(s.SoundData.(soundNames{1}).UniqueFreqs,1),2);
    rasterFreqLines(:,1) = s.SoundData.(soundNames{1}).ToneReps*...
        size(s.SoundData.(soundNames{1}).UniqueDBs,1):s.SoundData.(soundNames{1}).ToneReps*...
        size(s.SoundData.(soundNames{1}).UniqueDBs,1):length(s.SoundData.(soundNames{1}).Frequencies);
    rasterFreqLines(:,2) = s.SoundData.(soundNames{1}).UniqueFreqs;
    for k = 1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)
        plot(params.rasterWindow*s.SoundData.(soundNames{1}).ToneDur,[rasterFreqLines(k,1) rasterFreqLines(k,1)],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1));
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
    set(gca,'YTick',rasterFreqLines(:,1));
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
end
end