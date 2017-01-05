%plotting function that plots out information for simple tuning curves.

function [s] = functionOneTonePairingMasterPlot(numUnits,s,...
    desigNames,params,fileName,soundNames);

histBinVector = [(params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur) + ...
    params.histBin/2:params.histBin:(params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur)-params.histBin/2];

for i=1:numUnits
    hFig = figure;
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    set(hFig, 'Position', [5 5 1280 1000])
    %% plots average waveform
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    title({fileName;desigNames{i};strcat('AverageFiringRate:',num2str(s.(desigNames{i}).BaselineFiringRate))});
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
    
    %draws in tone!
    plot([0 0],[ylim],'r');
    plot([s.SoundData.(soundNames{1}).ToneDur ...
        s.SoundData.(soundNames{1}).ToneDur],[ylim],'r');
    %sets xlimits to avoid awkward graphs
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur...
        params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    title('Histogram K before B after')
    %% plot rasters organized by time
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
    set(gca,'YTick',rasterFreqLines1(2:end,1)-s.SoundData.(soundNames{1}).ToneReps/2*size(s.SoundData.(soundNames{1}).UniqueDBs,1));
    set(gca,'YTickLabel',rasterFreqLines1(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 length(s.SoundData.(soundNames{1}).Frequencies)])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    title({'Pre-Pairing Rasters';'dB Increase in Descending Order'})
    %% plot rasters organized by frequency
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
    set(gca,'YTick',rasterFreqLines2(2:end,1)-s.SoundData.(soundNames{2}).ToneReps/2*size(s.SoundData.(soundNames{2}).UniqueDBs,1));
    set(gca,'YTickLabel',rasterFreqLines2(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 length(s.SoundData.(soundNames{2}).Frequencies)])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{2}).ToneDur params.rasterWindow(2)*s.SoundData.(soundNames{2}).ToneDur])
    title({'Post-Pairing Rasters';'dB Increase in Descending Order'})
    
    %% plots heatmap by frequencies and time. Plots the first tuning.
    subplot(4,3,5)
    dataPrep1 = squeeze(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).FreqDBHistograms);
    imagesc(dataPrep1)
    colorbar
    cMinMax = [0 0];
    %sets limits so that next graph is displayed with the same color
    %settings. 
    cMinMax(1) = min(min(dataPrep1));
    cMinMax(2) = max(max(dataPrep1));
    set(gca,'YTick',s.SoundData.(soundNames{1}).OctaveRange(:,2));
    set(gca,'YTickLabel',s.SoundData.(soundNames{1}).OctaveRange(:,1));
    set(gca,'XTick',[1:20:size(histBinVector,2)]);
    set(gca,'XTickLabel',histBinVector(1:20:end));
    histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
    histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),s.SoundData.(soundNames{1}).ToneDur);
    line([histBinZero histBinZero],[-1 size(s.SoundData.(soundNames{1}).UniqueFreqs,1)],'LineWidth',3,'Color','green')
    line([histBinZero histBinZero],[-1 size(s.SoundData.(soundNames{1}).UniqueFreqs,1)],'LineWidth',2,'Color','black')
    line([histBinTone histBinTone],[-1 size(s.SoundData.(soundNames{1}).UniqueFreqs,1)],'LineWidth',3,'Color','green')
    line([histBinTone histBinTone],[-1 size(s.SoundData.(soundNames{1}).UniqueFreqs,1)],'LineWidth',2,'Color','black')
    
     %this will plot a colored box around the target frequency
    %target!
    if ~isempty(find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5)
        rectangle('Position',[0.5 find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5 size(dataPrep1,2) 1],'EdgeColor','g','LineWidth',2)
    else
        findClosest = round(s.SoundData.(soundNames{2}).UniqueFreqs) - s.SoundData.(soundNames{3}).Frequency;
        smallestDiff = min(abs(findClosest));
        findClosest = find(abs(findClosest) == smallestDiff);
        rectangle('Position',[0.5 findClosest-0.5 size(dataPrep1,2) 1],'EdgeColor','g','LineWidth',2)
    end
    title('Pre-Pairing Heatmap by Frequency and Time')
    
    %% plots heatmap by frequencies and time. Plots the second tuning.
    subplot(4,3,8)
    dataPrep2 = squeeze(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).FreqDBHistograms);
    if cMinMax(1) == cMinMax(2)
        imagesc(dataPrep2)
    else
        imagesc(dataPrep2, cMinMax)
    end
    colorbar
    set(gca,'YTick',s.SoundData.(soundNames{2}).OctaveRange(:,2));
    set(gca,'YTickLabel',s.SoundData.(soundNames{2}).OctaveRange(:,1));
    set(gca,'XTick',[1:20:size(histBinVector,2)]);
    set(gca,'XTickLabel',histBinVector(1:20:end));
    histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
    histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),s.SoundData.(soundNames{1}).ToneDur);
    line([histBinZero histBinZero],[-1 size(s.SoundData.(soundNames{2}).UniqueFreqs,1)],'LineWidth',3,'Color','green')
    line([histBinZero histBinZero],[-1 size(s.SoundData.(soundNames{2}).UniqueFreqs,1)],'LineWidth',2,'Color','black')
    line([histBinTone histBinTone],[-1 size(s.SoundData.(soundNames{2}).UniqueFreqs,1)],'LineWidth',3,'Color','green')
    line([histBinTone histBinTone],[-1 size(s.SoundData.(soundNames{2}).UniqueFreqs,1)],'LineWidth',2,'Color','black')
    
     %this will plot a colored box around the target frequency
    %target!
    if ~isempty(find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5)
        rectangle('Position',[0.5 find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5 size(dataPrep1,2) 1],'EdgeColor','g','LineWidth',2)
    else
        findClosest = round(s.SoundData.(soundNames{2}).UniqueFreqs) - s.SoundData.(soundNames{3}).Frequency;
        smallestDiff = min(abs(findClosest));
        findClosest = find(abs(findClosest) == smallestDiff);
        rectangle('Position',[0.5 findClosest-0.5 size(dataPrep1,2) 1],'EdgeColor','g','LineWidth',2)
    end
    title('Post-Pairing Heatmap by Frequency and Time') 
    
    %% plots difference in heatmaps by frequency. Raw subtraction
    subplot(4,3,11)
    imagesc(dataPrep2-dataPrep1)
    colorbar
    set(gca,'YTick',s.SoundData.(soundNames{2}).OctaveRange(:,2));
    set(gca,'YTickLabel',s.SoundData.(soundNames{2}).OctaveRange(:,1));
    set(gca,'XTick',[1:20:size(histBinVector,2)]);
    set(gca,'XTickLabel',histBinVector(1:20:end));
    histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
    histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),s.SoundData.(soundNames{1}).ToneDur);
    line([histBinZero histBinZero],[-1 size(s.SoundData.(soundNames{2}).UniqueFreqs,1)],'LineWidth',3,'Color','green')
    line([histBinZero histBinZero],[-1 size(s.SoundData.(soundNames{2}).UniqueFreqs,1)],'LineWidth',2,'Color','black')
    line([histBinTone histBinTone],[-1 size(s.SoundData.(soundNames{2}).UniqueFreqs,1)],'LineWidth',3,'Color','green')
    line([histBinTone histBinTone],[-1 size(s.SoundData.(soundNames{2}).UniqueFreqs,1)],'LineWidth',2,'Color','black')
    
     %this will plot a colored box around the target frequency
    %target!
    if ~isempty(find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5)
        rectangle('Position',[0.5 find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5 size(dataPrep1,2) 1],'EdgeColor','g','LineWidth',2)
    else
        findClosest = round(s.SoundData.(soundNames{2}).UniqueFreqs) - s.SoundData.(soundNames{3}).Frequency;
        smallestDiff = min(abs(findClosest));
        findClosest = find(abs(findClosest) == smallestDiff);
        rectangle('Position',[0.5 findClosest-0.5 size(dataPrep1,2) 1],'EdgeColor','g','LineWidth',2)
    end
    title('Difference Heatmap by Frequency and Time') 
    
    
    
    %% plot rasters from pairing
    subplot(4,3,3)
    plot(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).AllRasters(:,2),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{3}).ToneDuration s.SoundData.(soundNames{3}).ToneDuration],[ylim],'b');
    ylim([0 s.SoundData.(soundNames{3}).ToneRepetitions])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
    title('Pairing Rasters')
    
    %% plots binned responses relative to pairing presentation
    subplot(4,3,6)
    plot(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).BinSpikeTimes(:,s.Parameters.chosenSpikeBin))
    xlim([1 size(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).BinSpikeTimes,1)])
    try
        ylim([min(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).BinSpikeTimes)...
        max(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).BinSpikeTimes)])
    catch
    end
    ylabel('Binned Spikes')
    xlabel('Pairing Trial Number')
    title('Binned Responses During Pairing')
    
    %% Plot tuning curve! Non background subtracted
    subplot(4,3,9)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinSpikeStats(:,1,1,s.Parameters.chosenSpikeBin),'k');
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinSpikeStats(:,1,1,s.Parameters.chosenSpikeBin),'r');
    set(gca,'XTick',[1:1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs);
    title('RAW Tuning Curves K before R after')
    
    %% Plot tuning curve , with background subtracted
    subplot(4,3,12)
    plot(s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinSpikeStats(:,1,1,s.Parameters.chosenSpikeBin)-...
        s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinSpikeStats(:,1,1,1),'k');
    hold on
    plot(s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinSpikeStats(:,1,1,s.Parameters.chosenSpikeBin)-...
        s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinSpikeStats(:,1,1,1),'r');
    set(gca,'XTick',[1:1:size(s.SoundData.(soundNames{1}).UniqueFreqs,1)]);
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).UniqueFreqs);
    title('Back Subtract Tuning Curves K before R after')
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