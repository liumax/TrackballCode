%plotting function that plots out information for simple tuning curves.

function [s] = functionNewPairingMasterPlot(numUnits,s,...
    desigNames,params,fileName,soundNames);

histBinVector = [(params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur) + ...
    params.histBin/2:params.histBin:(params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur)-params.histBin/2];

for i=1:numUnits
    hFig = figure;
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    set(hFig, 'Position', [5 5 1280 1000])
    %plots average waveform
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms(:,2),'LineWidth',2)
    plot(s.(desigNames{i}).AverageWaveForms(:,1),'r','LineWidth',1)
    plot(s.(desigNames{i}).AverageWaveForms(:,3),'r','LineWidth',1)
    title({spikeGraphName;strcat('AverageFiringRate:',num2str(s.(desigNames{i}).BaselineFiringRate))});
    set(0, 'DefaulttextInterpreter', 'none')
    %plots ISI
    subplot(4,6,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([params.rpvTime params.rpvTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(params.clusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))});
    
    %plots histogram
    subplot(4,3,2)
    %plots histogram from first tuning.
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).AllHistograms,'k','LineWidth',2)
    hold on
    plot(histBinVector,s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).AllHistograms,'b','LineWidth',2)
    %draws in tone!
    plot([0 0],[ylim],'r');
    plot([s.SoundData.(soundNames{1}).ToneDur ...
        s.SoundData.(soundNames{1}).ToneDur],[ylim],'r');
    %sets xlimits to avoid awkward graphs
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur...
        params.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur])
    title('Histogram K before B after')
    
    %plots heatmap
    subplot(4,3,4)
    %plots heatmap for first tuning curve.
    dataPrep1 = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).BinSpikeStats(:,:,params.chosenSpikeBin);
    imagesc(dataPrep1')
    colormap hot
    colorbar
    cMinMax = [0 0];
    %sets limits so that next graph is displayed with the same color
    %settings. 
    cMinMax(1) = min(min(dataPrep1));
    cMinMax(2) = max(max(dataPrep1));
    set(gca,'XTick',s.SoundData.(soundNames{1}).OctaveRange(:,2));
    set(gca,'XTickLabel',s.SoundData.(soundNames{1}).OctaveRange(:,1));
    set(gca,'YTick',s.SoundData.(soundNames{1}).dBRange(:,2));
    set(gca,'YTickLabel',s.SoundData.(soundNames{1}).dBRange(:,1));
    title('Pre-Pairing Binned Response to Tone');
    %plots rectangle on target 
    if ~isempty(find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5)
        rectangle('Position',[find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5...
            find(round(s.SoundData.(soundNames{2}).UniqueDBs) == s.SoundData.(soundNames{3}).dB)-0.5 1 1],'EdgeColor','g','LineWidth',5)
    else
        findClosest = round(s.SoundData.(soundNames{2}).UniqueFreqs) - s.SoundData.(soundNames{3}).Frequency;
        smallestDiff = min(abs(findClosest));
        findClosest = find(abs(findClosest) == smallestDiff);
        rectangle('Position',[findClosest-0.5 find(round(s.SoundData.(soundNames{2}).UniqueDBs) == s.SoundData.(soundNames{3}).dB)-0.5 1 1],'EdgeColor','g','LineWidth',5)
    end
    
    %plots heatmap for second tuning curve.
    subplot(4,3,7)
    dataPrep2 = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).BinSpikeStats(:,:,params.chosenSpikeBin);
    if cMinMax(1) == cMinMax(2)
        imagesc(dataPrep2')
    else
        imagesc(dataPrep2', cMinMax)
    end
    colormap hot
    colorbar
    set(gca,'XTick',s.SoundData.(soundNames{2}).OctaveRange(:,2));
    set(gca,'XTickLabel',s.SoundData.(soundNames{2}).OctaveRange(:,1));
    set(gca,'YTick',s.SoundData.(soundNames{2}).dBRange(:,2));
    set(gca,'YTickLabel',s.SoundData.(soundNames{2}).dBRange(:,1));
    title('Post-Pairing Binned Response to Tone')
    %plots rectangle on target
    if ~isempty(find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5)
        rectangle('Position',[find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5...
            find(round(s.SoundData.(soundNames{2}).UniqueDBs) == s.SoundData.(soundNames{3}).dB)-0.5 1 1],'EdgeColor','g','LineWidth',5)
    else
        findClosest = round(s.SoundData.(soundNames{2}).UniqueFreqs) - s.SoundData.(soundNames{3}).Frequency;
        smallestDiff = min(abs(findClosest));
        findClosest = find(abs(findClosest) == smallestDiff);
        rectangle('Position',[findClosest-0.5 find(round(s.SoundData.(soundNames{2}).UniqueDBs) == s.SoundData.(soundNames{3}).dB)-0.5 1 1],'EdgeColor','g','LineWidth',5)
    end
    %plots difference of two heatmaps. This is raw subtraction of total
    %spikes.
    subplot(4,3,10)
    imagesc((dataPrep2-dataPrep1)')
    colormap hot
    colorbar
    set(gca,'XTick',s.SoundData.(soundNames{2}).OctaveRange(:,2));
    set(gca,'XTickLabel',s.SoundData.(soundNames{2}).OctaveRange(:,1));
    set(gca,'YTick',s.SoundData.(soundNames{2}).dBRange(:,2));
    set(gca,'YTickLabel',s.SoundData.(soundNames{2}).dBRange(:,1));
    title('Difference in Binned Responses')
    %plots rectangle on target.
    if ~isempty(find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5)
        rectangle('Position',[find(round(s.SoundData.(soundNames{2}).UniqueFreqs) == s.SoundData.(soundNames{3}).Frequency)-0.5...
            find(round(s.SoundData.(soundNames{2}).UniqueDBs) == s.SoundData.(soundNames{3}).dB)-0.5 1 1],'EdgeColor','g','LineWidth',5)
    else
        findClosest = round(s.SoundData.(soundNames{2}).UniqueFreqs) - s.SoundData.(soundNames{3}).Frequency;
        smallestDiff = min(abs(findClosest));
        findClosest = find(abs(findClosest) == smallestDiff);
        rectangle('Position',[findClosest-0.5 find(round(s.SoundData.(soundNames{2}).UniqueDBs) == s.SoundData.(soundNames{3}).dB)-0.5 1 1],'EdgeColor','g','LineWidth',5)
    end
    
    %plots heatmap by frequencies. Plots the first tuning.
    subplot(4,3,5)
    dataPrep1 = s.(desigNames{i}).(strcat(soundNames{1},'Analysis')).FrequencyHistograms;
    imagesc(dataPrep1)
    colorbar
    cMinMax = [0 0];
    %sets limits so that next graph is displayed with the same color
    %settings. 
    cMinMax(1) = min(min(dataPrep1));
    cMinMax(2) = max(max(dataPrep1));
    colormap hot
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
    
    %plots heatmap by frequencies. Plots the second tuning.
    subplot(4,3,8)
    dataPrep2 = s.(desigNames{i}).(strcat(soundNames{2},'Analysis')).FrequencyHistograms;
    if cMinMax(1) == cMinMax(2)
        imagesc(dataPrep2)
    else
        imagesc(dataPrep2, cMinMax)
    end
    colorbar
    colormap hot
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
    
    %plots difference in heatmaps by frequency. Raw subtraction
    subplot(4,3,11)
    imagesc(dataPrep2-dataPrep1)
    colormap hot
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
    
    %plot rasters organized by frequency and dB
    subplot(4,3,3)
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
    title({'Pre-Pairing Rasters';'dB Increase in Descending Order'})
    
    subplot(4,3,6)
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
    title({'Post-Pairing Rasters';'dB Increase in Descending Order'})
    
    %plot rasters from pairing
    subplot(4,3,9)
    plot(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).AllRasters(:,1),...
        s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).AllRasters(:,2),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([s.SoundData.(soundNames{3}).ToneDuration s.SoundData.(soundNames{3}).ToneDuration],[ylim],'b');
    ylim([0 s.SoundData.(soundNames{3}).ToneRepetitions])
    xlim([params.rasterWindow(1)*s.SoundData.(soundNames{3}).ToneDuration params.rasterWindow(2)*s.SoundData.(soundNames{3}).ToneDuration])
    title('Pairing Rasters')
    
    %plots binned responses relative to pairing presentation
    subplot(4,3,12)
    plot(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).BinSpikeTimes)
    xlim([1 size(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).BinSpikeTimes,1)])
    try
        ylim([min(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).BinSpikeTimes)...
        max(s.(desigNames{i}).(strcat(soundNames{3},'Analysis')).BinSpikeTimes)])
    catch
    end
    ylabel('Binned Spikes')
    xlabel('Pairing Trial Number')
    title('Binned Responses During Pairing')
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