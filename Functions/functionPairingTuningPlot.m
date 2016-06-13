%plotting function that plots out information for simple tuning curves.



function [masterStruct] = functionPairingTuningPlot(i,...
    truncatedNames,rpvTime,clusterWindow,rasterWindow,histBinVector,...
    clims1,numFreqs,numDBs,uniqueFreqs,uniqueDBs,...
    trodesDesignation,masterStruct,fileName,clusters,soundName);

numTrials = size(masterStruct.SoundData.(soundName).Frequencies,1);
octaveRange = masterStruct.SoundData.(soundName).OctaveRange;
dbRange = masterStruct.SoundData.(soundName).dBRange;

for j = 1:clusters
    hFig = figure;
    set(hFig, 'Position', [10 10 1280 1000])
    %plots average waveform
    subplot(4,3,1)
    hold on
    plot(masterStruct.(truncatedNames{i}).AverageWaveForms(:,j,2),'LineWidth',2)
    plot(masterStruct.(truncatedNames{i}).AverageWaveForms(:,j,1),'r','LineWidth',1)
    plot(masterStruct.(truncatedNames{i}).AverageWaveForms(:,j,3),'r','LineWidth',1)
    title(strcat('AverageFiringRate:',num2str...
        (masterStruct.(truncatedNames{i}).AverageFiringRates(j))))
    %plots ISI
    subplot(4,3,4)
    hist(masterStruct.(truncatedNames{i}).ISIData{j},1000)
    histMax = max(hist(masterStruct.(truncatedNames{i}).ISIData{j},1000));
    line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(clusterWindow)
    title(strcat('ISI RPV %: ',num2str(masterStruct.(truncatedNames{i}).RPVs(j))))

     %plots histogram
    subplot(2,3,4)
    plot(masterStruct.(truncatedNames{i}).(soundName).Histogram{j}(:,2),...
        masterStruct.(truncatedNames{i}).(soundName).Histogram{j}(:,1),'k','LineWidth',2)
    hold on
    plot(masterStruct.(truncatedNames{i}).(soundName).Histogram{j}(:,2),...
        masterStruct.(truncatedNames{i}).(soundName).StandardErrorPlotting(:,j,1),'b')
    plot(masterStruct.(truncatedNames{i}).(soundName).Histogram{j}(:,2),...
        masterStruct.(truncatedNames{i}).(soundName).StandardErrorPlotting(:,j,2),'b')
    plot([0 0],[ylim],'r');
    plot([masterStruct.SoundData.(soundName).ToneDur ...
        masterStruct.SoundData.(soundName).ToneDur],[ylim],'r');
    xlim([rasterWindow(1) rasterWindow(2)])
    title('Histogram')

    %plots simple rasters
    subplot(2,3,2)
    plot(masterStruct.(truncatedNames{i}).(soundName).Rasters{j}(:,2),...
        masterStruct.(truncatedNames{i}).(soundName).Rasters{j}(:,1),'k.','markersize',4)
    hold on
    ylim([0 numTrials])
    xlim([rasterWindow(1) rasterWindow(2)])
    plot([0 0],[ylim],'r');
    plot([masterStruct.SoundData.(soundName).ToneDur...
        masterStruct.SoundData.(soundName).ToneDur],[ylim],'r');

    title(strcat(fileName,truncatedNames{i},' Cluster ',num2str(j)))

    %plots rasters organized by frequency and amp
    subplot(2,3,5)
    plot(masterStruct.(truncatedNames{i}).(soundName).Rasters{j}(:,2),...
        masterStruct.(truncatedNames{i}).(soundName).Rasters{j}(:,5),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'r');
    plot([masterStruct.SoundData.(soundName).ToneDur...
        masterStruct.SoundData.(soundName).ToneDur],[ylim],'r');
    %removed blue lines to be able to see rasters better
%         for k = 1:size(matclustStruct.UniqueFreqs,1)*size(matclustStruct.UniqueDBs,1)
%             plot(matclustStruct.RasterLimits,...
%                 [soundFile.soundData.ToneRepetitions*k soundFile.soundData.ToneRepetitions*k])
%         end
    rasterFreqLines = zeros(numFreqs,2);
    rasterFreqLines(:,1) = masterStruct.SoundData.(soundName).ToneReps*...
        numDBs/2:masterStruct.SoundData.(soundName).ToneReps*numDBs:...
        numTrials;
    rasterFreqLines(:,2) = uniqueFreqs;
    %this generates green lines separating by Frequency
    for k = 1:numFreqs
        plot(rasterWindow,...
            [masterStruct.SoundData.(soundName).ToneReps*numDBs*k ...
            masterStruct.SoundData.(soundName).ToneReps*numDBs*k],...
            'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    ylim([0 numTrials])
    xlim(rasterWindow)
    title('Sorted Ascending')

    %plots heatmap. This uses a log10 scaling for change in firing
    %rates. This way, no change is essentially zero on the imagesc
    %scale, rather than having a linear scale where green actually
    %represents a fairly large increase in response, and there is
    %little room for inhibition
    subplot(4,3,3)
    imagesc(log10(masterStruct.(truncatedNames{i}).(soundName).FrequencyResponse{j}...
        /(masterStruct.SoundData.(soundName).ToneReps*...
        masterStruct.SoundData.(soundName).ToneDur*...
        masterStruct.(truncatedNames{i}).AverageFiringRates(j))),clims1)
    colormap hot
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title(strcat('Normalized Frequency Response.Max',...
        num2str(max(max(masterStruct.(truncatedNames{i}).(soundName).FrequencyResponse{j}))...
        /masterStruct.(truncatedNames{i}).AverageFiringRates(j)),...
        'Min',num2str(min(min(masterStruct.(truncatedNames{i}).(soundName).FrequencyResponse{j}))...
        /masterStruct.(truncatedNames{i}).AverageFiringRates(j))))

    %plots reliability of response in heat map
    subplot(4,3,6)
    clims = [0,1];
    imagesc(masterStruct.SoundData.(soundName).UniqueFreqs,...
        masterStruct.SoundData.(soundName).UniqueDBs,...
        masterStruct.(truncatedNames{i}).(soundName).ResponseReliability{j}',clims)
    colormap hot
%         colormap hot
    set(gca,'XTick',octaveRange(:,2));
    set(gca,'XTickLabel',octaveRange(:,1));
    set(gca,'YTick',dbRange(:,2));
    set(gca,'YTickLabel',dbRange(:,1));
    title(strcat('ResponseReliability.Max',...
        num2str(max(max(masterStruct.(truncatedNames{i}).(soundName).ResponseReliability{j}))),...
        'Min',num2str(min(min(masterStruct.(truncatedNames{i}).(soundName).ResponseReliability{j})))))

    %plots heatmap by frequencies
    subplot(2,3,6)
    x = masterStruct.(truncatedNames{i}).(soundName).AverageFrequencyHistogram{j}...
        /masterStruct.(truncatedNames{i}).AverageFiringRates(j);
    imagesc(log10(x'),clims1)
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    set(gca,'XTick',[1:10:size(histBinVector,2)]);
    set(gca,'XTickLabel',histBinVector(1:10:end));
    histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
    histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),...
        masterStruct.SoundData.(soundName).ToneDur);
    line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
    line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
    line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
    line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
%         title('Heatmap by Frequency and Time Max')
    title(strcat('Normalized Heatmap by F and T.Max',...
        num2str(max(max(masterStruct.(truncatedNames{i}).(soundName).AverageFrequencyHistogram{j}))...
        /masterStruct.(truncatedNames{i}).AverageFiringRates(j)),...
        'Min',num2str(min(min(masterStruct.(truncatedNames{i}).(soundName).AverageFrequencyHistogram{j}))...
        /masterStruct.(truncatedNames{i}).AverageFiringRates(j))))
    hold off
    %save as matlab figure with correct name (fileName+LFP)
    spikeGraphName = strcat(fileName,trodesDesignation{i},' Cluster ',num2str(j),'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end

end