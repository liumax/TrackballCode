%This function is meant to go through tuning data and plot out full
%figures, to allow for selection of tuned (1) and untuned (0) units.
%Designations are stored in decisionTuning

%Inputs: 
%s: the structured output of functionBasicTuning or analysisBasicTuning.
%fileName: the name of the file in question. should not have file
%extensions in the name.

%Outputs: 
%decisionTuning: an n x 1 vector with n = numUnits. 1s represent tuned
%units, 0s represent untuned units. 


function [decisionTuning] = functionTuningSelectionTool(s,fileName);

numUnits = length(s.DesignationName);
desigNames = s.DesignationName;

toneDur = s.SoundData.ToneDuration;
toneReps = s.SoundData.ToneRepetitions;
numDBs = s.SoundData.NumDBs;
numFreqs = s.SoundData.NumFreqs;
uniqueFreqs = s.SoundData.UniqueFrequencies;
uniqueDBs = s.SoundData.UniqueDBs;
totalTrialNum = length(s.SoundData.Frequencies);

histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2];

%set aside array for information about tuning. y means good, n means no
%tuning
decisionTuning = zeros(numUnits,1);

hFig = figure;
for i = 1:numUnits
    set(hFig, 'Position', [700 100 1280 1000])
    %plots average waveform
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms(:,2),'LineWidth',2)
    plot(s.(desigNames{i}).AverageWaveForms(:,1),'r','LineWidth',1)
    plot(s.(desigNames{i}).AverageWaveForms(:,3),'r','LineWidth',1)
    title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    %plots ISI
    subplot(4,6,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    %plots first spike latency
    subplot(4,3,4)
    imagesc(s.(desigNames{i}).FirstSpikeStats(:,:,1,s.Parameters.ChosenSpikeBin)')
    colormap hot
    colorbar
    set(gca,'XTick',s.Parameters.OctaveRange(:,2));
    set(gca,'XTickLabel',s.Parameters.OctaveRange(:,1));
    set(gca,'YTick',s.Parameters.DBRange(:,2));
    set(gca,'YTickLabel',s.Parameters.DBRange(:,1));
    title('Mean First Spike Latency')
    %plots heatmap of binned spikes to the chosen spike timing window.
    subplot(4,3,7)
    imagesc(squeeze(s.(desigNames{i}).BinSpikeStats(:,:,1,s.Parameters.ChosenSpikeBin))')
    colormap hot
    colorbar
    set(gca,'XTick',s.Parameters.OctaveRange(:,2));
    set(gca,'XTickLabel',s.Parameters.OctaveRange(:,1));
    set(gca,'YTick',s.Parameters.DBRange(:,2));
    set(gca,'YTickLabel',s.Parameters.DBRange(:,1));
    title('Binned Response')
    %plots heatmaps of response reliability in chosen bin 
    subplot(4,3,10)
    imagesc(squeeze(s.(desigNames{i}).FirstSpikeStats(:,:,3,s.Parameters.ChosenSpikeBin))')
    colormap hot
    colorbar
    set(gca,'XTick',s.Parameters.OctaveRange(:,2));
    set(gca,'XTickLabel',s.Parameters.OctaveRange(:,1));
    set(gca,'YTick',s.Parameters.DBRange(:,2));
    set(gca,'YTickLabel',s.Parameters.DBRange(:,1));
    title('Probability of Response')
    %plots rasters (chronological)
    subplot(3,3,2)
    plot(s.(desigNames{i}).AllRasters(:,1),...
        s.(desigNames{i}).AllRasters(:,2),'k.','markersize',4)
    hold on
    ylim([0 totalTrialNum])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    title({fileName;desigNames{i}})
    set(0, 'DefaulttextInterpreter', 'none')
    %plots rasters (frequency and amplitude organized)
    subplot(3,3,5)
    plot(s.(desigNames{i}).AllRasters(:,1),...
        s.(desigNames{i}).AllRasters(:,3),'k.','markersize',4)
    hold on
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    rasterFreqLines = zeros(numFreqs,2);
    rasterFreqLines(:,1) = toneReps*size(uniqueDBs,1)/2:toneReps*size(uniqueDBs,1):totalTrialNum;
    rasterFreqLines(:,2) = uniqueFreqs;
    %this generates green lines separating by Frequency
    for k = 1:size(uniqueFreqs,1)
        plot(s.Parameters.RasterWindow,[toneReps*numDBs*k toneReps*numDBs*k],'g','LineWidth',1)
    end
    set(gca,'YTick',rasterFreqLines(:,1));
    set(gca,'YTickLabel',rasterFreqLines(:,2));
    set(gca,'Ydir','reverse')
    ylim([0 totalTrialNum])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    title('Descending = increase in amplitude and freq')
    %plot heatmap organized by frequency
    subplot(3,3,8)
    imagesc(s.(desigNames{i}).FrequencyHistograms(:,:))
    colorbar
    set(gca,'YTick',s.Parameters.OctaveRange(:,2));
    set(gca,'YTickLabel',s.Parameters.OctaveRange(:,1));
    set(gca,'XTick',[1:10:size(histBinVector,2)]);
    set(gca,'XTickLabel',histBinVector(1:20:end));
    histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
    histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),toneDur);
    line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
    line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
    line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
    line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
%         title('Heatmap by Frequency and Time Max')
    title('Frequency Arranged Heatmap')
    % plot histogram.
    subplot(4,3,3)
    plot(histBinVector,s.(desigNames{i}).AllHistograms,'k','LineWidth',2)
    hold on
    plot(histBinVector,s.(desigNames{i}).AllHistograms - s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
    plot(histBinVector,s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
    plot([0 0],[ylim],'b');
    plot([toneDur toneDur],[ylim],'b');
    %plot significant values
    plot(s.(desigNames{i}).AllHistogramSig.Centers(...
        s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1),...
        s.(desigNames{i}).AllHistogramSig.Histogram(...
        s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1,1),...
        'b*')
    plot(s.(desigNames{i}).AllHistogramSig.Centers(...
        s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1),...
        s.(desigNames{i}).AllHistogramSig.Histogram(...
        s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1,1),...
        'c*')
    %plot negative values for first tuning
    plot(s.(desigNames{i}).AllHistogramSig.Centers(...
        s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
        s.(desigNames{i}).AllHistogramSig.Histogram(...
        s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
        'k*')
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    title('Histogram')
    
    hold off
    
    %ask for input! 
    promptCounter = 1; %This is used to run the while loop.
    whileCounter = 0; %this is the counter that gets updated to exit the loop
    
    while whileCounter ~= promptCounter
        try
            prompt = 'Is this unit tuned in some way? (y/n)';
            str = input(prompt,'s');
            if str~='n' & str~='y'
                error
            else
                whileCounter = 1;
            end
        catch
        end
    end
    if strfind(str,'y')
        decisionTuning(i) = 1;
    elseif strfind(str,'n')
        decisionTuning(i) = 0;
    end
    
    %clear figure.
    clf
end

close

end