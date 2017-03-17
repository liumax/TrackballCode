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


function [decisionTuning,tuningType] = functionTuningSelectionTool(s,fileName);

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

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
tuningType = cell(numUnits,1);

hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    for i = 1:numUnits
        %plots average waveform
        subplot(4,6,1)
        hold on
        plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
        title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
        %plots ISI
        subplot(4,6,2)
        hist(s.(desigNames{i}).ISIGraph,1000)
        histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
        line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
        xlim(s.Parameters.ClusterWindow)
        title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
            strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
        %Plot binned response during tone period
        subplot(4,3,4)
        imagesc(s.(desigNames{i}).BinTone')
        colorbar
        set(gca,'XTick',s.Parameters.OctaveRange(:,2));
        set(gca,'XTickLabel',s.Parameters.OctaveRange(:,1));
        set(gca,'YTick',s.Parameters.DBRange(:,2));
        set(gca,'YTickLabel',s.Parameters.DBRange(:,1));
        title('Mean Binned Response (tone)')
        %Plot binned response during general period
        subplot(4,3,7)
        imagesc(s.(desigNames{i}).BinGen')
        colorbar
        set(gca,'XTick',s.Parameters.OctaveRange(:,2));
        set(gca,'XTickLabel',s.Parameters.OctaveRange(:,1));
        set(gca,'YTick',s.Parameters.DBRange(:,2));
        set(gca,'YTickLabel',s.Parameters.DBRange(:,1));
        title('Mean Binned Response (general)')
        %Plot peak response during general period
        subplot(4,3,10)
        imagesc(s.(desigNames{i}).PeakMap')
        colorbar
        set(gca,'XTick',s.Parameters.OctaveRange(:,2));
        set(gca,'XTickLabel',s.Parameters.OctaveRange(:,1));
        set(gca,'YTick',s.Parameters.DBRange(:,2));
        set(gca,'YTickLabel',s.Parameters.DBRange(:,1));
        title('Peak Response (general)')
        %plot velocity data
        subplot(4,3,6)
        hold on
        plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'b')
        plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
        xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
        ylim([-0.1,1])
        if isfield(s.(desigNames{i}),'TrueAUC')
            title(strcat('Vel & Firing Rate. AUC:',num2str(s.(desigNames{i}).TrueAUC),'99%Range',num2str(prctile(s.(desigNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(desigNames{i}).ShuffleAUC,1))))
        else
            title('Vel & Firing Rate')
        end
        %plot probability of response (tone)
        subplot(4,3,9)
        imagesc(s.(desigNames{i}).ProbTone')
        colorbar
        set(gca,'XTick',s.Parameters.OctaveRange(:,2));
        set(gca,'XTickLabel',s.Parameters.OctaveRange(:,1));
        set(gca,'YTick',s.Parameters.DBRange(:,2));
        set(gca,'YTickLabel',s.Parameters.DBRange(:,1));
        title('Probability of Response (tone)')
        %plot probability of response (gen)
        subplot(4,3,12)
        imagesc(s.(desigNames{i}).ProbGen')
        colorbar
        set(gca,'XTick',s.Parameters.OctaveRange(:,2));
        set(gca,'XTickLabel',s.Parameters.OctaveRange(:,1));
        set(gca,'YTick',s.Parameters.DBRange(:,2));
        set(gca,'YTickLabel',s.Parameters.DBRange(:,1));
        title('Probability of Response (general)')

        %plots rasters (chronological)
        subplot(3,3,2)
        plot(s.(desigNames{i}).AllRasters(:,1),...
            s.(desigNames{i}).AllRasters(:,2),'k.','markersize',4)
        hold on
        ylim([0 totalTrialNum])
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        title({fileName;desigNames{i}},'fontweight','bold')
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
            'bo')
        %plot negative values for first tuning
        plot(s.(desigNames{i}).AllHistogramSig.Centers(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
            s.(desigNames{i}).AllHistogramSig.Histogram(...
            s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
            'k*')
        plot([0 0],[ylim],'b');
        plot([toneDur toneDur],[ylim],'b');
        xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
        title('Histogram')

        hold off
        spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')

        %ask for input! 
        promptCounter = 1; %This is used to run the while loop.
        whileCounter = 0; %this is the counter that gets updated to exit the loop

        while whileCounter ~= promptCounter
            try
                prompt = 'How is this unit tuned? (excite(e)/inhib(i)/both(b)/none(n))';
                str = input(prompt,'s');
                if str~='n' && str~='e' && str~='i' && str~='b'
                    error
                else
                    whileCounter = 1;
                end
            catch
            end
        end
        if str=='e' || str=='i' || str=='b'
            decisionTuning(i) = 1;
        elseif strfind(str,'n')
            decisionTuning(i) = 0;
        end
        
        tuningType{i} = str;

        %clear figure.
        clf

    end
    s.TuningDecision = decisionTuning;
    s.TuningType = tuningType;

close

end