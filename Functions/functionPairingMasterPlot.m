%plotting function that plots out information for simple tuning curves.

function [masterStruct] = functionPairingMasterPlot(numTrodes,masterStruct,...
    truncatedNames,rpvTime,clusterWindow,rasterWindow,histBin,...
    clims1,fileName,trodesDesignation,names);

histBinVector = [(rasterWindow(1)*masterStruct.SoundData.Tuning1.ToneDur) + ...
    histBin/2:histBin:(rasterWindow(2)*masterStruct.SoundData.Tuning1.ToneDur)-histBin/2];

master = masterStruct.SoundData.(names{1}).MasterArray;
numFreqs = size(unique(master(:,2)),1);

for i = 1:numTrodes
    for j = 1:masterStruct.(truncatedNames{i}).Clusters
        hFig = figure;
        spikeGraphName = strcat(fileName,trodesDesignation{i},' Cluster ',num2str(j),'SpikeAnalysis');
        set(hFig, 'Position', [5 5 1280 1000])
        %plots average waveform
        subplot(4,3,1)
        hold on
        plot(masterStruct.(truncatedNames{i}).AverageWaveForms(:,j,2),'LineWidth',2)
        plot(masterStruct.(truncatedNames{i}).AverageWaveForms(:,j,1),'r','LineWidth',1)
        plot(masterStruct.(truncatedNames{i}).AverageWaveForms(:,j,3),'r','LineWidth',1)
        title(strcat(spikeGraphName,'AverageFiringRate:',num2str(masterStruct.(truncatedNames{i}).AverageFiringRates(j))))
        %plots ISI
        subplot(4,3,4)
        hist(masterStruct.(truncatedNames{i}).ISIData{j},1000)
        histMax = max(hist(masterStruct.(truncatedNames{i}).ISIData{j},1000));
        line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
        xlim(clusterWindow)
        title(strcat('ISI RPV %: ',num2str(masterStruct.(truncatedNames{i}).RPVs(j))))
        
        %plots histogram
        subplot(2,3,4)
        %plots histogram from first tuning.
        plot(masterStruct.(truncatedNames{i}).(names{1}).Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{1}).Histogram{j}(:,1),'k','LineWidth',2)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{1}).Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{1}).StandardErrorPlotting(:,j,1),'g')
        plot(masterStruct.(truncatedNames{i}).(names{1}).Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{1}).StandardErrorPlotting(:,j,2),'g')
        %plots histogram from second tuning.
        plot(masterStruct.(truncatedNames{i}).(names{2}).Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{2}).Histogram{j}(:,1),'b','LineWidth',2)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{2}).Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{2}).StandardErrorPlotting(:,j,1),'c')
        plot(masterStruct.(truncatedNames{i}).(names{2}).Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{2}).StandardErrorPlotting(:,j,2),'c')
        %draws in tone!
        plot([0 0],[ylim],'r');
        plot([masterStruct.SoundData.(names{1}).ToneDur ...
            masterStruct.SoundData.(names{1}).ToneDur],[ylim],'r');
        %sets xlimits to avoid awkward graphs
        xlim([rasterWindow(1)*masterStruct.SoundData.(names{1}).ToneDur...
            rasterWindow(2)*masterStruct.SoundData.(names{1}).ToneDur])
        title('Histogram K before B after')
        
        %plots heatmap. This uses a log10 scaling for change in firing
        %rates. This way, no change is essentially zero on the imagesc
        %scale, rather than having a linear scale where green actually
        %represents a fairly large increase in response, and there is
        %little room for inhibition
        subplot(3,3,2)
        %plots heatmap for first tuning curve.
        imagesc(masterStruct.(truncatedNames{i}).(names{1}).FrequencyResponse{j})
        colormap hot
        set(gca,'XTick',masterStruct.SoundData.(names{1}).OctaveRange(:,2));
        set(gca,'XTickLabel',masterStruct.SoundData.(names{1}).OctaveRange(:,1));
        set(gca,'YTick',masterStruct.SoundData.(names{1}).dBRange(:,2));
        set(gca,'YTickLabel',masterStruct.SoundData.(names{1}).dBRange(:,1));
        title(strcat('Normalized Frequency Response.Max',...
            num2str(max(max(masterStruct.(truncatedNames{i}).(names{1}).FrequencyResponse{j}))....
            /masterStruct.(truncatedNames{i}).AverageFiringRates(j)),...
            'Min',num2str(min(min(masterStruct.(truncatedNames{i}).(names{1}).FrequencyResponse{j}))/...
            masterStruct.(truncatedNames{i}).AverageFiringRates(j))))
        %plots heatmap for second tuning curve.
        subplot(3,3,5)
        imagesc(masterStruct.(truncatedNames{i}).(names{2}).FrequencyResponse{j})
        colormap hot
        set(gca,'XTick',masterStruct.SoundData.(names{2}).OctaveRange(:,2));
        set(gca,'XTickLabel',masterStruct.SoundData.(names{2}).OctaveRange(:,1));
        set(gca,'YTick',masterStruct.SoundData.(names{2}).dBRange(:,2));
        set(gca,'YTickLabel',masterStruct.SoundData.(names{2}).dBRange(:,1));
        title(strcat('Normalized Frequency Response.Max',...
            num2str(max(max(masterStruct.(truncatedNames{i}).(names{2}).FrequencyResponse{j}))....
            /masterStruct.(truncatedNames{i}).AverageFiringRates(j)),...
            'Min',num2str(min(min(masterStruct.(truncatedNames{i}).(names{2}).FrequencyResponse{j}))/...
            masterStruct.(truncatedNames{i}).AverageFiringRates(j))))
        %plots difference of two heatmaps. This is raw subtraction of total
        %spikes.
        subplot(3,3,8)
        imagesc(masterStruct.(truncatedNames{i}).(names{2}).FrequencyResponse{j} -...
            masterStruct.(truncatedNames{i}).(names{1}).FrequencyResponse{j})
        colormap hot
        set(gca,'XTick',masterStruct.SoundData.(names{2}).OctaveRange(:,2));
        set(gca,'XTickLabel',masterStruct.SoundData.(names{2}).OctaveRange(:,1));
        set(gca,'YTick',masterStruct.SoundData.(names{2}).dBRange(:,2));
        set(gca,'YTickLabel',masterStruct.SoundData.(names{2}).dBRange(:,1));
        %this will plot a colored box around the target and control
        %frequencies
        %control
        hold on
        if ~isempty(find(round(masterStruct.SoundData.(names{2}).UniqueFreqs) == masterStruct.SoundData.(names{5}).ControlSet(1))-0.5)
            rectangle('Position',[find(round(masterStruct.SoundData.(names{2}).UniqueFreqs) == masterStruct.SoundData.(names{5}).ControlSet(1))-0.5...
                find(round(masterStruct.SoundData.(names{2}).UniqueDBs) == masterStruct.SoundData.(names{5}).ControlSet(2))-0.5 1 1],'EdgeColor','b','LineWidth',5)
        else
            findClosest = round(masterStruct.SoundData.(names{2}).UniqueFreqs) - masterStruct.SoundData.(names{5}).ControlSet(1);
            smallestDiff = min(abs(findClosest));
            findClosest = find(abs(findClosest) == smallestDiff);
            rectangle('Position',[findClosest-0.5 find(round(masterStruct.SoundData.(names{2}).UniqueDBs) == masterStruct.SoundData.(names{5}).ControlSet(2))-0.5 1 1],'EdgeColor','b','LineWidth',5)
        end
        %target!
        if ~isempty(find(round(masterStruct.SoundData.(names{2}).UniqueFreqs) == masterStruct.SoundData.(names{5}).TargetSet(1))-0.5)
            rectangle('Position',[find(round(masterStruct.SoundData.(names{2}).UniqueFreqs) == masterStruct.SoundData.(names{5}).TargetSet(1))-0.5...
                find(round(masterStruct.SoundData.(names{2}).UniqueDBs) == masterStruct.SoundData.(names{5}).TargetSet(2))-0.5 1 1],'EdgeColor','g','LineWidth',5)
        else
            findClosest = round(masterStruct.SoundData.(names{2}).UniqueFreqs) - masterStruct.SoundData.(names{5}).TargetSet(1);
            smallestDiff = min(abs(findClosest));
            findClosest = find(abs(findClosest) == smallestDiff);
            rectangle('Position',[findClosest-0.5 find(round(masterStruct.SoundData.(names{2}).UniqueDBs) == masterStruct.SoundData.(names{5}).TargetSet(2))-0.5 1 1],'EdgeColor','g','LineWidth',5)
        end
        
        title(strcat('Difference in Frequency Response.Max',...
            num2str(max(max(max(masterStruct.(truncatedNames{i}).(names{1}).FrequencyResponse{j}...
            -masterStruct.(truncatedNames{i}).(names{2}).FrequencyResponse{j})))),...
            'Min',num2str(min(min(min(masterStruct.(truncatedNames{i}).(names{1}).FrequencyResponse{j}...
            -masterStruct.(truncatedNames{i}).(names{2}).FrequencyResponse{j}))))))
        
        %plots heatmap by frequencies. Plots the first tuning.
        subplot(3,3,3)
        x = masterStruct.(truncatedNames{i}).(names{1}).AverageFrequencyHistogram{j};
        imagesc(x')
        colormap hot
        set(gca,'YTick',masterStruct.SoundData.(names{2}).OctaveRange(:,2));
        set(gca,'YTickLabel',masterStruct.SoundData.(names{2}).OctaveRange(:,1));
        set(gca,'XTick',[1:10:size(histBinVector,2)]);
        set(gca,'XTickLabel',histBinVector(1:10:end));
        histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
        histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),masterStruct.SoundData.(names{1}).ToneDur);
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
        %         title('Heatmap by Frequency and Time Max')
        title(strcat('Normalized Heatmap by F and T.Max',...
            num2str(max(max(masterStruct.(truncatedNames{i}).(names{1}).AverageFrequencyHistogram{j}))/masterStruct.(truncatedNames{i}).AverageFiringRates(j)),...
            'Min',num2str(min(min(masterStruct.(truncatedNames{i}).(names{1}).AverageFrequencyHistogram{j}))/masterStruct.(truncatedNames{i}).AverageFiringRates(j))))
        %plots second tuning
        subplot(3,3,6)
        x = masterStruct.(truncatedNames{i}).(names{2}).AverageFrequencyHistogram{j};
        imagesc(x')
        colormap hot
        set(gca,'YTick',masterStruct.SoundData.(names{2}).OctaveRange(:,2));
        set(gca,'YTickLabel',masterStruct.SoundData.(names{2}).OctaveRange(:,1));
        set(gca,'XTick',[1:10:size(histBinVector,2)]);
        set(gca,'XTickLabel',histBinVector(1:10:end));
        histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
        histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),masterStruct.SoundData.(names{2}).ToneDur);
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
        %         title('Heatmap by Frequency and Time Max')
        title(strcat('Normalized Heatmap by F and T.Max',...
            num2str(max(max(masterStruct.(truncatedNames{i}).(names{2}).AverageFrequencyHistogram{j}))/masterStruct.(truncatedNames{i}).AverageFiringRates(j)),...
            'Min',num2str(min(min(masterStruct.(truncatedNames{i}).(names{2}).AverageFrequencyHistogram{j}))/masterStruct.(truncatedNames{i}).AverageFiringRates(j))))
        %plots difference. raw subtraction
        subplot(3,3,9)
        x = masterStruct.(truncatedNames{i}).(names{2}).AverageFrequencyHistogram{j}-masterStruct.(truncatedNames{i}).(names{1}).AverageFrequencyHistogram{j};
        imagesc((x'))
        set(gca,'YTick',masterStruct.SoundData.(names{2}).OctaveRange(:,2));
        set(gca,'YTickLabel',masterStruct.SoundData.(names{2}).OctaveRange(:,1));
        set(gca,'XTick',[1:10:size(histBinVector,2)]);
        set(gca,'XTickLabel',histBinVector(1:10:end));
        histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
        histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),masterStruct.SoundData.(names{1}).ToneDur);
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
        line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
        title(strcat('Difference Heatmap by F and T. Max',...
            num2str(max(max(x))),...
            'Min',num2str(min(min(x)))))
        %% GENERATES SECOND FIGURE. This one will include information about
        %long tone presentations. EDITS TO HERE.
        hFig2 = figure;
        set(hFig2, 'Position', [5 5 1280 1000])
        subplot(2,2,1)
        %% plots histograms for control. for first presentation
        plot(masterStruct.(truncatedNames{i}).(names{3}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{3}).Control.Histogram{j}(:,1),'k','LineWidth',2)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{3}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{3}).Control.StandardErrorPlotting(:,j,1),'g')
        plot(masterStruct.(truncatedNames{i}).(names{3}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{3}).Control.StandardErrorPlotting(:,j,2),'g')
        
        %plots histograms for control. for second presentation
        plot(masterStruct.(truncatedNames{i}).(names{4}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{4}).Control.Histogram{j}(:,1),'b','LineWidth',2)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{4}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{4}).Control.StandardErrorPlotting(:,j,1),'c')
        plot(masterStruct.(truncatedNames{i}).(names{4}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{4}).Control.StandardErrorPlotting(:,j,2),'c')
        
        %plots histograms for control. for pairing presentation
        plot(masterStruct.(truncatedNames{i}).(names{5}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{5}).Control.Histogram{j}(:,1),'r','LineWidth',2)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{5}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{5}).Control.StandardErrorPlotting(:,j,1),'m')
        plot(masterStruct.(truncatedNames{i}).(names{5}).Control.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{5}).Control.StandardErrorPlotting(:,j,2),'m')
        
        %draws in tone!
        plot([0 0],[ylim],'r');
        plot([masterStruct.SoundData.(names{3}).ToneDur ...
            masterStruct.SoundData.(names{3}).ToneDur],[ylim],'r');
        %sets xlimits to avoid awkward graphs
        xlim([rasterWindow(1)*masterStruct.SoundData.(names{3}).ToneDur...
            rasterWindow(2)*masterStruct.SoundData.(names{3}).ToneDur])
        title('Histogram of Control Responses')
        %% plots histograms for Target. for first presentation
        subplot(2,2,3)
        plot(masterStruct.(truncatedNames{i}).(names{3}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{3}).Target.Histogram{j}(:,1),'k','LineWidth',2)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{3}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{3}).Target.StandardErrorPlotting(:,j,1),'g')
        plot(masterStruct.(truncatedNames{i}).(names{3}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{3}).Target.StandardErrorPlotting(:,j,2),'g')
        
        %plots histograms for Target. for second presentation
        plot(masterStruct.(truncatedNames{i}).(names{4}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{4}).Target.Histogram{j}(:,1),'b','LineWidth',2)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{4}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{4}).Target.StandardErrorPlotting(:,j,1),'c')
        plot(masterStruct.(truncatedNames{i}).(names{4}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{4}).Target.StandardErrorPlotting(:,j,2),'c')
        
        %plots histograms for Target. for pairing presentation
        plot(masterStruct.(truncatedNames{i}).(names{5}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{5}).Target.Histogram{j}(:,1),'r','LineWidth',2)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{5}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{5}).Target.StandardErrorPlotting(:,j,1),'m')
        plot(masterStruct.(truncatedNames{i}).(names{5}).Target.Histogram{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{5}).Target.StandardErrorPlotting(:,j,2),'m')
        
        %draws in tone!
        plot([0 0],[ylim],'r');
        plot([masterStruct.SoundData.(names{3}).ToneDur ...
            masterStruct.SoundData.(names{3}).ToneDur],[ylim],'r');
        %sets xlimits to avoid awkward graphs
        xlim([rasterWindow(1)*masterStruct.SoundData.(names{3}).ToneDur...
            rasterWindow(2)*masterStruct.SoundData.(names{3}).ToneDur])
        title('Histogram of Target Responses')
        %% Plot the rasters for controls
        %plots simple rasters for all presentations, with different parts
        %in different colors
        
        %plots simple rasters for first presentation
        subplot(2,2,2)
        plot(masterStruct.(truncatedNames{i}).(names{3}).Control.Rasters{j}(:,2),...%plots the initial presentations
            masterStruct.(truncatedNames{i}).(names{3}).Control.Rasters{j}(:,1),'k.','markersize',4)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{5}).Control.Rasters{j}(:,2),...%plots the pairing presentations
            masterStruct.(truncatedNames{i}).(names{5}).Control.Rasters{j}(:,1)+masterStruct.SoundData.(names{3}).ToneReps,'r.','markersize',4)
        plot(masterStruct.(truncatedNames{i}).(names{4}).Control.Rasters{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{4}).Control.Rasters{j}(:,1)+(masterStruct.SoundData.(names{3}).ToneReps + masterStruct.SoundData.(names{5}).ToneReps),'b.','markersize',4)
        ylim([0 masterStruct.SoundData.(names{3}).ToneReps + masterStruct.SoundData.(names{4}).ToneReps + masterStruct.SoundData.(names{5}).ToneReps])
        xlim([rasterWindow(1)*masterStruct.SoundData.(names{3}).ToneDur...
            rasterWindow(2)*masterStruct.SoundData.(names{3}).ToneDur])
        plot([0 0],[ylim],'r');
        plot([masterStruct.SoundData.(names{3}).ToneDur...
            masterStruct.SoundData.(names{3}).ToneDur],[ylim],'r');
        title('Control Rasters')
        
        %% Plots rasters for Targets.
        subplot(2,2,4)
        plot(masterStruct.(truncatedNames{i}).(names{3}).Target.Rasters{j}(:,2),...%plots the initial presentations
            masterStruct.(truncatedNames{i}).(names{3}).Target.Rasters{j}(:,1),'k.','markersize',4)
        hold on
        plot(masterStruct.(truncatedNames{i}).(names{5}).Target.Rasters{j}(:,2),...%plots the pairing presentations
            masterStruct.(truncatedNames{i}).(names{5}).Target.Rasters{j}(:,1)+masterStruct.SoundData.(names{3}).ToneReps,'r.','markersize',4)
        plot(masterStruct.(truncatedNames{i}).(names{4}).Target.Rasters{j}(:,2),...
            masterStruct.(truncatedNames{i}).(names{4}).Target.Rasters{j}(:,1)+(masterStruct.SoundData.(names{3}).ToneReps + masterStruct.SoundData.(names{5}).ToneReps),'b.','markersize',4)
        ylim([0 masterStruct.SoundData.(names{3}).ToneReps + masterStruct.SoundData.(names{4}).ToneReps + masterStruct.SoundData.(names{5}).ToneReps])
        xlim([rasterWindow(1)*masterStruct.SoundData.(names{3}).ToneDur...
            rasterWindow(2)*masterStruct.SoundData.(names{3}).ToneDur])
        plot([0 0],[ylim],'r');
        plot([masterStruct.SoundData.(names{3}).ToneDur...
            masterStruct.SoundData.(names{3}).ToneDur],[ylim],'r');
        title('Target Rasters')
        
        %%
        hold off
                %save as matlab figure with correct name (fileName+LFP)
                spikeGraphName = strcat(fileName,trodesDesignation{i},'Cluster',num2str(j),'SpikeAnalysis');
%                 savefig(hFig,spikeGraphName);
                spikeGraphName2 = strcat(fileName,trodesDesignation{i},'Cluster',num2str(j),'SpikeAnalysis2');
%                 savefig(hFig2,spikeGraphName2);
        
                %save as PDF with correct name
                set(hFig,'Units','Inches');
                pos = get(hFig,'Position');
                set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
                print(hFig,spikeGraphName,'-dpdf','-r0')
                set(hFig2,'Units','Inches');
                pos = get(hFig2,'Position');
                set(hFig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
                print(hFig2,spikeGraphName2,'-dpdf','-r0')
    end
end
end