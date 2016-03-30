function [statsGraph] = auditoryStatGraph(batchDataStructure,x);

soundWindow = [0,0.500];

namer = x;
windows = batchDataStructure.(namer).AnalysisVariables.StatsWindow;
bins = batchDataStructure.(namer).AnalysisVariables.StatsBins;

alphas = batchDataStructure.(namer).Statistics.AlphaValues;
tPoints = batchDataStructure.(namer).Statistics.TimePoints;

averageRate = batchDataStructure.(namer).Plotting.AverageFiringRate;
histograms = batchDataStructure.(namer).Plotting.HistogramData;
rasters = batchDataStructure.(namer).Plotting.RasterData;
rasterAxis = batchDataStructure.(namer).Plotting.RasterAxis;
sError = batchDataStructure.(namer).Plotting.HistogramErrorData;

unitWaves = batchDataStructure.(namer).Units.Waves;
unitNames = batchDataStructure.(namer).Units.Names;
unitIndex = batchDataStructure.(namer).Units.UnitIndex;


trialReps = batchDataStructure.(namer).Events.TrialRepetitions;

responseRasters = cell(length(rasters),1);
%finds all responses in the stimulus window
for i = 1:length(rasters)
    responseRasters{i} = rasters{i}(rasters{i}(:,2) > soundWindow(1) &...
        rasters{i}(:,2) < soundWindow(2),:);
end

uniqueResponses = cell(length(rasters),1);
percentResponses = zeros(length(rasters),1);
for i = 1:length(rasters)
    uniqueResponses{i} = unique(responseRasters{i}(:,1));
    percentResponses(i) = length(uniqueResponses{i})/trialReps;
end

%eliminates data outside of the stats window
for i=1:length(histograms)
    stError(:,i,:)= sError(histograms{i}(:,2) > windows(1) &...
        histograms{i}(:,2) < windows(2),i,:);
    histograms{i}= histograms{i}(histograms{i}(:,2) > windows(1) &...
        histograms{i}(:,2) < windows(2),:);
end

% %actually performs plotting
% set(0, 'DefaulttextInterpreter', 'none')
% 
% for i = 1:length(unitIndex)
%     hFig = figure;
%     set(hFig,'Units','inches');
%     set(hFig,'Position',[1 1 6 8]);
% 
%     subplot(3,2,1)
%     plot(unitWaves{i})
%     set(gca, 'Units', 'inches');
%     set(gca,'OuterPosition',[0 5.5 3 2.7]);
%     title('Average Waveform')
% 
%     mTextBox = uicontrol('style','text');
%     descr = {'File:';
%         namer;
%         'Unit:';
%         unitNames{unitIndex(i,2)}
%         'Average Firing Rate (Hz):';
%         averageRate(i)};
%     set(mTextBox,'String',descr);
%     set(mTextBox,'Units','inches');
%     set(mTextBox,'Position',[3,5.5,3,2.5])
% 
%     subplot(3,2,3)
%     plot(histograms{unitIndex(i,3)}(:,2),histograms{unitIndex(i,3)}(:,1),...
%         histograms{unitIndex(i,3)}(:,2),stError(:,unitIndex(i,3),1),'b--',...
%         histograms{unitIndex(i,3)}(:,2),stError(:,unitIndex(i,3),2),'b--')
%     xlim(windows);
%     set(gca, 'Units', 'inches');
%     set(gca,'OuterPosition',[0.1 3 3 2.5]);
%     title('Multiunit Histogram Relative to Tone')
% 
%     subplot(3,2,4)
%     plot(histograms{unitIndex(i,2)}(:,2),histograms{unitIndex(i,2)}(:,1),...
%         histograms{unitIndex(i,2)}(:,2),stError(:,unitIndex(i,2),1),'b--',...
%         histograms{unitIndex(i,2)}(:,2),stError(:,unitIndex(i,2),2),'b--');
%     xlim(windows);
%     set(gca, 'Units', 'inches');
%     set(gca,'OuterPosition',[3.1 3 3 2.5]);
%     title('Unit Histogram Relative to Tone') 
% 
%     subplot(3,2,5)
%     semilogy(tPoints{1},alphas{unitIndex(i,3),1},'r',...
%         tPoints{2},alphas{unitIndex(i,3),2},'b',...
%         tPoints{3},alphas{unitIndex(i,3),3},'k');
%     xlim([0,windows(2)]);
%     grid on
%     set(gca, 'Units', 'inches');
%     set(gca,'OuterPosition',[0.1 0.2 3 2.5]);
%     title('Multiunit Alpha Values (Rank Sum)')
% 
%     subplot(3,2,6)
%     semilogy(tPoints{1},alphas{unitIndex(i,2),1},'r',...
%         tPoints{2},alphas{unitIndex(i,2),2},'b',...
%         tPoints{3},alphas{unitIndex(i,2),3},'k');
%     xlim([0,windows(2)]);
%     grid on
%     set(gca, 'Units', 'inches');
%     set(gca,'OuterPosition',[3.1 0.2 3 2.5]);
%     title('Multiunit Alpha Values (Rank Sum)')     
% end

statsGraph.ResponseTrials = uniqueResponses;
statsGraph.ResponsePercentage = percentResponses;
statsGraph.ResponseRasters = responseRasters;

end

