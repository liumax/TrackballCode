%This function is meant to go through laser data and plot out full
%figures, to allow for selection of laser (1) and non-laser (0) units.
%Designations are stored in decisionTuning

%Inputs: 
%s: the structured output of functionBasicTuning or analysisBasicTuning.
%fileName: the name of the file in question. should not have file
%extensions in the name.

%Outputs: 
%decisionTuning: an n x 1 vector with n = numUnits. 1s represent laser
%responsive units, while 0s represent non-responsive units.


function [decisionLaser] = functionLaserSelectionTool(s,fileName);

numUnits = length(s.DesignationName);
desigNames = s.DesignationName;

%set aside array for information about tuning. y means good, n means no
%tuning
decisionLaser = zeros(numUnits,1);

hFig = figure;
for i = 1:numUnits
    set(hFig, 'Position', [700 100 1280 1000])
    %plots average waveform
    subplot(4,2,1)
    hold on
    plot(s.(desigNames{i}).LaserRelated.AverageNormalWave,'k','LineWidth',2)
    plot(s.(desigNames{i}).LaserRelated.AverageLaserWave,'c','LineWidth',2)
    waveCorrelation = corrcoef(s.(desigNames{i}).LaserRelated.AverageLaserWave,s.(desigNames{i}).LaserRelated.AverageNormalWave);

    title({strcat('OverallFiringRate:',num2str(s.(desigNames{i}).OverallFiringRate)),strcat('WaveCorrelation:',num2str(waveCorrelation(2))),...
        strcat('LaserResponsePercentage:',num2str(s.(desigNames{i}).LaserRelated.PercentLaserResponse))})
    %plots ISI
    subplot(4,2,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Params.RPVTime s.Params.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim([s.Params.ClusterWindow(1) s.Params.ClusterWindow(2)])
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    %plots rasters
    subplot(4,1,2)
    plot(s.(desigNames{i}).LaserRelated.LaserRaster(:,1),...
        s.(desigNames{i}).LaserRelated.LaserRaster(:,2),'k.')
    line([0 0],[0 size(s.LaserData.LaserTimes,1)],'LineWidth',1,'Color','blue')
    line([s.LaserData.LaserDuration/1000 s.LaserData.LaserDuration/1000],[0 size(s.LaserData.LaserTimes,1)],'LineWidth',1,'Color','blue')
    line([s.Params.PairingCutoff/1000 s.Params.PairingCutoff/1000],[0 size(s.LaserData.LaserTimes,1)],'LineWidth',2,'Color','red')
    ylim([0 size(s.LaserData.LaserTimes,1)]);
    xlim([s.Params.RasterWindow(1) s.Params.RasterWindow(2)]);
    h = title(strcat(fileName,desigNames{i}));
    set(h,'interpreter','none') 
    %plots histogram
    subplot(4,1,3)
    plot(s.(desigNames{i}).LaserRelated.LaserHistogram(:,2),...
        s.(desigNames{i}).LaserRelated.LaserHistogram(:,1),'k')
    line([0 0],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',1,'Color','blue')
    line([s.LaserData.LaserDuration/1000 s.LaserData.LaserDuration/1000],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',1,'Color','blue')
    line([s.Params.PairingCutoff/1000 s.Params.PairingCutoff/1000],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',2,'Color','red')
    title(strcat('Mean Laser Dur:',num2str(s.LaserData.LaserDuration),'ms'))
    %plots zScore
    subplot(4,1,4)
    plot(s.(desigNames{i}).LaserRelated.LaserHistogram(:,2),...
        s.(desigNames{i}).LaserRelated.zScoreHistogram,'k')
    line([0 0],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','blue')
    line([s.LaserData.LaserDuration/1000 s.LaserData.LaserDuration/1000],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','blue')
    line([s.Params.PairingCutoff/1000 s.Params.PairingCutoff/1000],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',2,'Color','red')
    %if there is a threshold crossing, plots it
    if ~isempty(s.(desigNames{i}).LaserRelated.FirstZCrossing)
        line([s.(desigNames{i}).LaserRelated.FirstZCrossing s.(desigNames{i}).LaserRelated.FirstZCrossing]...
            ,[0 max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','green')
        title(strcat('CELL IS LASER RESPONSIVE WITH ',num2str(s.(desigNames{i}).LaserRelated.FirstZCrossing*1000),' ms DELAY'))
    end
    if min(s.(desigNames{i}).LaserRelated.zScoreHistogram) ~= max(s.(desigNames{i}).LaserRelated.zScoreHistogram)
        ylim([min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)])
    end
    hold off
    
    %ask for input! 
    promptCounter = 1; %This is used to run the while loop.
    whileCounter = 0; %this is the counter that gets updated to exit the loop
    
    while whileCounter ~= promptCounter
        try
            prompt = 'Is this unit laser responsive? (y/n)';
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
        decisionLaser(i) = 1;
    elseif strfind(str,'n')
        decisionLaser(i) = 0;
    end
    
    %clear figure.
    clf
end

close

end