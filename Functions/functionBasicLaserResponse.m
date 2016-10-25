%This function analyzes data to evaluate recordings in which laser ID was
%performed. The code identifies laser epochs based on digital input 2
%(DIO2), and makes PSTH and raster plots. It then detects whether the
%response (if it exists) exceeds the threshold that is set at the
%beginning. 

%% Inputs
%fileName: name of the file with no extension. This should be the name of
%the sound file. 

%% Outputs

%s: this is a structured array organized as follows: 

%n units as fields, with the name ntXclusterY

%DesignationArray and DesignationName: designation array provides an array
%that decodes which probe index (column 1) and which cluster index (column
%2). The designation name is a cell array with the names of all units.

%LaserData: contains laser times and laser duration.

function [s] = functionBasicLaserResponse(fileName);

%generate placeholder structure
s = struct;

%set parameters.
%store in s.
s.Params.SampleRate = 30000;
s.Params.RasterWindow = [-0.1,0.1];
s.Params.ClusterWindow =[-0.01,0.05];
s.Params.RPVTime = 0.001; %time limit in seconds for consideration as an RPV
s.Params.RasterAxis = [s.Params.RasterWindow(1):0.001:s.Params.RasterWindow(2)-0.001];
s.Params.PairingCutoff = 15; %ms cutoff. Significant responses past this point will not be considered. 
s.Params.EarlyCutoff = 2; %minimum number of ms after a laser before a spike is considered laser related
s.Params.zScoreCutoff = 4; %zscore above which cell is considered IDed
s.Params.HistBin = 0.001; %bin size in seconds


s.Params.HistBin = 0.001; %bin size in seconds
s.Params.HistBinVector = [s.Params.RasterWindow(1)+s.Params.HistBin/2:s.Params.HistBin:s.Params.RasterWindow(2)-s.Params.HistBin/2]; %this is vector with midpoints of all histogram bins
%s.Params.HistBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes.
%% 
%Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';
%% Find and ID Matclust Files for Subsequent Analysis. Generates Structured Array for Data Storage
%pull matclust file names
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');

%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Params.RPVTime,...
    matclustFiles,s,s.Params.ClusterWindow);
%pulls number of units, designation names and array for making
%designations.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%% Extract DIO information from DIO2.
%finds DIO folder, extracts D2 specifically
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};

%extracts DIO stuffs!
[DIOData] = readTrodesExtractedDataFile(D2FileName);

%extracts port states and times of changes
dioState = double(DIOData.fields(2).data);
dioTime = double(DIOData.fields(1).data);

%pulls the times when state goes up!
s.LaserData.LaserTimes = dioTime(dioState == 1)/s.Params.SampleRate;


%pulls the duration of the pulse (this should give laser pulse duration
dioDiff = diff(dioTime);
dioDiff = dioDiff(dioDiff<5000);
s.LaserData.LaserDuration = mean(dioDiff)/30;

%% do actual extraction of data!
for i = 1:numUnits
    %first thing I want to do is separate out spikes that are within the
    %laser period and those that are out of the laser period. Spikes within
    %the laser period are treated as laser evoked spikes, while spikes
    %outside the period are considered baseline spikes. 
    
    %first I want to generate a matrix of n x m, where n is the number of
    %spikes and m is the number of laser pulses. each column will have the
    %spike times minus the time of the laser pulse. 
    
    %allocates zero matrix
    laserMatrix = zeros(length(s.(desigNames{i}).SpikeTimes),length(s.LaserData.LaserTimes));
    %performs subtraction by laser pulse
    for j = 1:length(s.LaserData.LaserTimes)
        laserMatrix(:,j) = s.(desigNames{i}).SpikeTimes - s.LaserData.LaserTimes(j);
    end
    %determines which waveforms are laser related. 
    [rows,cols] = find(laserMatrix > s.Params.EarlyCutoff/1000 & laserMatrix< s.Params.PairingCutoff/1000);
    %safeguard against no response
    if ~isempty(rows)
        %calculate percentage of trials with a response
        laserRespPercent = length(unique(cols))/length(s.LaserData.LaserTimes);
        %calculate average spiking response to laser
        laserRespSize = length(cols)/length(unique(cols));
        %identify laser waveforms, use these to identify non-laser
        %waveforms
        laserWaveIndex = rows;
        fullIndex = 1:1:length(s.(desigNames{i}).SpikeTimes);
        nonLaserIndex = setdiff(fullIndex,laserWaveIndex);
    else
        laserRespPercent = 0;
        laserRespSize = 0;
        laserWaveIndex = [];
        nonLaserIndex = 1:1:length(s.(desigNames{i}).SpikeTimes);
    end
    %calculates average waves derived from the above indices.
    averageWaveLaser = mean(s.(desigNames{i}).AllWaves(:,laserWaveIndex),2);
    waveSTDLaser = std(s.(desigNames{i}).AllWaves(:,laserWaveIndex),0,2);
    averageWaveNon = mean(s.(desigNames{i}).AllWaves(:,nonLaserIndex),2);
    waveSTDNon = std(s.(desigNames{i}).AllWaves(:,nonLaserIndex),0,2);
    
    %save these values
    s.(desigNames{i}).LaserRelated.PercentLaserResponse = laserRespPercent;
    s.(desigNames{i}).LaserRelated.SpikesPerLaser = laserRespSize;
    s.(desigNames{i}).LaserRelated.AverageLaserWave = averageWaveLaser;
    s.(desigNames{i}).LaserRelated.LaserWaveSTD = waveSTDLaser;
    s.(desigNames{i}).LaserRelated.AverageNormalWave = averageWaveNon;
    s.(desigNames{i}).LaserRelated.NormalWaveSTD = waveSTDNon;
    
    %move on to next step, aligning actual spikes. This involves generating
    %simple raster as well as histogram.
    %here, the raster will be easier to generate, since I already have the
    %massive subtraction array. 
    
    %find all events in which the values are within the range of the
    %raster window. Only save columns (which indicate laser pulse number
    [~,cols] = find(laserMatrix > s.Params.RasterWindow(1) & laserMatrix < s.Params.RasterWindow(2));
    %extract times of all these points as well. 
    rasterTimes = laserMatrix(laserMatrix > s.Params.RasterWindow(1) & laserMatrix < s.Params.RasterWindow(2));
    %insert pulse number in to save information. 
    rasterTimes(:,2) = cols;
    
    %generate histogram.
    [counts,centers] = hist(rasterTimes(:,1),s.Params.HistBinVector);
    countSize = size(counts);
    centerSize = size(centers);
    if countSize(1)>countSize(2)
        counts = counts';
    end
    if centerSize(1)>centerSize(2)
        centers = centers';
    end
    laserHist = [counts'*(1/s.Params.HistBin)/length(s.LaserData.LaserTimes),centers'];
    
    %calculate a z score (this is dirty, since it includes the laser
    %period)
    zScoreFiring = zscore(laserHist(:,1));
    acceptablePeriod = [find(s.Params.HistBinVector>0,1,'first') find(s.Params.HistBinVector>s.Params.PairingCutoff/1000,1,'first')]; %figures out the bins from laser onset to end of acceptable laser period (ex. 15 ms)
    firstZCrossing = laserHist(find(zScoreFiring(acceptablePeriod(1):acceptablePeriod(2))>s.Params.zScoreCutoff,1,'first')+acceptablePeriod(1)-1,2); %calculates the first time bin at which the firing rate crosses threshold.
    %note the -1 fudge factor is because I am adding the index of
    %acceptablePeriod(1), which produces an indexing error
    
    %store these values as well.
    s.(desigNames{i}).LaserRelated.LaserRaster = rasterTimes;
    s.(desigNames{i}).LaserRelated.LaserHistogram = laserHist;
    s.(desigNames{i}).LaserRelated.zScoreHistogram = zScoreFiring;
    s.(desigNames{i}).LaserRelated.FirstZCrossing = firstZCrossing;
end
%% Plot things out

% for i = 1:numUnits
%     hFig = figure;
%     set(hFig,'Position',[40 80 600 1000])
%     %plots average waveform
%     subplot(4,2,1)
%     hold on
%     plot(s.(desigNames{i}).LaserRelated.AverageNormalWave,'k','LineWidth',2)
%     plot(s.(desigNames{i}).LaserRelated.AverageLaserWave,'c','LineWidth',2)
%     waveCorrelation = corrcoef(s.(desigNames{i}).LaserRelated.AverageLaserWave,s.(desigNames{i}).LaserRelated.AverageNormalWave);
% 
%     title({strcat('OverallFiringRate:',num2str(s.(desigNames{i}).OverallFiringRate)),strcat('WaveCorrelation:',num2str(waveCorrelation(2))),...
%         strcat('LaserResponsePercentage:',num2str(s.(desigNames{i}).LaserRelated.PercentLaserResponse))})
%     %plots ISI
%     subplot(4,2,2)
%     hist(s.(desigNames{i}).ISIGraph,1000)
%     histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
%     line([s.Params.RPVTime s.Params.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
%     xlim([s.Params.ClusterWindow(1) s.Params.ClusterWindow(2)])
%     title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
%         strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
%     %plots rasters
%     subplot(4,1,2)
%     plot(s.(desigNames{i}).LaserRelated.LaserRaster(:,1),...
%         s.(desigNames{i}).LaserRelated.LaserRaster(:,2),'k.')
%     line([0 0],[0 size(s.LaserData.LaserTimes,1)],'LineWidth',1,'Color','blue')
%     line([s.LaserData.LaserDuration/1000 s.LaserData.LaserDuration/1000],[0 size(s.LaserData.LaserTimes,1)],'LineWidth',1,'Color','blue')
%     line([s.Params.PairingCutoff/1000 s.Params.PairingCutoff/1000],[0 size(s.LaserData.LaserTimes,1)],'LineWidth',2,'Color','red')
%     ylim([0 size(s.LaserData.LaserTimes,1)]);
%     xlim([s.Params.RasterWindow(1) s.Params.RasterWindow(2)]);
%     h = title(strcat(fileName,desigNames{i}));
%     set(h,'interpreter','none') 
%     %plots histogram
%     subplot(4,1,3)
%     plot(s.(desigNames{i}).LaserRelated.LaserHistogram(:,2),...
%         s.(desigNames{i}).LaserRelated.LaserHistogram(:,1),'k')
%     line([0 0],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',1,'Color','blue')
%     line([s.LaserData.LaserDuration/1000 s.LaserData.LaserDuration/1000],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',1,'Color','blue')
%     line([s.Params.PairingCutoff/1000 s.Params.PairingCutoff/1000],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',2,'Color','red')
%     title(strcat('Mean Laser Dur:',num2str(s.LaserData.LaserDuration),'ms'))
%     %plots zScore
%     subplot(4,1,4)
%     plot(s.(desigNames{i}).LaserRelated.LaserHistogram(:,2),...
%         s.(desigNames{i}).LaserRelated.zScoreHistogram,'k')
%     line([0 0],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','blue')
%     line([s.LaserData.LaserDuration/1000 s.LaserData.LaserDuration/1000],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','blue')
%     line([s.Params.PairingCutoff/1000 s.Params.PairingCutoff/1000],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',2,'Color','red')
%     %if there is a threshold crossing, plots it
%     if ~isempty(s.(desigNames{i}).LaserRelated.FirstZCrossing)
%         line([s.(desigNames{i}).LaserRelated.FirstZCrossing s.(desigNames{i}).LaserRelated.FirstZCrossing]...
%             ,[0 max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','green')
%         title(strcat('CELL IS LASER RESPONSIVE WITH ',num2str(s.(desigNames{i}).LaserRelated.FirstZCrossing*1000),' ms DELAY'))
%     end
%     if min(s.(desigNames{i}).LaserRelated.zScoreHistogram) ~= max(s.(desigNames{i}).LaserRelated.zScoreHistogram)
%         ylim([min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)])
%     end
%     %save as figure and PDF
%     spikeGraphName = strcat(fileName,desigNames{i},'LaserResponse');
%     savefig(hFig,spikeGraphName);
%     set(hFig,'Units','Inches');
%     pos = get(hFig,'Position');
%     set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(hFig,spikeGraphName,'-dpdf','-r0')
% end

% clearvars -except s fname pname
%saves s
pname = pwd;
fname = strcat(fileName,'LaserIDAnalysis');
save(fullfile(pname,fname),'s');

end















