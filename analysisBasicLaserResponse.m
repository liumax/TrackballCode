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

function [s] = analysisBasicLaserResponse(fileName);

%set parameters.
sampleRate = 30000;
rasterWindow = [-0.1,0.1];
clusterWindow = [-0.01,0.05];
rpvTime = 0.001; %time limit in seconds for consideration as an RPV
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
pairingCutoff = 15; %ms cutoff. Significant responses past this point will not be considered. 
pairingEarlyCutoff = 2; %minimum number of ms after a laser before a spike is considered laser related
zScoreCutoff = 4; %zscore above which cell is considered IDed

histBin = 0.001; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
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
%generate placeholder structure
s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(rpvTime,...
    matclustFiles,s,clusterWindow);
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
inTimes = dioTime(dioState == 1)/sampleRate;


%pulls the duration of the pulse (this should give laser pulse duration
dioDiff = diff(dioTime);
dioDiff = dioDiff(dioDiff<5000);
meanDiff = mean(dioDiff)/30;

s.LaserData.LaserTimes = inTimes;
s.LaserData.LaserDuration = meanDiff;
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
    laserMatrix = zeros(length(s.(desigNames{i}).SpikeTimes),length(inTimes));
    %performs subtraction by laser pulse
    for j = 1:length(inTimes)
        laserMatrix(:,j) = s.(desigNames{i}).SpikeTimes - inTimes(j);
    end
    %determines which waveforms are laser related. 
    [rows,cols] = find(laserMatrix > pairingEarlyCutoff/1000 & laserMatrix< pairingCutoff/1000);
    %safeguard against no response
    if ~isempty(rows)
        %calculate percentage of trials with a response
        laserRespPercent = length(unique(cols))/length(inTimes);
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
    averageWaveLaser = mean(s.(desigNames{i}).AllWaves(:,:,laserWaveIndex),3);
    averageWaveNon = mean(s.(desigNames{i}).AllWaves(:,:,nonLaserIndex),3);
    
    %save these values
    s.(desigNames{i}).LaserRelated.PercentLaserResponse = laserRespPercent;
    s.(desigNames{i}).LaserRelated.SpikesPerLaser = laserRespSize;
    s.(desigNames{i}).LaserRelated.AverageLaserWave = averageWaveLaser;
    s.(desigNames{i}).LaserRelated.AverageNormalWave = averageWaveNon;
    
    %move on to next step, aligning actual spikes. This involves generating
    %simple raster as well as histogram.
    %here, the raster will be easier to generate, since I already have the
    %massive subtraction array. 
    
    %find all events in which the values are within the range of the
    %raster window. Only save columns (which indicate laser pulse number
    [~,cols] = find(laserMatrix > rasterWindow(1) & laserMatrix < rasterWindow(2));
    %extract times of all these points as well. 
    rasterTimes = laserMatrix(laserMatrix > rasterWindow(1) & laserMatrix < rasterWindow(2));
    %insert pulse number in to save information. 
    rasterTimes(:,2) = cols;
    
    %generate histogram.
    [counts,centers] = hist(rasterTimes(:,1),histBinVector);
    countSize = size(counts);
    centerSize = size(centers);
    if countSize(1)>countSize(2)
        counts = counts';
    end
    if centerSize(1)>centerSize(2)
        centers = centers';
    end
    laserHist = [counts'*(1/histBin)/length(inTimes),centers'];
    
    %calculate a z score (this is dirty, since it includes the laser
    %period)
    zScoreFiring = zscore(laserHist(:,1));
    acceptablePeriod = [find(histBinVector>0,1,'first') find(histBinVector>pairingCutoff/1000,1,'first')]; %figures out the bins from laser onset to end of acceptable laser period (ex. 15 ms)
    firstZCrossing = laserHist(find(zScoreFiring(acceptablePeriod(1):acceptablePeriod(2))>zScoreCutoff,1,'first')+acceptablePeriod(1)-1,2); %calculates the first time bin at which the firing rate crosses threshold.
    %note the -1 fudge factor is because I am adding the index of
    %acceptablePeriod(1), which produces an indexing error
    
    %store these values as well.
    s.(desigNames{i}).LaserRelated.LaserRaster = rasterTimes;
    s.(desigNames{i}).LaserRelated.LaserHistogram = laserHist;
    s.(desigNames{i}).LaserRelated.zScoreHistogram = zScoreFiring;
    s.(desigNames{i}).LaserRelated.FirstZCrossing = firstZCrossing;
end
%% Plot things out

for i = 1:numUnits
    hFig = figure;
    set(hFig,'Position',[40 80 600 1000])
    %plots average waveform
    subplot(4,2,1)
    hold on
    plot(s.(desigNames{i}).LaserRelated.AverageNormalWave,'k','LineWidth',2)
    plot(s.(desigNames{i}).LaserRelated.AverageLaserWave,'c','LineWidth',2)
    %turns out corrcoef can function across multiple columns! Produces a
    %maybe inflated value?
    
    waveCorrelation = corrcoef(s.(desigNames{i}).LaserRelated.AverageLaserWave,s.(desigNames{i}).LaserRelated.AverageNormalWave);
% 
    title({strcat('OverallFiringRate:',num2str(s.(desigNames{i}).OverallFiringRate)),strcat('WaveCorrelation:',num2str(waveCorrelation(2))),...
        strcat('LaserResponsePercentage:',num2str(s.(desigNames{i}).LaserRelated.PercentLaserResponse))})
    %plots ISI
    subplot(4,2,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim([clusterWindow(1) clusterWindow(2)])
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    %plots rasters
    subplot(4,1,2)
    plot(s.(desigNames{i}).LaserRelated.LaserRaster(:,1),...
        s.(desigNames{i}).LaserRelated.LaserRaster(:,2),'k.')
    line([0 0],[0 size(inTimes,1)],'LineWidth',1,'Color','blue')
    line([meanDiff/1000 meanDiff/1000],[0 size(inTimes,1)],'LineWidth',1,'Color','blue')
    line([pairingCutoff/1000 pairingCutoff/1000],[0 size(inTimes,1)],'LineWidth',2,'Color','red')
    ylim([0 size(inTimes,1)]);
    xlim([rasterWindow(1) rasterWindow(2)]);
    h = title(strcat(fileName,desigNames{i}));
    set(h,'interpreter','none') 
    %plots histogram
    subplot(4,1,3)
    plot(s.(desigNames{i}).LaserRelated.LaserHistogram(:,2),...
        s.(desigNames{i}).LaserRelated.LaserHistogram(:,1),'k')
    line([0 0],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',1,'Color','blue')
    line([meanDiff/1000 meanDiff/1000],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',1,'Color','blue')
    line([pairingCutoff/1000 pairingCutoff/1000],[0 max(s.(desigNames{i}).LaserRelated.LaserHistogram(:,1))],'LineWidth',2,'Color','red')
    title(strcat('Mean Laser Dur:',num2str(meanDiff),'ms'))
    %plots zScore
    subplot(4,1,4)
    plot(s.(desigNames{i}).LaserRelated.LaserHistogram(:,2),...
        s.(desigNames{i}).LaserRelated.zScoreHistogram,'k')
    line([0 0],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','blue')
    line([meanDiff/1000 meanDiff/1000],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','blue')
    line([pairingCutoff/1000 pairingCutoff/1000],[min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',2,'Color','red')
    %if there is a threshold crossing, plots it
    if ~isempty(s.(desigNames{i}).LaserRelated.FirstZCrossing)
        line([s.(desigNames{i}).LaserRelated.FirstZCrossing s.(desigNames{i}).LaserRelated.FirstZCrossing]...
            ,[0 max(s.(desigNames{i}).LaserRelated.zScoreHistogram)],'LineWidth',1,'Color','green')
        title(strcat('CELL IS LASER RESPONSIVE WITH ',num2str(s.(desigNames{i}).LaserRelated.FirstZCrossing*1000),' ms DELAY'))
    end
    if min(s.(desigNames{i}).LaserRelated.zScoreHistogram) ~= max(s.(desigNames{i}).LaserRelated.zScoreHistogram)
        ylim([min(s.(desigNames{i}).LaserRelated.zScoreHistogram) max(s.(desigNames{i}).LaserRelated.zScoreHistogram)])
    end
    %save as figure and PDF
    spikeGraphName = strcat(fileName,desigNames{i},'LaserResponse');
    savefig(hFig,spikeGraphName);
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end

% clearvars -except s fname pname
%saves s
pname = pwd;
fname = strcat(fileName,'LaserIDAnalysis');
save(fullfile(pname,fname),'s');

end















