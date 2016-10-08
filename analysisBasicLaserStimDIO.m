
function [matclustStruct] = analysisBasicLaserStimDIO(fileName);
%This function is meant to perform the basic analysis of spike responses in
%response to laser stimulation for ID purposes. This will produce an
%average waveform, a graph of ISIs, and a raster and histogram of IDed
%responses. fileName should be just the name of the file without
%extensions.


rasterWindow = [-0.1,0.1];
clusterWindow = [-0.01,0.05];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
pairingCutoff = 15; %ms cutoff beyond which cell is not considered to be IDed 
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
%%
%finds matclust folder, extracts matclust file names
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');

%removes periods which allow structured array formation.
truncatedNames = matclustFiles;
numTrodes = length(truncatedNames);
for i = 1:length(truncatedNames)
    truncatedNames{i} = truncatedNames{i}(16:find(truncatedNames{i} == '.')-1);
end

%generates structured array for storage of data
matclustStruct = struct;
for i = 1:length(truncatedNames);
    matclustStruct.(truncatedNames{i}) = [];
end

%%
%finds DIO folder, extracts D2 specifically
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};

%extracts DIO stuffs!
[DIOData] = readTrodesExtractedDataFile(D2FileName);
%extracts port states and times of changes
dioState = double(DIOData.fields(2).data);
dioTime = double(DIOData.fields(1).data);

%pulls the times when state goes up!
inTimes = dioTime(dioState == 1)/30000;
master(:,1) = inTimes;

%pulls the duration of the pulse (this should give laser pulse duration
dioDiff = diff(dioTime);
dioDiff = dioDiff(dioDiff<5000);
meanDiff = mean(dioDiff)/30;
%%


matclustStruct.SoundTimes = master(:,1);

%%
%extracts times of clustered spikes.
for i = 1:numTrodes
    matclustFile = open(matclustFiles{i});
    %this extracts indexes for cluster components
    clusterSizer= size(matclustFile.clustattrib.clustersOn,1);
    matclustStruct.(truncatedNames{i}).ClusterNumber = clusterSizer;
    %this will find which clusters actually have info

    clusterIndex = cell(clusterSizer,1);

    for j = 1:clusterSizer
        clusterIndex{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
    end
    %replaces indices with real times
    for j = 1:clusterSizer
        clusterIndex{j} = matclustFile.clustdata.params(clusterIndex{j},1);
    end
    
    %puts clusterIndex into structured array
    matclustStruct.(truncatedNames{i}).SpikeIndices = clusterIndex;
    
    %preps series of arrays
    clusterSpikes = cell(clusterSizer,1);
    clusterWaves = cell(clusterSizer,1);
    clearedSpikes = cell(clusterSizer,1);
    %calculates inter-spike interval, saves this information
    for j = 1:clusterSizer
        %this line pulls the actual indices of spikes for the cluster
        clusterSpikes{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
        %this subtracts all adjacent spikes
        diffSpikes = diff(matclustFile.clustdata.params(clusterSpikes{j},1));
        %this then windows that out, to only plot the small ISIs
        clearedSpikes{j} = diffSpikes(diffSpikes<clusterWindow(2)); %removes long pauses
        diffSpikes = [];
    end
    matclustStruct.(truncatedNames{i}).ISIData = clearedSpikes;
    clearedSpikes = [];
    
    %prepares cluster indices. Finds actual points of index per cluster.
    %Then replaces with real times. 
    clusterIndex = cell(clusterSizer,1);
    for j = 1:clusterSizer
        clusterIndex{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
        clusterIndex{j} = matclustFile.clustdata.params(clusterIndex{j},1);
    end
    %puts clusterIndex into structured array
    matclustStruct.(truncatedNames{i}).SpikeTimes = clusterIndex;
    
    %calculates average firing rate
    totalTime = matclustFile.clustdata.datarange(2,1)-matclustFile.clustdata.datarange(1,1);
    averageFiring = zeros(clusterSizer,1);
    for j = 1:clusterSizer
        averageFiring(j) = size(clusterIndex{j},1)/totalTime;
    end
    matclustStruct.(truncatedNames{i}).AverageFiringRate = averageFiring;
    
    
    targetWaveName = strcat('waves_',truncatedNames{i},'.mat');
    waveLoader = open(targetWaveName);
    waveLoader =squeeze(waveLoader.waves);
    waveHolder = cell(clusterSizer,1);
    averageWaveHolder = zeros(size(waveLoader,1),clusterSizer,3);
    for j = 1:clusterSizer
        waveHolder{j} = waveLoader(:,clusterSpikes{j});
        averageWaveHolder(:,j,2) = mean(waveHolder{j},2);
        averageWaveHolder(:,j,1) = mean(waveHolder{j},2)-std(waveHolder{j},0,2)/sqrt(size(waveHolder{j},2));
        averageWaveHolder(:,j,3) = mean(waveHolder{j},2)+std(waveHolder{j},0,2)/sqrt(size(waveHolder{j},2));
    end
    matclustStruct.(truncatedNames{i}).AverageWaveForms = averageWaveHolder;
    matclustStruct.(truncatedNames{i}).AllWaveForms = waveHolder;
    
   
    %find laser evoked and not laser evoked spikes
    laserWaveIndices = cell(clusterSizer,2); %will hold indices for laser vs non-laser waveforms
    averageWaveHolderLaser = zeros(size(waveLoader,1),clusterSizer,3);
    averageWaveHolderNonLaser = zeros(size(waveLoader,1),clusterSizer,3);
    laserResponseCounter = zeros(clusterSizer,1);
    for j = 1:clusterSizer
        %first, have to find the indices of all spikes that come within the
        %acceptable period. We will consider all of these to be light
        %evoked.
        laserSpikeHolder = zeros(size(waveHolder{j},2),1);
        laserSpikeCounter = 1;
        for k = 1:size(inTimes,1)
            laserSpikeSubtractor = clusterIndex{j} - inTimes(k);
            laserSpikeSubtractor = find(laserSpikeSubtractor > (pairingEarlyCutoff/1000) & laserSpikeSubtractor < pairingCutoff/1000);
            if ~isempty(laserSpikeSubtractor) %safeguard against trials with no response
                laserResponseCounter(j) = laserResponseCounter(j) + 1/size(inTimes,1);
                laserSpikeHolder(laserSpikeCounter:size(laserSpikeSubtractor,1)+laserSpikeCounter-1) = laserSpikeSubtractor; %populates vector with correct values from subtraction
                laserSpikeCounter = laserSpikeCounter + size(laserSpikeSubtractor,1); %advances counter.
            end
        end
        laserSpikeHolder(laserSpikeHolder == 0) = []; %eliminates extra zeros
        laserWaveIndices{j,1} = laserSpikeHolder; %stores laser waveform indices in the first column
        nonLaserSpikeHolder = 1:1:size(clusterIndex{j},1);
        nonLaserSpikeHolder(laserSpikeHolder) = []; %stores all other waveforms
        laserWaveIndices{j,2} = nonLaserSpikeHolder;
        %computes average laser waveforms.
        averageWaveHolderLaser(:,j,2) = mean(waveHolder{j}(:,laserSpikeHolder),2);
        averageWaveHolderLaser(:,j,1) = mean(waveHolder{j}(:,laserSpikeHolder),2)-std(waveHolder{j}(:,laserSpikeHolder),0,2)/sqrt(size(waveHolder{j}(:,laserSpikeHolder),2));
        averageWaveHolderLaser(:,j,3) = mean(waveHolder{j}(:,laserSpikeHolder),2)+std(waveHolder{j}(:,laserSpikeHolder),0,2)/sqrt(size(waveHolder{j}(:,laserSpikeHolder),2));
        %computers average non-laser waveform.
        averageWaveHolderNonLaser(:,j,2) = mean(waveHolder{j}(:,nonLaserSpikeHolder),2);
        averageWaveHolderNonLaser(:,j,1) = mean(waveHolder{j}(:,nonLaserSpikeHolder),2)-std(waveHolder{j}(:,nonLaserSpikeHolder),0,2)/sqrt(size(waveHolder{j}(:,nonLaserSpikeHolder),2));
        averageWaveHolderNonLaser(:,j,3) = mean(waveHolder{j}(:,nonLaserSpikeHolder),2)+std(waveHolder{j}(:,nonLaserSpikeHolder),0,2)/sqrt(size(waveHolder{j}(:,nonLaserSpikeHolder),2));
    end
    
    matclustStruct.(truncatedNames{i}).LaserSpikeIndices = laserWaveIndices;
    matclustStruct.(truncatedNames{i}).AverageWaveFormLaser = averageWaveHolderLaser;
    matclustStruct.(truncatedNames{i}).AverageWaveFormNonLaser = averageWaveHolderNonLaser;
    matclustStruct.(truncatedNames{i}).LaserResponseProbability = laserResponseCounter;
    
    masterToneRaster = [];
    masterToneHist = [];
    zScoreFiring = cell(clusterSizer,1);
    firstZCrossing = cell(clusterSizer,1);
    indivToneHist = cell(clusterSizer,length(inTimes));
    %generate rasters and histograms, ignores frequency information
    for j = 1:clusterSizer
        rasterHolder = 1;
        for k = 1:length(inTimes)
            toneHolder = clusterIndex{j}(clusterIndex{j}>master(k,1)+rasterWindow(1) & clusterIndex{j}<master(k,1)+rasterWindow(2));
            toneHolder = toneHolder - master(k,1);
            %This will generate a histogram based on each individual trace,
            %which can be used to generate standard deviations.
            [counts,centers] = hist(toneHolder,histBinVector);
            %below code necessary to prevent bugs with row vectors
            countSize = size(counts);
            centerSize = size(centers);
            if countSize(1)>countSize(2)
                counts = counts';
            end
            if centerSize(1)>centerSize(2)
                centers = centers';
            end
            indivToneHist{j,k} = [counts'*(1/histBin),centers'];
            counts = [];
            centers = [];
            %fills in large raster plot. holder updated position.
            toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,1) = k;
            toneRaster(rasterHolder:rasterHolder + length(toneHolder)-1,2) = toneHolder;
            rasterHolder = rasterHolder + length(toneHolder);
            toneRaster(toneRaster(:,1) == 0,:) = [];
            toneHolder = [];
        end
        masterToneRaster{j} = toneRaster;
        toneRaster = zeros(100000,2);
        [counts,centers] = hist(masterToneRaster{j}(:,2),histBinVector);
        countSize = size(counts);
        centerSize = size(centers);
        if countSize(1)>countSize(2)
            counts = counts';
        end
        if centerSize(1)>centerSize(2)
            centers = centers';
        end
        masterToneHist{j} = [counts'*(1/histBin)/length(master(:,1)),centers'];
        counts = [];
        centers = [];
        zScoreFiring{j} = zscore(masterToneHist{j}(:,1)); %calculates the z score. This is technically dirty since it includes the laser period.
        acceptablePeriod = [find(histBinVector>0,1,'first') find(histBinVector>pairingCutoff/1000,1,'first')]; %figures out the bins from laser onset to end of acceptable laser period (ex. 15 ms)
        firstZCrossing{j} = masterToneHist{j}(find(zScoreFiring{j}(acceptablePeriod(1):acceptablePeriod(2))>zScoreCutoff,1,'first')+acceptablePeriod(1)-1,2); %calculates the first time bin at which the firing rate crosses threshold.
        %note the -1 fudge factor is because I am adding the index of
        %acceptablePeriod(1), which produces an indexing error
    end
    
    matclustStruct.(truncatedNames{i}).Rasters = masterToneRaster;
    matclustStruct.(truncatedNames{i}).Histogram = masterToneHist;
    matclustStruct.(truncatedNames{i}).zScore = zScoreFiring;
    matclustStruct.(truncatedNames{i}).zScoreCrossing = firstZCrossing;

    
    stdHolder = zeros(length(histBinVector),length(master(:,1)));
    steHolder = zeros(length(histBinVector),clusterSizer);

    for j = 1:clusterSizer
        for k = 1:length(master(:,1))
            stdHolder(:,k) = indivToneHist{j,k}(:,1);
        end
        steHolder(:,j) = std(stdHolder,0,2)/sqrt(length(master(:,1)));
    end

    %generates lines representing standard error.
    stePlotter = zeros(length(histBinVector),clusterSizer,2);
    for j = 1:clusterSizer
        stePlotter(:,j,1) = masterToneHist{j}(:,1)-steHolder(:,j);
        stePlotter(:,j,2) = masterToneHist{j}(:,1)+steHolder(:,j);
    end
    matclustStruct.(truncatedNames{i}).StandardError = steHolder;
    matclustStruct.(truncatedNames{i}).StandardErrorPlotting = stePlotter;
    
    masterToneRaster = [];
    masterToneHist = [];
    steHolder = [];
    stePlotter = [];
    indivToneHist = [];
end

%%Graphing!

for i = 1:numTrodes
    for j = 1:matclustStruct.(truncatedNames{i}).ClusterNumber
        hFig = figure;
        set(hFig,'Position',[40 80 600 1000])
        %plots average waveform
        subplot(4,2,1)
        hold on
        plot(matclustStruct.(truncatedNames{i}).AverageWaveFormNonLaser(:,j,2),'k','LineWidth',2)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveFormNonLaser(:,j,1),'r','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveFormNonLaser(:,j,3),'r','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveFormLaser(:,j,2),'b','LineWidth',2)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveFormLaser(:,j,1),'c','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveFormLaser(:,j,3),'c','LineWidth',1)
        waveCorrelation = corrcoef(matclustStruct.(truncatedNames{i}).AverageWaveFormLaser(:,j,2),matclustStruct.(truncatedNames{i}).AverageWaveFormNonLaser(:,j,2));

        title({strcat('AverageFiringRate:',num2str(matclustStruct.(truncatedNames{i}).AverageFiringRate(j))),strcat('WaveCorrelation:',num2str(waveCorrelation(2))),...
            strcat('LaserResponsePercentage:',num2str(matclustStruct.(truncatedNames{i}).LaserResponseProbability(j)))})
        %plots ISI
        subplot(4,2,2)
        hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000)
        xlim([clusterWindow(1) clusterWindow(2)])
        title('ISI')
        %plots rasters
        subplot(4,1,2)
        plot(matclustStruct.(truncatedNames{i}).Rasters{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Rasters{j}(:,1),'k.')
        line([0 0],[0 size(inTimes,1)],'LineWidth',1,'Color','blue')
        line([meanDiff/1000 meanDiff/1000],[0 size(inTimes,1)],'LineWidth',1,'Color','blue')
        line([pairingCutoff/1000 pairingCutoff/1000],[0 size(inTimes,1)],'LineWidth',2,'Color','red')
        ylim([0 size(inTimes,1)]);
        xlim([rasterWindow(1) rasterWindow(2)]);
        h = title(strcat(fileName,truncatedNames{i},' Cluster ',num2str(j)));
        set(h,'interpreter','none') 
        %plots histogram
        subplot(4,1,3)
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).Histogram{j}(:,1),'k')
        line([0 0],[0 max(matclustStruct.(truncatedNames{i}).Histogram{j}(:,1))],'LineWidth',1,'Color','blue')
        line([meanDiff/1000 meanDiff/1000],[0 max(matclustStruct.(truncatedNames{i}).Histogram{j}(:,1))],'LineWidth',1,'Color','blue')
        line([pairingCutoff/1000 pairingCutoff/1000],[0 max(matclustStruct.(truncatedNames{i}).Histogram{j}(:,1))],'LineWidth',2,'Color','red')
        title(strcat('Mean Laser Dur:',num2str(meanDiff),'ms'))
        %plots zScore
        subplot(4,1,4)
        plot(matclustStruct.(truncatedNames{i}).Histogram{j}(:,2),...
            matclustStruct.(truncatedNames{i}).zScore{j},'k')
        line([0 0],[min(matclustStruct.(truncatedNames{i}).zScore{j}) max(matclustStruct.(truncatedNames{i}).zScore{j})],'LineWidth',1,'Color','blue')
        line([meanDiff/1000 meanDiff/1000],[min(matclustStruct.(truncatedNames{i}).zScore{j}) max(matclustStruct.(truncatedNames{i}).zScore{j})],'LineWidth',1,'Color','blue')
        line([pairingCutoff/1000 pairingCutoff/1000],[min(matclustStruct.(truncatedNames{i}).zScore{j}) max(matclustStruct.(truncatedNames{i}).zScore{j})],'LineWidth',2,'Color','red')
        %if there is a threshold crossing, plots it
        if ~isempty(matclustStruct.(truncatedNames{i}).zScoreCrossing{j})
            line([matclustStruct.(truncatedNames{i}).zScoreCrossing{j} matclustStruct.(truncatedNames{i}).zScoreCrossing{j}]...
                ,[0 max(matclustStruct.(truncatedNames{i}).zScore{j})],'LineWidth',1,'Color','green')
            title(strcat('CELL IS LASER RESPONSIVE WITH ',num2str(matclustStruct.(truncatedNames{i}).zScoreCrossing{j}*1000),' ms DELAY'))
        end
        if min(matclustStruct.(truncatedNames{i}).zScore{j}) ~= max(matclustStruct.(truncatedNames{i}).zScore{j})
            ylim([min(matclustStruct.(truncatedNames{i}).zScore{j}) max(matclustStruct.(truncatedNames{i}).zScore{j})])
        end
        %save as figure and PDF
        spikeGraphName = strcat(fileName,truncatedNames{i},' Cluster ',num2str(j),'LaserResponse');
        savefig(hFig,spikeGraphName);
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end



% clearvars -except matclustStruct fname pname
%saves matclustStruct
pname = pwd;
fname = strcat(fileName,'LaserIDAnalysis');
save(fullfile(pname,fname),'matclustStruct');

end















