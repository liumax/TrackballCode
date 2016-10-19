%This is meant to be a function built code that analyzes the experiment in
%which the animal is delivered pulsed sounds. 

%This needs the following in the same folder: matclust file of picked
%spikes and matlab file with audio order

function [] = analysisPulsedSoundFunctionsLaser(fileName);


%% constants
sampleRate = 30000;

%% Variables
rpvTime = 0.001; %time limit in seconds for consideration as an RPV
clusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info

defaultBins = 0.001;% bin size for calculating significant responses
smoothingBins = [0.01 0.001];%bins for smoothing
calcWindow = [0 2]; %window for calculating significant responses
zLimit = 3; %zlimit for calculating significant responses
% firstSpikeWindow = [0 0.5 1 1.5]; %ratios! need to be multiplied by tone duration.
firstSpikeWindow = [0 1];
chosenSpikeBin = 1; %delineates which spike window I will graph.
baselineBin = [-1 0]; %ratio for bin from which baseline firing rate will be calculated
rasterWindow = [-1 3]; %ratio for raster window. will be multiplied by toneDur
histBin = 0.005; %histogram bin size in seconds

bigRasterWindow = [-0.5 1.5]; %ratio for big raster window, which will be the raster across the entire pip set.
bigHistBin = 0.05; %hist bin size in seconds for overall histograms
%% sets up file saving stuff
saveName = strcat(fileName,'PulsedSoundAnalysis','.mat');
fname = saveName;
pname = pwd;

%% Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%% Find and ID Matclust Files for Subsequent Analysis. Generates Structured Array for Data Storage
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');
%generate placeholder structure
s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(rpvTime,...
    matclustFiles,s,clusterWindow);

%% Extracts Sound Data from soundFile
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData;
%extract the data for the sounds
pipITI = soundData.PipITI;
pipReps = soundData.PipRepetitions;
pipPos = 1:1:pipReps;
blockNum = soundData.BlockRepetitions;
laserITI = soundData.LaserTTLITI;
laserSwitch = soundData.LaserOn;
toneDur = soundData.ToneDuration;
blockDur = pipITI*pipReps;
controlArray = soundData.ControlMatrix;

if laserSwitch == 1
    expectedTrialNum = blockNum * length(pipITI)*2;
else
    expectedTrialNum = blockNum * length(pipITI);
end

%set up raster windows
rasterWindow = rasterWindow*toneDur; %for plotting individual pips
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
bigRasterAxis = cell(length(blockDur),1);
for i = 1:length(blockDur)
    holder(:,i) = bigRasterWindow*blockDur(i);
    bigRasterAxis{i} = [holder(1,i):0.01:holder(2,i)-0.01];
end
bigRasterWindow = holder;

%set up histogram bin vectors
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2];
bigHistVector = cell(length(blockDur),1);
for i = 1:length(blockDur)
    bigHistVector{i} = [bigRasterWindow(1,i)+bigHistBin/2:bigHistBin:bigRasterWindow(2,i)-bigHistBin/2];
end

%set up baseline bin
baselineBin = baselineBin * toneDur;

%% DIO extractions
%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIO1Data] = readTrodesExtractedDataFile(D1FileName);
%extracts only points where going up to high state.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIO1Data,sampleRate);
DIO1True = dioTimes;
DIO1TrueDiff = dioTimeDiff;
dioTimes = [];
dioTimeDiff = [];

%find DIO folder and D2 file for analysis. Resistant to having no inputs
%come in (reports empty arrays)
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIO2Data,sampleRate);
DIO2True = dioTimes;
DIO2TrueDiff = dioTimeDiff;
dioTimes = [];
dioTimeDiff = [];

%%%%% ABOVE THIS POINT, THIS COULD BE LASER OR NON-LASER.

%% analyze the differences in the DIO data. 
%DIO 1 should be direct inputs from the audio card. DIO 2 is laser.

%analyze the DIO 1 data to find separate blocks

laserTiming = soundData.LaserTiming; %pulls the pip the laser is aligned to.
ttlInd = cumsum(controlArray(:,6) + 5); %generates indices of the end of each pip set
startTimesBlock = DIO1True([1;(ttlInd(1:end-1)+1)]); %pulls the starts of each set of pips.
laserTrials = find(controlArray(:,6) == 1); %determines which trials are laser trials based on column 6
laserInd = ttlInd(laserTrials)-(pipReps - laserTiming); %pulls indexes of laser TTLs in DIO1! 
pipInd = [1:length(DIO1True)]'; %maps indices from 1:end of DIO1True
pipInd(laserInd) = []; %removes laser indices
startTimesPip = DIO1True(pipInd); %converts to times
startTimesLaser = DIO1True(laserInd); %converts to times

%DIO2: pull first laser times!
firstLaserDIO2 = find(DIO2TrueDiff > 0.5);
firstLaserDIO2 = [1;firstLaserDIO2+1];
firstLaserDIO2 = DIO2True(firstLaserDIO2);


%% check that everything is kosher

if length(pipInd) + length(laserInd) ~= length(DIO1True)
    error('PROBLEM: Pip and Laser Indices Dont Add Up Correctly')
end

if length(startTimesBlock) ~= expectedTrialNum
    error('Incorrect number of pip blocks!')
end

if length(laserInd) ~= blockNum * length(pipITI)
    error('Incorrect number of laser TTLs calculated')
end

if length(firstLaserDIO2) ~= length(laserInd)
    error('Mismatch between DIO1 and DIO2 measurement of laser')
end

%% generate matrix for all pulses and blocks

allPulseMatrix = zeros(length(startTimesPip),7);
allPulseMatrix(:,1) = startTimesPip;
counter = 1;
for i = 1:length(startTimesBlock)
    allPulseMatrix(counter:counter+pipReps-1,2) = controlArray(i,1); %stores frequency
    allPulseMatrix(counter:counter+pipReps-1,3) = controlArray(i,2); %stores dB
    allPulseMatrix(counter:counter+pipReps-1,4) = i; %stores block number
    allPulseMatrix(counter:counter+pipReps-1,5) = controlArray(i,4);%stores pipITIs
    allPulseMatrix(counter:counter+pipReps-1,6) = controlArray(i,6); %stores laser status
    allPulseMatrix(counter:counter+pipReps-1,7) = 1:1:pipReps; %stores number of pip within block
    counter = counter + pipReps;
end

blockMatrix = zeros(length(startTimesBlock),6);
blockMatrix(:,1) = startTimesBlock; %stores TTL time.
blockMatrix(:,2) = controlArray(:,1);%stores frequency
blockMatrix(:,3) = controlArray(:,2);%stores dB
blockMatrix(:,4) = 1:1:length(startTimesBlock);%stores block number
blockMatrix(:,5) = controlArray(:,4);%stores pipITIs
blockMatrix(:,6) = controlArray(:,6);%stores laser status


%EDITS TO HERE. 

%convert indices to time
toneTTLFinder = DIO1True(toneTTLFinder);
firstTTLDIO1 = DIO1True(firstTTLDIO1);

allPulseMatrix = zeros(pipReps,length(firstTTLDIO1));
pulseCounter = 1;
for i = 1:length(firstTTLDIO1)
    allPulseMatrix(:,i) = toneTTLFinder(pulseCounter:pulseCounter+pipReps-1);
    pulseCounter = pulseCounter + pipReps;
end
allPulseDesig = zeros(pipReps,length(firstTTLDIO1));
for i = 1:size(allPulseMatrix,1)*size(allPulseMatrix,2)
    allPulseDesig(i) = i;
end

master = zeros(size(toneTTLFinder,1),5);
master(:,1) = toneTTLFinder;
if ~ischar(soundData.TargetFrequency)
    master(:,2) = soundData.ControlFrequency;
end
master(:,3) = 1:1:size(toneTTLFinder,1);
master(:,4) = fix((master(:,3)-1)/pipReps)+1;
master(:,5) = rem(master(:,3),pipReps);
master(master(:,5) == 0,5) = pipReps; %fixed zeros from remainder calculation
%split into control and laser trial sets for analysis
controlCounter = 1;
laserCounter = 1;
controlMaster = zeros(size(toneTTLFinder,1)/2,5);
laserMaster = zeros(size(toneTTLFinder,1)/2,5);


for counter1 = 1: size(master,1)
    if ismember(master(counter1,4),controlBlockFinder)
        controlMaster(controlCounter,:) = master(counter1,:);
        controlCounter = controlCounter + 1;
    elseif ismember(master(counter1,4),laserBlockFinder)
        laserMaster(laserCounter,:) = master(counter1,:);
        laserCounter = laserCounter + 1;
    end
end


bigMaster(:,1) = firstTTLDIO1;
bigMaster(:,2) = 0;
if ~ischar(soundData.TargetFrequency)
    bigMaster(laserBlockFinder,2) = soundData.TargetFrequency;
end
if ~ischar(soundData.ControlFrequency)
    bigMaster(controlBlockFinder,2) = soundData.ControlFrequency;
end

%split into control and laser trial sets for analysis
controlCounter = 1;
laserCounter = 1;
controlBigMaster = zeros(size(bigMaster,1)/2,size(bigMaster,2));
laserBigMaster = zeros(size(bigMaster,1)/2,size(bigMaster,2));


for counter1 = 1: size(bigMaster,1)
    if ismember(counter1,controlBlockFinder)
        controlBigMaster(controlCounter,:) = bigMaster(counter1,:);
        controlCounter = controlCounter + 1;
    elseif ismember(counter1,laserBlockFinder)
        laserBigMaster(laserCounter,:) = bigMaster(counter1,:);
        laserCounter = laserCounter + 1;
    end
end


pipSpikeStorage = cell(numTrodes,1);

for i = 1:numTrodes
    %% extract spikes
    [matclustStruct, clusterSizer] = functionSpikeWaveExtraction(rpvTime,...
        i,matclustFiles,matclustStruct,truncatedNames,clusterWindow);
    %% make rasters/histograms based on all tone presentations (small rasters)
    [matclustStruct] = functionRasterHistExtraction(i,clusterSizer,...
    master,baselineBins,matclustStruct,truncatedNames,...
    rasterWindow,histBin,histBinVector);
    
    %% make rasters/histograms of big picture (big rasters)
    [matclustStruct] = functionPulseBigRasterHistExtraction(i,clusterSizer,...
    bigMaster,bigBaselineBins,matclustStruct,truncatedNames,...
    bigRasterWindow,bigHistBin,bigHistBinVector,laserSwitch);

    %% make calculations for rasters per tone, concatenate for graphing.
    
    indivRasterControl = cell(clusterSizer,1);
    indivRasterLaser = cell(clusterSizer,1);
    for rasterCounter =1:clusterSizer
        rasterInfo = matclustStruct.(truncatedNames{i}).Rasters{rasterCounter};
        controlHolder = zeros(round(histBinNum)*pipReps,1);
        laserHolder = zeros(round(histBinNum)*pipReps,1);
        controlCounter = 1;
        laserCounter = 1;
        %do separate out trials
        controlRasters = rasterInfo(ismember(rasterInfo(:,6),controlBlockFinder),:);
        for k = 1:pipReps
            controlHolder(controlCounter:controlCounter + round(histBinNum)-1,1) = hist(controlRasters(controlRasters(:,5) == k,2),histBinVector);
            controlCounter = controlCounter + round(histBinNum);
        end
        laserRasters = rasterInfo(ismember(rasterInfo(:,6),laserBlockFinder),:);
        for k = 1:pipReps
            laserHolder(laserCounter:laserCounter + round(histBinNum)-1,1) = hist(laserRasters(laserRasters(:,5) == k,2),histBinVector);
            laserCounter = laserCounter + round(histBinNum);
        end
        indivRasterControl{rasterCounter} = controlHolder;
        indivRasterLaser{rasterCounter} = laserHolder;
    end
    matclustStruct.(truncatedNames{i}).IndivRastersControl = indivRasterControl;
    matclustStruct.(truncatedNames{i}).IndivRastersLaser = indivRasterLaser;
%     

%make zero points for graphing of individual hists.
zeroVector = zeros(pipReps,1);
zeroCounter = 0;
endVector = zeros(pipReps,1);
toneEndVector = zeros(pipReps,1);
for k = 1:pipReps
    zeroPoint = find(histBinVector > 0,1,'first')-0.5;
    tonePoint = find(histBinVector > toneDur,1,'first')-0.5;
    zeroVector(k) = zeroPoint + zeroCounter*length(histBinVector);
    toneEndVector(k) = tonePoint + zeroCounter*length(histBinVector);
    endVector(k) = (zeroCounter + 1)*length(histBinVector);
    zeroCounter = zeroCounter + 1;
end
    
    %% calculate first spike timing stuffs
%     for j = 1:clusterSizer
%         rasterData = matclustStruct.(truncatedNames{i}).Rasters{1,j};
%         [s] = functionFirstSpikeTiming(totalTrialNum,rasterData,uniqueFreqs,uniqueDBs,toneDur,toneReps);
%         matclustStruct.(truncatedNames{i}).Stats{j} = s;
%     end
% by design, the big raster function will detect average firing rate and
% should overwrite that of the functionRasterHistExtraction. 
    %% do simple binning procedure to look at changes to responses
    for rasterCounter = 1:clusterSizer
        pipSpikeMatrix = zeros(pipReps,length(firstTTLDIO1));
        for k = 1:size(allPulseMatrix,1)*size(allPulseMatrix,2)
            pipSpikeMatrix(k) = size(find(matclustStruct.(truncatedNames{i}).Rasters{rasterCounter}(:,1) == k &...
                matclustStruct.(truncatedNames{i}).Rasters{rasterCounter}(:,2) > 0 &...
                matclustStruct.(truncatedNames{i}).Rasters{rasterCounter}(:,2) < toneDur),1);
        end
        pipSpikeStorage{i,rasterCounter} = pipSpikeMatrix;
    end
end

%% Plotting
for i = 1:numTrodes
    for rasterCounter = 1:matclustStruct.(truncatedNames{i}).Clusters
        hFig = figure;
        set(hFig, 'Position', [10 10 1280 1000])
        
        %plots average waveform
        subplot(4,6,1)
        hold on
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,rasterCounter,2),'LineWidth',2)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,rasterCounter,1),'r','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,rasterCounter,3),'r','LineWidth',1)
        title(strcat('AverageFiringRate:',num2str(matclustStruct.(truncatedNames{i}).AverageFiringRate(rasterCounter))))
        
        %plots ISI
        subplot(4,6,2)
        hist(matclustStruct.(truncatedNames{i}).ISIData{rasterCounter},1000)
        histMax = max(hist(matclustStruct.(truncatedNames{i}).ISIData{rasterCounter},1000));
        line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
        xlim(clusterWindow)
        title({strcat('ISI RPV %: ',num2str(matclustStruct.(truncatedNames{i}).RPVs(rasterCounter)));...
            strcat(num2str(matclustStruct.(truncatedNames{i}).RPVNumber(rasterCounter)),'/',num2str(matclustStruct.(truncatedNames{i}).TotalSpikeNumber(rasterCounter)))})
        
        %plots histograms

        subplot(4,1,2)
        plot(matclustStruct.(truncatedNames{i}).BigHistogram{rasterCounter}{1}(:,2),...
            matclustStruct.(truncatedNames{i}).BigHistogram{rasterCounter}{1}(:,1),'k','LineWidth',2)
        hold on
        plot(matclustStruct.(truncatedNames{i}).BigHistogram{rasterCounter}{1}(:,2),...
            matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,rasterCounter,1,1),'k')
        plot(matclustStruct.(truncatedNames{i}).BigHistogram{rasterCounter}{1}(:,2),...
            matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,rasterCounter,2,1),'k')

        plot(matclustStruct.(truncatedNames{i}).BigHistogram{rasterCounter}{2}(:,2),...
            matclustStruct.(truncatedNames{i}).BigHistogram{rasterCounter}{2}(:,1),'r','LineWidth',2)
        hold on
        plot(matclustStruct.(truncatedNames{i}).BigHistogram{rasterCounter}{2}(:,2),...
            matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,rasterCounter,1,2),'r')
        plot(matclustStruct.(truncatedNames{i}).BigHistogram{rasterCounter}{2}(:,2),...
            matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,rasterCounter,2,2),'r')

        plot([0 0],[ylim],'r');
        plot([blockDur blockDur],[ylim],'r');
        for k = 1:pipReps
            plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
        end
        xlim([bigRasterWindow(1) bigRasterWindow(2)])
        title({fileName;strcat(truncatedNames{i},' Cluster ',num2str(rasterCounter))})
        set(0, 'DefaulttextInterpreter', 'none')

        %plots concatenated histograms
        subplot(4,1,3)
        hold on
        plot(matclustStruct.(truncatedNames{i}).IndivRastersControl{rasterCounter},'k')
        plot(matclustStruct.(truncatedNames{i}).IndivRastersLaser{rasterCounter},'r')
        for k = 1:pipReps
            plot([zeroVector(k) zeroVector(k)],[ylim],'b')
            plot([toneEndVector(k) toneEndVector(k)],[ylim],'b')
            plot([endVector(k) endVector(k)],[ylim],'k','LineWidth',2)
        end
        
        %plots concatenated histograms
        subplot(4,1,4)
        hold on
        plot(smooth(matclustStruct.(truncatedNames{i}).IndivRastersControl{rasterCounter},25),'k')
        plot(smooth(matclustStruct.(truncatedNames{i}).IndivRastersLaser{rasterCounter},25),'r')
        for k = 1:pipReps
            plot([zeroVector(k) zeroVector(k)],[ylim],'b')
            plot([toneEndVector(k) toneEndVector(k)],[ylim],'b')
            plot([endVector(k) endVector(k)],[ylim],'k','LineWidth',2)
        end
        
%         %plots z-scored histogram
%         subplot(4,3,3)
% 
%         plot(matclustStruct.(truncatedNames{i}).BigZScore{j}{1}(:,2),...
%             matclustStruct.(truncatedNames{i}).BigZScore{j}{1}(:,1),'k','LineWidth',2)
%         hold on
%         plot(matclustStruct.(truncatedNames{i}).BigZScore{j}{2}(:,2),...
%             matclustStruct.(truncatedNames{i}).BigZScore{j}{2}(:,1),'b','LineWidth',2)
% 
%         plot([0 0],[ylim],'r');
%         plot([blockDur blockDur],[ylim],'r');
%         for k = 1:pipReps
%             plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
%         end
%         xlim([bigRasterWindow(1) bigRasterWindow(2)])
%         title('Z Scored Histogram')

        
%         %plots simple rasters
%         subplot(2,3,2)
% 
%         plot(matclustStruct.(truncatedNames{i}).BigRasters{j}(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,4) == 0,2),...
%             matclustStruct.(truncatedNames{i}).BigRasters{j}(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,4) == 0,1),'k.','markersize',5)
%         hold on
%         plot(matclustStruct.(truncatedNames{i}).BigRasters{j}(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,4) == 1,2),...
%             matclustStruct.(truncatedNames{i}).BigRasters{j}(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,4) == 1,1),'b.','markersize',5)
%         ylim([0 length(firstTTLDIO1)])
%         xlim([bigRasterWindow(1) bigRasterWindow(2)])
%         set(gca,'YDir','Reverse')
%         plot([0 0],[ylim],'r','LineWidth',2);
%         plot([blockDur blockDur],[ylim],'r','LineWidth',2);
%         for k = 1:pipReps
%             plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
%         end
% 
%         title({fileName;strcat(truncatedNames{i},' Cluster ',num2str(j))})
%         set(0, 'DefaulttextInterpreter', 'none')

        
        % plot heatmap of binned responses
%         subplot(2,3,5)
%         imagesc(pipSpikeStorage{i,j}')
%         colormap hot
%         
%         % plot plots of binned responses
%         subplot(2,3,3)
%         plot(mean(pipSpikeStorage{i,j},1))
%         
%         subplot(2,3,6)
%         plot(mean(pipSpikeStorage{i,j},2))
        

        spikeGraphName = strcat(fileName,truncatedNames{i},' Cluster ',num2str(j),'PulseAnalysis');
        savefig(hFig,spikeGraphName);
        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end


save(fullfile(pname,fname),'matclustStruct');



end
