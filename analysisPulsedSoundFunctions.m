%This is meant to be a function built code that analyzes the experiment in
%which the animal is delivered pulsed sounds. 

%This needs the following in the same folder: matclust file of picked
%spikes and matlab file with audio order

function [] = analysisPulsedSoundFunctions(fileName);

%% constants
sampleRate = 30000;

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

%extracts matclust file names and removes periods which allow structured array formation.
truncatedNames = matclustFiles;
numTrodes = length(truncatedNames);
for i = 1:length(truncatedNames)
    truncatedNames{i} = truncatedNames{i}(1:find(truncatedNames{i} == '.')-1);
end

trodesDesignation = cell(size(truncatedNames));
for i = 1:length(truncatedNames)
    trodesDesignation{i} = truncatedNames{i}(17:end);
end

%generates structured array for storage of data
matclustStruct = struct;
for i = 1:length(truncatedNames);
    matclustStruct.(truncatedNames{i}) = [];
end
matclustStruct.NumberTrodes = numTrodes;

%% Extracts Sound Data from soundFile
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData;
%extract the data for the sounds
pipITI = soundData.PipITI;
pipReps = soundData.PipRepetitions;
blockITI = soundData.BlockITI;
blockNum = soundData.BlockRepetitions;
laserITI = soundData.LaserTTLITI;
laserSwitch = soundData.LaserOn;
toneDur = soundData.ToneDuration;

blockDur = pipITI*pipReps;

%set up raster windows
rasterWindow = [-soundData.ToneDuration 2*soundData.ToneDuration]; %for plotting individual pips
bigRasterWindow = [-blockDur/2 1.5*blockDur]; %for plotting entire block!
bigRasterAxis=[bigRasterWindow(1):0.01:bigRasterWindow(2)-0.01];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
%These are for refractory period violations and looking at spike ITIs
rpvTime = 0.0013; %limit to be considered an RPV.
clusterWindow = [-0.01,0.03]; %this is hardcoded for consistency
%Use raster parameters to set histogram settings.
histBin = 0.005; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes
baselineBins = [1,find(histBinVector<0,1,'last')];

bigHistBin = 0.05; %bin size in seconds
bigHistBinNum = (bigRasterWindow(2)-bigRasterWindow(1))/bigHistBin;
bigHistBinVector = [bigRasterWindow(1)+bigHistBin/2:bigHistBin:bigRasterWindow(2)-bigHistBin/2]; %this is vector with midpoints of all histogram bins
%bigHistBinVector is for the purposes of graphing. This provides a nice axis
bigBaselineBins = [1,find(bigHistBinVector<0,1,'last')];


%% DIO extractions
%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIO1Data] = readTrodesExtractedDataFile(D1FileName);
%extracts port states and times of changes. First column is time, second is
%state.
dio1Data(:,1) = double(DIO1Data.fields(1).data)/sampleRate;
dio1Data(:,2) = double(DIO1Data.fields(2).data);
%pull all differences in state for DIO data (diff subtracts x(2)-x(1)). +1
%fudge factor adjusts for indexing. MUST BE DOUBLE!!
DIO1Diff = find(diff(dio1Data(:,2))==1)+1;
DIO1High = find(dio1Data(:,2) == 1);
%finds all points which are down to up states, extracts these times.
DIO1True = intersect(DIO1Diff,DIO1High);
DIO1True = dio1Data(DIO1True,1);
%finds differences between time points
DIO1TrueDiff = diff(DIO1True);

if laserSwitch == 1; %only activates if laser is used.
    %find DIO folder and D2 file for analysis
    [D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
    D2FileName = D2FileName{1};
    [DIO2Data] = readTrodesExtractedDataFile(D2FileName);
    dio2Data(:,1) = double(DIO2Data.fields(1).data)/sampleRate;
    dio2Data(:,2) = double(DIO2Data.fields(2).data);
    %pull all differences in state for DIO data (diff subtracts x(2)-x(1)). +1
    %fudge factor adjusts for indexing. MUST BE DOUBLE!!
    DIO2Diff = find(diff(dio2Data(:,2))==1)+1;
    DIO2High = find(dio2Data(:,2) == 1);
    %finds all points which are down to up states, extracts these times.
    DIO2True = intersect(DIO2Diff,DIO2High);
    DIO2True = dio2Data(DIO2True,1);
    %finds differences between time points
    DIO2TrueDiff = diff(DIO2True);
end
%% analyze the differences in the DIO data. 
%DIO 1 should be direct inputs from the audio card. DIO 2 is laser.

%give enough of a fudge factor to account for some degree of wobbling
blockLimit = pipITI*2; %determines what classifies/separates blocks
pipITILimit = pipITI*0.9; %determines if something counts as a laser

%finds tone TTLs (excludes laser TTLs)
toneTTLFinder = find(DIO1TrueDiff > pipITILimit)+1;
toneTTLFinder = [1;toneTTLFinder];


%finds the first TTLs of a block of sounds
firstTTLDIO1 = find(DIO1TrueDiff > blockLimit)+1;
firstTTLDIO1 = [1;firstTTLDIO1];


%check to make sure things line up with what was supposed to be delivered
if laserSwitch == 1
    if length(firstTTLDIO1) == 2* blockNum
        disp('Correct Number of Blocks');
        disp(num2str(length(firstTTLDIO1)));
        laserBlockFinder = find(diff(firstTTLDIO1)>pipReps);
        laserBlockFinder = [laserBlockFinder;laserBlockFinder(end)+2];
        controlBlockFinder = find(diff(firstTTLDIO1)==pipReps);
        if length(laserBlockFinder) ~= blockNum
            error('INCORRECT IDENTIFICATION OF LASER BLOCKS')
        elseif length(laserBlockFinder) ~= length(controlBlockFinder)
            error('Mismatch between laser and control blocks')
        end
    else
        error('INCORRECT BLOCK NUMBER')
    end
else
    if length(firstTTLDIO1) == blockNum
        disp('Correct Number of Blocks');
        disp(num2str(length(firstTTLDIO1)));
    else
        error('INCORRECT BLOCK NUMBER')
    end
end

%check to make sure things line up with what was supposed to be delivered
if laserSwitch == 1
    %finds the first TTLs of laser stim
    firstTTLDIO2 = find(DIO2TrueDiff > blockLimit)+1;
    firstTTLDIO2 = [1;firstTTLDIO2];
    if length(firstTTLDIO2) == blockNum
        disp('Correct Number of Laser Blocks');
        disp(num2str(length(firstTTLDIO2)));
    else
        error('INCORRECT BLOCK NUMBER')
    end
end
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

master(:,1) = toneTTLFinder;
master(:,2:5) = 0;

bigMaster(:,1) = firstTTLDIO1;
bigMaster(:,2) = 0;
if ~ischar(soundData.TargetFrequency)
    bigMaster(laserBlockFinder,2) = soundData.TargetFrequency;
end
if ~ischar(soundData.ControlFrequency)
    bigMaster(controlBlockFinder,2) = soundData.ControlFrequency;
end

pipSpikeStorage = cell(numTrodes,1);

for i = 1:numTrodes
    %% extract spikes
    [matclustStruct, clusterSizer] = functionSpikeWaveExtraction(rpvTime,...
        i,matclustFiles,matclustStruct,truncatedNames,clusterWindow);
    %% make rasters/histograms of all tone presentations (small rasters)
    [matclustStruct] = functionRasterHistExtraction(i,clusterSizer,...
    master,baselineBins,matclustStruct,truncatedNames,...
    rasterWindow,histBin,histBinVector);
    
    %% make rasters/histograms of big picture (big rasters)
    [matclustStruct] = functionPulseBigRasterHistExtraction(i,clusterSizer,...
    bigMaster,bigBaselineBins,matclustStruct,truncatedNames,...
    bigRasterWindow,bigHistBin,bigHistBinVector,laserSwitch);
% by design, the big raster function will detect average firing rate and
% should overwrite that of the functionRasterHistExtraction. 
    %% do simple binning procedure to look at changes to responses
    for j = 1:clusterSizer
        pipSpikeMatrix = zeros(pipReps,length(firstTTLDIO1));
        for k = 1:size(allPulseMatrix,1)*size(allPulseMatrix,2)
            pipSpikeMatrix(k) = size(find(matclustStruct.(truncatedNames{i}).Rasters{j}(:,1) == k &...
                matclustStruct.(truncatedNames{i}).Rasters{j}(:,2) > 0 &...
                matclustStruct.(truncatedNames{i}).Rasters{j}(:,2) < toneDur),1);
        end
        pipSpikeStorage{i,j} = pipSpikeMatrix;
    end
end

%% Plotting
for i = 1:numTrodes
    for j = 1:matclustStruct.(truncatedNames{i}).Clusters
        hFig = figure;
        set(hFig, 'Position', [10 10 1280 1000])
        
        %plots average waveform
        subplot(4,3,1)
        hold on
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,2),'LineWidth',2)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,1),'r','LineWidth',1)
        plot(matclustStruct.(truncatedNames{i}).AverageWaveForms(:,j,3),'r','LineWidth',1)
        title(strcat('AverageFiringRate:',num2str(matclustStruct.(truncatedNames{i}).AverageFiringRate(j))))
        
        %plots ISI
        subplot(4,3,4)
        hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000)
        histMax = max(hist(matclustStruct.(truncatedNames{i}).ISIData{j},1000));
        line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
        xlim(clusterWindow)
        title({strcat('ISI RPV %: ',num2str(matclustStruct.(truncatedNames{i}).RPVs(j)));...
            strcat(num2str(matclustStruct.(truncatedNames{i}).RPVNumber(j)),'/',num2str(matclustStruct.(truncatedNames{i}).TotalSpikeNumber(j)))})
        
        %plots histograms
        if laserSwitch == 1
            subplot(4,3,7)
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}{1}(:,2),...
                matclustStruct.(truncatedNames{i}).BigHistogram{j}{1}(:,1),'k','LineWidth',2)
            hold on
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}{1}(:,2),...
                matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,j,1,1),'r')
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}{1}(:,2),...
                matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,j,2,1),'r')
            
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}{2}(:,2),...
                matclustStruct.(truncatedNames{i}).BigHistogram{j}{2}(:,1),'b','LineWidth',2)
            hold on
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}{2}(:,2),...
                matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,j,1,2),'c')
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}{2}(:,2),...
                matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,j,2,2),'c')

            plot([0 0],[ylim],'r');
            plot([blockDur blockDur],[ylim],'r');
            for k = 1:pipReps
                plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
            end
            xlim([bigRasterWindow(1) bigRasterWindow(2)])
            title('Overall Histogram')
        else
            subplot(4,3,7)
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}(:,2),...
                matclustStruct.(truncatedNames{i}).BigHistogram{j}(:,1),'k','LineWidth',2)
            hold on
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}(:,2),...
                matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,j,1),'b')
            plot(matclustStruct.(truncatedNames{i}).BigHistogram{j}(:,2),...
                matclustStruct.(truncatedNames{i}).BigStandardErrorPlotting(:,j,2),'b')

            plot([0 0],[ylim],'r');
            plot([blockDur blockDur],[ylim],'r');
            for k = 1:pipReps
                plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
            end
            xlim([bigRasterWindow(1) bigRasterWindow(2)])
            title('Overall Histogram')
        end
        
        
        %plots z-scored histogram
        subplot(4,3,10)
        if laserSwitch == 1
            plot(matclustStruct.(truncatedNames{i}).BigZScore{j}{1}(:,2),...
                matclustStruct.(truncatedNames{i}).BigZScore{j}{1}(:,1),'k','LineWidth',2)
            hold on
            plot(matclustStruct.(truncatedNames{i}).BigZScore{j}{2}(:,2),...
                matclustStruct.(truncatedNames{i}).BigZScore{j}{2}(:,1),'b','LineWidth',2)
            
            plot([0 0],[ylim],'r');
            plot([blockDur blockDur],[ylim],'r');
            for k = 1:pipReps
                plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
            end
            xlim([bigRasterWindow(1) bigRasterWindow(2)])
            title('Z Scored Histogram')
        else
            plot(matclustStruct.(truncatedNames{i}).BigZScore{j}(:,2),...
                matclustStruct.(truncatedNames{i}).BigZScore{j}(:,1),'k','LineWidth',2)
            hold on
            plot([0 0],[ylim],'r');
            plot([blockDur blockDur],[ylim],'r');
            for k = 1:pipReps
                plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
            end
            xlim([bigRasterWindow(1) bigRasterWindow(2)])
            title('Z Scored Histogram')
        end
        
        %plots simple rasters
        subplot(2,3,2)
        if laserSwitch == 1
            plot(matclustStruct.(truncatedNames{i}).BigRasters{j}(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,4) == 0,2),...
                matclustStruct.(truncatedNames{i}).BigRasters{j}(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,4) == 0,1),'k.','markersize',5)
            hold on
            plot(matclustStruct.(truncatedNames{i}).BigRasters{j}(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,4) == 1,2),...
                matclustStruct.(truncatedNames{i}).BigRasters{j}(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,4) == 1,1),'b.','markersize',5)
            ylim([0 length(firstTTLDIO1)])
            xlim([bigRasterWindow(1) bigRasterWindow(2)])
            set(gca,'YDir','Reverse')
            plot([0 0],[ylim],'r','LineWidth',2);
            plot([blockDur blockDur],[ylim],'r','LineWidth',2);
            for k = 1:pipReps
                plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
            end

            title({fileName;strcat(truncatedNames{i},' Cluster ',num2str(j))})
            set(0, 'DefaulttextInterpreter', 'none')
        else
            plot(matclustStruct.(truncatedNames{i}).BigRasters{j}(:,2),...
                matclustStruct.(truncatedNames{i}).BigRasters{j}(:,1),'k.','markersize',5)
            hold on
            ylim([0 length(firstTTLDIO1)])
            xlim([bigRasterWindow(1) bigRasterWindow(2)])
            set(gca,'YDir','Reverse')
            plot([0 0],[ylim],'r','LineWidth',2);
            plot([blockDur blockDur],[ylim],'r','LineWidth',2);
            for k = 1:pipReps
                plot([pipITI*k pipITI*k],[ylim],'r','LineWidth',1)
            end

            title({fileName;strcat(truncatedNames{i},' Cluster ',num2str(j))})
            set(0, 'DefaulttextInterpreter', 'none')
        end
        spikeGraphName = strcat(fileName,trodesDesignation{i},' Cluster ',num2str(j),'PulseAnalysis');
        savefig(hFig,spikeGraphName);
        
        % plot heatmap of binned responses
        subplot(2,3,5)
        imagesc(pipSpikeStorage{i,j}')
        colormap hot
        
        % plot plots of binned responses
        subplot(2,3,3)
        plot(mean(pipSpikeStorage{i,j},1))
        
        subplot(2,3,6)
        plot(mean(pipSpikeStorage{i,j},2))
        

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
end






end
