%This is meant to be a function built code that analyzes the experiment in
%which the animal is delivered pulsed sounds. 

%This needs the following in the same folder: matclust file of picked
%spikes and matlab file with audio order

function [] = testNewPulseLaserAnalysis(fileName);


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
numUnits = length(s.DesignationName);
desigName = s.DesignationName;

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

for i = 1:numUnits
    %first, I want to extract information about entire tone pip set. This
    %means the big rasters/histograms
    blockRasters = cell(length(pipITI),2);
    blockHist = cell(length(pipITI),2);
    for k = 1:2
        for j = 1:length(pipITI)
            matrixFinder = find(blockMatrix(:,5) == pipITI(j) & blockMatrix(:,6) == (k-1));
            [rasters] = functionBasicRaster(s.(desigName{i}).SpikeTimes,blockMatrix(matrixFinder,1),bigRasterWindow(:,j));
            blockRasters{j,k} = rasters;
            blockHist{j,k} = hist(rasters(:,1),bigHistVector{j});
        end
    end
    %second, extract all data from every single pip. Save this with data
    %about that pip
    indivRasters = cell(1,1);
    indivHist = zeros(length(allPulseMatrix),length(histBinVector));
    [rasters] = functionBasicRaster(s.(desigName{i}).SpikeTimes,allPulseMatrix(:,1),rasterWindow);
    indivRasters = rasters;
    for j = 1:length(allPulseMatrix)
        indivHist(j,:) = hist(rasters(rasters(:,2) == j,1),histBinVector);
    end
    s.(desigName{i}).ControlRasters = blockRasters(:,1);
    s.(desigName{i}).ControlHist = blockHist(:,1);
    s.(desigName{i}).LaserRasters = blockRasters(:,2);
    s.(desigName{i}).LaserHist = blockHist(:,2);
    s.(desigName{i}).IndivRasters = indivRasters;
    s.(desigName{i}).IndivHist = indivHist;
end


%Plot!
for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 10 1280 1000])
    subplot(4,6,1)
    hold on
    plot(s.(desigName{i}).AverageWaveForms(:,2),'LineWidth',2)
    plot(s.(desigName{i}).AverageWaveForms(:,1),'r','LineWidth',1)
    plot(s.(desigName{i}).AverageWaveForms(:,3),'r','LineWidth',1)
    title(strcat('OverallFiringRate:',num2str(s.(desigName{i}).OverallFiringRate)))

    %plots ISI
    subplot(4,6,2)
    hist(s.(desigName{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigName{i}).ISIGraph,1000));
    line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(clusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigName{i}).RPVPercent));...
        strcat(num2str(s.(desigName{i}).RPVNumber),'/',num2str(s.(desigName{i}).TotalSpikeNumber))})
    
    %plot out big rasters (for control)
    for j = 1:length(pipITI)
        subplot(length(pipITI)*2,3,(2+(j-1)*3))
        plot(s.(desigName{i}).ControlRasters{j}(:,1),s.(desigName{i}).ControlRasters{j}(:,2),'k.')
        xlabel('Time (seconds)')
        ylabel('Trial Number')
        xlim(bigRasterWindow(:,j))
        ylim([min(s.(desigName{i}).ControlRasters{j}(:,2)) max(s.(desigName{i}).ControlRasters{j}(:,2))])
        title(strcat('CONTROL RASTER PIP ITI:',num2str(pipITI(j))))
    end
    for j = 1:length(pipITI)
        subplot(length(pipITI)*2,3,(3+(j-1)*3))
        plot(bigHistVector{j},s.(desigName{i}).ControlHist{j},'k')
        xlabel('Time (seconds)')
        ylabel('Trial Number')
        xlim(bigRasterWindow(:,j))
        title(strcat('CONTROL HISTOGRAM PIP ITI:',num2str(pipITI(j))))
    end
    newStart1 = 2 +(length(pipITI))*3;
    newStart2 =  3 +(length(pipITI))*3;
    for j = 1:length(pipITI)
        subplot(length(pipITI)*2,3,(newStart1+(j-1)*3))
        plot(s.(desigName{i}).LaserRasters{j}(:,1),s.(desigName{i}).LaserRasters{j}(:,2),'k.')
        xlabel('Time (seconds)')
        ylabel('Trial Number')
        xlim(bigRasterWindow(:,j))
        ylim([min(s.(desigName{i}).LaserRasters{j}(:,2)) max(s.(desigName{i}).LaserRasters{j}(:,2))])
        title(strcat('LASER RASTER PIP ITI:',num2str(pipITI(j))))
    end
    for j = 1:length(pipITI)
        subplot(length(pipITI)*2,3,(newStart2+(j-1)*3))
        plot(bigHistVector{j},s.(desigName{i}).LaserHist{j},'k')
        xlabel('Time (seconds)')
        ylabel('Trial Number')
        xlim(bigRasterWindow(:,j))
        title(strcat('LASER HISTOGRAM PIP ITI:',num2str(pipITI(j))))
    end
    spikeGraphName = strcat(fileName,desigName{i},'PulseAnalysis');
    savefig(hFig,spikeGraphName);
    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end


save(fullfile(pname,fname),'s');



end
