%This code is meant for the basic analysis of tuning curve data, examining
%and plotting basic properties of the tuning curve and displaying them as a
%figure. 

%%Inputs
%fileName: name of the target file, without file extension. 

%%Outputs
%s: structured array that has all the data from the analysis. s will
%contain the following: 
%
%n number of units, which are independent of clusters/trodes. This is to
%say that each unit will have its own named field. For example, 2 clusters
%on channel 10 will have two separate fields.

%DesignationArray/DesignationName: The array indicates which probe sites
%and clusters were used, designation names should be the names of all
%units.

%SoundData: all the sound data from the tuning curve, with addition of
%unique frequencies/dbs and number of frequencies/dbs.

function [s] = analysisVarDur(fileName);
%% Constants and things you might want to tweak
s.Parameters.RasterWindow = [-1 1.5]; %ratio for raster window. will be multiplied by toneDur
s.Parameters.ToneWindow = [0 1];
s.Parameters.GenWindow = [0 1.5];
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate
s.Parameters.zLimit = 3; %zlimit for calculating significant responses
s.Parameters.FirstSpikeWindow = [0 1];
s.Parameters.BaselineBin = [-1 0]; %ratio for bin from which baseline firing rate will be calculated
s.Parameters.LFPWindow = [-1 2];

%stuff for significance
s.Parameters.calcWindow = [0 1.5]; %defines period for looking for responses, based on toneDur
s.Parameters.zLimit = [0.05 0.01 0.001];
s.Parameters.numShuffle = 1000;
s.Parameters.minSpikes = 100; %minimum number of spikes to do spike shuffling
s.Parameters.minSigSpikes = 2; %minimum number of significant points to record a significant response.
s.Parameters.BaselineWindow = [-1 0]; %window for counting baseline spikes, in ratio of toneDur
s.Parameters.BaselineCalcBins = 1; %bin size in seconds if there are insufficient baseline spikes to calculate a baseline rate.
s.Parameters.PercentCutoff = 99.9;
s.Parameters.BaselineCutoff = 95;
s.Parameters.latBin = 0.001;
s.Parameters.ThresholdHz = 4; %minimum response in Hz to be counted as significant.

%for duplicate elimination
s.Parameters.DownSampFactor = 3; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.01; % percentage overlap to trigger xcorr


%% sets up file saving stuff
saveName = strcat(fileName,'FullTuningAnalysis','.mat');
fname = saveName;
pname = pwd;

%% Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%% Find and ID Matclust Files for Subsequent Analysis. Generates Structured Array for Data Storage
%pull matclust file names
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');
[paramFiles] = functionFileFinder(subFoldersCell,'matclust','param');
s.NumberTrodes = length(paramFiles)-length(matclustFiles);

%generate placeholder structure
% s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Parameters.RPVTime,...
    matclustFiles,s,s.Parameters.ClusterWindow);

if length(s.DesignationName) > 1
    disp('Now Selecting Based on xCORR')
    [s] = functionDuplicateElimination(s);
end

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%% Extracts Sound Data from soundFile
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;

toneDurs = soundData.ToneDurations;
numDurs = length(toneDurs);
toneAmp = soundData.ToneAmplitude;
toneFreq = soundData.ToneFreq;
toneOrder = soundData.ToneOrder;
toneReps = soundData.ToneReps;

% Recalculate raster windows and others based on tone durations
%here, I just want ghetto code that tells me if a cell responds to long
%tone durations. Therefore, all I really need is the rasterWindow defined. 


s.Parameters.RasterWindow = s.Parameters.RasterWindow * max(toneDurs);
rasterAxis=[s.Parameters.RasterWindow(1):0.001:s.Parameters.RasterWindow(2)-0.001];
histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
s.Parameters.LFPWindow = s.Parameters.LFPWindow * max(toneDurs);


% Pull unique tone durations
uniqueDurs = unique(toneOrder);

%Identify the trials belonging to each tone duration
toneOrder(:,2) = 0;
sortCounter = 1;

for i = 1:numDurs
    sortFinder = find(toneOrder(:,1) == uniqueDurs(i));
    toneOrder(sortFinder,2) = sortCounter:1:sortCounter+size(sortFinder,1) -1;
    sortCounter = sortCounter + size(sortFinder,1);
end

%% Extract DIO information. Tuning curve should rely on just one DIO output, DIO1.

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);

%check number of dioTimes
if length(dioTimes) ~= (length(toneOrder))
    error('Bad number of TTLs')
end



for i = 1:numUnits
    respRasters = cell(numDurs,1);
    respHists = cell(numDurs,1);
    for j=1:numDurs
        spikeTimes = s.(desigNames{i}).SpikeTimes;
        alignTimes = dioTimes(toneOrder(:,1) == uniqueDurs(j));
        [rasters] = functionBasicRaster(spikeTimes,alignTimes,s.Parameters.RasterWindow);
        respRasters{j} = rasters;
        [counts centers] = hist(rasters,histBinVector);
        respHists{j} = counts;
    end
    s.(desigNames{i}).AllRasters = respRasters;
    s.(desigNames{i}).AllHistograms = respHists;
end


%% Plotting
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    %plots average waveform
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
%     title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    %plots ISI
    subplot(4,6,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    %plot rasters to tones of different durations
    subplot(3,1,2)
    hold on
    for j = 1:numDurs
        plot(s.(desigNames{i}).AllRasters{j}(:,1),s.(desigNames{i}).AllRasters{j}(:,2)+(j-1)*toneReps,'.','Color',[(j-1)/numDurs 0 0])
    end
    
    xlim(s.Parameters.RasterWindow)
    ylim([0 length(toneOrder)])
    plot([0 0],[ylim],'g')
    for j = 1:numDurs
        plot([toneDurs(j) toneDurs(j)],[ylim],'g')
    end
    if isnumeric(toneFreq)
        title({strcat('Frequency:',num2str(toneFreq),' Amp:',num2str(toneAmp));...
        strcat('Rasters for:',num2str(toneDurs),' in ascending duration, red longer')})
    else
        title({strcat('Frequency:',(toneFreq),' Amp:',num2str(toneAmp));...
        strcat('Rasters for:',num2str(toneDurs),' in ascending duration, red longer')})
    end
    
    subplot(3,1,3)
    hold on
    for j = 1:numDurs
        plot(histBinVector,s.(desigNames{i}).AllHistograms{j}(:,1),'Color',[(j-1)/numDurs 0 0])
    end
    xlim(s.Parameters.RasterWindow)
    plot([0 0],[ylim],'g')
    for j = 1:numDurs
        plot([toneDurs(j) toneDurs(j)],[ylim],'g')
    end
    title(strcat('Histograms in ascending duration, red longer'))
    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end

%% Saving
save(fullfile(pname,fname),'s');

end