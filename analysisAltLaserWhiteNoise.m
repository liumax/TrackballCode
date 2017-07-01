

function [s] = analysisAltLaserWhiteNoise(fileName);
%% Constants and things you might want to tweak
%lets set some switches to toggle things on and off.
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.

s.Parameters.RasterWindow = [-4 5]; %seconds for raster window. 
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.05; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate

%for duplicate elimination
s.Parameters.DownSampFactor = 10; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr

%for rotary encoder:
s.Parameters.InterpolationStepRotary = 0.01;
laserPeriod = [0 1]; %time window in seconds around laser onset that I want to analyze. This is if I want to select for trials with locomotion at a specifc time
locoThresh = 0.9; %threshold for time points in which the locomotion is active for trial to be considered locomotion trial
windowPref = [0 1.5]; %preference window in seconds. This is the period over which I determine whether locomotor starts are more or less common than expected by random chance
prefReps = 1000; %repetitions for bootstrapping to determine if locomotor starts are more common during stimulation
avLocoThresh = 0.1; %minimum speed for a trial to be considered a locomotion trial.

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%for selecting MSNs vs PV
s.Parameters.PVLim = 0.0005;

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
    matclustFiles,s,s.Parameters.ClusterWindow,subFoldersCell);

if toggleDuplicateElimination ==1
    if length(s.DesignationName) > 1
        disp('Now Selecting Based on xCORR')
        [s] = functionDuplicateElimination(s,s.Parameters.DownSampFactor,...
            s.Parameters.corrSlide,s.Parameters.ThresholdComparison,s.Parameters.trodesFS,s.Parameters.RPVTime,s.Parameters.ClusterWindow,...
            s.ShankDesignation,s.ShankMap,s.Shanks);
    end
else
    disp('NOT EXECUTING DUPLICATE ELIMINATION')
end

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%calculate window sizes and binning

histBinNum = (s.Parameters.RasterWindow(2)-s.Parameters.RasterWindow(1))/s.Parameters.histBin;
histBinVector = [s.Parameters.RasterWindow(1)+s.Parameters.histBin/2:s.Parameters.histBin:s.Parameters.RasterWindow(2)-s.Parameters.histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes

%% Pull data from sound file
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundData = soundFile.soundData; %pulls sound data info
s.SoundData = soundData;
toneDur = soundData.ToneDuration;
toneReps = soundData.ToneRepetitions;
laserLag = min(soundData.LaserLag);

%% Extract DIO information. Tuning curve should rely on just one DIO output, DIO1.

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);

%need to separate by trial type!
shortDIO = find(dioTimeDiff < laserLag*0.8);
shortDIODel = unique([shortDIO,shortDIO+1]);
dioTone = dioTimes;
dioTone(shortDIODel) = [];
%this isolates tones!
%now figure out which ones are laser
dioLaser = dioTimes(shortDIO);

%now find the tone laser DIO
medDIO = find(dioTimeDiff > laserLag * 0.8 & dioTimeDiff < laserLag * 1.2);
dioToneLaser1 = dioTimes(medDIO+1);
dioToneLaser2 = dioTimes(medDIO-1);

%now lets isolate tone only dio.
[C,ia,ib] = intersect(dioTone,dioToneLaser1);
dioToneOnly = dioTone;
dioToneOnly(ia) = [];

[C,ia,ib] = intersect(dioLaser,dioToneLaser2);
dioLaserOnly = dioLaser;
dioLaserOnly(ia) = [];

dioLaserOnlyFake = dioLaserOnly + laserLag; %This is a fake time point to generate rasters that are aligned to same times as tones.

%get a common time for all starts


totalTrialNum = 3*toneReps;


%% Extract data from rotary encoder.
[s] = functionRotaryExtraction(s,s.Parameters.trodesFS,s.Parameters.InterpolationStepRotary,subFoldersCell);
% 
% %rasterize this data
% jumpsBack = round(s.Parameters.RasterWindow(1)/s.Parameters.InterpolationStepRotary);
% jumpsForward = round(s.Parameters.RasterWindow(2)/s.Parameters.InterpolationStepRotary);
% velRaster = zeros(jumpsForward-jumpsBack+1,totalTrialNum);
% for i = 1:totalTrialNum
%     %find the time from the velocity trace closest to the actual stim time
%     targetInd = find(s.RotaryData.Velocity(:,1)-dioTimes(i) > 0,1,'first');
%     %pull appropriate velocity data
%     velRaster(:,i) = s.RotaryData.BinaryLocomotion([targetInd+jumpsBack:targetInd+jumpsForward]);
% end
% 
% %make average trace:
% averageVel = mean(velRaster,2);
% velSTD = std(velRaster,0,2);
% velVector = [s.Parameters.RasterWindow(1):s.Parameters.InterpolationStepRotary:s.Parameters.RasterWindow(2)];
% velDispVector = [s.Parameters.RasterWindow(1):1:s.Parameters.RasterWindow(2)];
% velDispIndex = [1:round(1/s.Parameters.InterpolationStepRotary):(jumpsForward-jumpsBack+1)];
% velZero = find(velVector >= 0,1,'first');
% 
% 
% %now lets find locomotion based on average from the trial
% avLoco = mean(velRaster);
% avLocoPos = find(avLoco>avLocoThresh);
% avLocoNeg = find(avLoco<=avLocoThresh);

%% Process spiking information: extract rasters and histograms, both general and specific to frequency/db
for i = 1:numUnits
    %allocate a bunch of empty arrays.
    fullHistTone = zeros(histBinNum,1);
    fullHistLaserTone = zeros(histBinNum,1);
    fullHistLaser = zeros(histBinNum,1);
    
    averageSpikeHolder = zeros(totalTrialNum,1);
    
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    
    %now calculate rasters for each case
    [rastersTone] = functionBasicRaster(spikeTimes,dioToneOnly,s.Parameters.RasterWindow);
    [rastersLaser] = functionBasicRaster(spikeTimes,dioLaserOnlyFake,s.Parameters.RasterWindow);
    [rastersToneLaser] = functionBasicRaster(spikeTimes,dioToneLaser1,s.Parameters.RasterWindow);
    
    %store data into structured array
    s.(desigNames{i}).RasterToneOnly = rastersTone;
    s.(desigNames{i}).RasterLaserOnly = rastersLaser;
    s.(desigNames{i}).RasterToneLaser = rastersToneLaser;
    
    %generate histograms and store.
    [counts centers] = hist(rastersTone(:,1),histBinVector);
    s.(desigNames{i}).HistogramToneOnly = counts/length(dioToneOnly)/s.Parameters.histBin;
    
    [counts centers] = hist(rastersLaser(:,1),histBinVector);
    s.(desigNames{i}).HistogramLaserOnly = counts/length(dioLaserOnlyFake)/s.Parameters.histBin;
    
    [counts centers] = hist(rastersToneLaser(:,1),histBinVector);
    s.(desigNames{i}).HistogramToneLaser = counts/length(dioToneLaser1)/s.Parameters.histBin;
    
end

%calculate and plot LFP information
% [lfpStruct] = functionLFPaverage(master, s.Parameters.LFPWindow, s,homeFolder,fileName, 1, 1, 1, 1);
% s.LFP = lfpStruct;

%% Plotting
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    %plots average waveform
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    title(strcat('OverallRate:',num2str(s.(desigNames{i}).OverallFiringRate)))
    %plots ISI
    subplot(4,6,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    %plot histograms
    subplot(3,3,2)
    hold on
    plot(histBinVector,s.(desigNames{i}).HistogramToneOnly,'k','LineWidth',2)
%     hold on
%     plot(histBinVector,s.(desigNames{i}).HistogramToneOnly - s.(desigNames{i}).HistogramStandardDeviationTone,'k','LineWidth',1)
%     plot(histBinVector,s.(desigNames{i}).HistogramToneOnly + s.(desigNames{i}).HistogramStandardDeviationTone,'k','LineWidth',1)
    
    plot(histBinVector,s.(desigNames{i}).HistogramToneLaser,'g','LineWidth',2)
    plot(histBinVector,s.(desigNames{i}).HistogramLaserOnly,'b','LineWidth',2)
%     plot(histBinVector,s.(desigNames{i}).OverallHistToneLaser - s.(desigNames{i}).HistogramStandardDeviationToneLaser,'g','LineWidth',1)
%     plot(histBinVector,s.(desigNames{i}).OverallHistToneLaser + s.(desigNames{i}).HistogramStandardDeviationToneLaser,'g','LineWidth',1)
   
    plot([0 0],[ylim],'b');
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    title({fileName;desigNames{i}})
    set(0, 'DefaulttextInterpreter', 'none')
    
    %plots rasters (chronological)
    subplot(3,3,3)
    plot(s.(desigNames{i}).RasterLaserOnly(:,1),...
        s.(desigNames{i}).RasterLaserOnly(:,2),'k.','markersize',5)
    hold on
    ylim([0 toneReps])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Laser Response')
    
    
    subplot(3,3,6)
    plot(s.(desigNames{i}).RasterToneOnly(:,1),...
        s.(desigNames{i}).RasterToneOnly(:,2),'k.','markersize',5)
    hold on
    ylim([0 toneReps])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Tone Response')
    
    subplot(3,3,9)
    plot(s.(desigNames{i}).RasterToneLaser(:,1),...
        s.(desigNames{i}).RasterToneLaser(:,2),'k.','markersize',5)
    hold on
    ylim([0 toneReps])
    xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
    plot([0 0],[ylim],'b');
    title('Tone Laser Response')
    
    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'threepeatAnalysis');
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