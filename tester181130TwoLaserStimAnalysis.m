
%TOGGLES FOR ENABLING/DISABLING FEATURES
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.
toggleROC = 0; %toggle for tuning on/off ROC analysis

%PARAMETERS FOR BASIC ARRANGEMENT OF DATA
s.Parameters.RasterWindow = [-6 4]; %ratio for raster window. will be multiplied by toneDur
s.Parameters.BaseWindow = [-4 0];
s.Parameters.ToneWindow = [0 1];
s.Parameters.GenWindow = [0 2];
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate
s.Parameters.LFPWindow = [-1 3];

%stuff for significance
s.Parameters.calcWindow = [0 2]; %defines period for looking for responses, based on toneDur
s.Parameters.zLimit = [0.05 0.01 0.001];
s.Parameters.minSpikes = 100; %minimum number of spikes in baseline, if lower, triggers a warning
s.Parameters.minSigSpikes = 5; %minimum number of significant points to record a significant response.
s.Parameters.PercentCutoff = 99.9; %for significance in latency calculations
s.Parameters.BaselineCutoff = 95; %for the onset in latency calculations
s.Parameters.latBin = 0.001; %histogram bins for latency and significance calculations
s.Parameters.SigSmoothWindow = 11; %window of smoothing for calculations of significance

%for duplicate elimination
s.Parameters.DownSampFactor = 10; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr

%for rotary encoder:
s.Parameters.InterpolationStepRotary = 0.01; %interpolation steps in seconds.

%for edr
s.Parameters.EDRdownsamp = 20; %number of samples to downsample by. Smoothing is likely unnecessary
s.Parameters.EDRTimeCol = 1;
s.Parameters.EDRTTLCol = 3;
s.Parameters.EDRPiezoCol = 2;

%for cell type
s.Parameters.PVLim = [0.0004 0.0005];
s.Parameters.ChatLim = 1.1;

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%for plotting the laser related stuff
s.Parameters.LaserWindow = [-0.3 0.4]; %plotting window for rasters
s.Parameters.LaserBin = 0.01; %histogram bin size
s.Parameters.LaserAnalysis = [-0.2,0;0.1,0.3];
% s.Parameters.LaserLim = 0.015; %maximum lag value for calculation.

%set other things
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
% format short
disp('Parameters Set')
%% sets up file saving stuff
saveName = strcat(fileName,'FullTuningAnalysis','.mat');
fname = saveName;
pname = pwd;
%check to see if analysis file already exists.
searchNames = what;
searchNames = searchNames.mat;
searchResult = strcmp(searchNames,saveName);
if find(searchResult == 1) 
    prevFile = 1;
    disp('Previous Analysis Found!!')
else
    prevFile = 0;
end

%% Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
if ismac
    subFoldersCell = strsplit(subFolders,':')';
elseif ispc
    subFoldersCell = strsplit(subFolders,';')';
end

%% Find and ID Matclust Files for Subsequent Analysis. Generates Structured Array for Data Storage
%pull matclust file names
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');
%pull lfp file names. this allows detection of the number of trodes. 
try
    [lfpFiles] = functionFileFinder(subFoldersCell,'LFP','LFP');
    s.NumberTrodes = length(lfpFiles);
catch
    [paramFiles] = functionFileFinder(subFoldersCell,'matclust','param');
    s.NumberTrodes = length(paramFiles) - length(matclustFiles);
end
disp('Matclust Files Extracted')
%generate placeholder structure
% s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Parameters.RPVTime,...
    matclustFiles,s,s.Parameters.ClusterWindow,subFoldersCell);

%if a previous analysis file exists, I want to default to that selection of
%duplicates, to avoid repeatedly choosing which units to save. 
if prevFile == 0
disp('No Previous Analysis, Performing Duplicate Elimination')
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
elseif prevFile == 1
    disp('Using Existing Analysis for Duplicate Elimination')
    n = load(saveName);
    %load deleted units
    if isfield(n.s,'DeletedUnits')
        s.DeletedUnits = n.s.DeletedUnits;
        %now delete bad units
        delUnits = fields(s.DeletedUnits);
        for indCount = 1:length(delUnits)
            s = rmfield(s,delUnits{indCount});
        end
        %replace designation array and names. 
        s.DesignationArray = n.s.DesignationArray;
        s.DesignationName = n.s.DesignationName;
    else
        %replace designation array and names. 
        s.DesignationArray = n.s.DesignationArray;
        s.DesignationName = n.s.DesignationName;
    end
    
    disp('Finished! Closing original file')
    clear n
end
    
%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;


%generate master array for 2-d storage of important values
masterData = zeros(numUnits,10);
masterHeader = cell(10,1);
masterInd = 1;

%% Produce Indices for sorting by position on shank
for reInd = 1:length(s.DesignationName) %go through every unit
    %extract information about tetrode
    testTet = s.DesignationArray(reInd,1);
    %extract biggest waveform
    [bigVal bigInd] = max(max(s.(s.DesignationName{reInd}).AverageWaveForms));
    %calculate overall index
    overInd(reInd) = (testTet - 1)*4 + bigInd;
end

s.PeakWaveIndex = overInd;
[s.SortedPeakWaveIndex s.SortedPeakWaveOrder] = sort(overInd);
s.ShankLength = max(s.ShankMap(:,1)/s.Shanks);
%now lets try and map this onto real space. We will have to do something to
%allow things that are on the same electrode to map out separately. 

%Go through each electrode and find matching sites. if exist, then store.
%if extras, then figure out spacing. 
cellCount = 1;
posArray = zeros(numUnits,2);
for i = 1:length(s.ShankMap)
    %see if there are cells centered on given electrode
    cellFind = find(s.SortedPeakWaveIndex == s.ShankMap(i,1));
    if length(cellFind) == 0
        disp(strcat('No Cells For Electrode',num2str(s.ShankMap(i,1))))
    elseif length(cellFind) == 1
        disp(strcat('Single Cell For Electrode',num2str(s.ShankMap(i,1))))
        posArray(cellCount,1) = s.ShankMap(i,1);
        posArray(cellCount,2) = s.ShankMap(i,3);
        cellCount = cellCount + 1;
    elseif length(cellFind) > 1
        disp(strcat('Multiple Cells For Electrode',num2str(s.ShankMap(i,1))))
        incAmount = 1/length(cellFind); %amount to increment different cells by
        for incInd = 1:length(cellFind)
            posArray(cellCount,1) = s.ShankMap(i,1) + (incInd-1)*incAmount;
            posArray(cellCount,2) = s.ShankMap(i,3);
            cellCount = cellCount + 1;
        end
    end
end
%make things negative so they plot out correctly
posArray(:,1) = -1*posArray(:,1);
posArray(posArray(:,2) == 2,1) = posArray(posArray(:,2) == 2,1)+s.ShankLength;

s.ElectrodePositionDisplayArray = posArray;

%store this into master
masterData(:,masterInd) = posArray(s.SortedPeakWaveOrder,1); masterHeader{masterInd} = 'Distance from Top of Shank'; masterInd = masterInd + 1;
masterData(:,masterInd) = posArray(s.SortedPeakWaveOrder,2); masterHeader{masterInd} = 'Shank Designation'; masterInd = masterInd + 1;
masterData(:,masterInd) = s.SortedPeakWaveOrder; masterHeader{masterInd} = 'SimpleShankOrder'; masterInd = masterInd + 1;

disp('Beginning DIO Extraction...')
%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
if length(D1FileName) == 0
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','Din2');
end
D1FileName = D1FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);
%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);
blueLaserData = dioTimes;
blueLaserDiff = dioTimeDiff;

[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D6');
if length(D1FileName) == 0
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','Din6');
end
D1FileName = D1FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);
%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);
redLaserData = dioTimes;
redLaserDiff = dioTimeDiff;


rasterWindowBlue = [-0.1 0.2];
rasterWindowRed = [-0.5 1];

binBlue = 0.01;
binRed = 0.05;
histVectBlue = [rasterWindowBlue(1):binBlue:rasterWindowBlue(2)];
histVectRed = [rasterWindowRed(1):binRed:rasterWindowRed(2)];
blueLaserHistStore = [];
redLaserHistStore = [];
for i = 1:numUnits
    %determine cell type
    waveForms = s.(desigNames{i}).AverageWaveForms;
    
    %find size of waveforms
    waveSize = size(waveForms);

    %find maximum waveform by max peak size
    maxWave = max(waveForms);
    [maxVal maxInd] = max(maxWave);

    %chose the big wave, interpolate to fine degree
    chosenWave = waveForms(:,maxInd);
    interpVect = [1:0.1:40];
    interpWave = interp1(1:40,chosenWave,interpVect,'spline');
    
    
    %now we need to find the peak. Find this starting at point 10. 

    [pkVal pkInd] = max(interpWave(100:end));
    pkInd = pkInd + 100 - 1;
    %now we need to find the minimum following the peak

    [troughVal troughInd] = min(interpWave(pkInd:end));
    troughInd = troughInd + pkInd - 1;

    peakTrough(i) = (troughInd - pkInd)/300000;
    
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    
    %determine ISI covariance 
    isiTimes = diff(spikeTimes);
    isiCov(i) = std(isiTimes)/mean(isiTimes);
    
    
   alignTimes = blueLaserData;
   [rasters] = functionBasicRaster(spikeTimes,alignTimes,rasterWindowBlue);
   blueLaserRasters = rasters;
   blueLaserHistStore(:,i) = hist(blueLaserRasters(:,1),histVectBlue)/length(blueLaserData)/binBlue;
   alignTimes = redLaserData;
   [rasters] = functionBasicRaster(spikeTimes,alignTimes,rasterWindowRed);
   redLaserRasters = rasters;
   redLaserHistStore(:,i) = hist(redLaserRasters(:,1),histVectRed)/length(redLaserData)/binRed;
   s.(desigNames{i}).BlueRaster = blueLaserRasters;
   s.(desigNames{i}).RedRaster = redLaserRasters;
   
   %find average rate
   alignTimes =[blueLaserData;redLaserData];
   [rasters] = functionBasicRaster(spikeTimes,alignTimes,rasterWindowBlue);
   averageSpikeHolder = zeros(length(alignTimes),1);
   for k = 1:length(alignTimes)
        averageSpikeHolder(k) = size(find(rasters(:,2) == k & rasters(:,1) > rasterWindowBlue(1) & rasters(:,1) < 0),1);
   end
   averageRate(i) = mean(averageSpikeHolder/(-rasterWindowBlue(1)));
    
end

%determine cell type
cellTypes = zeros(numUnits,1);
for i = 1:numUnits
    if peakTrough(i) < s.Parameters.PVLim(1) & isiCov(i) > s.Parameters.ChatLim
        cellTypes(i) = 1; %pv cell
    elseif peakTrough(i) > s.Parameters.PVLim(2) & isiCov(i) < s.Parameters.ChatLim & averageRate(i) > 2
        cellTypes(i) = 2; %ChAT Cell
    elseif peakTrough(i) > s.Parameters.PVLim(2) & isiCov(i) > s.Parameters.ChatLim
        cellTypes(i) = 0; %MSN
    else
        cellTypes(i) = NaN; %label as unknown
    end
end
findPVs = find(cellTypes == 1);
findMSNs = find(cellTypes == 0);
findCHATs = find(cellTypes == 2);


%calculate modulation indices. 
%want to select first 50 ms for ChR2 stim, whole period for halo. 
blueMod = (mean(blueLaserHistStore(12:16,:)) - mean(blueLaserHistStore(6:10,:)))./(mean(blueLaserHistStore(12:16,:)) + mean(blueLaserHistStore(6:10,:)));
redMod = (mean(redLaserHistStore(12:21,:)) - mean(redLaserHistStore(1:10,:)))./(mean(redLaserHistStore(12:21,:)) + mean(redLaserHistStore(1:10,:)));

hFig = figure
plot(blueMod,redMod,'r.')
xlim([-1 1])
ylim([-1 1])
xlabel('Modulation Index(by 405)')
ylabel('Modulation Index(by 638)')
title('Optrode Recordings in DLS, A2A Cre with DIO ChR2 and NpHR')
spikeGraphName = 'OverallModulationIndexChR2NpHRTest';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now plot by actual modulation

hFig = figure
plot(mean(blueLaserHistStore(12:16,:)) - mean(blueLaserHistStore(6:10,:)),redMod,'r.')
% xlim([-1 1])
ylim([-1 1])
xlabel('deltaFR (by 405)')
ylabel('Modulation Index(by 638)')
title('Optrode Recordings in DLS, A2A Cre with DIO ChR2 and NpHR')
spikeGraphName = 'OverallModulationFRChR2NpHRTest';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1900 1000])
    subplot(2,2,1)
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    if ismember(i,findPVs)
        title(strcat('PV AverageFiringRate:',num2str(averageRate(i))))
    elseif ismember(i,findMSNs)
        title(strcat('MSN AverageFiringRate:',num2str(averageRate(i))))
    elseif ismember(i,findCHATs)
        title(strcat('CHAT AverageFiringRate:',num2str(averageRate(i))))
    else
        title(strcat('UNK AverageFiringRate:',num2str(averageRate(i))))
    end
    subplot(2,2,3)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    subplot(2,2,2)
    plot(histVectBlue,blueLaserHistStore(:,i),'b')
    title(desigNames{i})
    xlim(rasterWindowBlue)
    subplot(2,2,4)
    plot(histVectRed,redLaserHistStore(:,i),'r')
    xlim(rasterWindowRed)
    
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')
end