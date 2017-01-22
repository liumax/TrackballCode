%testerABBAanalysis
%test code to analyze recording results of ABBA stimulation.

function [s] = analysisABBA(fileName);


%% Hardcoded Variables:
params = struct;
params.rpvTime = 0.002; %limit to be considered an RPV.
params.clusterWindow = [-0.01,0.03]; %this is hardcoded for consistency
params.trodesFS = 30000; %sampling rate of trodes box.
params.rasterWindow = [-4,3]; %duration of raster window. These are numbers that will
params.histBin = 0.005; %bin size in seconds

%variables for tuning analysis
params.baselineBin = [-4,0]; %defines duration of baseline period based on toneDur. 
params.calcWindow = [0 2]; %defines period for looking for responses, based on toneDur
params.zLimit = [0.05 0.01 0.001];
params.numShuffle = 1000;
params.firstSpikeWindow = [0 1];%defines period for looking for first spike, based on toneDur
params.chosenSpikeBin = 1; %spike bin selected in binSpike (in the event of multiple spike bins)
params.minSpikes = 100; %minimum number of spikes to do spike shuffling
params.minSigSpikes = 2; %minimum number of significant points to record a significant response.
params.BaselineWindow = [-0.4 0]; %window for counting baseline spikes, in SECONDS. NOTE THIS IS DIFFERENT FROM RASTER WINDOW
params.BaselineCalcBins = 1; %bin size in seconds if there are insufficient baseline spikes to calculate a baseline rate.

%variables for latBinPeak calculations
params.toneWindow = [0,1];
params.genWindow = [0,3];
params.latBin = 0.001;
params.percentCutoff = 99.9;
params.baseCutoff = 95;


disp('Parameters Set')
%% Establishes the folder and subfolders for analysis. Adds all folders to
%path for easy access to files.
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
[s, truncatedNames] = functionMatclustExtraction(params.rpvTime,...
    matclustFiles,s,params.clusterWindow);
disp('Structured Array Generated, Names Extracted')

%save params!
s.Parameters = params;

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;


%% Import Sound Data
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundFile = soundFile.soundData;

%extract names!
soundNames = soundFile.Names;
%generate spike names from these. These will be used as headers for storage
%of spikes from a particular time period.
spikeNames = soundNames;
for i = 1:size(soundNames,1)
    spikeNames{i} = strcat('spikes',spikeNames{i});
end
disp('Sound Data Imported')

%% This code extracts DIO timepoints and states, and uses that information to
%%extrapolate which timepoints represent upward changes in the DIO. These
%%represent true TTL pulses.

%extract DIO filenames
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};

%extract actual data
[DIO1Data] = readTrodesExtractedDataFile(D1FileName);
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%remove just time and state data. convert to double for easier processing
DIO1Data = double([DIO1Data.fields(1).data,DIO1Data.fields(2).data]);
DIO2Data = double([DIO2Data.fields(1).data,DIO2Data.fields(2).data]);

[dio1Times,dio1TimeDiff] = functionBasicDIOCheck(DIO1Data,params.trodesFS);
[dio2Times,dio2TimeDiff] = functionBasicDIOCheck(DIO2Data,params.trodesFS);




























end