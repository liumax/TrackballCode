%This is code to try and make STRFs!



%% Constants and things you might want to tweak
%TOGGLES FOR ENABLING/DISABLING FEATURES
s.Parameters.toggleRPV = 0; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 0; %1 means you want to eliminate duplicates.
toggleROC = 0; %toggle for tuning on/off ROC analysis

%PARAMETERS FOR BASIC ARRANGEMENT OF DATA
s.Parameters.STRFWin = [];
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate

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

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%for plotting the laser related stuff
s.Parameters.LaserWindow = [-0.1 0.1]; %plotting window for rasters
s.Parameters.LaserBin = 0.001; %histogram bin size
s.Parameters.LaserLim = 0.015; %maximum lag value for calculation.

%set other things
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
% format short
disp('Parameters Set')
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

%% Extract DIO information. Tuning curve should rely on just one DIO output, DIO1.
disp('Beginning DIO Extraction...')
%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
if length(D1FileName) == 0
    [D1FileName] = functionFileFinder(subFoldersCell,'DIO','Din1');
end
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);
s.SoundData.ToneTimes = dioTimes;
s.SoundData.ToneTimeDiff = dioTimeDiff;



%pull D2 info (laser ID)
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
if length(D2FileName) == 0
    [D2FileName] = functionFileFinder(subFoldersCell,'DIO','Din2');
end
D2FileName = D2FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%pulls out DIO up state onsets.
[dio2Times,dio2TimeDiff] = functionBasicDIOCheck(DIO2Data,s.Parameters.trodesFS);

%check for tons of DIO2. if alot this should be an ID tuning.
if length(dio2Times) > 10
    idToggle = 1;
    disp('Laser Pulses Detected: ID Session')
else
    idToggle = 0;
    disp('No Laser Pulses Detected')
end






