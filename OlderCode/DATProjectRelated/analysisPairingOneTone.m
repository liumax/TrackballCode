%This code is meant to analyze trials in which I perform pairing of one
%tone with dopamine stimulation. 

%Inputs:
% fileName: name of the sound file, without the .mat extension.
% 
% Outputs:
% s: structured array storing data from analysis. 


function [s] = analysisPairingOneTone(fileName);

%This is meant to be the analysis code for pairing experiments. This code
%should analyze the tuning curves before and after, the playing of long
%tones before and after, and the actual changes that may occur with laser
%light itself. 

% fileName = '160622_ML160410A_L_2504_reversePairingProtocol';

%% Hardcoded Variables:
params = struct;
params.rpvTime = 0.0013; %limit to be considered an RPV.
params.clusterWindow = [-0.01,0.03]; %this is hardcoded for consistency
params.trodesFS = 30000; %sampling rate of trodes box.
params.rasterWindow = [-4,3]; %duration of raster window. These are numbers that will
%be multiplied by the tone duration. EX: raster window for 0.1sec tone will
%be -100 to 300 ms.
params.pairingWindow = [-10,20];
params.histBin = 0.005; %bin size in seconds
params.clims1 = [-1 1]; %limits for the range of heatmaps for firing. Adjust if reach saturation. Currently based on log10

%variables for tuning analysis
params.baselineBin = [-4,0]; %defines duration of baseline period based on toneDur. 
params.calcWindow = [0 2]; %defines period for looking for responses, based on toneDur
params.zLimit = [0.05 0.01 0.001];
params.numShuffle = 1000;
params.firstSpikeWindow = [0 1];%defines period for looking for first spike, based on toneDur
params.chosenSpikeBin = 2; %spike bin selected in binSpike (in the event of multiple spike bins)
params.minSpikes = 80; %minimum number of spikes to do spike shuffling
params.minSigSpikes = 2; %minimum number of significant points to record a significant response.
params.BaselineWindow = [-0.4 0]; %window for counting baseline spikes, in SECONDS. NOTE THIS IS DIFFERENT FROM RASTER WINDOW
params.BaselineCalcBins = 1; %bin size in seconds if there are insufficient baseline spikes to calculate a baseline rate.
params.ThresholdHz = 4; %minimum response in Hz to be counted as significant.

%for latbinPeak
params.toneWindow = [0,1];
params.genWindow = [0,3];
params.latBin = 0.001;
params.percentCutoff = 99.9;
params.baseCutoff = 95;

%for duplicate elimination
params.DownSampFactor = 3; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
params.corrSlide = 0.05; % window in seconds for xcorr
params.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr

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
[paramFiles] = functionFileFinder(subFoldersCell,'matclust','param');

%generate placeholder structure
s = struct;

%put in parameters!
s.Parameters = params;
s.NumberTrodes = length(paramFiles)-length(matclustFiles);

%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(params.rpvTime,...
    matclustFiles,s,params.clusterWindow);
disp('Structured Array Generated, Names and Spikes Extracted')

if length(s.DesignationName) > 1
    disp('Now Selecting Based on xCORR')
    [s] = functionDuplicateElimination(s);
end

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;



%% Import Sound Data
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundFile = soundFile.fullData;
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

%pull all differences in state for DIO data (diff subtracts x(2)-x(1)). +1
%fudge factor adjusts for indexing. MUST BE DOUBLE!!
DIO1Diff = find(diff(DIO1Data(:,2))==1)+1;
DIO1High = find(DIO1Data(:,2) == 1);
%finds all points which are down to up states, extracts these times.
DIO1True = intersect(DIO1Diff,DIO1High);
DIO1True = DIO1Data(DIO1True,1);
%finds differences between time points
DIO1TrueDiff = diff(DIO1True);

disp('DIO Signals Processed')

%% pull variables from soundfile regarding signal TTLs to calculate acceptable range of signals for signal TTLs
signalPulseNum = soundFile.DividerTTLNumber;
signalITI = soundFile.DividerTTLiti; %signal TTL ITI in seconds
acceptRange = [0.9,1.1]; %range of values above or below ideal ITI that are acceptable
signalITI = signalITI*30000; %converts to trodes timestamps
signalRange = acceptRange*signalITI;

%% This next part finds all TTLs that are markers. This finds the ones that
%correspond with diff. This will correctly label the first three pulses of
%every sequence
findSignals = find(DIO1TrueDiff>signalRange(1) & DIO1TrueDiff<signalRange(2));
%generates array to hold first and last signal pulses, fills based on
%information about number of signal pulses
signalHolder = zeros(2,size(findSignals,1)/(signalPulseNum-1));
for i = 1:size(findSignals,1)/(signalPulseNum-1)
    signalHolder(1,i) = findSignals(1+(signalPulseNum-1)*(i-1));
    signalHolder(2,i) = signalHolder(1,i) + signalPulseNum - 1;
end

%% This now uses hardcoded values to extract the TTLs and time periods for each stage
%of the experiment. Currently, times are in Trodes samples, which are at 30
%kHz.
timesBaseline = [DIO1Data(1,1),DIO1True(signalHolder(1,1))];
timesTuningFirst = [DIO1True(signalHolder(2,1)),DIO1True(signalHolder(1,2))];
timesPairing = [DIO1True(signalHolder(2,2)),DIO1True(signalHolder(1,3))];
timesTuningSecond = [DIO1True(signalHolder(2,3)),DIO1True(end)];
%also using the same outputs from signalHolder, this determines which TTL
%pulses belong to which output. 
TTLsTuningFirst = DIO1True(signalHolder(2,1)+1:signalHolder(1,2)-1);
TTLsPairing = DIO1True(signalHolder(2,2)+1:signalHolder(1,3)-1);
TTLsTuningSecond = DIO1True(signalHolder(2,3)+1:end);
%stores these values into the master structure. Converts to seconds.
s.TTLs.(soundNames{1}) = TTLsTuningFirst/params.trodesFS;
s.TTLs.(soundNames{2}) = TTLsTuningSecond/params.trodesFS;
s.TTLs.(soundNames{3}) = TTLsPairing/params.trodesFS;
%Remember that this has converted time values to seconds.
s.TimePeriods.Baseline = timesBaseline/params.trodesFS;
s.TimePeriods.(soundNames{1}) = timesTuningFirst/params.trodesFS;
s.TimePeriods.(soundNames{2}) = timesTuningSecond/params.trodesFS;
s.TimePeriods.(soundNames{3}) = timesPairing/params.trodesFS;

if size(TTLsTuningFirst,1) == soundFile.TuningRepetitions;
    disp('Correct Number of Early Tuning Presentations')
else
    error('MISMATCHED EARLY TUNING PRESENTATIONS') 
end

if size(TTLsTuningSecond,1) == soundFile.TuningRepetitions;
    disp('Correct Number of Late Tuning Presentations')
else
    error('MISMATCHED LATE TUNING PRESENTATIONS') 
end


%% First thing I want to do is divvy up spikes to different time periods
[s] = functionNewPairingSpikeSeparator(s,...
   desigNames,spikeNames,soundNames);
disp('Spikes Separated')

%% calculates average firing rate for the initial period and overall firing rates across all periods
[s] = functionNewPairingAverageRate(s,desigNames,...
    spikeNames,soundNames);
disp('Average Rates Calculated')

%% Next thing is to analyze tuning curve chunks for differences.

%Analyze first tuning curve. 

%first pull out sound data for the relevant file
[s] = functionPairingSoundDataExtraction(s,...
    soundNames{1},soundFile); 
%next, pull tuning information. DOES NOT GRAPH
[s] = functionNewPairingTuning(s,desigNames,...
    spikeNames{1},soundNames{1},params); 
disp('First Tuning Curve Analyzed')
%now I analyze the second tuning curve.
%first pull out sound data for the relevant file
[s] = functionPairingSoundDataExtraction(s,...
    soundNames{2},soundFile); 
%next, pull tuning information. DOES NOT GRAPH
[s] = functionNewPairingTuning(s,desigNames,...
    spikeNames{2},soundNames{2},params); 
disp('Second Tuning Curve Analyzed')
%% next pull data from pairing session:

%extract sound data:
s.SoundData.(soundNames{3}) = soundFile.(soundNames{3});
%process TTLs (since these are not just unitary TTLs signaling tone onset).
%generates new TTL file under s.SoundData.soundNames{3}Original, which
%preserves the original information. 
[s] = functionNewPairingTTLAdjust(s,soundNames{3});

%now align spikes to tone presentations and store data as histogram.
[s] = functionNewPairingPairedToneHist(s,desigNames,...
    spikeNames{3},soundNames{3},params); 
disp('Pairing Analyzed')
%% NOW I NEED TO PLOT EVERYTHING IN A WAY THAT MAKES SENSE
[s] = functionOneTonePairingMasterPlot(numUnits,s,...
    desigNames,params,fileName,soundNames);



pname = pwd;
fname = strcat(fileName,'PairingAnalysis');
save(fullfile(pname,fname),'s');

end














