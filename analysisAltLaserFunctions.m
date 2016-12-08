function [] = analysisAltLaserFunctions(fileName);



% fileName = '160711TestingAltPairing';

%% Hardcoded Variables:
rpvTime = 0.0013; %limit to be considered an RPV.
clusterWindow = [-0.01,0.03]; %this is hardcoded for consistency
trodesFS = 30000; %sampling rate of trodes box.
rasterWindow = [-2,3]; %duration of raster window. These are numbers that will
%be multiplied by the tone duration. EX: raster window for 0.1sec tone will
%be -100 to 300 ms.
histBin = 0.005; %bin size in seconds
histBinLong = 0.05; %bin size for long presentations. 
clims1 = [-1 1]; %limits for the range of heatmaps for firing. Adjust if reach saturation. Currently based on log10

%% Establishes the folder and subfolders for analysis. Adds all folders to
%path for easy access to files.
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%% generate a master structure for storage of everything.
masterStruct = struct;

%% Find and ID Matclust Files for Subsequent Analysis.
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');

%extracts matclust file names and removes periods which allow structured
%array formation.
truncatedNames = matclustFiles;
numTrodes = length(truncatedNames);
for i = 1:length(truncatedNames)
    truncatedNames{i} = truncatedNames{i}(1:find(truncatedNames{i} == '.')-1);
    masterStruct.(truncatedNames{i}) = [];
end

%saves number of electodes that I'm reading from.
masterStruct.NumberTrodes = numTrodes;

%generates even shorter names for figure presentation. 
trodesDesignation = cell(size(truncatedNames));
for i = 1:length(truncatedNames)
    trodesDesignation{i} = truncatedNames{i}(17:end);
end

%% Goes and pulls spike times, overall firing rate, waveforms, rpvs from matclust files
[masterStruct] = functionPairingSpikeExtractor(truncatedNames,...
    matclustFiles,rpvTime,clusterWindow,masterStruct);

%saves overall firing rate as average firing rate
for i = 1:size(truncatedNames,2)
    masterStruct.(truncatedNames{i}).AverageFiringRates = masterStruct.(truncatedNames{i}).OverallFiringRate;
end

%% First thing is to import the sound file data.
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundFile = soundFile.soundData;
%

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

%% pull variables from soundfile regarding signal TTLs to calculate acceptable range of signals for laser paired trials
signalITI = -soundFile.LaserTiming; %signal TTL ITI in seconds
acceptRange = [0.9,1.1]; %range of values above or below ideal ITI that are acceptable
signalITI = signalITI*30000; %converts to trodes timestamps
signalRange = acceptRange*signalITI;

%% This next part finds all TTLs that are markers. This finds the ones that
%correspond with diff. This will correctly label the first three pulses of
%every sequence
findSignals = find(DIO1TrueDiff>signalRange(1) & DIO1TrueDiff<signalRange(2));
%Separates laser and non laser trials.
trialsLaser = DIO1True(findSignals+1)/trodesFS;
trialsNoLaser = DIO1True(findSignals+2)/trodesFS; %this is a cheating method: just assumes that pulse following paired pulse is unpaired. Will not work if not alternating.


%% This now saves these TTLs into the structured array.
masterStruct.TTLs.UnpairedStimuli = trialsNoLaser;
masterStruct.TTLs.PairedStimuli = trialsLaser;

%% check size!if the wrong size will throw error!
if size(trialsLaser,1) == size(soundFile.PairedStimuli.Frequencies,1);
    disp('Correct Number of Laser Trials')
else
    error('MISMATCHED LASER TRIALS') 
end

if size(trialsNoLaser,1) == size(soundFile.UnpairedStimuli.Frequencies,1);
    disp('Correct Number of No Laser Trials')
else
    error('MISMATCHED UNPAIRED TRIALS') 
end

%% Next thing is to analyze tuning curve chunks for differences.

%generate a names array for calling different trial types
names = cell(2,1);
names{1} = 'UnpairedStimuli';
names{2} = 'PairedStimuli';


%Analyze paired stimuli tuning curve. 

%first pull out sound data for the relevant file
[masterStruct] = functionPairingSoundDataExtraction(masterStruct,...
    names{1},soundFile); 
%next, pull tuning information. DOES NOT GRAPH
[masterStruct] = functionPairingTuning(masterStruct,truncatedNames,...
    'SpikeTimes',names{1},names{1},rasterWindow,histBin,clusterWindow,...
    clims1,rpvTime,trodesDesignation,fileName); 

%now I analyze the unpaired tuning curve.
%first pull out sound data for the relevant file
[masterStruct] = functionPairingSoundDataExtraction(masterStruct,...
    names{2},soundFile); 
%next, pull tuning information. DOES NOT GRAPH
[masterStruct] = functionPairingTuning(masterStruct,truncatedNames,...
    'SpikeTimes',names{2},names{2},rasterWindow,histBin,clusterWindow,...
    clims1,rpvTime,trodesDesignation,fileName); 

%NOW I NEED TO PLOT EVERYTHING IN A WAY THAT MAKES SENSE
[masterStruct] = functionAltMasterPlot(numTrodes,masterStruct,...
    truncatedNames,rpvTime,clusterWindow,rasterWindow,histBin,...
    clims1,fileName,trodesDesignation,names);

pname = pwd;
fname = strcat(fileName,'AltLaserAnalysis');
save(fullfile(pname,fname),'masterStruct');

end














