function [s] = analysisPairingFunctions(fileName);

%This is meant to be the analysis code for pairing experiments. This code
%should analyze the tuning curves before and after, the playing of long
%tones before and after, and the actual changes that may occur with laser
%light itself. 

% fileName = '160622_ML160410A_L_2504_reversePairingProtocol';


%% Hardcoded Variables:
%lets set some switches to toggle things on and off.
s.Parameters.toggleRPV = 1; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 1; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 1; %1 means you want to eliminate duplicates.
toggleROC = 0; %toggle for tuning on/off ROC analysis
toggleNumSpike = 0; %toggle for turning on/off selection based on number of total spikes. 

%parameters for data analysis
s.Parameters.rpvTime = 0.002; %limit to be considered an RPV.
s.Parameters.clusterWindow = [-0.01,0.03]; %this is hardcoded for consistency
s.Parameters.trodesFS = 30000; %sampling rate of trodes box.
s.Parameters.rasterWindow = [-4,3]; %duration of raster window. These are numbers that will
%be multiplied by the tone duration. EX: raster window for 0.1sec tone will
%be -100 to 300 ms.
s.Parameters.pairingWindow = [-10,20];
s.Parameters.histBin = 0.005; %bin size in seconds
s.Parameters.clims1 = [-1 1]; %limits for the range of heatmaps for firing. Adjust if reach saturation. Currently based on log10

%variables for tuning analysis
s.Parameters.baselineBin = [-4,0]; %defines duration of baseline period based on toneDur. 
s.Parameters.calcWindow = [0 2]; %defines period for looking for responses, based on toneDur
s.Parameters.zLimit = [0.05 0.01 0.001];
s.Parameters.numShuffle = 1000;
s.Parameters.firstSpikeWindow = [0 1];%defines period for looking for first spike, based on toneDur
s.Parameters.chosenSpikeBin = 2; %spike bin selected in binSpike (in the event of multiple spike bins)
s.Parameters.minSpikes = 80; %minimum number of spikes to do spike shuffling
s.Parameters.minSigSpikes = 2; %minimum number of significant points to record a significant response.
s.Parameters.BaselineWindow = [-0.4 0]; %window for counting baseline spikes, in SECONDS. NOTE THIS IS DIFFERENT FROM RASTER WINDOW
s.Parameters.BaselineCalcBins = 1; %bin size in seconds if there are insufficient baseline spikes to calculate a baseline rate.
s.Parameters.ThresholdHz = 4; %minimum response in Hz to be counted as significant.

%for latbinPeak
s.Parameters.toneWindow = [0,0.5];
s.Parameters.genWindow = [0,3];
s.Parameters.latBin = 0.001;
s.Parameters.percentCutoff = 99.9;
s.Parameters.baseCutoff = 95;

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

%for summary figure
s.Parameters.SumSmoothBin = 0.200; %time in sec for smoothing function

%for toggleNumSpike
numSpikeLim = 400; %limit for total number of spikes below which data will not be analyzed.

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
s.NumberTrodes = length(paramFiles)-length(matclustFiles);
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Parameters.rpvTime,...
    matclustFiles,s,s.Parameters.clusterWindow);
disp('Structured Array Generated, Names and Spikes Extracted')


if toggleDuplicateElimination ==1
    if length(s.DesignationName) > 1
        disp('Now Selecting Based on xCORR')
        [s] = functionDuplicateElimination(s,s.Parameters.DownSampFactor,...
            s.Parameters.corrSlide,s.Parameters.ThresholdComparison,s.Parameters.trodesFS,...
            s.Parameters.rpvTime,s.Parameters.clusterWindow);
    end
else
    disp('NOT EXECUTING DUPLICATE ELIMINATION')
end

sizeCount = 1;
if toggleNumSpike == 1
    for i = 1:size(s.DesignationName,2)
        sizeCheck = s.(s.DesignationName{i}).TotalSpikeNumber;
        if sizeCheck < numSpikeLim
            sizeFlag(sizeCount) = i;
            sizeCount = sizeCount + 1;
            disp(strcat('Eliminating',(s.DesignationName{i}),' with ',num2str(sizeCheck),' Spikes'))
            %removes the unit
            s = rmfield(s,(s.DesignationName{i}));
        end
    end
    %now lets amend designation array and designation names
    s.DesignationArray(sizeFlag,:) = [];
    s.DesignationName(sizeFlag) = [];
else
    disp('NOT SELECTING BY TOTAL SPIKE NUMBER')
end

%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%Extract data from rotary encoder.
[s] = functionRotaryExtraction(s,s.Parameters.trodesFS,s.Parameters.InterpolationStepRotary,subFoldersCell);

%% Import Sound Data
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);
soundFile = soundFile.fullData;
%extract names!
names = soundFile.Names;
%generate spike names from these. These will be used as headers for storage
%of spikes from a particular time period.
spikeNames = names;
for i = 1:size(names,1)
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
timesPresentationFirst = [DIO1True(signalHolder(2,2)),DIO1True(signalHolder(1,3))];
timesPairing = [DIO1True(signalHolder(2,3)),DIO1True(signalHolder(1,4))];
timesPresentationSecond = [DIO1True(signalHolder(2,4)),DIO1True(signalHolder(1,5))];
timesTuningSecond = [DIO1True(signalHolder(2,5)),DIO1Data(end,1)];
%also using the same outputs from signalHolder, this determines which TTL
%pulses belong to which output. 
TTLsTuningFirst = DIO1True(signalHolder(2,1)+1:signalHolder(1,2)-1);
TTLsPresentationFirst = DIO1True(signalHolder(2,2)+1:signalHolder(1,3)-1);
TTLsPairing = DIO1True(signalHolder(2,3)+1:signalHolder(1,4)-1);
TTLsPresentationSecond = DIO1True(signalHolder(2,4)+1:signalHolder(1,5)-1);
TTLsTuningSecond = DIO1True(signalHolder(2,5)+1:end);
%stores these values into the master structure. Converts to seconds.
s.TTLs.(names{1}) = TTLsTuningFirst/s.Parameters.trodesFS;
s.TTLs.(names{2}) = TTLsTuningSecond/s.Parameters.trodesFS;
s.TTLs.(names{3}) = TTLsPresentationFirst/s.Parameters.trodesFS;
s.TTLs.(names{4}) = TTLsPresentationSecond/s.Parameters.trodesFS;
s.TTLs.(names{5}) = TTLsPairing/s.Parameters.trodesFS;
%Remember that this has converted time values to seconds.
s.TimePeriods.Baseline = timesBaseline/s.Parameters.trodesFS;
s.TimePeriods.(names{1}) = timesTuningFirst/s.Parameters.trodesFS;
s.TimePeriods.(names{2}) = timesTuningSecond/s.Parameters.trodesFS;
s.TimePeriods.(names{3}) = timesPresentationFirst/s.Parameters.trodesFS;
s.TimePeriods.(names{4}) = timesPresentationSecond/s.Parameters.trodesFS;
s.TimePeriods.(names{5}) = timesPairing/s.Parameters.trodesFS;

%% check size!if the wrong size will throw error!


if size(TTLsPresentationFirst,1) == soundFile.PresentationRepetitions;
    disp('Correct Number of Early Long Presentations')
else
    disp('MISMATCHED EARLY LONG PRESENTATIONS') 
    disp('Initiating Repair Code')
    [s,repairedTTLs] = functionTTLRepairSystem(soundFile.PresentationRepetitions,soundFile.Tones1.Delays,s.TTLs.Tones1,0,0,0,s);
    s.TTLs.Tones1 = repairedTTLs;
end

if size(TTLsPresentationSecond,1) == soundFile.PresentationRepetitions;
    disp('Correct Number of Second Long Presentations')
else
    disp('MISMATCHED SECOND LONG PRESENTATIONS') 
    disp('Initiating Repair Code')
    [s,repairedTTLs] = functionTTLRepairSystem(soundFile.PresentationRepetitions,soundFile.Tones2.Delays,s.TTLs.Tones2,0,0,0,s);
    s.TTLs.Tones2 = repairedTTLs;
end

if size(TTLsTuningFirst,1) == soundFile.TuningRepetitions;
    disp('Correct Number of Early Tuning Presentations')
else
    disp('MISMATCHED EARLY TUNING PRESENTATIONS') 
    disp('Initiating Repair Code')
    [s,repairedTTLs] = functionTTLRepairSystem(soundFile.TuningRepetitions,soundFile.Tuning1.Delays,s.TTLs.Tuning1,0,0,0,s);
    s.TTLs.Tuning1 = repairedTTLs;
end

if size(TTLsTuningSecond,1) == soundFile.TuningRepetitions*soundFile.SecondTuningRatio;
    disp('Correct Number of Late Tuning Presentations')
else
    disp('MISMATCHED LATE TUNING PRESENTATIONS') 
    disp('Initiating Repair Code')
    [s,repairedTTLs] = functionTTLRepairSystem(soundFile.TuningRepetitions,soundFile.Tuning2.Delays,s.TTLs.Tuning2,0,laserSig,laserLag,s);
    s.TTLs.Tuning2 = repairedTTLs;
end

if size(TTLsPairing,1) == soundFile.PairingRepetitions;
    disp('Correct Number of Pairings')
else
    disp('MISMATCHED PAIRING PRESENTATIONS') 
    disp('Initiating Repair Code')
    [s,repairedTTLs] = functionTTLRepairSystem(soundFile.PairingRepetitions,soundFile.Pairing.ITI,s.TTLs.Pairing,1,soundFile.Pairing.LaserTriggerPulseITI,soundFile.Pairing.OptoStimDelay,s);
    s.TTLs.Pairing = repairedTTLs;
end


%% First thing I want to do is divvy up spikes to different time periods
[s] = functionPairingSpikeSeparator(s,...
   desigNames,spikeNames,names);
disp('Spikes Separated')

%% calculates average firing rate for the initial period and overall firing rates across all periods
[s] = functionPairingAverageRate(s,desigNames,...
    spikeNames,names);
disp('Overall Rates Calculated')

%% Next thing is to analyze tuning curve chunks for differences.

%Analyze first tuning curve. 

%first pull out sound data for the relevant file
[s] = functionPairingSoundDataExtraction(s,...
    names{1},soundFile); 
%next, pull tuning information. DOES NOT GRAPH
[s] = functionNewPairingTuning(s,desigNames,...
    spikeNames{1},names{1},s.Parameters); 
disp('First Tuning Curve Analyzed')
%now I analyze the second tuning curve.
%first pull out sound data for the relevant file
[s] = functionPairingSoundDataExtraction(s,...
    names{2},soundFile); 
%next, pull tuning information. DOES NOT GRAPH
[s] = functionNewPairingTuning(s,desigNames,...
    spikeNames{2},names{2},s.Parameters); 
disp('Second Tuning Curve Analyzed')

%%% EDITS TO HERE

%% Next I analyze the first and second long tone presentations

%extract sound data:

for i = 3:5
    s.SoundData.(names{i}) = soundFile.(names{i});
end

%next, pair this data with spiking data. Stores under "soundName3" divided
%into two structured arrays, one for target and one for control.
[s] = functionPairingToneAnalysis(s,desigNames,...
    spikeNames{3},names{3},s.Parameters); 

%pull data from the second set.
%next, pair this data with spiking data. Stores under "soundName3" divided
%into two structured arrays, one for target and one for control.
[s] = functionPairingToneAnalysis(s,desigNames,...
    spikeNames{4},names{4},s.Parameters);  

%% next pull data from pairing session:

%process TTLs (since these are not just unitary TTLs signaling tone onset.
[s] = functionPairingPairedTTLAdjust(s,names{5},...
    soundFile,names{5});
%next, pair this data with spiking data. Stores under "soundName3" divided
%into two structured arrays, one for target and one for control.
[s] = functionPairingToneAnalysis(s,desigNames,...
    spikeNames{5},names{5},s.Parameters);  

%NOW I NEED TO PLOT EVERYTHING IN A WAY THAT MAKES SENSE
%% get stuff for the summary figure
s.SumPlot = [];
timeWindow = [s.TimePeriods.Baseline(1),s.TimePeriods.Tuning2(2)];
s.SumPlot.TimeSpan = [timeWindow(1):10:timeWindow(2)];
smoothFact = round(s.Parameters.SumSmoothBin/s.Parameters.histBin);
for i = 1:numUnits
    %pull baseline firing rate
    s.SumPlot.TotalSpikes(i) = s.(desigNames{i}).TotalSpikeNumber;
    s.SumPlot.BaselineFiring(i) = s.(desigNames{i}).OverallFiringRate;
    %pull peak values for waveforms to plot against time. Find the biggest
    %waveform
    waveMax = max(s.(desigNames{i}).AverageWaveForms);
    [M maxInd] = max(waveMax);
    waveSize = max(squeeze(s.(desigNames{i}).AllWaves(:,maxInd,:)));
    s.SumPlot.WaveSize(i,:) = interp1(s.(desigNames{i}).SpikeTimes,waveSize,[timeWindow(1):10:timeWindow(2)]);
    
%     s.SumPlot.CtrlPre(i,:) = s.(desigNames{i}).Tones1ControlAnalysis.Histograms;
%     s.SumPlot.CtrlPost(i,:) = s.(desigNames{i}).Tones2ControlAnalysis.Histograms;
%     s.SumPlot.TarPre(i,:) = s.(desigNames{i}).Tones1TargetAnalysis.Histograms;
%     s.SumPlot.TarPost(i,:) = s.(desigNames{i}).Tones2TargetAnalysis.Histograms;
%     s.SumPlot.CtrlPair(i,:) = s.(desigNames{i}).PairingControlAnalysis.Histograms;
%     s.SumPlot.TarPair(i,:) = s.(desigNames{i}).PairingTargetAnalysis.Histograms;
%     s.SumPlot.CtrlSub(i,:) = s.SumPlot.CtrlPost(i,:) - s.SumPlot.TarPost(i,:);
%     s.SumPlot.TarSub(i,:) = s.SumPlot.TarPost(i,:) - s.SumPlot.TarPre(i,:);
%     
    s.SumPlot.CtrlPre(i,:) = interp1([1:1:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)],smooth(s.(desigNames{i}).Tones1ControlAnalysis.Histograms,smoothFact),[1:10:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)]);
    s.SumPlot.CtrlPost(i,:) = interp1([1:1:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)],smooth(s.(desigNames{i}).Tones2ControlAnalysis.Histograms,smoothFact),[1:10:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)]);
    s.SumPlot.TarPre(i,:) = interp1([1:1:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)],smooth(s.(desigNames{i}).Tones1TargetAnalysis.Histograms,smoothFact),[1:10:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)]);
    s.SumPlot.TarPost(i,:) = interp1([1:1:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)],smooth(s.(desigNames{i}).Tones2TargetAnalysis.Histograms,smoothFact),[1:10:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)]);
    s.SumPlot.CtrlPair(i,:) = interp1([1:1:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)],smooth(s.(desigNames{i}).PairingControlAnalysis.Histograms,smoothFact),[1:10:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)]);
    s.SumPlot.TarPair(i,:) = interp1([1:1:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)],smooth(s.(desigNames{i}).PairingTargetAnalysis.Histograms,smoothFact),[1:10:length(s.(desigNames{i}).Tones1ControlAnalysis.Histograms)]);
    s.SumPlot.CtrlSub(i,:) = s.SumPlot.CtrlPost(i,:) - s.SumPlot.TarPost(i,:);
    s.SumPlot.TarSub(i,:) = s.SumPlot.TarPost(i,:) - s.SumPlot.TarPre(i,:);
    %do peak normalization
    s.SumPlot.NormCtrlPre(i,:) = (s.SumPlot.CtrlPre(i,:)-min(s.SumPlot.CtrlPre(i,:)))/(max(s.SumPlot.CtrlPre(i,:))-min(s.SumPlot.CtrlPre(i,:)));
    s.SumPlot.NormCtrlPost(i,:) = (s.SumPlot.CtrlPost(i,:)-min(s.SumPlot.CtrlPost(i,:)))/(max(s.SumPlot.CtrlPost(i,:))-min(s.SumPlot.CtrlPost(i,:)));
    s.SumPlot.NormTarPre(i,:) = (s.SumPlot.TarPre(i,:)-min(s.SumPlot.TarPre(i,:)))/(max(s.SumPlot.TarPre(i,:))-min(s.SumPlot.TarPre(i,:)));
    s.SumPlot.NormTarPost(i,:) = (s.SumPlot.TarPost(i,:)-min(s.SumPlot.TarPost(i,:)))/(max(s.SumPlot.TarPost(i,:))-min(s.SumPlot.TarPost(i,:)));
    s.SumPlot.NormCtrlPair(i,:) = (s.SumPlot.CtrlPair(i,:)-min(s.SumPlot.CtrlPair(i,:)))/(max(s.SumPlot.CtrlPair(i,:))-min(s.SumPlot.CtrlPair(i,:)));
    s.SumPlot.NormTarPair(i,:) = (s.SumPlot.TarPair(i,:)-min(s.SumPlot.TarPair(i,:)))/(max(s.SumPlot.TarPair(i,:))-min(s.SumPlot.TarPair(i,:)));
    s.SumPlot.NormCtrlSub(i,:) = (s.SumPlot.CtrlSub(i,:)-min(s.SumPlot.CtrlSub(i,:)))/(max(s.SumPlot.CtrlSub(i,:))-min(s.SumPlot.CtrlSub(i,:)));
    s.SumPlot.NormTarSub(i,:) = (s.SumPlot.TarSub(i,:)-min(s.SumPlot.TarSub(i,:)))/(max(s.SumPlot.TarSub(i,:))-min(s.SumPlot.TarSub(i,:)));
    %do z scoring
    s.SumPlot.ZCtrlPre(i,:) = zscore(s.SumPlot.CtrlPre(i,:));
    s.SumPlot.ZCtrlPost(i,:) = zscore(s.SumPlot.CtrlPost(i,:));
    s.SumPlot.ZTarPre(i,:) = zscore(s.SumPlot.TarPre(i,:));
    s.SumPlot.ZTarPost(i,:) = zscore(s.SumPlot.TarPost(i,:));
    s.SumPlot.ZCtrlPair(i,:) = zscore(s.SumPlot.CtrlPair(i,:));
    s.SumPlot.ZTarPair(i,:) = zscore(s.SumPlot.TarPair(i,:));
    s.SumPlot.ZCtrlSub(i,:) = zscore(s.SumPlot.CtrlSub(i,:));
    s.SumPlot.ZTarSub(i,:) = zscore(s.SumPlot.TarSub(i,:));
end

[B,I] = sort(s.SumPlot.BaselineFiring);
s.SumPlot.AxisVector = [s.Parameters.rasterWindow(1)*s.SoundData.Tones1.ToneDuration:1:s.Parameters.rasterWindow(2)*s.SoundData.Tones1.ToneDuration];
s.SumPlot.AxisDesig = [0:length(s.SumPlot.CtrlPre(i,:))/(length(s.SumPlot.AxisVector)-1):length(s.SumPlot.CtrlPre(i,:))];
s.SumPlot.AxisDesig(1) = 1;
%% do ROC calculation for velocity
if toggleROC == 1
    for i = 1:numUnits
        targetName = desigNames{i}
        [s] = functionLocomotionROC(s,targetName);
    end
else
    for i = 1:numUnits
        s.(desigNames{i}).TrueAUC = [];
        s.(desigNames{i}).ShuffleAUC = [];
    end
end

%% Plotting
%first plot summary figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure;
set(hFig, 'Position', [10 80 1240 850])
%plot histogram of firing rates
subplot(4,3,1)
hist(s.SumPlot.BaselineFiring,100)
title('Histogram of Firing Rates')
%now plot peak waveform voltage over time
subplot(4,3,4)
plot(s.SumPlot.WaveSize')
xlim([0 length(s.SumPlot.WaveSize)])
title('Peak Waveform over Time')
%now plot heatmaps of neural responses!
subplot(4,3,2)
imagesc(s.SumPlot.NormCtrlPre(I,:))
title({fileName;'Norm Control Pre'})
set(gca,'XTick',s.SumPlot.AxisDesig)
set(gca,'XTickLabel',s.SumPlot.AxisVector)
subplot(4,3,3)
imagesc(s.SumPlot.NormTarPre(I,:))
title('Target Pre')
set(gca,'XTick',s.SumPlot.AxisDesig)
set(gca,'XTickLabel',s.SumPlot.AxisVector)
subplot(4,3,5)
imagesc(s.SumPlot.NormCtrlPair(I,:))
title('Ctrl Pair')
set(gca,'XTick',s.SumPlot.AxisDesig)
set(gca,'XTickLabel',s.SumPlot.AxisVector)
subplot(4,3,6)
imagesc(s.SumPlot.NormTarPair(I,:))
title('Target Pair')
set(gca,'XTick',s.SumPlot.AxisDesig)
set(gca,'XTickLabel',s.SumPlot.AxisVector)
subplot(4,3,8)
imagesc(s.SumPlot.NormCtrlPost(I,:))
title('Ctrl Post')
set(gca,'XTick',s.SumPlot.AxisDesig)
set(gca,'XTickLabel',s.SumPlot.AxisVector)
subplot(4,3,9)
imagesc(s.SumPlot.NormTarPost(I,:))
title('Target Post')
set(gca,'XTick',s.SumPlot.AxisDesig)
set(gca,'XTickLabel',s.SumPlot.AxisVector)
subplot(4,3,11)
imagesc(s.SumPlot.NormCtrlSub(I,:))
title('Ctrl Sub')
set(gca,'XTick',s.SumPlot.AxisDesig)
set(gca,'XTickLabel',s.SumPlot.AxisVector)
subplot(4,3,12)
imagesc(s.SumPlot.NormTarSub(I,:))
title('Target Sub')
set(gca,'XTick',s.SumPlot.AxisDesig)
set(gca,'XTickLabel',s.SumPlot.AxisVector)
spikeGraphName = strcat(fileName,'SummaryFigure');
savefig(hFig,spikeGraphName);
%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

if toggleTuneSelect == 1 %if you want tuning selection...
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    decisionTuning = zeros(numUnits,1);
    for i = 1:numUnits
        [s] = functionPairingMasterPlot(i,s,...
        desigNames,s.Parameters,fileName,names,hFig);
        %ask for input! 
        promptCounter = 1; %This is used to run the while loop.
        whileCounter = 0; %this is the counter that gets updated to exit the loop

        while whileCounter ~= promptCounter
            try
                prompt = 'Is this unit tuned? (y/n)';
                str = input(prompt,'s');
                if str~='n' & str~='y'
                    error
                else
                    whileCounter = 1;
                end
            catch
            end
        end
        if strfind(str,'y')
            decisionTuning(i) = 1;
        elseif strfind(str,'n')
            decisionTuning(i) = 0;
        end

        %clear figure.
        clf
    end
    s.TuningDecision = decisionTuning;
else
    for i = 1:numUnits
        hFig = figure;
        set(hFig, 'Position', [10 80 1240 850])
        [s] = functionPairingMasterPlot(i,s,...
        desigNames,s.Parameters,fileName,names,hFig);
    end
end

    


pname = pwd;
fname = strcat(fileName,'Analysis');
save(fullfile(pname,fname),'s');

end














