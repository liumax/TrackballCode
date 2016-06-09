%This is meant to be a function built code that replicates the functions of
%pre-existing tuning curve analysis

%This needs the following in the same folder: matclust file of picked
%spikes, matlab file with audio order
%and log file from MBED. 

function [] = functionTuningCruveFunctions(fileName);

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

%% Extracts Sound Data from soundFile, including freq, db, reps.
soundName = strcat(fileName,'.mat');
soundFile = open(soundName);

%extracts frequency information.
master = zeros(size(soundFile.soundData.Frequencies,1),5);
master(:,2) = soundFile.soundData.Frequencies;
uniqueFreqs = unique(master(:,2));
master(:,3) = soundFile.soundData.dBs;
uniqueDBs = unique(master(:,3));

numFreqs = size(uniqueFreqs,1);
numDBs = size(uniqueDBs,1);

%% Store Info into structured array.
matclustStruct.UniqueFreqs = uniqueFreqs;
matclustStruct.UniqueDBs = uniqueDBs;
matclustStruct.SoundTimes = master(:,1);
matclustStruct.Frequencies = master(:,2);
matclustStruct.dBs = master(:,3);
%also stores parameters for rep number and tone duration.
matclustStruct.ToneReps = soundFile.soundData.ToneRepetitions;
matclustStruct.ToneDur = soundFile.soundData.ToneDuration;

%% Set and/or Generate Raster and Histogram Parameters, store in Structured Array
rasterWindow = [-matclustStruct.ToneDur,matclustStruct.ToneDur*3];
lfpWindow = [-matclustStruct.ToneDur,matclustStruct.ToneDur*3];
rasterAxis=[rasterWindow(1):0.001:rasterWindow(2)-0.001];
%These are for refractory period violations and looking at spike ITIs
rpvTime = 0.0013; %limit to be considered an RPV.
clusterWindow = [-0.01,0.03]; %this is hardcoded for consistency
%clims1 is a value that sets the limits for heatmaps for displaying firing.
clims1 = [-1 1]; %limits for the range of heatmaps for firing. Adjust if reach saturation. Currently based on log10
%Use raster parameters to set histogram settings.
histBin = 0.005; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins
%histBinVector is for the purposes of graphing. This provides a nice axis
%for graphing purposes
baselineBins = [1,find(histBinVector<0,1,'last')];

%store values in structured array.
matclustStruct.RasterAxis = rasterAxis;
matclustStruct.HistogramAxis = histBinVector;
matclustStruct.RasterLimits = rasterWindow;
%%

%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%extracts port states and times of changes
dioState = double(DIOData.fields(2).data);
dioTime = double(DIOData.fields(1).data);

inTimes = dioTime(dioState == 1)/30000;
inTimesSpacing = diff(inTimes);
dioState = [];
dioTime = [];
if size(inTimes,1) ~= size(master,1)
    if size(inTimes,1) == size(master,1) + 1
        inTimes(end) = [];
        master(:,1) = inTimes;
    elseif inTimesSpacing(size(master,1)) > 3*mean(inTimesSpacing)
        master(:,1) = inTimes(1:size(master,1));
        disp('HOLY SHIT YOUR TTL PULSES TO DIO 1 ARE FUCKED, BUT ADJUSTED')
    else
        disp('I DONT KNOW WHATS GOING ON')
        pause
    end
elseif size(inTimes,1) == size(master,1)
    master(:,1) = inTimes;
end

inTimes = [];

%% Finds the number of octaves, makes array of octave steps. This will be used for imagesc graphing applications
totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));
%This then makes an array of the full octave steps I've made
octaveRange = zeros(totalOctaves + 1,2);
octaveRange(1,1) = uniqueFreqs(1);
for i = 1:totalOctaves
    octaveRange (i+1,1) = octaveRange(i,1)*2;
end
%next, I find the positions from uniqueFreqs that match octaveRange
for i = 1:size(octaveRange,1);
    octaveRange(i,2) = find(round(uniqueFreqs) == octaveRange(i,1));
end

%% Does the same for dBs. 
dbSteps = uniqueDBs(2) - uniqueDBs(1);
totalDBs = (uniqueDBs(end) - uniqueDBs(1))/dbSteps;
dbRange = zeros(totalDBs + 1,2);
dbRange(:,1) = uniqueDBs(1):dbSteps:uniqueDBs(end);
for i = 1:size(dbRange,1)
    dbRange(i,2) = find(uniqueDBs == dbRange(i,1));
end


%% Generates index that goes from low freq to high freq, and within each freq, goes from low to high amplitude
master(:,4) = 1:1:size(master,1);
master(:,5) = zeros;

sortingCounter = 1;

for i = 1:numFreqs
    for j = 1:numDBs
        sortingFinder = find(master(:,2) == uniqueFreqs(i) & master(:,3) == uniqueDBs(j));
        master(sortingFinder,5) = sortingCounter:1:sortingCounter + size(sortingFinder,1) - 1;
        sortingCounter = sortingCounter + size(sortingFinder,1);
    end
end
%now, master(:,5) is the index if I want to sort rasters by frequency and
%amplitude

%% Generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for i =1:numFreqs
    freqNameHolder{i} = num2str(uniqueFreqs(i));
end

%% big for loop that goes through every trode available and performs analysis
for i = 1:numTrodes
    %this should extract spikes, put them into matcluststruct. Will also
    %calculate RPVs and overall firing rate, as well as average waveform.
    [matclustStruct, clusterSizer] = functionSpikeWaveExtraction(rpvTime,...
        i,matclustFiles,matclustStruct,truncatedNames,clusterWindow);
    %this code should pull out rasters and histograms of all trials,
    %compute average firing rate from pre-trigger period, and generate a
    %z-scored histogram
    [matclustStruct] = functionRasterHistExtraction(i,clusterSizer,...
    master,baselineBins,numDBs,numFreqs,matclustStruct,truncatedNames,...
    rasterWindow,histBin,histBinVector);
    %This code will calculate histograms and rasters for every frequency
    %amplitude pair. It will also extract information about response
    %reliability. Also gets frequency specific histograms for plotting
    %purposes.
    [matclustStruct] = functionFreqAmpRasterHist(i, clusterSizer,...
    matclustStruct,histBinNum,histBinVector,histBin,truncatedNames);
end

%% graphing!
[matclustStruct] = functionTuningPlot(numTrodes,matclustStruct,...
    truncatedNames,rpvTime,clusterWindow,rasterWindow,histBinVector,...
    clims1,numFreqs,numDBs,uniqueFreqs,uniqueDBs,soundFile,octaveRange,...
    dbRange,fileName,trodesDesignation);

save(fullfile(pname,fname),'matclustStruct');

end




