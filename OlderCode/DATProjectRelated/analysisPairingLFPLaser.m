function [] = analysisPairingLFPLaser(fileName);

lfpWindow = [-1 2];
pulseTiming = 0.05;

%% sets up file saving stuff
saveName = strcat(fileName,'LFPAnalysis','.mat');
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
matclustStruct.NumberTrodes = numTrodes;

%%

%find DIO folder and D1 file for analysis
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};

[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D2FileName);
[DIOData1] = readTrodesExtractedDataFile(D1FileName);

%extracts port states and times of changes
dioState = double(DIOData.fields(2).data);
dioTime = double(DIOData.fields(1).data);

inTimes = dioTime(dioState == 1)/30000;
inTimesSpacing = diff(inTimes);
finder = find(inTimesSpacing > pulseTiming);
finder = [1;finder+1];
inTimes = inTimes(finder);


dioState = [];
dioTime = [];  

%%Calculate LFP
disp('Starting LFP Analysis')
%this code picks out the LFP folder and moves to that folder
lfpFinder =dir;
lfpFinder = {lfpFinder.name};
lfpIndex = strfind(lfpFinder,'LFP');
lfpIndex = find(not(cellfun('isempty', lfpIndex)));
lfpFolder = lfpFinder{lfpIndex};
newDir = strcat(pwd,'\',lfpFolder);
cd(newDir)
%this code pulls the file names of the LFP files, and reorders them in
%natural order (1,2,10, not 1,10,2)
lfpFinder = dir;
lfpFinder = {lfpFinder.name};
lfpIndex = strfind(lfpFinder,'LFP');
lfpIndex = find(not(cellfun('isempty', lfpIndex)));
lfpFiles = cell(size(lfpIndex,2),1);
lfpFiles= lfpFinder(lfpIndex);
[cs index] = sort_nat(lfpFiles);
lfpFiles = cs;
numLFPs = size(lfpFiles,2);

%opens LFP file
lfp = readTrodesExtractedDataFile(lfpFiles{1});
viewSamples = round((lfpWindow(2)-lfpWindow(1))*lfp.clockrate/lfp.decimation);

lfpAverages = zeros(numLFPs,viewSamples);

for i = 1:numLFPs
    %opens LFP file
    lfp = readTrodesExtractedDataFile(lfpFiles{i});
    %counts number of LFP samples
    lfpSamples = size(lfp.fields.data,1);
    %makes LFP time points based on # of samples and decimation
    %adjust time to actual time (to match master)
    lfpTimes = (((0:1:lfpSamples-1)*lfp.decimation)+lfp.first_timestamp)'/30000;
    lfpSignals = lfp.fields.data;
    
    %generates a temporary holder for LFP traces
    tempHolder = zeros(length(inTimes),viewSamples);
    for j = 1:1:length(inTimes)
        finder = find(lfpTimes > inTimes(j)+lfpWindow(1),1);
        tempHolder(j,:) = lfpSignals(finder:finder + viewSamples - 1);
    end
    lfpAverages(i,:) = mean(tempHolder);
    disp('Electrode')
    disp(i)
end

minLFP = min(min(lfpAverages));
maxLFP = max(max(lfpAverages));

%calculates the zero point based on number of samples. This is to generate
%lines that indicate tone duration
totalLFPWindow = lfpWindow(2)-lfpWindow(1);
lfpZero = abs(lfpWindow(1))*viewSamples/totalLFPWindow;

%save everything to a dummy structured array
s.Means = lfpAverages;
s.AnalysisWindow = lfpWindow;
s.ToneStartTime = lfpZero;

h = figure;
set(h, 'Position', [100 100 1000 1000])
%generates text box with information!
mTextBox = uicontrol('style','text');
descr = {'LFP By Frequency and DB';
    'File:';
    fileName;
    'Window (s):';
    lfpWindow};
set(mTextBox,'String',descr);
set(mTextBox,'Units','normalized');
set(mTextBox,'Position',[0.36,0.74,0.3,0.2])


for i = 1:numLFPs
    subplot(numLFPs,3,1+(3*(i-1)))
    plot([lfpZero lfpZero],[minLFP maxLFP],'k');
    hold on
    plot(squeeze(lfpAverages(j,:)));
    
    xlim([0 viewSamples])
    ylim([minLFP maxLFP])
    set(gca, 'XTick', [], 'YTick', [])
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1 1.4]) % stretch its width and height
end









