%This code will be for extraction of spikes, calculation of RPVs, and
%calculation of overall firing rate. This code will also extract the
%average waveform of the cluster.
%Addendum 170105 Adding code to eliminate duplicate units.

function [s, truncatedNames] = functionMatclustExtraction(rpvTime,...
    matclustFiles,s,clusterWindow,subFoldersCell);

%Parameters for selection of duplicates
dupSelectLim = 0.20;
spikeHolder = cell(0,0);
expPower = 4;
rpvPercentCutoff = 5; %percent cutoff for RPVs. in %, so this is actually like 0.05
% rpvSpikeCutoff = 100; %raw spike number cutoff for RPVs

artifactJitter = 0.0005; %in seconds!
artifactCutoff = 0.7; %ratio cutoff for laser aligned spikes. 

%check if there is a toggleRPV variable. If present, respect it. If not,
%then set to 1 (will use RPV to eliminate units)
toggleCheck = isfield(s.Parameters,'toggleRPV');

if toggleCheck == 0
    toggleRPV = 1;
    disp('NO TOGGLE DETECTED DEFAULTING TO RPV ELIMINATION')
else
    toggleRPV = s.Parameters.toggleRPV;
    if toggleRPV == 1
        disp('TOGGLE DETECTED, ELIMINATING BY RPV')
    else
        disp('TOGGLE DETECTED, NOT ELIMINATING UNITS')
    end
end

%pull names for parameter files
[paramFiles] = functionFileFinder(subFoldersCell,'matclust','param');

%first, lets see if there was any laser playback. If there was, this
%initiates code for looking at laser alignment. 

homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
if ismac
    subFoldersCell = strsplit(subFolders,':')';
elseif ispc
    subFoldersCell = strsplit(subFolders,';')';
end

%finds DIO folder, extracts D2 specifically

[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
if length(D2FileName) == 0
    [D2FileName] = functionFileFinder(subFoldersCell,'DIO','Din2');
end

D2FileName = D2FileName{1};

%extracts DIO stuffs!
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%extracts port states and times of changes
dio2State = double(DIO2Data.fields(2).data);
dio2Time = double(DIO2Data.fields(1).data)/s.Parameters.trodesFS;

%if there isnt anything coming through, this will only have one value.
laserSwitch = 0;
if length(dio2State) > 1
    laserSwitch = 1;
    %pulls the times when state goes up/down
    upTimes = dio2Time(dio2State == 1);
    downTimes = dio2Time(dio2State == 0);

    %we want the duration of the pulse. Therefore, we need to determine if
    %upTimes == downTimes. Because the program always begins on a zero, lets
    %use this a system to eliminate extras. We'll find the first up, and remove
    %all downs before it
    firstUp = upTimes(1);
    downTimes(downTimes<firstUp) = [];

    %check!
    if length(upTimes) == length(downTimes)
        disp('Up Times and Down Times Matched')
    else
        error(strcat('Mismatched Pulses Up (',num2str(length(upTimes)),') vs Down(',num2str(length(downTimes)),')'))
    end

    %pulls the duration of the pulse (this should give laser pulse duration
    dioDiff = downTimes - upTimes;
    meanDiff = mean(dioDiff);
    windowArtifact = [-0.1,meanDiff*2];
    s.LaserData.LaserStartTimes = upTimes;
    s.LaserData.LaserEndTimes = downTimes;
    s.LaserData.LaserDuration = meanDiff;
end

%The code below extracts out the names of the nTrodes and generates both
%names for nTrode-Cluster combinations, as well as an array for indexing
%things later on. This also generates a holder inside the s
%structured array for storage of data relating to that structure. 

%first, resort file names into normal order
matclustFiles = sort_nat(matclustFiles);
truncatedNames = matclustFiles;
numTrodes = length(truncatedNames);
%sets asside arrays and counters
desigNames = cell(0,0);
desigArray = [];
desigCounter = 1;
timeFilterHolder = []; %this will store first and last time values

for clusterCount = 1:numTrodes
    %extracts the "ntX" part of the nTrode name.
    truncatedNames{clusterCount} = truncatedNames{clusterCount}(16:find(truncatedNames{clusterCount} == '.')-1);
    %opens matclust file and extracts cluster numbers. 
    matclustFile = open(matclustFiles{clusterCount});
    timeFilterSize = size(matclustFile.clustdata.timefilterranges,1);
    if timeFilterSize > 1;
        timeFilterHolder(clusterCount,:) = matclustFile.clustdata.timefilterranges(end,:);
    else
       timeFilterHolder(clusterCount,:) = matclustFile.clustdata.timefilterranges;
    end
    clusterLength = length(matclustFile.clustattrib.clustersOn);
    %based on number of clusters, generates for loop to fit generate names
    %of nTrode-cluster pairs, and designates space in s. Also
    %fills out array, which will be used for calling data from matclust
    %file.
    for nameCount = 1:clusterLength
        desigNames{desigCounter} = strcat(truncatedNames{clusterCount},'cluster',num2str(matclustFile.clustattrib.clustersOn(nameCount)));
        
        desigArray(desigCounter,1) = str2double(truncatedNames{clusterCount}(3:end));
        desigArray(desigCounter,2) = nameCount;
        %now I do my conventional stuff for analyzing spike timing
        %pull all spike indices
        clusterSpikes = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(nameCount)}.index;
        %find spike times.
        spikeTimes = matclustFile.clustdata.params(clusterSpikes,1);
        spikeHolder{desigCounter} = round(spikeTimes*(10^expPower))/(10^expPower); %stores spikes for separating out duplicate units.
        %find difference in time between spikes. 
        diffSpikes = diff(spikeTimes);
        %based on differences, calculate RPV number and percentage
        rpvNumber = size(diffSpikes(diffSpikes<rpvTime),1);
        totalNumber = size(clusterSpikes,1);
        rpvPercent = rpvNumber/totalNumber*100;
        %for graphing purposes, save a small subset of spikes to show RPV
        selectedSpikes = diffSpikes(diffSpikes<clusterWindow(2));
        %now lets remove things with shitty RPVs and/or laser response
        if laserSwitch == 0
            if rpvPercent > rpvPercentCutoff && toggleRPV == 1
                disp(strcat('Unit',desigNames{desigCounter},'Fails RPVs with:',num2str(rpvPercent),'% violations out of',num2str(totalNumber),' Spikes'))
                desigNames(desigCounter) = [];
                desigArray(desigCounter,:) = [];
            else
                %calculate overall firing rate
                totalTime = matclustFile.clustdata.datarange(2,1)-matclustFile.clustdata.datarange(1,1);
                overallFiring = totalNumber/totalTime;
                %calculate average waveform with standard error.
                %pull waves and squeeze to a 40 x X array
                targetWaveName = strcat('waves_',truncatedNames{clusterCount},'.mat');
                waveLoader = open(targetWaveName);
                %extract relevant waves
                waveHolder = waveLoader.waves(:,:,clusterSpikes);
                %pull out average and standard error
                averageWaveHolder = mean(waveHolder,3);
                nTrodeWidth(desigCounter) = size(averageWaveHolder,2);
                %store into structured array!
                s.(desigNames{desigCounter}) = [];
                s.(desigNames{desigCounter}).ISIGraph = selectedSpikes;
                s.(desigNames{desigCounter}).RPVNumber = rpvNumber;
                s.(desigNames{desigCounter}).RPVPercent = rpvPercent;
                s.(desigNames{desigCounter}).TotalSpikeNumber = totalNumber;
                s.(desigNames{desigCounter}).SpikeTimes = spikeTimes;
                s.(desigNames{desigCounter}).OverallFiringRate = overallFiring;
                s.(desigNames{desigCounter}).TotalTimeRecording = totalTime;

                s.(desigNames{desigCounter}).AverageWaveForms = averageWaveHolder;
                s.(desigNames{desigCounter}).AllWaves = waveHolder;
                desigCounter = desigCounter + 1;
            end
        elseif laserSwitch == 1
            %find laser aligned spikes
            [rasters] = functionBasicRaster(spikeTimes,upTimes,windowArtifact);
            %find all spikes within 0.5 ms window of onset, offset
            spikesOnset = rasters(rasters(:,1) < 0+artifactJitter & rasters(:,1) > 0-artifactJitter,1);
            spikesOffset = rasters(rasters(:,1) < meanDiff+artifactJitter & rasters(:,1) > meanDiff-artifactJitter,1);

            numSpikes = length(spikeTimes);
            propOnset = length(spikesOnset)/numSpikes;
            propOffset = length(spikesOffset)/numSpikes;
            if rpvPercent > rpvPercentCutoff && toggleRPV == 1
                disp(strcat('Unit',desigNames{desigCounter},'Fails RPVs with:',num2str(rpvPercent),'% violations out of',num2str(totalNumber),' Spikes'))
                desigNames(desigCounter) = [];
                desigArray(desigCounter,:) = [];
            elseif propOnset > artifactCutoff 
                disp(strcat('Unit',desigNames{desigCounter},'Fails LaserOnset with:',num2str(propOnset),'laser aligned spikes'))
                desigNames(desigCounter) = [];
                desigArray(desigCounter,:) = [];
            elseif propOffset > artifactCutoff
                disp(strcat('Unit',desigNames{desigCounter},'Fails LaserOffset with:',num2str(propOffset),'laser aligned spikes'))
                desigNames(desigCounter) = [];
                desigArray(desigCounter,:) = [];
            else
                %calculate overall firing rate
                totalTime = matclustFile.clustdata.datarange(2,1)-matclustFile.clustdata.datarange(1,1);
                overallFiring = totalNumber/totalTime;
                %calculate average waveform with standard error.
                %pull waves and squeeze to a 40 x X array
                targetWaveName = strcat('waves_',truncatedNames{clusterCount},'.mat');
                waveLoader = open(targetWaveName);
                %extract relevant waves
                waveHolder = waveLoader.waves(:,:,clusterSpikes);
                %pull out average and standard error
                averageWaveHolder = mean(waveHolder,3);
                nTrodeWidth(desigCounter) = size(averageWaveHolder,2);
                %store into structured array!
                s.(desigNames{desigCounter}) = [];
                s.(desigNames{desigCounter}).ISIGraph = selectedSpikes;
                s.(desigNames{desigCounter}).RPVNumber = rpvNumber;
                s.(desigNames{desigCounter}).RPVPercent = rpvPercent;
                s.(desigNames{desigCounter}).TotalSpikeNumber = totalNumber;
                s.(desigNames{desigCounter}).SpikeTimes = spikeTimes;
                s.(desigNames{desigCounter}).OverallFiringRate = overallFiring;
                s.(desigNames{desigCounter}).TotalTimeRecording = totalTime;

                s.(desigNames{desigCounter}).AverageWaveForms = averageWaveHolder;
                s.(desigNames{desigCounter}).AllWaves = waveHolder;
                desigCounter = desigCounter + 1;
            end
        end
    end
end

%% Now lets put in some code to determine shanks etc. 
% first determine number of total nTrodes.
numberTrodes = length(paramFiles)-length(matclustFiles);

trueWidth = round(mean(nTrodeWidth));

if ~ismember(numberTrodes,[4,8,16,32]);
    disp(strcat('WARNING: UNCONVENTIONAL NUMBER OF nTRODES DETECTED',num2str(numberTrodes)))
    %extract params file names
    testParams = paramFiles(length(matclustFiles)+1:end);
    %extract numbers from these names
    for i = 1:length(testParams)
        nameHold = testParams{i};
        ntFind = strfind(nameHold,'nt');
        dotFind = strfind(nameHold,'.');
        numParam(i) = str2num(nameHold(ntFind+2:dotFind-1));
    end
    %find maximum value
    paramMax = max(numParam);
    if ismember(paramMax,[4,8,16,32]);
        disp('Exact Match for Parameter Files')
        numberTrodes = paramMax;
    else
        disp('No Exact Match for Parameter Files')
        numSub = [4,8,16,32];
        testSub = numSub - paramMax;
        findBig = find(testSub > 0,1,'first');
        numberTrodes = numSub(findBig);
    end
end


%since I can assume tetrode configurations, I want to generate a map of the
%electrodes. These should be standard.

%16 channel map
chanMap16 = zeros(16,3);
chanMap16(:,1) = [1:1:16];
chanMap16(:,2) = floor((chanMap16(:,1)-1)/4)+1;
chanMap16(:,3) = 1;

%32 channel map
chanMap32 = zeros(32,3);
chanMap32(:,1) = [1:1:32];
chanMap32(:,2) = floor((chanMap32(:,1)-1)/4)+1;
chanMap32(:,3) = floor((chanMap32(:,1)-1)/16)+1;

%64 channel map
chanMap64 = zeros(64,3);
chanMap64(:,1) = [1:1:64];
chanMap64(:,2) = floor((chanMap64(:,1)-1)/4)+1;
chanMap64(:,3) = floor((chanMap64(:,1)-1)/32)+1;

if trueWidth == 4 %tetrodes
    disp('TETRODE CONFIGURATION')
    if numberTrodes == 4;
        %this means i'm on the 16 channel single shank
        shanks = 1;
        probeMap = chanMap16;
    elseif numberTrodes == 8;
        %this means i'm on the 32 channel double shank
        shanks = 2;
        probeMap = chanMap32;
    elseif numberTrodes == 16;
        %this means i'm on the 64 channel double shank
        shanks = 2;
        probeMap = chanMap64;
    end
else trueWidth == 1 %singles
    disp('SINGLE ELECTRODE CONFIGURATION')
    if numberTrodes == 16;
        %this means i'm on the 16 channel single shank
        shanks = 1;
        probeMap = chanMap16;
        probeMap(:,2) = [1:1:16];
    elseif numberTrodes == 32;
        %this means i'm on the 32 channel double shank
        shanks = 2;
        probeMap = chanMap32;
        probeMap(:,2) = [1:1:32];
    elseif numberTrodes == 64;
        %this means i'm on the 64 channel double shank
        shanks = 2;
        probeMap = chanMap64;
        probeMap(:,2) = [1:1:64];
    end
end

%generate shank designation for duplicate elimination
shankDesig = zeros(numberTrodes/shanks,shanks);
holderShank = 1;
for shnkInd = 1:shanks
    shankDesig(:,shnkInd) = [holderShank:holderShank+numberTrodes/shanks-1];
    holderShank = holderShank + numberTrodes/shanks;
end


s.DesignationArray = desigArray;
s.DesignationName = desigNames;
firstPoint = min(timeFilterHolder(:,1));
lastPoint = max(timeFilterHolder(:,2));
s.TimeFilterRange = [firstPoint lastPoint];
s.ShankMap = probeMap;
s.ShankDesignation = shankDesig;
s.Shanks = shanks;



end

