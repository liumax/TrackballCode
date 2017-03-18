%This code will be for extraction of spikes, calculation of RPVs, and
%calculation of overall firing rate. This code will also extract the
%average waveform of the cluster.
%Addendum 170105 Adding code to eliminate duplicate units.

function [s, truncatedNames] = functionMatclustExtraction(rpvTime,...
    matclustFiles,s,clusterWindow);

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

%first, lets see if there was any laser playback. If there was, this
%initiates code for looking at laser alignment. 

homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%finds DIO folder, extracts D2 specifically
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
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
        
        desigArray(desigCounter,1) = clusterCount;
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

% %now that things are all squared away, lets eliminate duplicates!
% c = nchoosek([1:length(desigNames)],2);
% compHolder = zeros(length(c),3);
% 
% for clusterCounter = 1:length(c)
%     compHolder(clusterCounter,1) = length(intersect(spikeHolder{c(clusterCounter,1)},spikeHolder{c(clusterCounter,2)}));
%     compHolder(clusterCounter,2) = compHolder(clusterCounter,1)/length(spikeHolder{c(clusterCounter,1)});
%     compHolder(clusterCounter,3) = compHolder(clusterCounter,1)/length(spikeHolder{c(clusterCounter,2)});
% end
% 
% %find duplicates with cutoff set above. 
% dupTargets = find(compHolder(:,2) > dupSelectLim & compHolder(:,3) > dupSelectLim);
% 
% %find the true indices of the units, choose which to eliminate.
% while length(dupTargets)>0
%     targetInds = c(dupTargets(1),:);
%     waves1 = mean(max(s.(desigNames{targetInds(1)}).AverageWaveForms));
%     waves2 = mean(max(s.(desigNames{targetInds(2)}).AverageWaveForms));
%     if waves1 > waves2
%         disp(strcat(desigNames{targetInds(1)},'>',desigNames{targetInds(2)},'Saving Former'))
%         s = rmfield(s,desigNames{targetInds(2)});
%         desigNames(targetInds(2)) = [];
%         desigArray(targetInds(2),:) = [];
%     elseif waves2 > waves1
%         disp(strcat(desigNames{targetInds(1)},'<',desigNames{targetInds(2)},'Saving Latter'))
%         s = rmfield(s,desigNames{targetInds(1)});
%         desigNames(targetInds(1)) = [];
%         desigArray(targetInds(1),:) = [];
%     else
%         disp(strcat(desigNames{targetInds(1)},' and_',desigNames{targetInds(2)},'Equal, Saving Both'))
%     end
%     
%     c = nchoosek([1:length(desigNames)],2);
%     compHolder = zeros(length(c),3);
%     for clusterCounter = 1:length(c)
%         compHolder(clusterCounter,1) = length(intersect(spikeHolder{c(clusterCounter,1)},spikeHolder{c(clusterCounter,2)}));
%         compHolder(clusterCounter,2) = compHolder(clusterCounter,1)/length(spikeHolder{c(clusterCounter,1)});
%         compHolder(clusterCounter,3) = compHolder(clusterCounter,1)/length(spikeHolder{c(clusterCounter,2)});
%     end
%     %find duplicates with cutoff set above. 
%     dupTargets = find(compHolder(:,2) > dupSelectLim & compHolder(:,3) > dupSelectLim);
%     
% end



s.DesignationArray = desigArray;
s.DesignationName = desigNames;
firstPoint = min(timeFilterHolder(:,1));
lastPoint = max(timeFilterHolder(:,2));
s.TimeFilterRange = [firstPoint lastPoint];


end

