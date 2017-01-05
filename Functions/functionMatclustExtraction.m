%This code will be for extraction of spikes, calculation of RPVs, and
%calculation of overall firing rate. This code will also extract the
%average waveform of the cluster.

function [s, truncatedNames] = functionMatclustExtraction(rpvTime,...
    matclustFiles,s,clusterWindow);

%The code below extracts out the names of the nTrodes and generates both
%names for nTrode-Cluster combinations, as well as an array for indexing
%things later on. This also generates a holder inside the s
%structured array for storage of data relating to that structure. 
truncatedNames = matclustFiles;
numTrodes = length(truncatedNames);
%sets asside arrays and counters
desigNames = cell(0,0);
desigArray = [];
desigCounter = 1;
for clusterCount = 1:numTrodes
    %extracts the "ntX" part of the nTrode name.
    truncatedNames{clusterCount} = truncatedNames{clusterCount}(16:find(truncatedNames{clusterCount} == '.')-1);
    %opens matclust file and extracts cluster numbers. 
    matclustFile = open(matclustFiles{clusterCount});
    clusterLength = length(matclustFile.clustattrib.clustersOn);
    %based on number of clusters, generates for loop to fit generate names
    %of nTrode-cluster pairs, and designates space in s. Also
    %fills out array, which will be used for calling data from matclust
    %file.
    for nameCount = 1:clusterLength
        desigNames{desigCounter} = strcat(truncatedNames{clusterCount},'cluster',num2str(matclustFile.clustattrib.clustersOn(nameCount)));
        s.(desigNames{desigCounter}) = [];
        desigArray(desigCounter,1) = clusterCount;
        desigArray(desigCounter,2) = nameCount;
        %now I do my conventional stuff for analyzing spike timing
        %pull all spike indices
        clusterSpikes = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(nameCount)}.index;
        %find spike times.
        spikeTimes = matclustFile.clustdata.params(clusterSpikes,1);
        %find difference in time between spikes. 
        diffSpikes = diff(spikeTimes);
        %based on differences, calculate RPV number and percentage
        rpvNumber = size(diffSpikes(diffSpikes<rpvTime),1);
        totalNumber = size(clusterSpikes,1);
        rpvPercent = rpvNumber/totalNumber*100;
        %for graphing purposes, save a small subset of spikes to show RPV
        selectedSpikes = diffSpikes(diffSpikes<clusterWindow(2));
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

s.DesignationArray = desigArray;
s.DesignationName = desigNames;





end

